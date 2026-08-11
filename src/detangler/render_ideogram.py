"""The chromosome panel, as SVG and as interactive HTML."""
from __future__ import annotations

import argparse
import glob
import gzip
import json
import math
import os
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field, asdict
from typing import Dict, Iterable, List, Optional, Sequence as Seq, Set, Tuple

try:  # optional, only affects config file format
    import yaml  # type: ignore

    HAVE_YAML = True
except Exception:  # pragma: no cover
    yaml = None  # type: ignore
    HAVE_YAML = False
from .common import (
    figure_bp,
    PALETTE,
    TANGLE_LABEL,
    TANGLE_STYLE,
    esc,
    human_bp,
    median,
    wrap_text,
)
from .palette import (
    SEGMENT_COLOURS,
    _segment_number,
    _text_on,
)
from .model import (
    Model,
    _confidence,
)
from .render_common import (
    GAP,
    drawn_length_px,
    BAR_W,
    COV_W,
    FS_ANNOT,
    FS_HEADING,
    FS_PRIMARY,
    KEY_W,
    Layout,
    MARGIN_L,
    MAX_BAR_H,
    _annotation_colour,
    _bar_path,
    MIN_BLOCK_PX,
    PANEL_W,
)



def _block_bounds(lay, name, blocks, floor):
    """
    Drawn y-boundaries for one molecule's blocks.

    A block that is a real part of a chromosome must be VISIBLE, or the figure
    asserts a join and then hides the thing that makes it one chromosome - a
    51 kb centromere inside a 12 Mb chromosome is 0.4% of the bar, about two
    pixels. So any block below `floor` is raised to it.

    The pixels are stolen from the larger blocks in the same bar, in proportion
    to their size, so the molecule's OVERALL height does not change. That
    matters: the bars are sorted and read against each other by height, and
    inflating one bar to fit its small pieces would corrupt that ranking - the
    same failure that moving end-caps above the bar caused earlier.
    """
    y0 = lay.y(name, blocks[0][0])
    nat = [lay.y(name, b[1]) - lay.y(name, b[0]) for b in blocks]
    if sum(nat) <= 0:
        return [(y0, y0) for _ in blocks]
    small = [i for i, h in enumerate(nat) if h < floor]
    if not small:
        return _stack(y0, nat)
    # Floors cannot cost more than the bar has to give. If they would, every
    # floor shrinks together rather than one block swallowing the molecule.
    want = sum(floor - nat[i] for i in small)
    spare = sum(nat[i] for i in range(len(nat)) if i not in small)
    budget = min(want, spare * 0.6)
    scale = (budget / want) if want else 0.0
    out = list(nat)
    for i in small:
        out[i] = nat[i] + (floor - nat[i]) * scale
    if spare > 0 and budget > 0:
        for i in range(len(nat)):
            if i not in small:
                out[i] = nat[i] - budget * (nat[i] / spare)
    return _stack(y0, out)


def _stack(y0, heights):
    bounds, y = [], y0
    for h in heights:
        bounds.append((y, y + h))
        y += h
    return bounds


def ideogram_geometry(model: Model) -> Tuple[Layout, List[str], bool, float]:
    """
    The ideogram's layout, header text and header height. A thin view over
    _ideogram_frame so there is only ever one layout calculation: if the
    renderer and the paired figure disagreed by a pixel, every leader line
    would point to the wrong place.
    """
    f = _ideogram_frame(model)
    return f["lay"], f["head_lines"], f["show_cov"], f["header_h"]  # type: ignore


def ideogram_block_anchors(model: Model) -> Dict[str, Tuple[float, float]]:
    """Left edge and vertical centre of each segment block, in ideogram coordinates."""
    lay, _, _, _ = ideogram_geometry(model)
    out: Dict[str, Tuple[float, float]] = {}
    for s in lay.order:
        for b_start, b_end, seg, _colour in s.blocks:
            y = (lay.y(s.name, b_start) + lay.y(s.name, b_end)) / 2.0
            out.setdefault(seg, (lay.x[s.name], y))
    return out


def _reach_above_baseline(s, lay: Layout) -> float:
    """
    How far above the shared baseline a single molecule puts ink.

    A circular molecule is drawn as a ring rather than as its notional bar, and
    the ring is wider than the tiny bar an organelle would otherwise get, so
    asking the bar how tall it is would understate it and let it collide with
    the heading.
    """
    if getattr(s, "circular", False):
        return BAR_W * 1.9 + FS_ANNOT + 14
    # Caps are drawn inside the bar, so a molecule reaches exactly its own
    # height - nothing sits above it.
    return lay.height[s.name]


def _headline(model: Model) -> str:
    """
    The answer, in one line, under the panel heading.

    Reading the count off the bars means counting bars, which is exactly the
    thing the tool is supposed to have already done - and it silently drops the
    uncertainty, because a figure of five bars cannot show that three would also
    fit the evidence.
    """
    drawable = model.drawable()
    circular = [s for s in drawable if getattr(s, "circular", False)]
    rng = model.count_range()
    if rng:
        best, low, high = rng
        text = f"{best} linear molecule{'s' if best != 1 else ''}"
        # A collapsed range is not a range; printing "(range 5-5)" would make the
        # figure look hedged where the evidence actually pins the number down.
        if low != high:
            text += f" (range {low}-{high})"
    else:
        n = len(drawable) - len(circular)
        text = f"{n} linear molecule{'s' if n != 1 else ''}"
    if circular:
        text += f" + {len(circular)} circular"
    unassigned = model.unassigned()
    if unassigned:
        text += f"; {human_bp(sum(s.length for s in unassigned))} unassigned"
    return text


def _doubt_badge(cx: float, cy: float) -> str:
    """
    A neutral marker for "the depth does not agree with what is drawn here".

    It used to be carried by the cap itself, drawn faint and dashed. That made
    one contig two different colours depending on which claim it supported, and
    the figure's whole promise is that a contig looks the same everywhere. The
    doubt is real and has to stay on the page, so it moves off the colour and
    onto a badge that is deliberately not in any contig's palette.
    """
    r = 11.0
    return (
        f'<g class="doubt">'
        f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r:.1f}" fill="{PALETTE["muted"]}" '
        f'stroke="#ffffff" stroke-width="1.6"/>'
        f'<text x="{cx:.1f}" y="{cy + r * 0.43:.1f}" font-size="{r * 1.3:.0f}" '
        f'text-anchor="middle" fill="#ffffff" font-weight="700">!</text>'
        f'</g>'
    )


def _ideogram_frame(model: Model) -> Dict[str, object]:
    """
    Geometry for the chromosome figure, computed once and shared by the renderer
    and by anything that needs to point at a block from outside - the paired
    figure draws leader lines to these exact coordinates.
    """
    show_cov = bool(model.coverage) and model.settings.get("coverage", True)
    bar_h = float(getattr(model, 'max_bar_h', 0) or MAX_BAR_H)
    probe = Layout(model, show_cov, max_bar_h=bar_h)
    drawn = {s.name for s in probe.order}

    # v9: no summary paragraph under the panel title. The reasoning belongs in the
    # report; the figure carries only what points at something it draws.
    head_lines: List[str] = []
    # Space above the baseline is claimed by whichever object reaches highest,
    # not by a global allowance. Once the molecules stand on a baseline a short
    # one with a tall cap stack no longer approaches the heading, so reserving
    # the worst cap stack for everybody would open a band of white above every
    # bar that nothing ever occupies.
    reach = max(
        (_reach_above_baseline(s, probe) for s in probe.order),
        default=bar_h,
    )
    # The unassigned column shares the same baseline, so it can crowd the heading
    # in its own right - it is floored at a minimum drawn length, not scaled down.
    reach = max(
        [reach] + [drawn_length_px(s.length, probe.scale) for s in model.unassigned()]
    )
    # 54 clears the panel title. There is no headline line under it any more.
    header_h = 54 + 15 * (probe.max_label_lines - 1) + max(0.0, reach - bar_h)
    # Caps are side tabs now, so nothing extends below the baseline except the
    # size label itself.
    bottom_caps = 0.0

    lay = Layout(model, show_cov, header_h, max_bar_h=bar_h)
    legend_svg, legend_bottom = _legend_svg(model, lay)
    total_h = max(legend_bottom + 26, lay.baseline + bottom_caps + 90)
    return {
        "lay": lay,
        "head_lines": head_lines,
        "drawn": drawn,
        "show_cov": show_cov,
        "header_h": header_h,
        "legend_svg": legend_svg,
        "total_h": total_h,
    }


def render_svg(model: Model, interactive: bool = False) -> str:
    frame = _ideogram_frame(model)
    lay, head_lines, drawn = frame["lay"], frame["head_lines"], frame["drawn"]
    show_cov, header_h = frame["show_cov"], frame["header_h"]
    legend_svg, total_h = frame["legend_svg"], frame["total_h"]
    P: List[str] = []
    add = P.append

    add(
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{lay.width:.0f}" '
        f'height="{total_h:.0f}" viewBox="0 0 {lay.width:.0f} {total_h:.0f}" '
        f'font-family="Helvetica, Arial, sans-serif">'
    )
    add(
        '<defs>'
        '<pattern id="nts" width="6" height="6" patternUnits="userSpaceOnUse" '
        'patternTransform="rotate(45)">'
        '<rect width="6" height="6" fill="#ffffff" fill-opacity="0"/>'
        '<line x1="0" y1="0" x2="0" y2="6" stroke="#ffffff" stroke-width="2" stroke-opacity="0.55"/>'
        '</pattern>'
        '</defs>'
    )
    add(f'<rect width="100%" height="100%" fill="{PALETTE["bg"]}"/>')

    # ---- title ----
    add(
        f'<text x="{MARGIN_L}" y="40" font-size="{FS_HEADING}" font-weight="600" '
        f'fill="{PALETTE["text"]}">{esc(model.title)}</text>'
    )

    # The result as a sentence fragment, sitting with the heading. In the paired
    # figure the heading is drawn by the combined canvas and this panel's own
    # title is blanked, so the line moves up into the space the title vacated
    # rather than floating a heading's height below nothing.

    # No scale ruler. Each molecule carries its own size label, and a bare
    # axis down the left invited reading positions off it that a hypothesis
    # does not actually claim.

    # ---- tangle arcs (behind the bars) ----
    add('<g id="layer-tangles">')
    # v9 drops the arcs joining one chromosome to another. They were the biggest
    # source of visual noise and duplicated what the shared segment colours
    # already say: a repeat linking two chains is drawn as a copy on each. The
    # single-point markers that survived that cut are gone too: a nine-pixel
    # triangle beside a bar had no legend and no label, so it asserted that
    # something was there without ever saying what, and the reader had no way to
    # find out. The layer is kept because the interactive view toggles it.
    add("</g>")

    # ---- chromosome bars ----
    # Numbered by SIZE RANK, not by draw order. The bars now run smallest to
    # largest left to right, but chromosome 1 is the largest by near-universal
    # convention, so tying the number to position would have inverted it. The
    # numbers count down towards the right-hand end of the panel.
    _rank = {
        q.name: i + 1
        for i, q in enumerate(
            sorted((q for q in lay.order if q.role == "chromosome"),
                   key=lambda q: -q.length)
        )
    }
    add('<g id="layer-chromosomes">')
    for s in lay.order:
        x, top, h = lay.x[s.name], lay.top[s.name], lay.height[s.name]
        fill = PALETTE.get(s.role, PALETTE["chromosome"])
        rx = BAR_W / 2.0
        battrs = ""
        if interactive:
            battrs = (
                f' class="chrom" data-name="{esc(s.name)}" data-role="{esc(s.role)}"'
                f' data-length="{s.length}" data-depth="{s.depth if s.depth is not None else ""}"'
            )

        # v9: a circular molecule is drawn as a ring, not as a bar with a little
        # circle underneath it. Nothing circular is to scale against the nuclear
        # chromosomes anyway, so a ring is both truer and less misleading.
        if s.circular:
            # the segment's OWN colour, the one it has in the graph panel. Falling
            # back to the role colour here broke the figure's one promise: edge_11
            # came out cyan on the left and orange on the right.
            seg_colour = (
                s.blocks[0][3] if s.blocks
                else model.segment_colours.get(s.name, fill)
            )
            # an organelle is not on the nuclear scale, so it is drawn thinner
            # than a chromosome bar as well as round: nothing about it should
            # invite being read off the Mb axis
            ring_w = BAR_W * 0.45
            r = BAR_W * 0.95
            # The ring rests on the shared baseline like everything else. Hung
            # from its notional bar top instead, it would float at a height that
            # varies with a scale it is explicitly not drawn to.
            ccx, ccy = x + BAR_W / 2.0, lay.baseline - r

            add(
                f'<circle cx="{ccx:.1f}" cy="{ccy:.1f}" r="{r:.1f}" fill="none" '
                f'stroke="{seg_colour}" stroke-width="{ring_w:.1f}"/>'
            )
            # The segment NUMBER, on the ring itself. Every other molecule in the
            # panel carries the numbers of the contigs that compose it; without
            # one here the organelle was the only thing in the figure with no
            # way back to the graph panel except its colour.
            ring_seg = s.blocks[0][2] if s.blocks else s.name
            add(
                f'<text x="{ccx:.1f}" y="{ccy - r + FS_ANNOT * 0.36:.1f}" '
                f'font-size="{FS_ANNOT}" font-weight="700" text-anchor="middle" '
                f'fill="{_text_on(seg_colour)}">{esc(_segment_number(ring_seg))}</text>'
            )
            add(
                f'<text x="{ccx:.1f}" y="{lay.baseline + 34:.1f}" '
                f'font-size="{FS_ANNOT}" font-weight="700" text-anchor="middle" '
                f'fill="{PALETTE["text"]}">'
                f'{esc("mito" if s.role == "mitochondrion" else s.role)}</text>'
            )
            add(
                f'<text x="{ccx:.1f}" y="{lay.baseline + 34 + FS_ANNOT + 6:.1f}" '
                f'font-size="{FS_ANNOT}" text-anchor="middle" '
                f'fill="{PALETTE["muted"]}">{figure_bp(s.length)}</text>'
            )
            continue

        # No outline-and-fill. An outline only means something if it is drawn
        # from something OTHER than the contigs inside it - an expected size -
        # and the tool does not have that. Drawing the contigs' own extent as an
        # outline and then filling it makes a figure that is full by
        # construction, which looks like a finding and is not one.
        add(
            f'<path d="{_bar_path(x, top, BAR_W, h, rx, rx)}" '
            f'fill="{fill}" fill-opacity="0.82" stroke="none" '
            f'{battrs}/>'
        )

        # Segment blocks: the same colour the segment has in the graph figure, so
        # a node over there can be found on a chromosome over here.
        if s.blocks:
            add('<g class="blocks">')
            last = len(s.blocks) - 1
            inset = 0.0 if s.blocks_tile else 4.0
            bw = BAR_W - 2 * inset
            bounds = _block_bounds(lay, s.name, s.blocks, MIN_BLOCK_PX)
            for bi, (b_start, b_end, seg, colour) in enumerate(s.blocks):
                y1, y2 = bounds[bi]
                y2 = max(y2, y1 + 1.6)
                bl = ""
                if interactive:
                    bl = (
                        f' class="block" data-desc="{esc(seg)}, {human_bp(b_end - b_start)}, '
                        f'on {esc(s.display)}"'
                    )
                # Tiled blocks ARE the bar, so the outer ones keep its rounded
                # ends - UNLESS a cap sits against that end, in which case the
                # molecule continues and the corner must be square. Only the
                # outermost piece of the whole molecule is rounded.
                rt = rx if (s.blocks_tile and bi == 0) else 0
                rb = rx if (s.blocks_tile and bi == last) else 0
                add(
                    f'<path d="{_bar_path(x + inset, y1, bw, y2 - y1, rt, rb)}" '
                    f'fill="{colour}" fill-opacity="{0.95 if s.blocks_tile else 0.92}" '
                    f'stroke="none"{bl}/>'
                )
                # v9: the segment NUMBER, set inside the block, inked white or
                # dark for contrast against that block's own colour
                if y2 - y1 >= FS_PRIMARY + 4:
                    add(
                        f'<text x="{x + BAR_W / 2:.1f}" y="{(y1 + y2) / 2 + FS_PRIMARY * 0.35:.1f}" '
                        f'font-size="{FS_PRIMARY}" text-anchor="middle" fill="{_text_on(colour)}" '
                        f'font-weight="700">{esc(_segment_number(seg))}</text>'
                    )
            add("</g>")
        if s.name in lay.not_to_scale:
            add(
                f'<rect x="{x:.1f}" y="{top:.1f}" width="{BAR_W}" height="{h:.1f}" rx="{rx:.1f}" '
                f'ry="{rx:.1f}" fill="url(#nts)" stroke="none"/>'
            )

        # annotation bands
        for feat in model.annotations:
            if feat["seqname"] != s.name:
                continue
            y1 = lay.y(s.name, feat["start"])
            y2 = max(lay.y(s.name, feat["end"]), y1 + 1.2)
            c = _annotation_colour(str(feat.get("kind", "feature")))
            fattrs = (
                f' class="annot" data-desc="{esc(feat.get("name") or feat.get("kind"))} '
                f'{esc(s.name)}:{feat["start"]:,}-{feat["end"]:,}"'
                if interactive
                else ""
            )
            add(
                f'<rect x="{x:.1f}" y="{y1:.1f}" width="{BAR_W}" height="{y2 - y1:.1f}" '
                f'fill="{c}" fill-opacity="0.85" stroke="none"{fattrs}/>'
            )

        # coverage anomaly stripes
        for an in model.coverage_anomalies:
            if an["seqname"] != s.name:
                continue
            y1 = lay.y(s.name, an["start"])
            y2 = max(lay.y(s.name, an["end"]), y1 + 1.5)
            c = "#d62728" if an["kind"] == "high" else "#1f77b4"
            add(
                f'<rect x="{x - 5:.1f}" y="{y1:.1f}" width="4" height="{y2 - y1:.1f}" '
                f'fill="{c}" fill-opacity="0.9"/>'
            )

        # No second outline pass here. It existed to re-stroke the bar so the
        # annotation bands could not spill past its rounded ends, but the bars
        # lost their outline in v9, so it was a path with neither fill nor stroke
        # emitted once per molecule.

        # Repeats attached to a free end, drawn hanging OFF the bar rather than
        # inside it: they are not part of the molecule and not on the Mb scale,
        # but they are what tells you this end is a telomere or an rDNA block.
        # Flush against the bar, not floating beside it, so a molecule reads as
        # one object: cap, backbone, cap. Only the OUTER corner of the outermost
        # cap is rounded; every join between blocks is square so they abut.
        # Repeats are drawn IN LINE, at the end of the molecule they attach to,
        # and they EAT INTO the bar rather than extending it. Three positions
        # were tried: above the bar, which silently added ~2 Mb of apparent
        # length and made the tallest bar stop being the largest molecule; below
        # the baseline, which fixed the height but lost which end a repeat sits
        # on; and beside it, which kept both but read as unexplained chunks
        # floating next to the chromosome. In line is what a karyogram does, and
        # the molecule's total height stays exactly its true length. The repeat's
        # own share is overstated - it has to be, to be visible at all - so the
        # size it appears to occupy is not its size.
        # Repeats are drawn close to TRUE SCALE here, with a floor of just over
        # half the bar width - below that a rounded end cannot be drawn and the
        # molecule stops looking like a chromosome. So
        # a 15 kb repeat on a 9 Mb chromosome is a thin band rather than a fifth
        # of the molecule. The graph panel keeps its large floor, because there a
        # short contig has to be big enough to see which of its ends things join
        # to; here that question is already answered and the honest proportion
        # matters more. The two panels therefore agree on the big contigs, which
        # is what the eye tracks, and differ on the small ones by design.
        for side in ("top", "bottom"):
            entries = s.caps.get(side, [])
            stacked = 0.0
            for ci, (seg, colour, seg_len) in enumerate(entries):
                cap_h = max(seg_len * lay.scale, BAR_W * 0.62)
                cy = (top + stacked) if side == "top" else (
                    lay.baseline - stacked - cap_h
                )
                stacked += cap_h
                rt = rx if (side == "top" and ci == 0) else 0.0
                rb = rx if (side == "bottom" and ci == 0) else 0.0
                add(
                    f'<path d="{_bar_path(x, cy, BAR_W, cap_h, rt, rb)}" '
                    f'fill="{colour}" fill-opacity="0.95" stroke="none"/>'
                )
                add(
                    f'<text x="{x + BAR_W / 2:.1f}" y="{cy + cap_h / 2 + FS_ANNOT * 0.36:.1f}" '
                    f'font-size="{FS_ANNOT}" text-anchor="middle" fill="{_text_on(colour)}" '
                    f'font-weight="700">{esc(_segment_number(seg))}</text>'
                )
                # No badge here. The depth-versus-placement contradictions are
                # real and worth knowing, but a mark on nine of nine repeats
                # turned the chromosome figure into a warning display. It should
                # answer "does this assembly make sense" at a glance first; the
                # report states every contradiction in full, and the interactive
                # view is where they belong on the drawing.

        # size only. Chain headings are gone: which contigs belong together is
        # shown by the numbered blocks in the bar, not by a caption above it.
        n_top = len(s.caps.get("top", []))
        # Clear of the badge that may sit on the outer bottom corner: at the old
        # offset the badge landed on the last two characters of the size.
        # A name, then its size under it. The numbers inside the bar identify
        # CONTIGS; without a name for the molecule there was nothing to call a
        # chromosome by when talking about the figure.
        add(
            f'<text x="{x + BAR_W / 2:.1f}" y="{lay.baseline + 34:.1f}" '
            f'font-size="{FS_ANNOT}" font-weight="700" text-anchor="middle" '
            f'fill="{PALETTE["text"]}">chr {_rank.get(s.name, 0)}</text>'
        )
        add(
            f'<text x="{x + BAR_W / 2:.1f}" y="{lay.baseline + 34 + FS_ANNOT + 6:.1f}" '
            f'font-size="{FS_ANNOT}" text-anchor="middle" '
            f'fill="{PALETTE["muted"]}">{figure_bp(s.length)}</text>'
        )
    add("</g>")

    # ---- coverage track ----
    if show_cov:
        add('<g id="layer-coverage">')
        gm = median([v for v in model.coverage_median.values() if v]) or 1.0
        cap = gm * 2.0
        for s in lay.order:
            ws = model.coverage.get(s.name) or []
            if not ws:
                continue
            x0 = lay.x[s.name] + BAR_W + 16
            add(
                f'<line x1="{x0:.1f}" y1="{lay.top[s.name]:.1f}" x2="{x0:.1f}" '
                f'y2="{lay.top[s.name] + lay.height[s.name]:.1f}" stroke="{PALETTE["grid"]}"/>'
            )
            pts = []
            for w in ws:
                y = lay.y(s.name, (w.start + w.end) / 2.0)
                xx = x0 + min(w.depth / cap, 1.0) * COV_W
                pts.append(f"{xx:.1f},{y:.1f}")
            if pts:
                add(
                    f'<polyline points="{" ".join(pts)}" fill="none" stroke="#333333" '
                    f'stroke-width="1.3" stroke-opacity="0.9" stroke-linejoin="round"/>'
                )
            mid = x0 + min(gm / cap, 1.0) * COV_W
            add(
                f'<line x1="{mid:.1f}" y1="{lay.top[s.name]:.1f}" x2="{mid:.1f}" '
                f'y2="{lay.top[s.name] + lay.height[s.name]:.1f}" stroke="#999999" '
                f'stroke-dasharray="2 3"/>'
            )
        add(
            f'<text x="{MARGIN_L}" y="{lay.baseline + 52:.1f}" font-size="10" '
            f'fill="{PALETTE["muted"]}">Coverage track (right of each bar): 0 to 2x the genome '
            f'median ({gm:.0f}x); dashed line = median. Red/blue ticks left of a bar mark '
            f'depth outliers.</text>'
        )
        add("</g>")

    # ---- what was expected but not found ----
    # Told how many chromosomes to expect, the figure should show the shortfall
    # rather than leave the reader to count bars. An empty slot is the honest
    # shape for "this should be here and is not": it occupies the space a
    # chromosome would, at the median size of the ones that were found, and is
    # drawn as outline only so it cannot be mistaken for sequence.
    expected = getattr(model, "expected_chromosomes", None)
    if expected:
        found = [q for q in lay.order if q.role == "chromosome"]
        missing = expected - len(found)
        if missing > 0:
            heights = sorted(lay.height[q.name] for q in found)
            typical = heights[len(heights) // 2] if heights else lay.max_bar_h * 0.5
            for k in range(missing):
                gx = lay.x[found[-1].name] + (k + 1) * (BAR_W + GAP) if found else MARGIN_L
                gtop = lay.baseline - typical
                add(
                    f'<path d="{_bar_path(gx, gtop, BAR_W, typical, BAR_W / 2, BAR_W / 2)}" '
                    f'fill="none" stroke="{PALETTE["muted"]}" stroke-width="2.5" '
                    f'stroke-dasharray="7 6"/>'
                )
                add(
                    f'<text x="{gx + BAR_W / 2:.1f}" y="{gtop + typical / 2:.1f}" '
                    f'font-size="{FS_PRIMARY + 8}" text-anchor="middle" '
                    f'fill="{PALETTE["muted"]}" font-weight="700">?</text>'
                )
                add(
                    f'<text x="{gx + BAR_W / 2:.1f}" y="{lay.baseline + 34:.1f}" '
                    f'font-size="{FS_ANNOT}" font-weight="700" text-anchor="middle" '
                    f'fill="{PALETTE["muted"]}">chr {len(found) + k + 1}</text>'
                )
                add(
                    f'<text x="{gx + BAR_W / 2:.1f}" y="{lay.baseline + 34 + FS_ANNOT + 6:.1f}" '
                    f'font-size="{FS_ANNOT}" text-anchor="middle" '
                    f'fill="{PALETTE["muted"]}">not found</text>'
                )

    # ---- unassigned panel ----
    if lay.panel:
        add(_unassigned_panel_svg(model, lay, interactive))

    # ---- key ----
    add(_key_svg(model, lay))

    # ---- legend ----
    add(legend_svg)
    add("</svg>")
    return "\n".join(P)


def _unassigned_panel_svg(model: Model, lay: Layout, interactive: bool) -> str:
    """
    Sequences that fit no chromosome, kept visibly separate rather than being
    forced into the karyotype or dropped from the figure.

    Drawn as upright bars like the chromosomes, but deliberately narrower and on
    their own side of a divider, so they read as the same kind of object without
    implying they belong to the karyotype. Labels only - no sentences.
    """
    items = model.unassigned()
    x, y = lay.panel_x, lay.header_h
    out = ['<g id="layer-unassigned" font-family="Helvetica, Arial, sans-serif">']

    col_w = BAR_W * 2.4
    bw = float(BAR_W)
    shown = items[: max(int(PANEL_W / col_w), 1)]
    for i, s_ in enumerate(shown):
        cx = x + i * col_w
        # As short as it can be and still read as a bar. The shared
        # MIN_DRAWN_PX floor is set for the GRAPH panel, where a short contig
        # has to be big enough to show which of its ends things join to; here
        # that question is already answered and the honest thing is to draw an
        # unplaced fragment as the small thing it is. The floor is 1.35 bar
        # widths, which leaves both ends visibly rounded - at exactly one width
        # a stadium is a circle, and a circle already means something else in
        # this figure.
        h = max(s_.length * lay.scale, BAR_W * 1.35)
        top = lay.baseline - h
        colour = model.segment_colours.get(s_.name, PALETTE["unassigned"])
        attrs = ""
        if interactive:
            attrs = (
                f' class="chrom" data-name="{esc(s_.name)}" data-role="unassigned"'
                f' data-length="{s_.length}"'
                f' data-depth="{s_.depth if s_.depth is not None else ""}"'
            )
        # A question mark, directly above the contig. These are contigs the graph
        # gives no way to place - most often nothing links to them at all - and
        # the honest label for that is a question, not a category. It sits on the
        # contig it refers to, so it needs no heading and no divider to explain
        # which column it belongs to.
        out.append(
            f'<text x="{cx + bw / 2:.1f}" y="{top - 14:.1f}" font-size="{FS_PRIMARY + 6}" '
            f'font-weight="700" text-anchor="middle" fill="{PALETTE["muted"]}">?</text>'
        )
        out.append(
            f'<path d="{_bar_path(cx, top, bw, h, bw / 2, bw / 2)}" fill="{colour}" '
            f'fill-opacity="0.95" stroke="none"{attrs}/>'
        )
        if h >= FS_PRIMARY + 4:
            out.append(
                f'<text x="{cx + bw / 2:.1f}" y="{top + h / 2 + FS_PRIMARY * 0.35:.1f}" '
                f'font-size="{FS_PRIMARY}" text-anchor="middle" fill="{_text_on(colour)}" '
                f'font-weight="700">{esc(_segment_number(s_.name))}</text>'
            )
        out.append(
            f'<text x="{cx + bw / 2:.1f}" y="{lay.baseline + 34:.1f}" '
            f'font-size="{FS_ANNOT}" font-weight="700" text-anchor="middle" '
            f'fill="{PALETTE["text"]}">unplaced</text>'
        )
        out.append(
            f'<text x="{cx + bw / 2:.1f}" y="{lay.baseline + 34 + FS_ANNOT + 6:.1f}" '
            f'font-size="{FS_ANNOT}" text-anchor="middle" '
            f'fill="{PALETTE["muted"]}">{figure_bp(s_.length)}</text>'
        )
    out.append("</g>")
    return "\n".join(out)


def _key_svg(model: Model, lay: Layout) -> str:
    """
    No key. Everything on this panel is either labelled where it sits or is a
    convention the report explains at length; a legend in the corner was one
    more thing to read before the picture made sense, which is the opposite of
    what a look-and-see figure is for.
    """
    return ""

def _nice_tick(max_len: int) -> int:
    raw = max_len / 8.0
    mag = 10 ** int(math.floor(math.log10(max(raw, 1))))
    for m in (1, 2, 5, 10):
        if raw <= m * mag:
            return int(m * mag)
    return int(10 * mag)


def _tick_label(v: int) -> str:
    if v == 0:
        return "0"
    if v >= 1e6:
        return f"{v / 1e6:g} Mb"
    if v >= 1e3:
        return f"{v / 1e3:g} kb"
    return str(v)


def _legend_svg(model: Model, lay: Layout) -> Tuple[str, float]:
    """
    v9: no floating key, no footnote block. Everything the reader needs is a
    label attached to the thing it describes, so this now draws nothing. The
    function survives because the layout asks it where the figure ends.
    """
    return "", lay.baseline + 40




# ==========================================================================
# interactive HTML
# ==========================================================================
HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<title>__TITLE__</title>
<style>
  :root { --fg:#1a1a1a; --muted:#6b6b6b; --line:#e2e2e2; --panel:#fafafa; }
  * { box-sizing:border-box; }
  body { margin:0; font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif;
         color:var(--fg); background:#fff; }
  header { padding:16px 22px; border-bottom:1px solid var(--line); }
  h1 { font-size:18px; margin:0 0 4px; }
  .sub { color:var(--muted); font-size:13px; }
  .wrap { display:flex; align-items:flex-start; gap:0; }
  .canvas { flex:1 1 auto; overflow:auto; padding:10px 0 40px 0; }
  aside { width:360px; flex:0 0 360px; border-left:1px solid var(--line); height:calc(100vh - 78px);
          overflow:auto; background:var(--panel); padding:16px 18px; }
  .controls { padding:10px 22px; border-bottom:1px solid var(--line); display:flex; gap:18px;
              flex-wrap:wrap; align-items:center; font-size:13px; }
  label.chk { display:inline-flex; gap:6px; align-items:center; cursor:pointer; user-select:none; }
  h2 { font-size:13px; text-transform:uppercase; letter-spacing:.05em; color:var(--muted);
       margin:20px 0 8px; }
  h2:first-child { margin-top:0; }
  table { width:100%; border-collapse:collapse; font-size:12.5px; }
  th,td { text-align:left; padding:5px 6px; border-bottom:1px solid var(--line); vertical-align:top; }
  th { color:var(--muted); font-weight:600; }
  .card { border:1px solid var(--line); border-radius:7px; padding:9px 11px; margin-bottom:8px;
          background:#fff; cursor:pointer; }
  .card:hover, .card.on { border-color:#888; }
  .card .t { font-weight:600; font-size:12.5px; display:flex; align-items:center; gap:7px; }
  .swatch { width:11px; height:11px; border-radius:2px; flex:0 0 auto; }
  .card .d { color:var(--muted); font-size:12px; margin-top:4px; }
  .why { color:var(--muted); font-size:11.5px; margin-top:4px; font-style:italic; }
  .pill { display:inline-block; font-size:11px; padding:1px 6px; border-radius:9px;
          background:#eee; color:#444; margin-left:4px; }
  #tip { position:fixed; pointer-events:none; background:#111; color:#fff; padding:6px 9px;
         border-radius:5px; font-size:12px; max-width:340px; opacity:0; transition:opacity .1s;
         z-index:20; }
  .dim { opacity:.12 !important; }
  details { margin:6px 0; } summary { cursor:pointer; font-size:12.5px; }
  .warn { background:#fff5e6; border:1px solid #f0c987; border-radius:6px; padding:9px 11px;
          font-size:12.5px; margin-bottom:8px; }
  svg { display:block; margin:0 auto; }
</style></head><body>
<header><h1>__TITLE__</h1><div class="sub">__SUMMARY__</div></header>
<div class="controls">
  <label class="chk"><input type="checkbox" id="c-tangles" checked> Graph features</label>
  <label class="chk"><input type="checkbox" id="c-coverage" checked> Coverage</label>
  <label class="chk"><input type="checkbox" id="c-annot" checked> Annotations</label>
  <label class="chk"><input type="checkbox" id="c-legend" checked> Legend</label>
  <span style="margin-left:auto;color:var(--muted);font-size:12px">
    zoom <input type="range" id="zoom" min="50" max="220" value="100" style="vertical-align:middle">
    <span id="zv">100%</span></span>
</div>
<div class="wrap">
  <div class="canvas"><div id="svgbox">__SVG__</div></div>
  <aside>__SIDE__</aside>
</div>
<div id="tip"></div>
<script>
const tip = document.getElementById('tip');
function showTip(e, html){ tip.innerHTML = html; tip.style.opacity = 1;
  const x = Math.min(e.clientX + 14, window.innerWidth - 360);
  tip.style.left = x + 'px'; tip.style.top = (e.clientY + 16) + 'px'; }
function hideTip(){ tip.style.opacity = 0; }

document.querySelectorAll('.tangle').forEach(el => {
  el.style.cursor = 'pointer';
  el.addEventListener('mousemove', e => showTip(e,
    '<b>' + el.dataset.type.replace(/_/g,' ') + '</b><br>' + el.dataset.desc));
  el.addEventListener('mouseleave', hideTip);
  el.addEventListener('click', () => select(el.dataset.id));
});
document.querySelectorAll('.chrom').forEach(el => {
  el.addEventListener('mousemove', e => showTip(e, '<b>' + el.dataset.name + '</b><br>' +
    el.dataset.role + ', ' + Number(el.dataset.length).toLocaleString() + ' bp' +
    (el.dataset.depth ? '<br>depth ' + el.dataset.depth + 'x' : '')));
  el.addEventListener('mouseleave', hideTip);
});
document.querySelectorAll('.annot').forEach(el => {
  el.addEventListener('mousemove', e => showTip(e, el.dataset.desc));
  el.addEventListener('mouseleave', hideTip);
});

let current = null;
function select(id){
  current = (current === id) ? null : id;
  document.querySelectorAll('.tangle').forEach(el => {
    el.classList.toggle('dim', current !== null && el.dataset.id !== current); });
  document.querySelectorAll('.card').forEach(c => {
    c.classList.toggle('on', c.dataset.id === current); });
  if (current){ const c = document.querySelector('.card[data-id="'+current+'"]');
    if (c) c.scrollIntoView({block:'nearest', behavior:'smooth'}); }
}
document.querySelectorAll('.card').forEach(c =>
  c.addEventListener('click', () => select(c.dataset.id)));

function toggle(id, sel){ document.getElementById(id).addEventListener('change', e => {
  document.querySelectorAll(sel).forEach(el => el.style.display = e.target.checked ? '' : 'none');
}); }
toggle('c-tangles', '#layer-tangles');
toggle('c-coverage', '#layer-coverage');
toggle('c-annot', '.annot');
toggle('c-legend', '#legend');

const svg = document.querySelector('#svgbox svg');
const baseW = svg ? svg.getAttribute('width') : 0;
document.getElementById('zoom').addEventListener('input', e => {
  const z = e.target.value; document.getElementById('zv').textContent = z + '%';
  if (svg){ svg.style.width = (baseW * z / 100) + 'px'; svg.style.height = 'auto'; }
});
</script></body></html>
"""


def render_html(model: Model) -> str:
    svg = render_svg(model, interactive=True)
    side: List[str] = []

    if model.warnings:
        side.append("<h2>Warnings</h2>")
        for w in model.warnings:
            side.append(f'<div class="warn">{esc(w)}</div>')

    side.append("<h2>Karyotype calls</h2><table>")
    side.append("<tr><th>Sequence</th><th>Length</th><th>Call</th><th>Confidence</th></tr>")
    for s in model.drawable() + model.unplaced()[:15]:
        side.append(
            f"<tr><td>{esc(s.display)}</td><td>{figure_bp(s.length)}</td>"
            f"<td>{esc(s.role)}</td><td>{esc(_confidence(s))}</td></tr>"
        )
    side.append("</table>")

    for s in model.drawable():
        if not s.evidence:
            continue
        side.append(
            f"<details><summary>Why {esc(s.display)} was called {esc(s.role)}</summary><ul>"
            + "".join(f"<li>{esc(e.as_text())}</li>" for e in s.evidence)
            + "</ul></details>"
        )

    side.append(f"<h2>Graph features ({len(model.tangles)})</h2>")
    if not model.tangles:
        side.append('<div class="sub">No tangles detected, or no assembly graph supplied.</div>')
    for t in model.tangles:
        colour = TANGLE_STYLE.get(t.type, ("#888", ""))[0]
        mult = (
            f'<span class="pill">~{t.multiplicity:g} copies</span>'
            if t.multiplicity
            else ""
        )
        side.append(
            f'<div class="card" data-id="{esc(t.id)}">'
            f'<div class="t"><span class="swatch" style="background:{colour}"></span>'
            f"{esc(TANGLE_LABEL.get(t.type, t.type))}{mult}</div>"
            f'<div class="d">{esc(t.description)}</div>'
            f'<div class="why">on {esc(", ".join(t.sequences)) or "unplaced"}'
            + (f" &middot; {esc('; '.join(t.evidence))}" if t.evidence else "")
            + "</div></div>"
        )

    if model.coverage_anomalies:
        side.append(f"<h2>Coverage outliers ({len(model.coverage_anomalies)})</h2><table>")
        side.append("<tr><th>Region</th><th>Type</th><th>vs median</th></tr>")
        for a in sorted(model.coverage_anomalies, key=lambda a: -abs(a["peak"] - 1))[:40]:
            side.append(
                f"<tr><td>{esc(a['seqname'])}:{a['start']:,}-{a['end']:,}</td>"
                f"<td>{esc(a['kind'])}</td><td>{a['peak']:.1f}x</td></tr>"
            )
        side.append("</table>")

    return (
        HTML_TEMPLATE.replace("__TITLE__", esc(model.title))
        .replace("__SUMMARY__", esc(model.summary_sentence()))
        .replace("__SVG__", svg)
        .replace("__SIDE__", "\n".join(side))
    )
