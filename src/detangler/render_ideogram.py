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
    PALETTE,
    TANGLE_LABEL,
    TANGLE_STYLE,
    esc,
    human_bp,
    median,
    wrap_text,
)
from .palette import (
    _segment_number,
    _text_on,
)
from .model import (
    Model,
    _confidence,
)
from .render_common import (
    BAR_W,
    COV_W,
    FS_HEADING,
    FS_LABEL,
    FS_SUB,
    Layout,
    MARGIN_L,
    MAX_BAR_H,
    _annotation_colour,
    _bar_path,
    representative_anchors,
)



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


def _ideogram_frame(model: Model) -> Dict[str, object]:
    """
    Geometry for the chromosome figure, computed once and shared by the renderer
    and by anything that needs to point at a block from outside - the paired
    figure draws leader lines to these exact coordinates.
    """
    show_cov = bool(model.coverage) and model.settings.get("coverage", True)
    probe = Layout(model, show_cov)
    drawn = {s.name for s in probe.order}

    # v9: no summary paragraph under the panel title. The reasoning belongs in the
    # report; the figure carries only what points at something it draws.
    head_lines: List[str] = []
    header_h = 54 + 26 + 15 * (probe.max_label_lines - 1)

    lay = Layout(model, show_cov, header_h)
    legend_svg, legend_bottom = _legend_svg(model, lay)
    total_h = max(legend_bottom + 26, header_h + MAX_BAR_H + 90)
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

    # ---- scale ruler ----
    add(f'<g id="ruler" stroke="{PALETTE["grid"]}">')
    tick = _nice_tick(lay.max_len)
    v = 0
    while v <= lay.max_len:
        y = header_h + v * lay.scale
        add(f'<line x1="{MARGIN_L - 46}" y1="{y:.1f}" x2="{MARGIN_L - 12}" y2="{y:.1f}"/>')
        add(
            f'<text x="{MARGIN_L - 50}" y="{y + 4:.1f}" font-size="10.5" text-anchor="end" '
            f'fill="{PALETTE["muted"]}" stroke="none">{_tick_label(v)}</text>'
        )
        v += tick
    add("</g>")

    # ---- tangle arcs (behind the bars) ----
    add('<g id="layer-tangles">')
    # v9 drops the arcs joining one chromosome to another. They were the biggest
    # source of visual noise and duplicated what the shared segment colours
    # already say: a repeat linking two chains is drawn as a copy on each. Only
    # single-point features survive, as a marker beside the bar.
    for t in model.tangles:
        anchors = representative_anchors(t, drawn)
        if len(anchors) >= 2:
            continue
        colour, _dash = TANGLE_STYLE.get(t.type, ("#888888", ""))
        attrs = ""
        if interactive:
            attrs = (
                f' class="tangle" data-id="{esc(t.id)}" data-type="{esc(t.type)}"'
                f' data-desc="{esc(t.description)}"'
            )
        for a in anchors:
            y = lay.y(a.seqname, (a.start + a.end) / 2.0)
            x = lay.x[a.seqname] + BAR_W + 5
            add(
                f'<path d="M {x:.1f} {y:.1f} l 9 -5 l 0 10 z" fill="{colour}" '
                f'fill-opacity="0.9"{attrs}/>'
            )
    add("</g>")

    # ---- chromosome bars ----
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
            ccx, ccy = x + BAR_W / 2.0, top + r + 6
            add(
                f'<circle cx="{ccx:.1f}" cy="{ccy:.1f}" r="{r:.1f}" fill="none" '
                f'stroke="{PALETTE["bar_edge"]}" stroke-width="{ring_w + 2.0:.1f}"{battrs}/>'
            )
            add(
                f'<circle cx="{ccx:.1f}" cy="{ccy:.1f}" r="{r:.1f}" fill="none" '
                f'stroke="{seg_colour}" stroke-width="{ring_w:.1f}"/>'
            )
            add(
                f'<text x="{ccx:.1f}" y="{top - 10:.1f}" font-size="{FS_SUB + 2}" '
                f'text-anchor="middle" fill="{PALETTE["text"]}" font-weight="600">'
                f'{esc(s.role)}</text>'
            )
            add(
                f'<text x="{ccx:.1f}" y="{ccy + r * 2 + FS_SUB + 6:.1f}" '
                f'font-size="{FS_SUB}" text-anchor="middle" fill="{PALETTE["muted"]}">'
                f'{human_bp(s.length)}</text>'
            )
            continue

        add(
            f'<path d="{_bar_path(x, top, BAR_W, h, 0.0 if s.caps.get("top") else rx, 0.0 if s.caps.get("bottom") else rx)}" '
            f'fill="{fill}" fill-opacity="0.82" stroke="{PALETTE["bar_edge"]}" '
            f'stroke-width="0.9"{battrs}/>'
        )

        # Segment blocks: the same colour the segment has in the graph figure, so
        # a node over there can be found on a chromosome over here.
        if s.blocks:
            add('<g class="blocks">')
            last = len(s.blocks) - 1
            inset = 0.0 if s.blocks_tile else 3.5
            bw = BAR_W - 2 * inset
            for bi, (b_start, b_end, seg, colour) in enumerate(s.blocks):
                y1 = lay.y(s.name, b_start)
                y2 = max(lay.y(s.name, b_end), y1 + 1.6)
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
                rt = rx if (s.blocks_tile and bi == 0 and not s.caps.get("top")) else 0
                rb = rx if (s.blocks_tile and bi == last and not s.caps.get("bottom")) else 0
                add(
                    f'<path d="{_bar_path(x + inset, y1, bw, y2 - y1, rt, rb)}" '
                    f'fill="{colour}" fill-opacity="{0.95 if s.blocks_tile else 0.92}" '
                    f'stroke="#ffffff" stroke-width="0.6"{bl}/>'
                )
                # v9: the segment NUMBER, set inside the block, inked white or
                # dark for contrast against that block's own colour
                if y2 - y1 >= FS_LABEL + 4:
                    add(
                        f'<text x="{x + BAR_W / 2:.1f}" y="{(y1 + y2) / 2 + FS_LABEL * 0.35:.1f}" '
                        f'font-size="{FS_LABEL}" text-anchor="middle" fill="{_text_on(colour)}" '
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

        # re-stroke the outline so bands do not spill past the rounded ends
        add(
            f'<path d="{_bar_path(x, top, BAR_W, h, 0.0 if s.caps.get("top") else rx, 0.0 if s.caps.get("bottom") else rx)}" '
            f'fill="none" stroke="{PALETTE["bar_edge"]}" stroke-width="1.1"/>'
        )

        # Repeats attached to a free end, drawn hanging OFF the bar rather than
        # inside it: they are not part of the molecule and not on the Mb scale,
        # but they are what tells you this end is a telomere or an rDNA block.
        # Flush against the bar, not floating beside it, so a molecule reads as
        # one object: cap, backbone, cap. Only the OUTER corner of the outermost
        # cap is rounded; every join between blocks is square so they abut.
        cap_h = BAR_W * 0.62
        for side, entries in sorted(s.caps.items()):
            n_side = len(entries)
            for ci, (seg, colour) in enumerate(entries):
                if side == "top":
                    cy = top - cap_h * (ci + 1)
                    rt = rx if ci == n_side - 1 else 0.0
                    rb = 0.0
                else:
                    cy = top + h + cap_h * ci
                    rt = 0.0
                    rb = rx if ci == n_side - 1 else 0.0
                add(
                    f'<path d="{_bar_path(x, cy, BAR_W, cap_h, rt, rb)}" fill="{colour}" '
                    f'fill-opacity="0.95" stroke="{PALETTE["bar_edge"]}" stroke-width="1.1"/>'
                )
                add(
                    f'<text x="{x + BAR_W / 2:.1f}" y="{cy + cap_h / 2 + FS_SUB * 0.36:.1f}" '
                    f'font-size="{FS_SUB + 1}" text-anchor="middle" fill="{_text_on(colour)}" '
                    f'font-weight="700">{esc(_segment_number(seg))}</text>'
                )

        # size only. Chain headings are gone: which contigs belong together is
        # shown by the numbered blocks in the bar, not by a caption above it.
        n_top = len(s.caps.get("top", []))
        add(
            f'<text x="{x + BAR_W / 2:.1f}" y="{top + h + 20 + cap_h * len(s.caps.get("bottom", [])):.1f}" '
            f'font-size="{FS_SUB}" text-anchor="middle" fill="{PALETTE["muted"]}">'
            f'{human_bp(s.length)}</text>'
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
            f'<text x="{MARGIN_L}" y="{header_h + MAX_BAR_H + 52:.1f}" font-size="10" '
            f'fill="{PALETTE["muted"]}">Coverage track (right of each bar): 0 to 2x the genome '
            f'median ({gm:.0f}x); dashed line = median. Red/blue ticks left of a bar mark '
            f'depth outliers.</text>'
        )
        add("</g>")

    # ---- unassigned panel ----
    if lay.panel:
        add(_unassigned_panel_svg(model, lay, interactive))

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
    out.append(
        f'<line x1="{x - 14:.1f}" y1="{y - 26:.1f}" x2="{x - 14:.1f}" '
        f'y2="{y + MAX_BAR_H:.1f}" stroke="{PALETTE["grid"]}" stroke-dasharray="3 4"/>'
    )
    out.append(
        f'<text x="{x:.1f}" y="{y - 24:.1f}" font-size="{FS_SUB + 2}" font-weight="600" '
        f'fill="{PALETTE["text"]}">Not assigned</text>'
    )

    col_w = BAR_W * 2.4
    bw = BAR_W * 0.62
    max_len = max((s.length for s in items), default=1)
    top = y + 34.0
    max_h = 210.0
    shown = items[: max(int((lay.width - x) / col_w), 1)]
    for i, s in enumerate(shown):
        cx = x + i * col_w
        # log height: these span kb to tens of kb and would otherwise vanish
        h = max(
            18.0,
            max_h * math.log10(max(s.length, 10)) / math.log10(max(max_len, 100)),
        )
        colour = model.segment_colours.get(s.name, PALETTE["unassigned"])
        attrs = ""
        if interactive:
            attrs = (
                f' class="chrom" data-name="{esc(s.name)}" data-role="unassigned"'
                f' data-length="{s.length}"'
                f' data-depth="{s.depth if s.depth is not None else ""}"'
            )
        out.append(
            f'<rect x="{cx:.1f}" y="{top:.1f}" width="{bw:.1f}" height="{h:.1f}" '
            f'rx="{bw / 2:.1f}" ry="{bw / 2:.1f}" fill="{colour}" fill-opacity="0.9" '
            f'stroke="{PALETTE["bar_edge"]}" stroke-width="1.1"{attrs}/>'
        )
        if h >= FS_LABEL + 4:
            out.append(
                f'<text x="{cx + bw / 2:.1f}" y="{top + h / 2 + FS_LABEL * 0.35:.1f}" '
                f'font-size="{FS_LABEL}" text-anchor="middle" fill="{_text_on(colour)}" '
                f'font-weight="700">{esc(_segment_number(s.name))}</text>'
            )
        out.append(
            f'<text x="{cx + bw / 2:.1f}" y="{top + h + FS_SUB + 4:.1f}" '
            f'font-size="{FS_SUB}" text-anchor="middle" fill="{PALETTE["muted"]}">'
            f'{human_bp(s.length)}</text>'
        )
    out.append("</g>")
    return "\n".join(out)


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
    return "", lay.header_h + MAX_BAR_H + 40


def _legend_svg_unused(model: Model, lay: Layout) -> Tuple[str, float]:
    y0 = lay.header_h + MAX_BAR_H + 76
    out = [f'<g id="legend" font-size="11" fill="{PALETTE["text"]}">']
    x = MARGIN_L
    if any(s.blocks for s in model.sequences):
        out.append(
            f'<text x="{x}" y="{y0 + 1}" fill="{PALETTE["muted"]}">'
            f'Blocks within each bar are graph segments, coloured as in the assembly graph '
            f'figure so they can be traced between the two.</text>'
        )
        y0 += 20
    roles = [r for r in ("chromosome", "mitochondrion", "plastid") if any(s.role == r for s in model.sequences)]
    for r in roles:
        out.append(
            f'<rect x="{x}" y="{y0 - 9}" width="13" height="13" rx="6" fill="{PALETTE[r]}" '
            f'fill-opacity="0.82" stroke="{PALETTE["bar_edge"]}" stroke-width="0.8"/>'
        )
        out.append(f'<text x="{x + 19}" y="{y0 + 1}">{r}</text>')
        x += 24 + 7.2 * len(r)

    types = []
    for t in model.tangles:
        if t.type not in types:
            types.append(t.type)
    y = y0 + 24
    x = MARGIN_L
    for tt in types:
        colour, dash = TANGLE_STYLE.get(tt, ("#888888", ""))
        label = TANGLE_LABEL.get(tt, tt)
        dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
        out.append(
            f'<line x1="{x}" y1="{y - 4}" x2="{x + 26}" y2="{y - 4}" stroke="{colour}" '
            f'stroke-width="2.6"{dash_attr}/>'
        )
        out.append(f'<text x="{x + 32}" y="{y}">{esc(label)}</text>')
        x += 46 + 6.6 * len(label)
        if x > lay.width - 260:
            x = MARGIN_L
            y += 20

    notes: List[str] = []
    if lay.not_to_scale:
        notes.append(
            "* drawn at a minimum height so it stays visible: hatched bars are not to scale."
        )
    up = model.unplaced()
    if up:
        total = max(sum(s.length for s in model.sequences), 1)
        notes.append(
            f"{len(up)} unplaced sequence(s) are not drawn, totalling "
            f"{human_bp(sum(s.length for s in up))} "
            f"({100.0 * sum(s.length for s in up) / total:.1f}% of the assembly); "
            f"longest {human_bp(max(s.length for s in up))}."
        )
    for note in notes:
        for line in wrap_text(note, lay.text_cols):
            y += 18
            out.append(f'<text x="{MARGIN_L}" y="{y}" fill="{PALETTE["muted"]}">{esc(line)}</text>')
    out.append("</g>")
    return "\n".join(out), y


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
            f"<tr><td>{esc(s.display)}</td><td>{human_bp(s.length)}</td>"
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
