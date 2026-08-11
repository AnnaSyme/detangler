"""The paired figure: graph beside chromosomes."""
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
    Log,
    PALETTE,
    esc,
)
from .records import (
    GfaLink,
)
from .graph import (
    build_adjacency,
)
from .palette import (
    CLASS_LABEL,
)
from .calls import (
    SegmentCall,
)
from .model import (
    Model,
)
from .render_common import (
    BAR_W,
    GAP,
    FS_HEADING,
    _place_svg,
    _svg_height,
    _svg_width,
    embed_image,
    image_size,
)
from .render_ideogram import (
    ideogram_block_anchors,
    ideogram_geometry,
    render_svg,
)
from .render_graph import (
    segment_thickness,
    _graph_layout,
    graph_svg_for_style,
    render_graph_figure,
)



# ==========================================================================
# pipeline
# ==========================================================================
PAIR_GUTTER = 96.0


_NUM_RE = re.compile(r"-?\d+(?:\.\d+)?")
_STROKE_MARGIN = 40.0


def _ink_right(svg: str) -> Optional[float]:
    """
    Rightmost x at which the graph panel has ink. The panel's declared width
    includes padding sized for the worst case - a ring's radius, a label's
    reach - so butting the chromosome column against that width leaves a band
    of white that belongs to neither panel.
    """
    best: Optional[float] = None
    for m in re.finditer(r'\sd="([^"]+)"', svg):
        nums = [float(v) for v in _NUM_RE.findall(m.group(1))]
        for i in range(0, len(nums) - 1, 2):
            if best is None or nums[i] > best:
                best = nums[i]
    for m in re.finditer(
        r'<circle[^>]*\scx="([-\d.]+)"[^>]*\sr="([-\d.]+)"', svg
    ):
        v = float(m.group(1)) + float(m.group(2))
        if best is None or v > best:
            best = v
    return None if best is None else best + _STROKE_MARGIN


def _ink_bottom_in_column(svg: str, x0: float, x1: float) -> Optional[float]:
    """
    Lowest y at which the graph panel has ink, counting only points inside the
    x range [x0, x1]. Used to decide how far the chromosome panel can rise into
    an empty corner. Coarse on purpose - it reads path coordinates rather than
    rendering, so a curve that bulges between two control points is under-
    measured by at most the bulge, which the gap covers.
    """
    best: Optional[float] = None
    for m in re.finditer(r'\sd="([^"]+)"', svg):
        nums = [float(v) for v in _NUM_RE.findall(m.group(1))]
        for i in range(0, len(nums) - 1, 2):
            x, y = nums[i], nums[i + 1]
            if x0 <= x <= x1 and (best is None or y > best):
                best = y
    if best is not None:
        best += _STROKE_MARGIN   # ribbons are stroked, so ink sits below the spine
    for m in re.finditer(
        r'<circle[^>]*\scx="([-\d.]+)"[^>]*\scy="([-\d.]+)"[^>]*\sr="([-\d.]+)"', svg
    ):
        x, y, rad = float(m.group(1)), float(m.group(2)), float(m.group(3))
        # a ring is drawn as a STROKED circle, so its ink reaches r plus half the
        # stroke below the centre - taking cy alone put the bottom of every
        # circular contig a full radius higher than it really is
        sw = 0.0
        tail = svg[m.end(): m.end() + 200]
        sm = re.search(r'stroke-width="([-\d.]+)"', tail)
        if sm:
            sw = float(sm.group(1)) / 2.0
        low = y + rad + sw
        if x0 <= x <= x1 and (best is None or low > best):
            best = low
    return best


def _content_span(lay) -> Tuple[float, float]:
    """
    Left and right edge of the chromosome panel's actual INK, in the panel's own
    coordinates. `lay.width` is the panel's canvas, which carries trailing room
    for a key and a right margin that this figure does not draw - measuring
    against it left a band of white to the right of the tallest chromosome and
    scaled every bar down to pay for it.
    """
    xs = [lay.x[q.name] for q in lay.order]
    if not xs:
        return 0.0, float(lay.width)
    lo = min(xs)
    if getattr(lay, "panel", False):
        lo = min(lo, float(lay.panel_x))   # the unplaced column sits left of the bars
    hi = max(xs) + BAR_W
    # Slots for chromosomes the user says should exist but the graph did not
    # produce are drawn to the RIGHT of the largest one, and they are ink too.
    # Leaving them out of the span pushed them off the edge of the figure as
    # soon as the panel was right-justified.
    ghosts = int(getattr(lay, "ghost_slots", 0) or 0)
    if ghosts:
        hi += ghosts * (BAR_W + GAP)
    return lo, hi


def _tallest_bar_px(lay) -> float:
    """Height of the tallest drawn molecule, in the panel's own coordinates."""
    return max([0.0] + [lay.height.get(q.name, 0.0) for q in lay.order])


def render_paired_svg(
    model: Model,
    calls: List[SegmentCall],
    links: List[GfaLink],
    colours: Dict[str, str],
    args,
    log: Log,
) -> str:
    """
    One figure: the assembly graph on the left, the chromosomes it resolves into
    on the right, and faint leader lines joining a graph node to the block it
    became. Both panels are nested <svg> elements, which keeps each renderer
    independent and avoids the two fighting over coordinates.
    """
    adj = build_adjacency(links)
    pos, gw, gh, _ = _graph_layout(calls, adj)
    # the combined figure carries the title, so the panels must not repeat it
    real_title = model.title
    # Both panel headings are drawn by the COMBINED figure, not one by each
    # panel. Drawn inside the panel the right-hand heading inherits that panel's
    # offset and sits on a different baseline from the left one, which reads as
    # a size difference even though both are the same size.
    right_label = "Possible chromosomes"
    # One heading across the top, not one per panel. Read left to right it
    # states what the figure is FOR - this graph, resolved into those molecules -
    # which two separate captions never quite said.
    pair_label = ""
    model.title = ""
    try:
        ideo_svg = render_svg(model)
        lay, _, _, _ = ideogram_geometry(model)
        anchors = ideogram_block_anchors(model)
    finally:
        model.title = real_title

    # Build the chromosome panel FIRST so the graph panel can borrow its scale:
    # a contig should be the same size in both halves of the figure.
    args.graph_px_per_bp = lay.scale
    graph_svg = graph_svg_for_style(calls, links, "", colours, args, log)
    if args.graph_style == "bandage":
        gw, gh = _svg_width(graph_svg), _svg_height(graph_svg)

    iw, ih = lay.width, _svg_height(ideo_svg)

    # ---- left panel: a real Bandage export if given, otherwise our redraw ----
    external = args.bandage_image and os.path.exists(args.bandage_image)
    if args.bandage_image and not external:
        log.warn(f"--bandage-image {args.bandage_image} not found; using our own redraw instead")
    if external:
        src = image_size(args.bandage_image)
        if src is None:
            log.warn(
                f"could not read the dimensions of {args.bandage_image}; expected PNG, JPEG or "
                f"SVG. Using our own redraw instead."
            )
            external = False
    if external:
        # scale the Bandage export to the height of the chromosome panel
        panel_h = ih - 70
        panel_w = panel_h * (src[0] / src[1])
        if panel_w > args.bandage_max_width:
            panel_w = args.bandage_max_width
            panel_h = panel_w * (src[1] / src[0])
        gw_eff, gh_eff = panel_w, panel_h
        left_label = "Assembly graph, as drawn by Bandage"
    elif args.rotate_graph:
        # rotating our redraw a quarter turn trades a very wide figure for a
        # taller, narrower one that fits a page or a slide
        gw_eff, gh_eff = gh, gw
        left_label = "Assembly graph"
    else:
        gw_eff, gh_eff = gw, gh
        left_label = "Assembly graph"

    # STACKED, not side by side. The graph panel comes out wide and shallow and
    # the chromosome panel narrow and deep, so putting them side by side left
    # most of the canvas empty and forced everything to be shrunk to fit. One
    # above the other fills the space, which means both can be drawn larger.
    top = FS_HEADING * 1.5 + 34.0
    gap = PAIR_GUTTER * 0.35

    # SIDE BY SIDE: graph on the left, chromosomes on the right.
    #
    # Earlier versions stacked the panels, then overlapped them on a diagonal to
    # save the white space each one's empty corner cost. Both worked, and both
    # needed the two panels to know about each other's shape - a per-bar
    # clearance test, a keep-out in the graph's solver, a rotation that had to
    # aim somewhere. Two panels beside each other need none of that: each is a
    # rectangle, they cannot collide, and the reader gets two pictures rather
    # than one picture with two halves. The complexity was buying compactness,
    # which was never the thing the figure is for.
    #
    # The chromosome panel is scaled so its tallest bar matches the graph
    # panel's height, which is what keeps the two reading as equals.
    canvas_h = gh_eff
    if iw > 1.0:
        bar_target = canvas_h * 0.66
        model.max_bar_h = bar_target
        model.title = ""
        try:
            ideo_svg = render_svg(model)
            lay, _, _, _ = ideogram_geometry(model)
            anchors = ideogram_block_anchors(model)
        finally:
            model.title = real_title
        iw, ih = lay.width, _svg_height(ideo_svg)
    c_lo, c_hi = _content_span(lay)
    ideo_scale = 1.0
    iw_eff, ih_eff = (c_hi - c_lo) * ideo_scale, ih * ideo_scale

    HEAD_BAND = 24.0
    BORDER_PAD = 16.0
    pair_label = f"{left_label}  \u2192  {right_label}"
    # Butt the chromosome column against the graph's INK, not against its
    # declared width, which carries worst-case padding for rings and labels.
    right_edge = gw_eff
    if not external:
        ink_r = _ink_right(graph_svg)
        if ink_r is not None:
            right_edge = min(gw_eff, ink_r)
    ideo_x = right_edge + gap - c_lo * ideo_scale
    width = right_edge + gap + iw_eff + 60 + BORDER_PAD
    height = max(top + gh_eff, top + ih_eff) + 20 + HEAD_BAND + BORDER_PAD
    # Both panels stand on the same bottom line, so the chromosome baseline and
    # the foot of the graph agree.
    ideo_y = height - HEAD_BAND - 20 - ih_eff
    P = [
        f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
        f'width="{width:.0f}" height="{height:.0f}" '
        f'viewBox="0 0 {width:.0f} {height:.0f}" font-family="Helvetica, Arial, sans-serif">',
        f'<rect width="100%" height="100%" fill="{PALETTE["bg"]}"/>',
        # A dotted rule around everything, heading included. The two panels are
        # one figure making one claim; without a border they read as two images
        # that happen to have been saved together, and the heading looks like a
        # caption for whichever one it sits nearest.
        f'<rect x="{BORDER_PAD:.0f}" y="{BORDER_PAD:.0f}" '
        f'width="{width - 2 * BORDER_PAD:.0f}" height="{height - 2 * BORDER_PAD:.0f}" '
        f'rx="10" fill="none" stroke="{PALETTE["muted"]}" stroke-width="2" '
        f'stroke-dasharray="6 7" stroke-opacity="0.75"/>',
        f'<text x="{width / 2:.0f}" y="{top + 6:.0f}" text-anchor="middle" '
        f'font-size="{FS_HEADING * 1.5:.0f}" font-weight="700" '
        f'fill="{PALETTE["text"]}">{esc(pair_label)}</text>',
    ]

    if external:
        P.append(embed_image(args.bandage_image, 0, top, gw_eff, gh_eff))
    elif args.rotate_graph:
        P.append(_place_svg(graph_svg, 0, top, rotate=-90))
    else:
        P.append(_place_svg(graph_svg, 0, top))
    P.append(_place_svg(ideo_svg, ideo_x, ideo_y, scale=ideo_scale))
    P.append("</svg>")
    return "\n".join(P)


def write_bandage_colour_csv(
    calls: List[SegmentCall], colours: Dict[str, str], path: str, log: Log
) -> str:
    """
    A CSV Bandage can load to colour the graph the way we do, so a real Bandage
    render and our chromosome figure agree.

    Bandage keys rows by segment name in the first column and recognises a
    colour column; the remaining columns show up as node labels. If the colours
    do not take, check the column name against Bandage's own colour-schemes
    documentation for your version rather than assuming this header is right.
    """
    import csv

    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Name", "Colour", "Color", "Class", "CopyNumber", "Depth", "Length"])
        for c in sorted(calls, key=lambda c: -c.length):
            colour = colours.get(c.name, "#cfcfcf")
            w.writerow([
                c.name, colour, colour, CLASS_LABEL.get(c.cls, c.cls),
                f"{c.copy_number:.2f}" if c.copy_number is not None else "",
                f"{c.depth:.1f}" if c.depth is not None else "",
                c.length,
            ])
    log.info(f"wrote Bandage colour CSV for {len(calls)} segments to {path}")
    return path


def write_figures(
    model: Model,
    calls: List[SegmentCall],
    links: List[GfaLink],
    colours: Dict[str, str],
    base: str,
    args,
    log: Log,
) -> Dict[str, str]:
    """The graph panel, and the paired graph-plus-chromosomes figure."""
    out: Dict[str, str] = {}
    graph_path = render_graph_figure(
        calls, links, colours, model.title, base + "_graph.svg", args, log
    )
    if not graph_path:
        return out
    out["assembly graph figure"] = graph_path
    out["Bandage colour CSV (load this in Bandage)"] = write_bandage_colour_csv(
        calls, colours, base + "_bandage_colours.csv", log
    )
    paired = base + "_paired.svg"
    with open(paired, "w") as fh:
        fh.write(render_paired_svg(model, calls, links, colours, args, log))
    out["PAIRED FIGURE: graph and chromosomes"] = paired
    return out
