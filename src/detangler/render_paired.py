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
    FS_HEADING,
    FS_TITLE,
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
    _graph_layout,
    graph_svg_for_style,
    render_graph_figure,
)



# ==========================================================================
# pipeline
# ==========================================================================
PAIR_GUTTER = 96.0


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
    right_label = "Hypothesis of chromosome structure"
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
    top = 44.0
    # The graph panel carries its own bottom padding, so a symmetric gap put the
    # lower heading nearer the panel ABOVE it than the one it labels. The gap is
    # tightened and the heading sits just above its own panel.
    gap = PAIR_GUTTER * 0.55
    width = max(gw_eff, iw) + 60
    ideo_y = top + gh_eff + gap
    height = ideo_y + ih + 20

    P = [
        f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
        f'width="{width:.0f}" height="{height:.0f}" '
        f'viewBox="0 0 {width:.0f} {height:.0f}" font-family="Helvetica, Arial, sans-serif">',
        f'<rect width="100%" height="100%" fill="{PALETTE["bg"]}"/>',
        f'<text x="60" y="{top - 10:.0f}" font-size="{FS_HEADING}" font-weight="600" '
        f'fill="{PALETTE["text"]}">{esc(left_label)}</text>',
    ]

    if external:
        P.append(embed_image(args.bandage_image, 0, top, gw_eff, gh_eff))
    elif args.rotate_graph:
        P.append(_place_svg(graph_svg, 0, top, rotate=-90))
    else:
        P.append(_place_svg(graph_svg, 0, top))
    P.append(_place_svg(ideo_svg, 0, ideo_y))
    P.append(
        f'<text x="60" y="{ideo_y - 10:.0f}" font-size="{FS_HEADING}" '
        f'font-weight="600" fill="{PALETTE["text"]}">{esc(right_label)}</text>'
    )
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
