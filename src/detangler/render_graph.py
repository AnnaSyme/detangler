"""The assembly graph panel."""
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
    _Rand,
    esc,
    human_bp,
)
from .records import (
    GfaLink,
)
from .graph import (
    build_adjacency,
)
from .palette import (
    CLASS_COLOUR,
    CLASS_LABEL,
    _segment_number,
    _text_on,
    assign_segment_colours,
)
from .calls import (
    SegmentCall,
)
from .render_common import (
    MIN_DRAWN_PX,
    drawn_length_px,
    BAR_W,
    FS_HEADING,
    FS_LABEL,
    FS_SUB,
    _arc_path,
)



# --------------------------------------------------------------------------
# assembly graph figure (deterministic layout, class colours)
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# Bandage-style graph rendering
#
# Bandage draws each segment as a thick tapered path whose drawn length tracks
# its sequence length, arranged by a force-directed layout, so a repeat that
# joins several contigs shows up as a knot. This reproduces that look rather
# than that layout: it is our own spring model, seeded deterministically, so two
# runs give the same picture. It will not be pixel-identical to Bandage.
# --------------------------------------------------------------------------
def segment_draw_length(length: int, args) -> float:
    """Drawn length of a segment. Square-root scaled, as a compromise between
    Bandage's proportional default (which makes a 2 kb repeat invisible next to
    a 9 Mb contig) and a log scale (which makes them nearly equal)."""
    px_per_bp = getattr(args, "graph_px_per_bp", None)
    if px_per_bp:
        # Same bases-per-pixel AND the same floor as the chromosome panel, so a
        # contig is the same size in both halves of the figure. Anything shorter
        # than the floor is drawn LARGER than true - unavoidable if a 2 kb
        # segment is to be visible beside a 9 Mb one - but at least it is
        # overstated by the same amount in both places.
        return drawn_length_px(length, px_per_bp)
    return min(max(10.0 + args.graph_length_scale * math.sqrt(max(length, 1)),
                   MIN_DRAWN_PX),
               args.graph_max_segment_px)


def segment_thickness(depth: Optional[float] = None) -> float:
    """
    Uniform, and equal to the chromosome bar width (v9 design). Thickness used to
    track read depth, but that put a second variable into the ribbon width and
    made the two panels hard to match up; depth is now carried by the label only.
    """
    return float(BAR_W)


def bandage_layout(
    calls: List[SegmentCall], links: List[GfaLink], args, log: Log
) -> Tuple[Dict[str, List[Tuple[float, float]]], float, float]:
    """
    Lay the graph out as flexible polylines, one per contig.

    The approach is Bandage's, reimplemented: a contig is not a rigid stick
    between two endpoints but a CHAIN of beads spaced a fixed distance apart,
    joined by stiff springs. The force model then runs over every bead, so a
    9 Mb contig can bend around its neighbours instead of ploughing through
    them. Modelling each contig as a single spring - which is what this did
    before - puts a hard ceiling on how readable a busy graph can ever be, and
    no amount of repulsion tuning lifts it.

    Returns {segment: [(x, y), ...]} plus the canvas size. The polyline's total
    length is the contig's drawn length, so it still matches the chromosome
    panel.
    """
    by_name = {c.name: c for c in calls}
    names = sorted(by_name)
    if not names:
        return {}, 100.0, 100.0

    spacing = max(segment_thickness() * 0.9, 8.0)
    chain: Dict[str, List[int]] = {}
    pts: List[List[float]] = []
    springs: List[Tuple[int, int, float, float]] = []
    rnd = _Rand(20260810)

    radius = 60.0 + 14.0 * math.sqrt(len(names))
    for si, n in enumerate(names):
        drawn = segment_draw_length(by_name[n].length, args)
        beads = max(int(math.ceil(drawn / spacing)), 1) + 1
        rest = drawn / (beads - 1) if beads > 1 else drawn
        a = 2 * math.pi * si / len(names)
        ox = radius * math.cos(a) + rnd.uniform(-10, 10)
        oy = radius * math.sin(a) + rnd.uniform(-10, 10)
        # lay the chain out straight, pointing away from the centre, so it
        # starts untangled rather than folded on itself
        idxs = []
        for b in range(beads):
            pts.append([ox + rest * b * math.cos(a), oy + rest * b * math.sin(a)])
            idxs.append(len(pts) - 1)
        chain[n] = idxs
        for b in range(beads - 1):
            springs.append((idxs[b], idxs[b + 1], rest, 6.0))
        # a weak brace across every second bead keeps a contig from crumpling
        # into a ball while still letting it curve
        for b in range(beads - 2):
            springs.append((idxs[b], idxs[b + 2], rest * 1.94, 0.55))

    def terminal(seg: str, end: str) -> int:
        return chain[seg][0] if end == "s" else chain[seg][-1]

    link_pairs: Set[Tuple[int, int]] = set()
    for l in links:
        if l.a == l.b or l.a not in chain or l.b not in chain:
            continue
        a_i = terminal(l.a, "e" if l.a_orient == "+" else "s")
        b_i = terminal(l.b, "s" if l.b_orient == "+" else "e")
        if a_i == b_i:
            continue
        springs.append((a_i, b_i, spacing * 0.9, 5.0))
        link_pairs.add((min(a_i, b_i), max(a_i, b_i)))

    n_pts = len(pts)
    k = spacing * 1.6
    iters = int(min(600, max(120, 26000 / max(n_pts, 1))))
    log.info(f"graph layout: {len(names)} contigs as {n_pts} points, {iters} iterations")
    temp = spacing * 2.0

    for _step in range(iters):
        disp = [[0.0, 0.0] for _ in range(n_pts)]
        for i in range(n_pts):
            xi, yi = pts[i]
            for j in range(i + 1, n_pts):
                if (i, j) in link_pairs:
                    continue
                dx = xi - pts[j][0]
                dy = yi - pts[j][1]
                d2 = dx * dx + dy * dy
                if d2 < 1e-6:
                    dx, dy, d2 = rnd.uniform(-1, 1), rnd.uniform(-1, 1), 1.0
                d = math.sqrt(d2)
                # repulsion is capped in range: beyond a few bead-widths two
                # contigs do not need to push each other around, and letting
                # them do so inflates the canvas
                if d > k * 6.0:
                    continue
                f = (k * k) / d
                ux, uy = dx / d, dy / d
                disp[i][0] += ux * f
                disp[i][1] += uy * f
                disp[j][0] -= ux * f
                disp[j][1] -= uy * f
        for a, b, rest, strength in springs:
            dx = pts[a][0] - pts[b][0]
            dy = pts[a][1] - pts[b][1]
            d = math.hypot(dx, dy) or 1e-6
            f = strength * (d - rest) * 0.28
            ux, uy = dx / d, dy / d
            disp[a][0] -= ux * f
            disp[a][1] -= uy * f
            disp[b][0] += ux * f
            disp[b][1] += uy * f
        for i in range(n_pts):
            disp[i][0] -= pts[i][0] * 0.004
            disp[i][1] -= pts[i][1] * 0.004
        for i in range(n_pts):
            dx, dy = disp[i]
            d = math.hypot(dx, dy) or 1e-6
            lim = min(d, temp)
            pts[i][0] += dx / d * lim
            pts[i][1] += dy / d * lim
        temp = max(temp * 0.985, 0.4)

    poly = {n: [(pts[i][0], pts[i][1]) for i in chain[n]] for n in names}

    # ---- pack the connected components ----
    parent = {n: n for n in names}

    def find(a: str) -> str:
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for l in links:
        if l.a in parent and l.b in parent:
            ra, rb = find(l.a), find(l.b)
            if ra != rb:
                parent[ra] = rb
    comps: Dict[str, List[str]] = {}
    for n in names:
        comps.setdefault(find(n), []).append(n)

    def bbox(g: List[str]) -> Tuple[float, float, float, float]:
        xs = [p[0] for n in g for p in poly[n]]
        ys = [p[1] for n in g for p in poly[n]]
        return min(xs), min(ys), max(xs), max(ys)

    groups = sorted(comps.values(), key=lambda g: -sum(by_name[n].length for n in g))
    gap = segment_thickness() * 3.0
    _mx0, _my0, mx1, my1 = bbox(groups[0])
    cur_x, row_y = bbox(groups[0])[0], my1 + gap
    for g in groups[1:]:
        gx0, gy0, gx1, _gy1 = bbox(g)
        dx, dy = cur_x - gx0, row_y - gy0
        for n in g:
            poly[n] = [(x + dx, y + dy) for x, y in poly[n]]
        cur_x += (gx1 - gx0) + gap

    # landscape, so the figure is not forced tall by an arbitrary rotation
    xs = [p[0] for n in names for p in poly[n]]
    ys = [p[1] for n in names for p in poly[n]]
    if (max(ys) - min(ys)) > (max(xs) - min(xs)):
        poly = {n: [(y, -x) for x, y in v] for n, v in poly.items()}

    xs = [p[0] for n in names for p in poly[n]]
    ys = [p[1] for n in names for p in poly[n]]
    pad = 18.0 + segment_thickness() * 0.5
    minx, miny = min(xs), min(ys)
    width = (max(xs) - minx) + 2 * pad
    height = (max(ys) - miny) + 2 * pad
    out = {
        n: [(x - minx + pad, y - miny + pad) for x, y in v] for n, v in poly.items()
    }
    return out, width, height


def _trim_polyline(
    points: List[Tuple[float, float]], trim: float
) -> List[Tuple[float, float]]:
    """
    Shorten a polyline by `trim` at each end.

    Ribbons are stroked with round caps, which add half a stroke width beyond
    each endpoint. Left uncorrected a contig draws one ribbon-width longer than
    it is, which is invisible on a 9 Mb contig and enormous on a 2 kb one - it
    was why a short contig looked several times bigger in the graph than on the
    chromosome. Trimming by half a width first makes the cap land exactly on the
    true endpoint.
    """
    if len(points) < 2:
        return points
    segs = [
        math.hypot(points[i + 1][0] - points[i][0], points[i + 1][1] - points[i][1])
        for i in range(len(points) - 1)
    ]
    total = sum(segs)
    if total <= 2 * trim:
        # too short to trim: collapse to the midpoint and let the caps be it
        mid = len(points) // 2
        return [points[mid], points[mid]]

    def walk(pts, sl, want):
        acc = 0.0
        for i, seg in enumerate(sl):
            if acc + seg >= want:
                t = (want - acc) / seg if seg else 0.0
                x = pts[i][0] + (pts[i + 1][0] - pts[i][0]) * t
                y = pts[i][1] + (pts[i + 1][1] - pts[i][1]) * t
                return [(x, y)] + list(pts[i + 1:])
            acc += seg
        return list(pts[-1:])

    head = walk(points, segs, trim)
    rev = list(reversed(head))
    rsegs = [
        math.hypot(rev[i + 1][0] - rev[i][0], rev[i + 1][1] - rev[i][1])
        for i in range(len(rev) - 1)
    ]
    return list(reversed(walk(rev, rsegs, trim)))


def _smooth_path(points: List[Tuple[float, float]]) -> str:
    """A rounded path through every bead of a contig's polyline."""
    if len(points) == 1:
        x, y = points[0]
        return f"M {x:.1f} {y:.1f} L {x + 0.1:.1f} {y:.1f}"
    if len(points) == 2:
        return (f"M {points[0][0]:.1f} {points[0][1]:.1f} "
                f"L {points[1][0]:.1f} {points[1][1]:.1f}")
    # quadratic through midpoints: each bead becomes a control point, so the
    # curve passes smoothly along the chain instead of cornering at every bead
    d = [f"M {points[0][0]:.1f} {points[0][1]:.1f}"]
    for i in range(1, len(points) - 1):
        cx, cy = points[i]
        mx = (points[i][0] + points[i + 1][0]) / 2.0
        my = (points[i][1] + points[i + 1][1]) / 2.0
        d.append(f"Q {cx:.1f} {cy:.1f} {mx:.1f} {my:.1f}")
    d.append(f"L {points[-1][0]:.1f} {points[-1][1]:.1f}")
    return " ".join(d)


def render_bandage_style_svg(
    calls: List[SegmentCall],
    links: List[GfaLink],
    title: str,
    colours: Dict[str, str],
    args,
    log: Log,
) -> str:
    """The graph drawn Bandage-fashion: thick flexible contigs, force-directed."""
    by_name = {c.name: c for c in calls}
    geom, width, height = bandage_layout(calls, links, args, log)
    colours = colours or assign_segment_colours(calls)

    P = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width:.0f}" height="{height:.0f}" '
        f'viewBox="0 0 {width:.0f} {height:.0f}" font-family="Helvetica, Arial, sans-serif">',
        f'<rect width="100%" height="100%" fill="{PALETTE["bg"]}"/>',
    ]
    if title:
        P.append(
            f'<text x="40" y="40" font-size="{FS_HEADING}" font-weight="600" '
            f'fill="{PALETTE["text"]}">{esc(title)}</text>'
        )

    circular = {l.a for l in links if l.a == l.b}
    w = segment_thickness()

    def terminal(seg: str, end: str) -> Tuple[float, float]:
        return geom[seg][0] if end == "s" else geom[seg][-1]

    # junction connectors, behind the contigs
    P.append(
        f'<g id="layer-links" fill="none" stroke="{PALETTE["bar_edge"]}" '
        f'stroke-linecap="round">'
    )
    for l in links:
        if l.a not in geom or l.b not in geom or l.a == l.b:
            continue
        ax, ay = terminal(l.a, "e" if l.a_orient == "+" else "s")
        bx, by = terminal(l.b, "s" if l.b_orient == "+" else "e")
        P.append(
            f'<line x1="{ax:.1f}" y1="{ay:.1f}" x2="{bx:.1f}" y2="{by:.1f}" '
            f'stroke-width="{w * 0.26:.1f}"/>'
        )
    P.append("</g>")

    # contigs, drawn along their polyline
    P.append('<g id="layer-segments" fill="none">')
    rings: Dict[str, Tuple[float, float, float]] = {}
    for name in sorted(geom, key=lambda n: -by_name[n].length):
        c = by_name[name]
        colour = colours.get(name, "#cfcfcf")
        pointset = geom[name]
        if name in circular:
            seg_len = max(segment_draw_length(c.length, args), 30.0)
            r = max(seg_len / (2 * math.pi), 13.0)
            cx = sum(p[0] for p in pointset) / len(pointset)
            cy = sum(p[1] for p in pointset) / len(pointset)
            rings[name] = (cx, cy, r)
            P.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r:.1f}" '
                     f'stroke="{PALETTE["bar_edge"]}" stroke-width="{w + 2.0:.1f}"/>')
            P.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r:.1f}" '
                     f'stroke="{colour}" stroke-width="{w:.1f}"/>')
            continue
        d = _smooth_path(_trim_polyline(pointset, w / 2.0))
        P.append(f'<path d="{d}" stroke="{PALETTE["bar_edge"]}" stroke-width="{w + 2.0:.1f}" '
                 f'stroke-linecap="round" stroke-linejoin="round"/>')
        P.append(f'<path d="{d}" stroke="{colour}" stroke-width="{w:.1f}" '
                 f'stroke-linecap="round" stroke-linejoin="round"/>')
    P.append("</g>")

    # labels: the contig number inside the ribbon, coverage beside it
    label_all = len(calls) <= args.graph_label_limit
    P.append(f'<g id="layer-labels" fill="{PALETTE["text"]}">')
    for name in sorted(geom):
        c = by_name[name]
        if not label_all and c.cls == "backbone" and c.length < args.backbone_min_length:
            continue
        colour = colours.get(name, "#cfcfcf")
        pointset = geom[name]
        if name in rings:
            cx, cy, r = rings[name]
            mx, my = cx, cy - r
            nx, ny = 0.0, -1.0
        else:
            mid = len(pointset) // 2
            mx, my = pointset[mid]
            a = pointset[max(mid - 1, 0)]
            b = pointset[min(mid + 1, len(pointset) - 1)]
            tx, ty = b[0] - a[0], b[1] - a[1]
            tlen = math.hypot(tx, ty) or 1.0
            nx, ny = -ty / tlen, tx / tlen
        P.append(
            f'<text x="{mx:.1f}" y="{my + FS_LABEL * 0.35:.1f}" text-anchor="middle" '
            f'font-size="{FS_LABEL}" font-weight="700" fill="{_text_on(colour)}">'
            f'{esc(_segment_number(name))}</text>'
        )
        span = sum(
            math.hypot(pointset[i + 1][0] - pointset[i][0], pointset[i + 1][1] - pointset[i][1])
            for i in range(len(pointset) - 1)
        )
        if c.depth is not None and (name in rings or span >= w * 2.2):
            off = w / 2 + FS_SUB + 8
            lx, ly = mx + nx * off, my + ny * off
            P.append(
                f'<text x="{lx:.1f}" y="{ly:.1f}" text-anchor="middle" '
                f'font-size="{FS_SUB}" fill="{PALETTE["muted"]}">{c.depth:.0f}x</text>'
            )
    P.append("</g>")

    P.append("</svg>")
    return "\n".join(P)




GRAPH_COL_W, GRAPH_ROW_H, GRAPH_NODE_H = 210.0, 88.0, 26.0


def graph_node_width(length: int) -> float:
    return 46.0 + 26.0 * math.log10(max(length, 10))


def _graph_layout(
    calls: List[SegmentCall], adj: Dict[str, Set[str]]
) -> Tuple[Dict[str, Tuple[float, float]], float, float, float]:
    """
    Node positions for the graph figure: components stacked top to bottom, and
    within a component, BFS layers from the longest segment left to right.
    Deterministic, so two runs of the tool are directly comparable, and so the
    paired figure can draw leader lines to these exact positions.
    """
    by_name = {c.name: c for c in calls}
    comps: Dict[int, List[SegmentCall]] = defaultdict(list)
    for c in calls:
        comps[c.component].append(c)
    ordered = sorted(comps.values(), key=lambda cs: -sum(c.length for c in cs))

    pos: Dict[str, Tuple[float, float]] = {}
    y_cursor, max_x = 118.0, 0.0
    for cs in ordered:
        names = {c.name for c in cs}
        root = max(cs, key=lambda c: c.length).name
        layer: Dict[str, int] = {root: 0}
        queue = [root]
        while queue:
            cur = queue.pop(0)
            for nb in sorted(adj.get(cur, ())):
                if nb in names and nb not in layer:
                    layer[nb] = layer[cur] + 1
                    queue.append(nb)
        for c in cs:
            layer.setdefault(c.name, 0)
        rows: Dict[int, List[str]] = defaultdict(list)
        for name in sorted(layer, key=lambda n: (layer[n], -by_name[n].length, n)):
            rows[layer[name]].append(name)
        depth_rows = max((len(v) for v in rows.values()), default=1)
        comp_top = y_cursor
        for lidx, members in sorted(rows.items()):
            for k, name in enumerate(members):
                pos[name] = (90.0 + lidx * GRAPH_COL_W, comp_top + k * GRAPH_ROW_H)
                max_x = max(max_x, pos[name][0])
        y_cursor += depth_rows * GRAPH_ROW_H + 46

    return pos, max_x + 260, y_cursor + 150, y_cursor


def render_graph_svg(
    calls: List[SegmentCall],
    links: List[GfaLink],
    model_title: str,
    colours: Optional[Dict[str, str]] = None,
) -> str:
    """
    A redraw of the assembly graph with segments coloured by what we inferred
    them to be, and labels laid out in fixed slots.

    Not a claim that Bandage cannot colour a graph: it colours by depth, and it
    accepts a Color column in a CSV to set nodes explicitly. The difference is
    that these colours are derived from the copy-number classification rather
    than supplied, and the same map drives the chromosome figure, so a node here
    and a block there are the same colour by construction.
    """
    by_name = {c.name: c for c in calls}
    adj = build_adjacency(links)
    colours = colours or assign_segment_colours(calls)
    pos, width, height, y_cursor = _graph_layout(calls, adj)
    P = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width:.0f}" height="{height:.0f}" '
        f'viewBox="0 0 {width:.0f} {height:.0f}" font-family="Helvetica, Arial, sans-serif">',
        f'<rect width="100%" height="100%" fill="{PALETTE["bg"]}"/>',
        f'<text x="60" y="40" font-size="18" font-weight="600" fill="{PALETTE["text"]}">'
        + (
            f'{esc(model_title)} - assembly graph by inferred class'
            if model_title
            else "Assembly graph, as the assembler left it"
        )
        + "</text>",
        f'<text x="60" y="62" font-size="12" fill="{PALETTE["muted"]}">'
        f'Each backbone segment has its own colour, reused in the chromosome figure so a node '
        f'here and a block there can be matched. Colours follow the inferred class; '
        f'layout is deterministic.</text>',
    ]


    # edges first
    for l in links:
        if l.a not in pos or l.b not in pos:
            continue
        x1, y1 = pos[l.a]
        x2, y2 = pos[l.b]
        if l.a == l.b:
            w = graph_node_width(by_name[l.a].length)
            P.append(
                f'<path d="M {x1 + w * 0.35:.1f} {y1 - GRAPH_NODE_H / 2:.1f} '
                f'C {x1 + w * 0.2:.1f} {y1 - GRAPH_NODE_H - 26:.1f}, '
                f'{x1 + w * 0.8:.1f} {y1 - GRAPH_NODE_H - 26:.1f}, '
                f'{x1 + w * 0.65:.1f} {y1 - GRAPH_NODE_H / 2:.1f}" fill="none" stroke="#666" '
                f'stroke-width="1.6"/>'
            )
            continue
        ax = x1 + graph_node_width(by_name[l.a].length) if x2 >= x1 else x1
        bx = x2 if x2 >= x1 else x2 + graph_node_width(by_name[l.b].length)
        P.append(
            f'<path d="{_arc_path(ax, y1, bx, y2)}" fill="none" stroke="#8a8a8a" '
            f'stroke-width="1.4" stroke-opacity="0.8"/>'
        )

    for c in calls:
        if c.name not in pos:
            continue
        x, y = pos[c.name]
        w = graph_node_width(c.length)
        colour = colours.get(c.name, CLASS_COLOUR.get(c.cls, "#cfcfcf"))
        stroke, sw = PALETTE["bar_edge"], 0.9
        if c.at_rich:  # composition flag, shown without overriding the class colour
            stroke, sw = CLASS_COLOUR["at_rich"], 2.6
        P.append(
            f'<rect x="{x:.1f}" y="{y - GRAPH_NODE_H / 2:.1f}" width="{w:.1f}" height="{GRAPH_NODE_H}" '
            f'rx="5" fill="{colour}" fill-opacity="0.88" stroke="{stroke}" '
            f'stroke-width="{sw}"/>'
        )
        P.append(
            f'<text x="{x + w / 2:.1f}" y="{y + 4.5:.1f}" font-size="11.5" text-anchor="middle" '
            f'fill="#ffffff" font-weight="600">{esc(c.name)}</text>'
        )
        cn = f"{c.copy_number:.1f}x copies" if c.copy_number is not None else "copies unknown"
        dp = f"{c.depth:.0f}x depth" if c.depth is not None else "depth unknown"
        for i, line in enumerate((f"{human_bp(c.length)}, {dp}", f"{cn} - {CLASS_LABEL[c.cls]}")):
            P.append(
                f'<text x="{x:.1f}" y="{y + GRAPH_NODE_H / 2 + 13 + i * 12:.1f}" font-size="10" '
                f'fill="{PALETTE["muted"]}">{esc(line)}</text>'
            )

    # legend
    ly = y_cursor + 44
    n_bb = sum(1 for c in calls if c.cls == "backbone")
    P.append(
        f'<text x="60" y="{ly - 18:.1f}" font-size="11" fill="{PALETTE["muted"]}">'
        f'Backbone segments ({n_bb}) each have their own colour, repeated in the chromosome '
        f'figure. Remaining colours are by inferred class:</text>'
    )
    lx = 60.0
    for cls_name in [
        c for c in CLASS_COLOUR if c != "backbone" and any(x.cls == c for x in calls)
    ]:
        P.append(
            f'<rect x="{lx:.1f}" y="{ly - 9:.1f}" width="13" height="13" rx="3" '
            f'fill="{CLASS_COLOUR[cls_name]}"/>'
        )
        P.append(f'<text x="{lx + 19:.1f}" y="{ly + 2:.1f}" font-size="11">'
                 f'{esc(CLASS_LABEL[cls_name])}</text>')
        lx += 34 + 6.6 * len(CLASS_LABEL[cls_name])
        if lx > width - 220:
            lx, ly = 60.0, ly + 22
    P.append("</svg>")
    return "\n".join(P)


def render_graph_figure(
    calls: List[SegmentCall],
    links: List[GfaLink],
    colours: Dict[str, str],
    title: str,
    path: str,
    args,
    log: Log,
) -> Optional[str]:
    """The Bandage graph redrawn with our own, fixed colours."""
    if len(calls) > args.max_graph_nodes:
        log.warn(
            f"the graph has {len(calls)} segments, more than --max-graph-nodes "
            f"({args.max_graph_nodes}), so the graph figure was skipped. It would be unreadable "
            f"at that size; raise the limit if you want it anyway."
        )
        return None
    with open(path, "w") as fh:
        fh.write(graph_svg_for_style(calls, links, title, colours, args, log))
    return path


def graph_svg_for_style(calls, links, title, colours, args, log) -> str:
    if args.graph_style == "bandage":
        return render_bandage_style_svg(calls, links, title, colours, args, log)
    return render_graph_svg(calls, links, title, colours)
