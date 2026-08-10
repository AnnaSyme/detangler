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
    floor = segment_thickness() * 1.6
    px_per_bp = getattr(args, "graph_px_per_bp", None)
    if px_per_bp:
        # Same bases-per-pixel as the chromosome panel, so a contig is the same
        # size in both. Anything too short to draw at that scale is clamped to
        # the floor and is therefore drawn LARGER than true - unavoidable if a
        # 2 kb segment is to be visible beside a 9 Mb one, but it only ever
        # overstates the small ones.
        return max(length * px_per_bp, floor)
    return min(max(10.0 + args.graph_length_scale * math.sqrt(max(length, 1)), floor),
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
) -> Tuple[Dict[str, Tuple[float, float, float, float]], float, float]:
    """
    Force-directed layout over segment ENDS, not segment centres, which is what
    gives the Bandage look: each segment is a stiff spring of its own drawn
    length, and each link is a short spring tying one segment's end to another's.

    Returns {segment: (x1, y1, x2, y2)} plus the canvas size.
    """
    by_name = {c.name: c for c in calls}
    names = sorted(by_name)
    if not names:
        return {}, 100.0, 100.0

    # two point masses per segment: its start (+) and its end (-)
    pts: List[str] = []
    for n in names:
        pts += [n + "\x00s", n + "\x00e"]
        idx = {p: i for i, p in enumerate(pts)}

    n_pts = len(pts)
    rnd = _Rand(20260809)
    # deterministic ring start, jittered, so components unfold rather than
    # starting on top of one another
    radius = 40.0 + 9.0 * math.sqrt(n_pts)
    pos = []
    for i in range(n_pts):
        a = 2 * math.pi * i / n_pts
        pos.append([
            radius * math.cos(a) + rnd.uniform(-8, 8),
            radius * math.sin(a) + rnd.uniform(-8, 8),
        ])

    springs: List[Tuple[int, int, float, float]] = []  # a, b, rest, strength
    for n in names:
        springs.append((idx[n + "\x00s"], idx[n + "\x00e"],
                        segment_draw_length(by_name[n].length, args), 1.0))
    for l in links:
        if l.a not in by_name or l.b not in by_name:
            continue
        # a link leaves the end of a + oriented segment and enters the start of
        # the next; a - orientation flips which terminal is involved
        a_pt = l.a + ("\x00e" if l.a_orient == "+" else "\x00s")
        b_pt = l.b + ("\x00s" if l.b_orient == "+" else "\x00e")
        if a_pt == b_pt:
            continue
        # A junction is given a RADIUS rather than being a single point. Pulling
        # every linked end onto one coordinate made four contigs meeting at the
        # same hub (edge_9 has five ends on one side) collapse into a pile you
        # could not read. Held about a ribbon-width apart, they spread into an
        # arc and each join shows as its own short connector.
        springs.append((idx[a_pt], idx[b_pt], segment_thickness() * 1.15, 4.0))

    k = max(radius / max(math.sqrt(n_pts), 1.0), 135.0)
    iters = int(min(500, max(80, 9000 / max(n_pts, 1))))
    if n_pts > 400:
        log.info(f"graph layout: {n_pts} endpoints, {iters} iterations (this can take a moment)")
    temp = radius * 0.35

    # Endpoint pairs joined by a link must be free to touch. Left in the
    # all-pairs repulsion they settle at k^2/d against the spring, which for
    # k=46 parks them 40-135 px apart and turns every junction into a long bar.
    linked_pairs = {
        (min(a, b), max(a, b)) for a, b, rest, _s in springs
        if rest <= segment_thickness() * 1.2
    }

    for step in range(iters):
        disp = [[0.0, 0.0] for _ in range(n_pts)]
        # repulsion, all pairs except those a link is trying to hold together
        for i in range(n_pts):
            xi, yi = pos[i]
            for j in range(i + 1, n_pts):
                if (i, j) in linked_pairs:
                    continue
                dx = xi - pos[j][0]
                dy = yi - pos[j][1]
                d2 = dx * dx + dy * dy
                if d2 < 1e-6:
                    dx, dy, d2 = rnd.uniform(-1, 1), rnd.uniform(-1, 1), 1.0
                d = math.sqrt(d2)
                f = (k * k) / d
                ux, uy = dx / d, dy / d
                disp[i][0] += ux * f
                disp[i][1] += uy * f
                disp[j][0] -= ux * f
                disp[j][1] -= uy * f
        # springs
        for a, b, rest, strength in springs:
            dx = pos[a][0] - pos[b][0]
            dy = pos[a][1] - pos[b][1]
            d = math.hypot(dx, dy) or 1e-6
            f = strength * (d - rest) * 0.9
            ux, uy = dx / d, dy / d
            disp[a][0] -= ux * f
            disp[a][1] -= uy * f
            disp[b][0] += ux * f
            disp[b][1] += uy * f
        # pull to the centre so detached components (an unplaced contig, an
        # organelle) stay near the main mass instead of stranding themselves in a
        # far corner and stretching the canvas around a lot of white space
        for i in range(n_pts):
            disp[i][0] -= pos[i][0] * 0.038
            disp[i][1] -= pos[i][1] * 0.038
        # move, capped by the cooling temperature
        for i in range(n_pts):
            dx, dy = disp[i]
            d = math.hypot(dx, dy) or 1e-6
            lim = min(d, temp)
            pos[i][0] += dx / d * lim
            pos[i][1] += dy / d * lim
        temp = max(temp * 0.965, 0.6)

    raw = {
        n: (
            pos[idx[n + "\x00s"]][0], pos[idx[n + "\x00s"]][1],
            pos[idx[n + "\x00e"]][0], pos[idx[n + "\x00e"]][1],
        )
        for n in names
    }

    # ---- pack the connected components ----
    # A spring model left to itself flings a detached contig or an organelle into
    # a far corner, and the canvas then has to stretch around all that white
    # space. The components are laid out independently and then packed: the
    # largest keeps its position, the rest are set apart in a row beneath it.
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
        cx = [v for n in g for v in (raw[n][0], raw[n][2])]
        cy = [v for n in g for v in (raw[n][1], raw[n][3])]
        return min(cx), min(cy), max(cx), max(cy)

    thick = segment_thickness()
    groups = sorted(comps.values(), key=lambda g: -sum(by_name[n].length for n in g))
    placed: Dict[str, Tuple[float, float, float, float]] = {n: raw[n] for n in groups[0]}
    mx0, _my0, _mx1, my1 = bbox(groups[0])
    gap = thick * 3.0
    cur_x, row_y = mx0, my1 + gap
    for g in groups[1:]:
        gx0, gy0, gx1, _gy1 = bbox(g)
        dx, dy = cur_x - gx0, row_y - gy0
        for n in g:
            x1, y1, x2, y2 = raw[n]
            placed[n] = (x1 + dx, y1 + dy, x2 + dx, y2 + dy)
        cur_x += (gx1 - gx0) + gap

    # Orient the graph so its long axis is horizontal. A spring layout comes out
    # in an arbitrary rotation, and a tall one forces a tall figure: the
    # chromosome panel beside it is much shorter, so most of the canvas ends up
    # empty and everything has to be shrunk to fit a preview.
    xs0 = [v for t in placed.values() for v in (t[0], t[2])]
    ys0 = [v for t in placed.values() for v in (t[1], t[3])]
    if xs0 and (max(ys0) - min(ys0)) > (max(xs0) - min(xs0)):
        placed = {
            n: (y1, -x1, y2, -x2) for n, (x1, y1, x2, y2) in placed.items()
        }

    # normalise into a padded canvas
    xs = [v for t in placed.values() for v in (t[0], t[2])]
    ys = [v for t in placed.values() for v in (t[1], t[3])]
    pad = 18.0 + thick * 0.5
    minx, maxx, miny, maxy = min(xs), max(xs), min(ys), max(ys)
    width = (maxx - minx) + 2 * pad
    height = (maxy - miny) + 2 * pad
    out: Dict[str, Tuple[float, float, float, float]] = {}
    for n in names:
        x1, y1, x2, y2 = placed[n]
        out[n] = (
            x1 - minx + pad, y1 - miny + pad,
            x2 - minx + pad, y2 - miny + pad,
        )
    return out, width, height


def render_bandage_style_svg(
    calls: List[SegmentCall],
    links: List[GfaLink],
    title: str,
    colours: Dict[str, str],
    args,
    log: Log,
) -> str:
    """The graph drawn Bandage-fashion: thick tapered segments, force-directed."""
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

    # Parallel segments - two contigs running between the same pair of junctions -
    # land on almost the same chord and the second one disappears underneath the
    # first. Bucket by the endpoints they share and fan the bow out, alternating
    # sign, so each is visible. Without this, edge_7 hides entirely behind edge_2.
    # Detected from the GRAPH, not from the drawn coordinates: two segments are
    # parallel when they have the same set of neighbours. Bucketing by pixel
    # position looked simpler but is far too brittle - edge_2 and edge_7 land in
    # adjacent buckets and stack anyway.
    neighbours_of: Dict[str, Set[str]] = defaultdict(set)
    for l in links:
        if l.a == l.b:
            continue
        neighbours_of[l.a].add(l.b)
        neighbours_of[l.b].add(l.a)

    parallel: Dict[Tuple[str, ...], List[str]] = defaultdict(list)
    for name in sorted(geom):
        nb = neighbours_of.get(name, set())
        if len(nb) >= 2:
            parallel[tuple(sorted(nb))].append(name)
    bow_slot: Dict[str, int] = {}
    label_t: Dict[str, float] = {}
    for group in parallel.values():
        if len(group) < 2:
            continue
        # 0, +1, -1, +2, -2 ... so the bundle spreads either side of the chord
        for i, nm in enumerate(sorted(group)):
            bow_slot[nm] = ((i + 1) // 2) * (1 if i % 2 else -1)
            # and stagger the labels ALONG the ribbons. Fanning separates the
            # middles but the ends still converge, so labels placed at the same
            # fraction of two bundled segments collide however wide the fan.
            label_t[nm] = 0.30 + 0.40 * (i / max(len(group) - 1, 1))

    # Segments carrying a self-link are circular molecules and are drawn as rings
    # rather than as ribbons with a loop hanging off one end.
    circular = {l.a for l in links if l.a == l.b}
    w = segment_thickness()

    # Junction stubs, behind the segments. Linked ends already abut after layout,
    # so a link is a short dark connector rather than a long thin line.
    P.append(
        f'<g id="layer-links" fill="none" stroke="{PALETTE["bar_edge"]}" '
        f'stroke-linecap="round">'
    )
    for l in links:
        if l.a not in geom or l.b not in geom or l.a == l.b:
            continue
        ax, ay = (geom[l.a][2], geom[l.a][3]) if l.a_orient == "+" else (geom[l.a][0], geom[l.a][1])
        bx, by = (geom[l.b][0], geom[l.b][1]) if l.b_orient == "+" else (geom[l.b][2], geom[l.b][3])
        # Thin. The ribbons are drawn with ROUND ends, so two of them meeting at
        # an angle no longer leave a white wedge at the corner and the connector
        # does not have to be wide enough to cover one. A junction should read as
        # a join, not as another piece of sequence.
        P.append(
            f'<line x1="{ax:.1f}" y1="{ay:.1f}" x2="{bx:.1f}" y2="{by:.1f}" '
            f'stroke-width="{w * 0.26:.1f}"/>'
        )
    P.append("</g>")

    # segments as thick ribbons of uniform width
    P.append('<g id="layer-segments" fill="none">')
    rings: Dict[str, Tuple[float, float, float]] = {}
    for name, (x1, y1, x2, y2) in sorted(geom.items(), key=lambda kv: -by_name[kv[0]].length):
        c = by_name[name]
        colour = colours.get(name, "#cfcfcf")
        if name in circular:
            # a ring whose circumference matches the drawn length of the segment
            seg_len = max(segment_draw_length(c.length, args), 30.0)
            r = max(seg_len / (2 * math.pi), 13.0)
            cx, cy = (x1 + x2) / 2, (y1 + y2) / 2
            rings[name] = (cx, cy, r)
            P.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r:.1f}" '
                     f'stroke="{PALETTE["bar_edge"]}" stroke-width="{w + 2.0:.1f}"/>')
            P.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r:.1f}" '
                     f'stroke="{colour}" stroke-width="{w:.1f}"/>')
            continue
        # a gentle bow, so nothing looks like a ruler
        mx, my = (x1 + x2) / 2, (y1 + y2) / 2
        nx, ny = -(y2 - y1), (x2 - x1)
        nlen = math.hypot(nx, ny) or 1.0
        bow = min(26.0, math.hypot(x2 - x1, y2 - y1) * 0.14)
        slot = bow_slot.get(name, 0)
        if slot:
            # fan a bundle of parallel segments apart rather than stacking them
            bow = slot * max(
                abs(bow),
                min(math.hypot(x2 - x1, y2 - y1) * 0.42, segment_thickness() * 6.0),
            )
        cx, cy = mx + nx / nlen * bow, my + ny / nlen * bow
        d = f"M {x1:.1f} {y1:.1f} Q {cx:.1f} {cy:.1f} {x2:.1f} {y2:.1f}"
        # dark casing then the colour, flat butt caps so segments abut cleanly
        P.append(f'<path d="{d}" stroke="{PALETTE["bar_edge"]}" stroke-width="{w + 2.0:.1f}" '
                 f'stroke-linecap="round"/>')
        P.append(f'<path d="{d}" stroke="{colour}" stroke-width="{w:.1f}" '
                 f'stroke-linecap="round"/>')
    P.append("</g>")

    # labels: the segment number sits INSIDE the ribbon, coverage beside it
    label_all = len(calls) <= args.graph_label_limit
    P.append(f'<g id="layer-labels" fill="{PALETTE["text"]}">')
    for name, (x1, y1, x2, y2) in geom.items():
        c = by_name[name]
        if not label_all and c.cls == "backbone" and c.length < args.backbone_min_length:
            continue
        colour = colours.get(name, "#cfcfcf")
        if name in rings:
            cx, cy, r = rings[name]
            nx, ny, mx, my = 0.0, -1.0, cx, cy - r
        else:
            cxm, cym = (x1 + x2) / 2, (y1 + y2) / 2
            nx, ny = -(y2 - y1), (x2 - x1)
            nlen = math.hypot(nx, ny) or 1.0
            nx, ny = nx / nlen, ny / nlen
            bow = min(26.0, math.hypot(x2 - x1, y2 - y1) * 0.14)
            slot = bow_slot.get(name, 0)
            if slot:
                bow = slot * max(
                    abs(bow),
                    min(math.hypot(x2 - x1, y2 - y1) * 0.42, segment_thickness() * 6.0),
                )
            # Evaluate the drawn quadratic Bezier at this segment's own t, and
            # take the normal from the tangent there. Bundled segments get
            # different t values so their labels never stack.
            qx, qy = cxm + nx * bow, cym + ny * bow  # the control point
            t = label_t.get(name, 0.5)
            u = 1.0 - t
            mx = u * u * x1 + 2 * u * t * qx + t * t * x2
            my = u * u * y1 + 2 * u * t * qy + t * t * y2
            tx = 2 * u * (qx - x1) + 2 * t * (x2 - qx)
            ty = 2 * u * (qy - y1) + 2 * t * (y2 - qy)
            tlen = math.hypot(tx, ty) or 1.0
            nx, ny = -ty / tlen, tx / tlen
        # number inside the ribbon, inked for contrast against its own colour
        P.append(
            f'<text x="{mx:.1f}" y="{my + FS_LABEL * 0.35:.1f}" text-anchor="middle" '
            f'font-size="{FS_LABEL}" font-weight="700" fill="{_text_on(colour)}">'
            f'{esc(_segment_number(name))}</text>'
        )
        # coverage only, set clear of the ribbon. Skipped on very short segments,
        # where there is no room for it to sit anywhere it would not collide.
        drawn_len = math.hypot(x2 - x1, y2 - y1)
        if c.depth is not None and (name in rings or drawn_len >= w * 2.2):
            off = w / 2 + FS_SUB + 8
            lx, ly = mx + nx * off, my + ny * off
            P.append(
                f'<text x="{lx:.1f}" y="{ly:.1f}" text-anchor="middle" '
                f'font-size="{FS_SUB}" fill="{PALETTE["muted"]}">{c.depth:.0f}x</text>'
            )
    P.append("</g>")

    P.append("</svg>")
    return "\n".join(P)


def _graph_legend_svg(calls: List[SegmentCall], height: float) -> str:
    out = [f'<g id="legend" font-size="10" fill="{PALETTE["text"]}">']
    x, y = 40.0, height - 24.0
    n_bb = sum(1 for c in calls if c.cls == "backbone")
    out.append(
        f'<text x="{x}" y="{y - 16:.1f}" fill="{PALETTE["muted"]}" font-size="9.5">'
        f'{n_bb} backbone segment(s), each its own colour and reused in the chromosome figure. '
        f'Other colours are by inferred class:</text>'
    )
    for cls_name in [c for c in CLASS_COLOUR if c != "backbone" and any(x2.cls == c for x2 in calls)]:
        out.append(f'<rect x="{x:.1f}" y="{y - 8:.1f}" width="11" height="11" rx="2" '
                   f'fill="{CLASS_COLOUR[cls_name]}"/>')
        out.append(f'<text x="{x + 16:.1f}" y="{y + 1:.1f}">{esc(CLASS_LABEL[cls_name])}</text>')
        x += 30 + 6.0 * len(CLASS_LABEL[cls_name])
    out.append("</g>")
    return "\n".join(out)


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
