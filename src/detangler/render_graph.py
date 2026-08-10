"""The assembly graph panel."""
from __future__ import annotations

import argparse
import glob
import gzip
import json
import math
import os
import re
import subprocess
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


def _graphviz_positions(
    pts: List[List[float]],
    springs: List[Tuple[int, int, float, float]],
    spacing: float,
    log: Log,
) -> bool:
    """
    Ask neato for an initial layout, in place. Returns False if graphviz is not
    installed, in which case the caller falls back to its own model.

    Two things neato gives that a plain spring model does not: a desired length
    per edge, honoured properly, and `overlap=prism`, which removes node overlap
    as a separate optimisation rather than hoping repulsion sorts it out.
    """
    seen: Dict[Tuple[int, int], float] = {}
    for a, b, rest, _st in springs:
        key = (min(a, b), max(a, b))
        if key not in seen or rest < seen[key]:
            seen[key] = rest
    if not seen:
        return False

    scale = 72.0  # graphviz works in inches; len is in inches, output in points
    lines = ["graph g {", "  node [shape=point, width=0.01, label=\"\"];"]
    for (a, b), rest in sorted(seen.items()):
        lines.append(f"  n{a} -- n{b} [len={rest / scale:.4f}];")
    lines.append("}")
    dot = "\n".join(lines)

    try:
        proc = subprocess.run(
            ["neato", "-Tplain", "-Goverlap=prism", "-Gsep=+6", "-Gmodel=subset"],
            input=dot, capture_output=True, text=True, timeout=120,
        )
    except (OSError, subprocess.SubprocessError):
        return False
    if proc.returncode != 0 or not proc.stdout.strip():
        return False

    got = 0
    for line in proc.stdout.splitlines():
        f = line.split()
        if len(f) >= 4 and f[0] == "node" and f[1].startswith("n"):
            try:
                i = int(f[1][1:])
                pts[i][0] = float(f[2]) * scale
                pts[i][1] = float(f[3]) * scale
            except (ValueError, IndexError):
                continue
            got += 1
    if got < len(pts) * 0.9:
        return False
    log.info(f"graph layout: initial placement from graphviz neato ({got} points)")
    return True


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
            springs.append((idxs[b], idxs[b + 1], rest, 30.0))
        # a brace across every second bead keeps a contig from crumpling into a
        # ball while still letting it curve
        for b in range(beads - 2):
            springs.append((idxs[b], idxs[b + 2], rest * 1.94, 3.0))

    def terminal(seg: str, end: str) -> int:
        return chain[seg][0] if end == "s" else chain[seg][-1]

    # Beads that a spring is holding together must NOT also repel each other.
    # Leaving them in the all-pairs sum inflates every chain - a bead pair at
    # rest length still feels k^2/d of push, so the contig stretches until the
    # spring can match it, and its drawn length stops meaning anything.
    bonded: Set[Tuple[int, int]] = {
        (min(a, b), max(a, b)) for a, b, _rest, _st in springs
    }
    link_pairs: Set[Tuple[int, int]] = set()
    loose_links: Set[Tuple[int, int]] = set()
    for l in links:
        if l.a == l.b or l.a not in chain or l.b not in chain:
            continue
        a_i = terminal(l.a, "e" if l.a_orient == "+" else "s")
        b_i = terminal(l.b, "s" if l.b_orient == "+" else "e")
        if a_i == b_i:
            continue
        # A contig with BOTH ends on the same neighbour is a loop. Held at the
        # normal junction length it folds into a hairpin lying on top of itself;
        # it needs room to come back round.
        loop = sum(
            1 for o in links
            if o.a != o.b and {o.a, o.b} == {l.a, l.b}
        ) > 1
        rest_here = spacing * (3.2 if loop else 0.9)
        springs.append((a_i, b_i, rest_here, 6.0 if loop else 12.0))
        link_pairs.add((min(a_i, b_i), max(a_i, b_i)))
        bonded.add((min(a_i, b_i), max(a_i, b_i)))
        if loop:
            loose_links.add((min(a_i, b_i), max(a_i, b_i)))

    # ---- initial placement, by graphviz when it is available ----
    # neato is a mature implementation of exactly this problem: it honours a
    # desired length per edge (`len`), and `overlap=prism` removes node overlap
    # properly afterwards. Bandage leans on OGDF's FMMM for the same reason. The
    # hand-rolled model below still runs, but starting from a good global layout
    # is worth far more than any amount of tuning it.
    if _graphviz_positions(pts, springs, spacing, log):
        iters_scale = 0.25
    else:
        iters_scale = 1.0

    n_pts = len(pts)
    k = spacing * 1.15
    iters = int(min(600, max(120, 26000 / max(n_pts, 1))) * iters_scale)
    log.info(f"graph layout: {len(names)} contigs as {n_pts} points, {iters} iterations")
    temp = spacing * 2.0

    for _step in range(iters):
        disp = [[0.0, 0.0] for _ in range(n_pts)]
        for i in range(n_pts):
            xi, yi = pts[i]
            for j in range(i + 1, n_pts):
                if (i, j) in bonded:
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

    # ---- enforce the drawn lengths ----
    # A force model cannot guarantee a length: a bead pair sitting at its rest
    # length still gets pushed by every other bead, so chains stretch, and a
    # contig's drawn length stops meaning anything. This is a constraint
    # projection pass - walk each chain and move consecutive beads back to
    # exactly their rest separation, repeatedly. It converges in a few dozen
    # passes and barely disturbs the shape the force model found, but it makes
    # the drawn length exact, which is the whole point of drawing it to scale.
    chain_rest = {
        n: (segment_draw_length(by_name[n].length, args) / max(len(chain[n]) - 1, 1))
        for n in names
    }
    link_rest = spacing * 0.9
    w_sep = segment_thickness() * 1.12
    owner = {}
    for n in names:
        for i in chain[n]:
            owner[i] = n
    for _pass in range(120):
        worst = 0.0
        for n in names:
            idxs = chain[n]
            rest = chain_rest[n]
            for b in range(len(idxs) - 1):
                i, j = idxs[b], idxs[b + 1]
                dx = pts[j][0] - pts[i][0]
                dy = pts[j][1] - pts[i][1]
                d = math.hypot(dx, dy) or 1e-6
                corr = (d - rest) / d * 0.5
                worst = max(worst, abs(d - rest))
                pts[i][0] += dx * corr
                pts[i][1] += dy * corr
                pts[j][0] -= dx * corr
                pts[j][1] -= dy * corr
        # Only RUNAWAY junctions are reined in. A junction stretched right across
        # the drawing has to pass under whatever lies between its two ends, which
        # is how a connector ended up running beneath contig 2. Pulling every
        # junction shut instead collapses the whole graph onto its hubs and the
        # labels pile up, so this fires well above the resting length and stops
        # at a comfortable one.
        # Contigs must not lie on top of one another. Repulsion during the force
        # phase discourages it but cannot forbid it, and the length constraints
        # above can push two ribbons back into each other. This is a hard
        # separation: any two beads on DIFFERENT contigs closer than a ribbon
        # width get pushed apart, so no two ribbons can overlap.
        clear = w_sep
        for i in range(n_pts):
            for j in range(i + 1, n_pts):
                if (i, j) in bonded or owner[i] == owner[j]:
                    continue
                dx = pts[j][0] - pts[i][0]
                dy = pts[j][1] - pts[i][1]
                d2 = dx * dx + dy * dy
                if d2 >= clear * clear:
                    continue
                d = math.sqrt(d2) or 1e-6
                corr = (d - clear) / d * 0.5
                worst = max(worst, clear - d)
                pts[i][0] += dx * corr
                pts[i][1] += dy * corr
                pts[j][0] -= dx * corr
                pts[j][1] -= dy * corr

        far = link_rest * 2.6
        target = link_rest * 2.0
        for i, j in link_pairs:
            if (i, j) in loose_links:
                continue
            dx = pts[j][0] - pts[i][0]
            dy = pts[j][1] - pts[i][1]
            d = math.hypot(dx, dy) or 1e-6
            if d <= far:
                continue
            corr = (d - target) / d * 0.25
            worst = max(worst, d - target)
            pts[i][0] += dx * corr
            pts[i][1] += dy * corr
            pts[j][0] -= dx * corr
            pts[j][1] -= dy * corr
        if worst < 0.05:
            break

    # ---- smooth the chains ----
    # The constraint passes fix distances but not directions, so a chain settles
    # with a small kink at every bead and the contig looks wobbly rather than
    # drawn. Each interior bead is pulled hard towards the midpoint of its
    # neighbours and the spacing is then restored; the two ends are held, so
    # junctions do not move.
    # A contig laid out between two fixed ends settles as a straight run, which
    # is accurate but reads as a wire diagram. Each chain is bowed to one side
    # first - alternating, so neighbouring contigs curve away from each other -
    # and the smoothing below then turns that into a clean arc rather than a
    # kinked one.
    for si, n in enumerate(sorted(names)):
        idxs = chain[n]
        if len(idxs) < 3:
            continue
        ax_, ay_ = pts[idxs[0]]
        bx_, by_ = pts[idxs[-1]]
        vx_, vy_ = bx_ - ax_, by_ - ay_
        vlen_ = math.hypot(vx_, vy_) or 1.0
        nx_, ny_ = -vy_ / vlen_, vx_ / vlen_
        amp = vlen_ * 0.22 * (1 if si % 2 else -1)
        for b in range(1, len(idxs) - 1):
            t = b / (len(idxs) - 1)
            k_ = math.sin(math.pi * t)
            pts[idxs[b]][0] += nx_ * amp * k_
            pts[idxs[b]][1] += ny_ * amp * k_

    for _round in range(60):
        for n in names:
            idxs = chain[n]
            if len(idxs) < 3:
                continue
            new_xy = []
            for b in range(1, len(idxs) - 1):
                i, prv, nxt = idxs[b], idxs[b - 1], idxs[b + 1]
                mx_ = (pts[prv][0] + pts[nxt][0]) / 2.0
                my_ = (pts[prv][1] + pts[nxt][1]) / 2.0
                new_xy.append((i, pts[i][0] * 0.25 + mx_ * 0.75,
                                  pts[i][1] * 0.25 + my_ * 0.75))
            for i, nx_, ny_ in new_xy:
                pts[i][0], pts[i][1] = nx_, ny_
        for n in names:
            idxs = chain[n]
            rest = chain_rest[n]
            for _sub in range(3):
                for b in range(len(idxs) - 1):
                    i, j = idxs[b], idxs[b + 1]
                    dx = pts[j][0] - pts[i][0]
                    dy = pts[j][1] - pts[i][1]
                    d = math.hypot(dx, dy) or 1e-6
                    corr = (d - rest) / d * 0.5
                    pts[i][0] += dx * corr
                    pts[i][1] += dy * corr
                    pts[j][0] -= dx * corr
                    pts[j][1] -= dy * corr

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
    gap = segment_thickness() * 1.5
    mx0, my0, mx1, my1 = bbox(groups[0])
    # Stacked BESIDE the main component, not beneath it. The figure as a whole is
    # graph above chromosomes, so height is the scarce dimension - a row of
    # isolated contigs underneath made the whole thing too tall to view at once,
    # while there was empty space to the side.
    col_x = mx1 + gap
    cur_y = my0
    for g in groups[1:]:
        gx0, gy0, _gx1, gy1 = bbox(g)
        dx, dy = col_x - gx0, cur_y - gy0
        for n in g:
            poly[n] = [(x + dx, y + dy) for x, y in poly[n]]
        cur_y += (gy1 - gy0) + gap

    # Rotate so the layout's LONG axis lies horizontal. A spring layout comes out
    # at an arbitrary angle, and the figure stacks the graph above the
    # chromosomes, so height is the scarce dimension. A quarter turn only helps
    # when the spread happens to be axis-aligned; the principal axis is the
    # general answer.
    allp = [p for n in names for p in poly[n]]
    if len(allp) > 2:
        cx0 = sum(p[0] for p in allp) / len(allp)
        cy0 = sum(p[1] for p in allp) / len(allp)
        sxx = sum((p[0] - cx0) ** 2 for p in allp)
        syy = sum((p[1] - cy0) ** 2 for p in allp)
        sxy = sum((p[0] - cx0) * (p[1] - cy0) for p in allp)
        theta = 0.5 * math.atan2(2 * sxy, sxx - syy)
        ct, st = math.cos(-theta), math.sin(-theta)
        poly = {
            n: [
                ((x - cx0) * ct - (y - cy0) * st, (x - cx0) * st + (y - cy0) * ct)
                for x, y in v
            ]
            for n, v in poly.items()
        }
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
    colours = colours or assign_segment_colours(calls, links)

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

    def poly_len(pointset: List[Tuple[float, float]]) -> float:
        return sum(
            math.hypot(pointset[i + 1][0] - pointset[i][0],
                       pointset[i + 1][1] - pointset[i][1])
            for i in range(len(pointset) - 1)
        )

    def raw_end(seg: str, end: str) -> Tuple[float, float]:
        return geom[seg][0] if end == "s" else geom[seg][-1]

    # Which physical ends of each contig actually carry links, and what they
    # attach to. A contig with links on ONE end only is a tip - the neighbours
    # meet at that end, they do not pass through - and the drawing has to say so.
    ends_used: Dict[str, Set[str]] = defaultdict(set)
    partners: Dict[Tuple[str, str], List[Tuple[float, float]]] = defaultdict(list)
    for l in links:
        if l.a == l.b or l.a not in geom or l.b not in geom:
            continue
        ae = "e" if l.a_orient == "+" else "s"
        be = "s" if l.b_orient == "+" else "e"
        ends_used[l.a].add(ae)
        ends_used[l.b].add(be)
        partners[(l.a, ae)].append(raw_end(l.b, be))
        partners[(l.b, be)].append(raw_end(l.a, ae))

    # A contig shorter than one ribbon width can only be drawn as a dot, since a
    # round-capped stroke is never shorter than its own width. Where it is a tip,
    # the dot is pushed off the junction along the direction away from its
    # neighbours, so it reads as hanging OFF the join rather than sitting in the
    # middle of it. Centring it on the junction would draw the neighbours as
    # passing through - the very join this graph does not support.
    dot_centre: Dict[str, Tuple[float, float]] = {}
    for name, pointset in geom.items():
        if name in circular or len(pointset) < 2:
            continue
        # a shade above one width: the constraint pass lands chains within a
        # small tolerance of their target, so an exact comparison misses the
        # very segments this is for
        if poly_len(pointset) > w * 1.15:
            continue
        used = ends_used.get(name, set())
        mid = (
            (pointset[0][0] + pointset[-1][0]) / 2.0,
            (pointset[0][1] + pointset[-1][1]) / 2.0,
        )
        if len(used) == 1:
            live = next(iter(used))
            jx, jy = raw_end(name, live)
            near = partners.get((name, live), [])
            if near:
                cx = sum(p[0] for p in near) / len(near)
                cy = sum(p[1] for p in near) / len(near)
                dx, dy = jx - cx, jy - cy
                d = math.hypot(dx, dy) or 1.0
                dot_centre[name] = (jx + dx / d * w / 2.0, jy + dy / d * w / 2.0)
                continue
        dot_centre[name] = mid

    allpts = [p for v in geom.values() for p in v]
    gcx = sum(p[0] for p in allpts) / max(len(allpts), 1)
    gcy = sum(p[1] for p in allpts) / max(len(allpts), 1)

    def terminal(seg: str, end: str) -> Tuple[float, float]:
        # a tip is joined AT its live end, on the rim of the dot; a short contig
        # that really is a bridge keeps its centre so both sides meet in it
        if seg in dot_centre and len(ends_used.get(seg, set())) == 1:
            return raw_end(seg, end)
        if seg in dot_centre:
            return dot_centre[seg]
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
        # Bowed, not straight. A junction drawn as a straight line reads as a
        # ruled edge in a diagram; the contigs around it are all curves, and the
        # straight lines were the thing making the panel look mechanical. The bow
        # goes away from the centre of the drawing so connectors at a hub splay
        # rather than fold over one another.
        mxx, myy = (ax + bx) / 2.0, (ay + by) / 2.0
        vx, vy = mxx - gcx, myy - gcy
        vlen = math.hypot(vx, vy) or 1.0
        chord = math.hypot(bx - ax, by - ay)
        bow = min(chord * 0.55, w * 2.4)
        qx, qy = mxx + vx / vlen * bow, myy + vy / vlen * bow
        P.append(
            f'<path d="M {ax:.1f} {ay:.1f} Q {qx:.1f} {qy:.1f} {bx:.1f} {by:.1f}" '
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
            # the ring must have a visible hole, or it renders as a solid blob
            # and stops being distinguishable from a short linear contig
            seg_len = max(segment_draw_length(c.length, args), 30.0)
            ring_w = w * 0.55
            r = max(seg_len / (2 * math.pi), w * 0.85)
            cx = sum(p[0] for p in pointset) / len(pointset)
            cy = sum(p[1] for p in pointset) / len(pointset)
            rings[name] = (cx, cy, r)
            P.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r:.1f}" '
                     f'stroke="{colour}" stroke-width="{ring_w:.1f}"/>')
            continue
        if name in dot_centre:
            # Drawn as a short BAR with square ends, not a disc. A disc reads as
            # a circular contig - which is a real and different thing in a graph
            # (edge_11, a mitochondrion) - so a linear contig too short to draw
            # at scale must not borrow that shape.
            cx, cy = dot_centre[name]
            ex, ey = pointset[-1][0] - pointset[0][0], pointset[-1][1] - pointset[0][1]
            elen = math.hypot(ex, ey) or 1.0
            ux, uy = ex / elen, ey / elen
            half = max(poly_len(pointset), w) / 2.0
            d = (f"M {cx - ux * half:.1f} {cy - uy * half:.1f} "
                 f"L {cx + ux * half:.1f} {cy + uy * half:.1f}")
            P.append(f'<path d="{d}" stroke="{colour}" stroke-width="{w:.1f}" '
                     f'stroke-linecap="butt"/>')
            continue
        d = _smooth_path(_trim_polyline(pointset, w / 2.0))
        P.append(f'<path d="{d}" stroke="{colour}" stroke-width="{w:.1f}" '
                 f'stroke-linecap="round" stroke-linejoin="round"/>')
    P.append("</g>")

    # labels: the contig number inside the ribbon, coverage beside it
    label_all = len(calls) <= args.graph_label_limit
    placed_labels: List[Tuple[float, float]] = []
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
            mx, my = dot_centre.get(name, pointset[mid])
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
        # Coverage is shown for EVERY contig. It used to be suppressed on short
        # ones for want of room, which quietly hid it on exactly the segments
        # where depth decides the call - a repeat, a tip, an organelle.
        if c.depth is not None:
            # Coverage labels crowd badly around a hub, where several short
            # contigs meet. Each label tries a ring of candidate positions and
            # takes the one furthest from the labels already placed.
            off = w / 2 + FS_SUB + 12
            best_pt, best_score = None, -1.0
            for k in range(12):
                ang = 2 * math.pi * k / 12
                ox_, oy_ = math.cos(ang), math.sin(ang)
                # bias towards the ribbon's normal, so a label still reads as
                # belonging to its own contig
                bias = 1.0 + 0.6 * (ox_ * nx + oy_ * ny)
                cxp, cyp = mx + ox_ * off, my + oy_ * off
                near = min(
                    (math.hypot(cxp - px, cyp - py) for px, py in placed_labels),
                    default=1e6,
                )
                score = min(near, 400.0) * bias
                if score > best_score:
                    best_pt, best_score = (cxp, cyp), score
            lx, ly = best_pt  # type: ignore
            placed_labels.append((lx, ly))
            P.append(
                f'<text x="{lx:.1f}" y="{ly:.1f}" text-anchor="middle" '
                f'font-size="{FS_SUB + 5}" font-weight="600" '
                f'fill="{PALETTE["text"]}">{c.depth:.0f}x</text>'
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
    colours = colours or assign_segment_colours(calls, links)
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
