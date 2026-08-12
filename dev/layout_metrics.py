#!/usr/bin/env python3
"""
Measure how legible a detangler graph panel is.

Every number here is DIMENSIONLESS - distances in ribbon widths, angles in
degrees - so a figure of eleven contigs and a figure of six hundred are scored
on the same scale. That is the whole point: a layout constant tuned to look
nice on one fungal graph tells you nothing about a plant one.

Metrics
    clearance     smallest gap between the centrelines of two DIFFERENT contigs,
                  in ribbon widths. 1.0 means the ribbons just touch. Below 1.0
                  they overlap. The figure reads as crowded well before 1.0 -
                  around 1.1 two ribbons are separated by a hairline.
    r_cross       ribbon x ribbon centreline crossings. Must be 0.
    rc_cross      ribbon x connector crossings. Must be 0.
    winding       how much a ribbon curls: total absolute turning along its
                  centreline, in degrees. THE metric that matches human
                  judgement of "is this clear". A gentle arc is ~40 deg; a
                  contig folded back on itself is 200+. Reported as the mean
                  over ribbons and the worst single one.
    hub_gap       smallest angle between two connectors leaving the same hub,
                  in degrees. Reported, but NOT scored: it turned out to rank
                  a figure people call clear BELOW one they call crowded, so
                  it does not measure what it looks like it measures.
    fill          fraction of the canvas covered by ribbon ink. Guards against
                  the degenerate fix of spreading everything out until nothing
                  touches: that scores well on clearance and looks empty.
    area          canvas area in ribbon-widths squared. Same guard, absolute.

Usage
    python3 layout_metrics.py FIGURE.svg [FIGURE.svg ...]
    python3 layout_metrics.py --json FIGURE.svg

Reads the SVG that detangler already writes. Nothing to install.
"""

from __future__ import annotations

import json
import math
import re
import sys
from typing import Dict, List, Optional, Sequence, Tuple

Point = Tuple[float, float]

# ---------------------------------------------------------------------------
# reading the SVG
# ---------------------------------------------------------------------------

_NUM = re.compile(r"-?\d*\.?\d+(?:e-?\d+)?")


def _nums(s: str) -> List[float]:
    return [float(x) for x in _NUM.findall(s)]


def sample_path(d: str, per_curve: int = 12) -> List[Point]:
    """
    Turn an SVG path into a polyline.

    Handles the three commands detangler emits: M, L, and quadratic Q for the
    ribbons, cubic C for the connectors. Curves are sampled rather than solved
    because every measurement below is a distance, and a dense polyline gives
    those to well under a pixel.
    """
    pts: List[Point] = []
    cur: Point = (0.0, 0.0)
    for cmd, args in re.findall(r"([MLQCmlqc])([^MLQCmlqc]*)", d):
        v = _nums(args)
        up = cmd.upper()
        if up == "M" and len(v) >= 2:
            cur = (v[0], v[1])
            pts.append(cur)
        elif up == "L":
            for i in range(0, len(v) - 1, 2):
                cur = (v[i], v[i + 1])
                pts.append(cur)
        elif up == "Q":
            for i in range(0, len(v) - 3, 4):
                p0, p1, p2 = cur, (v[i], v[i + 1]), (v[i + 2], v[i + 3])
                for k in range(1, per_curve + 1):
                    t = k / per_curve
                    m = 1 - t
                    pts.append((m * m * p0[0] + 2 * m * t * p1[0] + t * t * p2[0],
                                m * m * p0[1] + 2 * m * t * p1[1] + t * t * p2[1]))
                cur = p2
        elif up == "C":
            for i in range(0, len(v) - 5, 6):
                p0 = cur
                p1, p2, p3 = (v[i], v[i+1]), (v[i+2], v[i+3]), (v[i+4], v[i+5])
                for k in range(1, per_curve + 1):
                    t = k / per_curve
                    m = 1 - t
                    pts.append((
                        m**3 * p0[0] + 3*m*m*t * p1[0] + 3*m*t*t * p2[0] + t**3 * p3[0],
                        m**3 * p0[1] + 3*m*m*t * p1[1] + 3*m*t*t * p2[1] + t**3 * p3[1]))
                cur = p3
    return pts


def read_figure(path: str) -> Dict:
    """
    Pull the ribbons and connectors out of a detangler SVG.

    Works on both `_graph.svg` and the paired figure, because it keys off the
    layer ids rather than the document structure.
    """
    svg = open(path).read()

    def layer(name: str) -> str:
        m = re.search(rf'<g id="{name}".*?(?=<g id="|</svg>)', svg, re.S)
        return m.group(0) if m else ""

    seg_layer, link_layer = layer("layer-segments"), layer("layer-links")

    ribbons, widths = [], []
    for m in re.finditer(r'<path d="([^"]+)"([^/>]*)/>', seg_layer):
        d, attrs = m.group(1), m.group(2)
        w = re.search(r'stroke-width="([0-9.]+)"', attrs)
        pts = sample_path(d)
        if len(pts) >= 2:
            ribbons.append(pts)
            widths.append(float(w.group(1)) if w else 1.0)

    connectors = []
    for m in re.finditer(r'<path d="([^"]+)"', link_layer):
        pts = sample_path(m.group(1))
        if len(pts) >= 2:
            connectors.append(pts)

    W = re.search(r'width="([0-9.]+)"', svg)
    H = re.search(r'height="([0-9.]+)"', svg)
    return {
        "ribbons": ribbons,
        "connectors": connectors,
        # the ribbon width is the natural unit; they are all drawn the same
        "unit": max(widths) if widths else 1.0,
        "canvas": (float(W.group(1)) if W else 0.0, float(H.group(1)) if H else 0.0),
    }


# ---------------------------------------------------------------------------
# geometry
# ---------------------------------------------------------------------------

def _seg_dist(a: Point, b: Point, c: Point, d: Point) -> float:
    """Shortest distance between two line segments."""
    def pt_seg(p, q, r):
        qx, qy = r[0] - q[0], r[1] - q[1]
        L = qx * qx + qy * qy
        if L < 1e-12:
            return math.hypot(p[0] - q[0], p[1] - q[1])
        t = max(0.0, min(1.0, ((p[0]-q[0])*qx + (p[1]-q[1])*qy) / L))
        return math.hypot(p[0] - (q[0] + t*qx), p[1] - (q[1] + t*qy))
    if _crosses(a, b, c, d):
        return 0.0
    return min(pt_seg(a, c, d), pt_seg(b, c, d), pt_seg(c, a, b), pt_seg(d, a, b))


def _side(p: Point, q: Point, r: Point) -> float:
    return (q[0]-p[0])*(r[1]-p[1]) - (q[1]-p[1])*(r[0]-p[0])


def _crosses(a: Point, b: Point, c: Point, d: Point) -> bool:
    d1, d2 = _side(c, d, a), _side(c, d, b)
    d3, d4 = _side(a, b, c), _side(a, b, d)
    return ((d1 > 0) != (d2 > 0)) and ((d3 > 0) != (d4 > 0))


def _grid(polys: List[List[Point]], cell: float):
    """Bucket every polyline segment so neighbour queries stay linear-ish."""
    g: Dict[Tuple[int, int], List[Tuple[int, int]]] = {}
    for pi, pts in enumerate(polys):
        for si in range(len(pts) - 1):
            (x1, y1), (x2, y2) = pts[si], pts[si+1]
            for gx in range(int(min(x1, x2)//cell), int(max(x1, x2)//cell) + 1):
                for gy in range(int(min(y1, y2)//cell), int(max(y1, y2)//cell) + 1):
                    g.setdefault((gx, gy), []).append((pi, si))
    return g


def _pairs_near(polys, cell):
    """Candidate segment pairs from different polylines, via the grid."""
    g = _grid(polys, cell)
    seen = set()
    for cellkey, items in g.items():
        gx, gy = cellkey
        near = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                near.extend(g.get((gx+dx, gy+dy), ()))
        for (pi, si) in items:
            for (pj, sj) in near:
                if pi >= pj:
                    continue
                key = (pi, si, pj, sj)
                if key in seen:
                    continue
                seen.add(key)
                yield pi, si, pj, sj


def winding(pts: List[Point]) -> float:
    """
    Total absolute turning along a centreline, in degrees.

    This is the metric that agrees with the eye. A contig drawn as a gentle arc
    turns through ~40 degrees end to end; one that folds back on itself turns
    through 200 or more, and reads as a knot however much clear space is around
    it. Crucially it moves OPPOSITE to clearance under naive tuning: pushing
    ribbons apart makes them curl to fit, so optimising separation alone makes
    the picture worse. Both have to be in the objective.
    """
    tot = 0.0
    for k in range(1, len(pts) - 1):
        a = math.atan2(pts[k][1]-pts[k-1][1], pts[k][0]-pts[k-1][0])
        b = math.atan2(pts[k+1][1]-pts[k][1], pts[k+1][0]-pts[k][0])
        tot += abs((b - a + math.pi) % (2*math.pi) - math.pi)
    return math.degrees(tot)


def measure(fig: Dict) -> Dict:
    unit = fig["unit"] or 1.0
    ribbons, conns = fig["ribbons"], fig["connectors"]
    cw, ch = fig["canvas"]

    # --- clearance and ribbon x ribbon crossings ---
    best = float("inf")
    r_cross = 0
    cell = unit * 2.5
    for pi, si, pj, sj in _pairs_near(ribbons, cell):
        a, b = ribbons[pi][si], ribbons[pi][si+1]
        c, d = ribbons[pj][sj], ribbons[pj][sj+1]
        dist = _seg_dist(a, b, c, d)
        if dist < best:
            best = dist
        if dist == 0.0:
            r_cross += 1
    clearance = (best / unit) if best < float("inf") else float("nan")

    # --- ribbon x connector crossings ---
    rc = 0
    allp = ribbons + conns
    nr = len(ribbons)
    for pi, si, pj, sj in _pairs_near(allp, cell):
        if (pi < nr) == (pj < nr):
            continue                      # both ribbons, or both connectors
        a, b = allp[pi][si], allp[pi][si+1]
        c, d = allp[pj][sj], allp[pj][sj+1]
        if _crosses(a, b, c, d):
            rc += 1

    # --- hub angles ---
    # A hub is where connector ENDS gather. Cluster the endpoints, then measure
    # the angles at which the connectors leave. This is what decides whether a
    # reader can see what attaches to a busy contig.
    ends: List[Tuple[Point, float]] = []
    for pts in conns:
        if len(pts) < 3:
            continue
        ends.append((pts[0], math.atan2(pts[2][1]-pts[0][1], pts[2][0]-pts[0][0])))
        ends.append((pts[-1], math.atan2(pts[-3][1]-pts[-1][1], pts[-3][0]-pts[-1][0])))
    hubs: List[List[float]] = []
    used = [False] * len(ends)
    for i, (p, th) in enumerate(ends):
        if used[i]:
            continue
        group = [th]
        used[i] = True
        for j in range(i+1, len(ends)):
            if used[j]:
                continue
            q, th2 = ends[j]
            if math.hypot(p[0]-q[0], p[1]-q[1]) <= unit * 1.2:
                group.append(th2)
                used[j] = True
        if len(group) >= 2:
            hubs.append(group)
    hub_gap = 180.0
    for group in hubs:
        a = sorted(group)
        gaps = [a[k+1]-a[k] for k in range(len(a)-1)] + [2*math.pi - (a[-1]-a[0])]
        hub_gap = min(hub_gap, math.degrees(min(gaps)))
    if not hubs:
        hub_gap = float("nan")

    # --- ink and canvas ---
    ink = 0.0
    for pts in ribbons:
        ink += sum(math.hypot(pts[k+1][0]-pts[k][0], pts[k+1][1]-pts[k][1])
                   for k in range(len(pts)-1)) * unit
    winds = sorted((winding(p) for p in ribbons), reverse=True) or [0.0]
    # The area that matters is the box the GRAPH occupies, not the canvas. A
    # paired figure carries the chromosome panel too, so canvas area would
    # compare a graph-only file against a two-panel one and call the graph
    # tighter for reasons that have nothing to do with the layout.
    xs = [q[0] for pts in ribbons for q in pts] or [0.0]
    ys = [q[1] for pts in ribbons for q in pts] or [0.0]
    area = (max(xs) - min(xs) + unit) * (max(ys) - min(ys) + unit)
    return {
        "wind_mean": round(sum(winds) / len(winds), 1),
        "wind_max": round(winds[0], 1),
        "clearance": round(clearance, 3),
        "r_cross": r_cross,
        "rc_cross": rc,
        "hub_gap": round(hub_gap, 1),
        "fill": round(ink / area, 4) if area else 0.0,
        "area_u2": round(area / (unit * unit), 1),
        "n_ribbons": len(ribbons),
        "n_connectors": len(conns),
    }


# ---------------------------------------------------------------------------
# scoring
# ---------------------------------------------------------------------------

# What "clear enough" means, in ribbon widths and degrees. Both are targets, not
# maxima: exceeding them earns nothing, so the optimiser cannot win by simply
# inflating the canvas.
# Calibrated against the figure Anna judges clear (the one in the README):
# clearance 1.56 ribbon widths, mean winding 39 deg. Everything she called
# crowded sits at 78-104 deg of winding. So the winding target is set at the
# good figure's value and the clearance target just below it - the good figure
# should score close to full marks, and nothing should be able to beat it by
# trading one against the other.
CLEAR_TARGET = 1.5
WIND_TARGET = 45.0        # mean degrees of turning per ribbon; lower is better
WIND_MAX_TARGET = 170.0   # worst single ribbon


def score(m: Dict, area_ref: Optional[float] = None) -> float:
    if m["r_cross"] or m["rc_cross"]:
        return -1.0                                  # crossings are disqualifying
    c = m["clearance"]
    s = 0.0
    # separation, capped: exceeding the target earns nothing
    s += min(c / CLEAR_TARGET, 1.0) if c == c else 0.0
    # straightness, and it is worth double. Curling is what makes a figure
    # unreadable, and it is the thing naive tuning makes worse.
    s += 2.0 * min(WIND_TARGET / max(m["wind_mean"], 1e-6), 1.0)
    s += 0.5 * min(WIND_MAX_TARGET / max(m["wind_max"], 1e-6), 1.0)
    # and pay for canvas: doubling the area costs a whole point, so spreading
    # everything out is never a free win
    if area_ref:
        s -= max(0.0, m["area_u2"] / area_ref - 1.0)
    return round(s, 3)


def main(argv: Sequence[str]) -> int:
    as_json = "--json" in argv
    files = [a for a in argv if not a.startswith("--")]
    if not files:
        print(__doc__)
        return 2
    rows = []
    for f in files:
        try:
            m = measure(read_figure(f))
        except Exception as exc:                      # noqa: BLE001
            print(f"{f}: could not measure ({exc})", file=sys.stderr)
            continue
        m["file"] = f
        rows.append(m)
    if as_json:
        print(json.dumps(rows, indent=2))
        return 0
    ref = min((r["area_u2"] for r in rows), default=None)
    print(f"{'figure':<40} {'wind':>6} {'wmax':>6} {'clear':>6} {'hub°':>6} "
          f"{'xr':>3} {'xrc':>4} {'area':>8} {'score':>6}")
    for r in sorted(rows, key=lambda r: -score(r, ref)):
        name = r["file"]
        if len(name) > 39:
            name = "..." + name[-36:]
        print(f"{name:<40} {r['wind_mean']:>6.0f} {r['wind_max']:>6.0f} "
              f"{r['clearance']:>6.2f} {r['hub_gap']:>6.1f} "
              f"{r['r_cross']:>3} {r['rc_cross']:>4} "
              f"{r['area_u2']:>8.0f} {score(r, ref):>6.2f}")
    print(f"\ntargets: mean winding <= {WIND_TARGET}deg, worst <= {WIND_MAX_TARGET}deg, "
          f"clearance >= {CLEAR_TARGET} widths, crossings 0")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
