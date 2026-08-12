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
    RIBBON_W,
    MIN_DRAWN_PX,
    drawn_length_px,
    BAR_W,
    FS_ANNOT,
    FS_HEADING,
    FS_PRIMARY,
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
    Uniform. Thickness used to track read depth, but that put a second variable
    into the ribbon width and made the two panels hard to match up; depth is
    carried by the label instead. Wider than a chromosome bar, because a ribbon
    in the graph has to hold a number and be followed through a tangle, whereas
    a chromosome reads better slim.
    """
    return float(RIBBON_W)


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

    # K is graphviz's natural spring length and it is in INCHES. It was 1.4,
    # i.e. 100.8 pt - about 3.7x the bead spacing everything downstream assumes,
    # so `-Gsep=+18` (a node-separation budget meant to be read against the
    # ribbon width) was being applied at a scale we had not chosen. Deriving K
    # from `spacing` makes the whole call dimensionally coherent.
    k_inches = max(spacing, 1.0) / 72.0

    try:
        # sfdp, not neato. Bandage leans on OGDF's FMMM, which is a MULTILEVEL
        # force method - it coarsens the graph, lays the coarse version out, then
        # refines - and that is where its clean global arrangement comes from.
        # sfdp is the multilevel engine in graphviz; neato is single-level.
        # sfdp implements Yifan Hu, "Efficient and high quality force-directed
        # graph drawing", Mathematica Journal 10(1), 2005.
        # http://yifanhu.net/PUB/graph_draw_small.pdf
        #
        # Re-measured 11 Aug 2026 after finding that sfdp SILENTLY IGNORES `len`
        # (the graphviz docs mark it neato/fdp only, and a direct test confirms
        # it: neato honours a len=3.0 edge among len=0.2 edges, sfdp returns
        # them uniform). The original crossing measurement was therefore taken
        # under that bug, so it was redone. It still holds: on this graph sfdp
        # gives 0 ribbon crossings against neato's 2, with or without the K fix.
        # Chain length survives sfdp through BEAD COUNT rather than through
        # `len`, and the constraint-projection pass restores it exactly anyway.
        proc = subprocess.run(
            ["sfdp", "-Tplain", "-Goverlap=prism", "-Gsep=+18", f"-GK={k_inches:.4f}"],
            input=dot, capture_output=True, text=True, timeout=120,
        )
        if proc.returncode != 0 or not proc.stdout.strip():
            proc = subprocess.run(
                ["neato", "-Tplain", "-Goverlap=prism", "-Gsep=+6", "-Gmodel=circuit"],
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
    log.info(f"graph layout: initial placement from graphviz sfdp ({got} points)")
    return True


def _fan_angles(theta: List[float], gap: float) -> List[float]:
    """
    Spread angles apart to a minimum gap, moving them AS LITTLE AS POSSIBLE.

    Solves  min sum (theta'_i - theta_i)^2  subject to  theta'_{i+1} - theta'_i
    >= gap, with the cyclic order preserved. Substituting phi_i = theta_i - i*gap
    turns the gap constraint into plain monotonicity, so the exact optimum is the
    non-decreasing L2 isotonic regression of phi - computed here by
    pool-adjacent-violators in O(k). The result is recentred on the original
    circular mean so the fan does not drift.

    Minimising displacement is the point. An even fan would splay the contigs
    with mechanical symmetry, which is a different kind of artificial; this keeps
    whatever asymmetry the layout found and only opens the gaps that are too
    tight to read.

    PAVA: Ayer, Brunk, Ewing, Reid & Silverman (1955), Ann. Math. Statist.
    26:641-647. https://projecteuclid.org/euclid.aoms/1177728423

    Returns angles in the ORDER GIVEN, not sorted.
    """
    k = len(theta)
    if k < 2:
        return list(theta)
    gap = min(gap, 2.0 * math.pi / k)      # k gaps cannot exceed a full turn
    order = sorted(range(k), key=lambda i: theta[i])
    base = theta[order[0]]
    # unroll onto a line starting at the first angle
    lin = [(theta[i] - base) % (2.0 * math.pi) for i in order]
    phi = [lin[j] - j * gap for j in range(k)]
    # PAVA: pool adjacent violators of non-decreasing order, weights all 1
    vals: List[float] = []
    wts: List[int] = []
    for v in phi:
        vals.append(v)
        wts.append(1)
        while len(vals) > 1 and vals[-2] > vals[-1]:
            v2, w2 = vals.pop(), wts.pop()
            v1, w1 = vals.pop(), wts.pop()
            vals.append((v1 * w1 + v2 * w2) / (w1 + w2))
            wts.append(w1 + w2)
    flat: List[float] = []
    for v, wt in zip(vals, wts):
        flat.extend([v] * wt)
    out_lin = [flat[j] + j * gap for j in range(k)]
    shift = (sum(lin) - sum(out_lin)) / k     # recentre
    res = [0.0] * k
    for j, i in enumerate(order):
        res[i] = base + out_lin[j] + shift
    return res


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

    # Bead spacing is derived from the LONGEST contig, not fixed.
    #
    # A contig is a chain of beads, and how far it can curl is set by how many
    # beads it has. With a fixed spacing a 12 Mb contig gets a long floppy chain
    # and folds back on itself, while a 500 kb one gets three beads and cannot
    # bend at all - so one figure contains both a knot and a stick, and it gets
    # worse the wider the contig-length spread is. Capping the bead count of the
    # longest contig makes every contig's stiffness comparable on any genome.
    #
    # Measured on the Fusarium graph: mean winding (total turning per contig)
    # fell from 78 to 25 degrees, and it is winding, not separation, that tracks
    # whether a reader calls the figure clear. Pushing ribbons further apart
    # without this makes winding WORSE, because they curl to fit.
    MAX_BEADS = 8
    _longest = max((segment_draw_length(by_name[n].length, args) for n in names),
                   default=0.0)
    spacing = max(segment_thickness() * 0.42, 8.0, _longest / MAX_BEADS)
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
        iters_scale = 0.75
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
                if d > k * 14.0:
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
        if args.graph_triangle:
            # Keep-out for the lower-right triangle, applied as a FORCE inside
            # the solver rather than as a projection after it. That is the whole
            # difference: here it negotiates with repulsion and the springs in
            # the same iteration, so the layout finds a shape satisfying all
            # three. Bolted on afterwards it simply fought the length pass and
            # squashed contigs. The chromosome row is a staircase rising to the
            # right and so fills the complementary triangle; the two figures
            # then interlock instead of each claiming a rectangle.
            #
            # Applied in the FORCE phase ONLY. Interleaving it with the
            # constraint projection as well was what wrecked the layout: it kept
            # shoving beads together after the separation pass had pulled them
            # apart, and ribbons ended up 50 px apart on a 64 px stroke. Left to
            # the force phase alone, the 120 pure length-and-separation passes
            # that follow have the last word - measured 89 px apart, with the
            # ink past the diagonal still down from 0.071 to 0.045.
            xs_ = [q[0] for q in pts]
            ys_ = [q[1] for q in pts]
            X0_, Y0_ = min(xs_), min(ys_)
            W_ = (max(xs_) - X0_) or 1.0
            H_ = (max(ys_) - Y0_) or 1.0
            nx_, ny_ = 1.0 / W_, 1.0 / H_
            nl_ = math.hypot(nx_, ny_) or 1.0
            nx_, ny_ = nx_ / nl_, ny_ / nl_
            for i in range(n_pts):
                t_ = (pts[i][0] - X0_) / W_ + (pts[i][1] - Y0_) / H_ - 1.0
                if t_ > 0:
                    push = (t_ / nl_) * 0.6
                    disp[i][0] -= nx_ * push
                    disp[i][1] -= ny_ * push
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
    # Connector length is judged against the RIBBON WIDTH, not the bead spacing.
    # Bead spacing is a mechanical quantity - how finely a contig is
    # subdivided - and tying connectors to it means that making contigs stiffer
    # also lengthens every connector, which has nothing to do with it. It also
    # raises the threshold of the reel-in pass below (link_rest * 2.6), so a
    # contig the seed happened to drop far from its neighbour was never pulled
    # back. Measured on the Fusarium graph: longest connector 10.5 -> 7.2 ribbon
    # widths, and the figure 45% smaller in area.
    link_rest = segment_thickness() * 0.38
    # Centre-to-centre, in ribbon widths, so the VISIBLE gap is (w_sep - 1)
    # widths. At 1.12 that gap was 0.12 of a ribbon - technically not touching,
    # and reading as touching. On the Fusarium graph the old value actually
    # produced a clearance of 0.99, i.e. overlapping ribbons, which breaks the
    # invariant in detangler_graph-layout-research_v1.md.
    w_sep = segment_thickness() * 1.8
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
    # No artificial bow. A sinusoidal displacement added to every chain, with
    # its direction alternating by index, produced curves that were identical in
    # shape and unrelated to anything in the graph - decoration, and it read as
    # decoration. A contig should bend because its neighbours push it, which is
    # where the curve in a Bandage layout actually comes from.
    for _round in range(18):
        for n in names:
            idxs = chain[n]
            if len(idxs) < 3:
                continue
            new_xy = []
            for b in range(1, len(idxs) - 1):
                i, prv, nxt = idxs[b], idxs[b - 1], idxs[b + 1]
                mx_ = (pts[prv][0] + pts[nxt][0]) / 2.0
                my_ = (pts[prv][1] + pts[nxt][1]) / 2.0
                new_xy.append((i, pts[i][0] * 0.62 + mx_ * 0.38,
                                  pts[i][1] * 0.62 + my_ * 0.38))
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

    # ---- fan the hubs and clamp the chain end tangents ----
    # Two problems, one fix.
    #
    # (1) The smoothing above holds each chain's end POSITIONS but leaves its end
    # TANGENTS free. A chain with pinned endpoints and free tangents, smoothed to
    # convergence, is a CIRCULAR ARC - constant curvature, which is the most
    # machine-looking shape there is, and measurably what three of the four long
    # contigs had become. Pinning the tangents as well makes the minimum-bending
    # curve an Euler elastica instead, whose curvature VARIES along the stroke.
    # That variation is what reads as a confident drawn line, and it is derived
    # entirely from real geometry - the layout's endpoints and the graph's joins.
    # On minimum-bending-energy curves being piecewise elastica: Levien & Sequin,
    # "Interpolating splines: which is the fairest of them all?", CAD &
    # Applications 6 (2009), sec. 2.1.
    # https://people.eecs.berkeley.edu/~sequin/PAPERS/2009_CAD_Levien_Sequin.pdf
    #
    # (2) Where several contigs meet at one end, they leave it at whatever angles
    # the physics happened to give, often nearly on top of each other.
    #
    # So: compute a target outward direction per end, spread the ones that share
    # a hub, then rotate each chain's first few beads to match. A rotation about
    # the pinned end bead is an isometry, so drawn length is preserved exactly.
    def _end_idx(seg: str, end: str) -> int:
        return chain[seg][0] if end == "s" else chain[seg][-1]

    def _tangent_out(seg: str, end: str, span: float) -> Tuple[float, float]:
        idxs = chain[seg] if end == "e" else list(reversed(chain[seg]))
        tip = pts[idxs[-1]]
        ref = pts[idxs[0]]
        acc = 0.0
        for k2 in range(len(idxs) - 1, 0, -1):
            a_, b_ = pts[idxs[k2]], pts[idxs[k2 - 1]]
            acc += math.hypot(a_[0] - b_[0], a_[1] - b_[1])
            if acc >= span:
                ref = b_
                break
        dx, dy = tip[0] - ref[0], tip[1] - ref[1]
        d = math.hypot(dx, dy) or 1.0
        return dx / d, dy / d

    partners: Dict[Tuple[str, str], List[Tuple[str, str]]] = {}
    for l in links:
        if l.a not in chain or l.b not in chain or l.a == l.b:
            continue
        ae = "e" if l.a_orient == "+" else "s"
        be = "s" if l.b_orient == "+" else "e"
        partners.setdefault((l.a, ae), []).append((l.b, be))
        partners.setdefault((l.b, be), []).append((l.a, ae))

    w_rib = segment_thickness()
    # Two ribbons of width w are visibly clear of one another at radial distance
    # d from a hub once their angular gap exceeds ~1.15*w/d. Taking d = 2w - the
    # distance at which a reader wants them already separate - gives ~33 degrees.
    min_gap = 1.15 / 2.0
    target: Dict[Tuple[str, str], Tuple[float, float]] = {}
    for hub, ps in partners.items():
        if len(ps) < 2:
            continue
        hx, hy = pts[_end_idx(*hub)]
        incident = [hub] + [q for q in ps if q in partners]
        seen_i: List[Tuple[str, str]] = []
        for q in incident:
            if q not in seen_i:
                seen_i.append(q)
        if len(seen_i) < 2:
            continue
        # angle of the direction each contig DEPARTS the hub (= -outward)
        th = []
        for seg, end in seen_i:
            ux, uy = _tangent_out(seg, end, w_rib * 0.75)
            th.append(math.atan2(-uy, -ux))
        for (seg, end), a in zip(seen_i, _fan_angles(th, min_gap)):
            target[(seg, end)] = (-math.cos(a), -math.sin(a))

    if target:
        for _pass in range(6):
            for (seg, end), (tx, ty) in target.items():
                idxs = chain[seg] if end == "s" else list(reversed(chain[seg]))
                if len(idxs) < 3:
                    continue
                m = max(3, min(len(idxs) - 1,
                               int(round(w_rib * 2.5 / max(chain_rest[seg], 1e-6)))))
                if m >= len(idxs):
                    continue
                # How far the FIRST segment is from where it should point.
                ax_, ay_ = pts[idxs[1]][0] - pts[idxs[0]][0], pts[idxs[1]][1] - pts[idxs[0]][1]
                al = math.hypot(ax_, ay_) or 1e-6
                cur = math.atan2(ay_ / al, ax_ / al)
                dth = math.atan2(-ty, -tx) - cur
                dth = (dth + math.pi) % (2.0 * math.pi) - math.pi
                dth *= 0.5   # damped, so the chain eases round instead of snapping
                # Spread the correction over m joints with a smoothly DECAYING
                # profile, applied as nested suffix rotations: rotate beads
                # k..end about bead k-1 by the step in the profile, so segment
                # k's direction ends up turned by phi_k and segment m by nothing.
                #
                # What this replaces: rotating the first m beads RIGIDLY about
                # the pinned end. That put the whole correction into a single
                # joint and left everything before it perfectly straight, which
                # is why a hairpin read as two straight legs and a bend rather
                # than as one curve. Every step here is still a rigid rotation
                # of a suffix, so segment lengths are preserved exactly.
                prev_phi = 0.0
                for k in range(m):
                    phi = dth * (1.0 - k / float(m)) ** 2
                    step = phi - prev_phi
                    prev_phi = phi
                    if abs(step) < 1e-12:
                        continue
                    px, py = pts[idxs[k]]
                    cs, sn = math.cos(step), math.sin(step)
                    for b in idxs[k + 1:]:
                        rx, ry = pts[b][0] - px, pts[b][1] - py
                        pts[b][0] = px + rx * cs - ry * sn
                        pts[b][1] = py + rx * sn + ry * cs
            for n in names:
                idxs = chain[n]
                rest = chain_rest[n]
                for _sub in range(4):
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

    # Rotate to fill a WIDE band rather than to align a principal axis.
    #
    # This used to take the PCA principal axis and lay it horizontal. That gets
    # the long direction right but says nothing about the shape of what is left
    # over: a sprawling component can be principal-axis-aligned and still leave a
    # whole corner of the panel empty, which is what happened once the fanning
    # spread the components out. Searching rotations and scoring the resulting
    # BOUNDING BOX directly is the fix, and it is what Bandage does with
    # stepsForRotatingComponents (program/graphlayoutworker.cpp).
    #
    # The score wants a box near ASPECT_TARGET wide-to-tall, because the figure
    # stacks the graph above the chromosomes and so height is the scarce
    # dimension, with a light preference for smaller area to break ties. Applied
    # per component first - each one is packed separately, so each gets to
    # choose its own angle - and then once more to the assembled panel.
    ASPECT_TARGET = 1.9

    def _rot(v, ang, ox, oy):
        ca, sa = math.cos(ang), math.sin(ang)
        return [((x - ox) * ca - (y - oy) * sa, (x - ox) * sa + (y - oy) * ca)
                for x, y in v]

    def _best_angle(pointsets) -> float:
        pts = [p for v in pointsets for p in v]
        if len(pts) < 3:
            return 0.0
        ox = sum(p[0] for p in pts) / len(pts)
        oy = sum(p[1] for p in pts) / len(pts)
        best, best_sc = 0.0, None
        for step in range(90):                     # 2 degree steps over a half turn
            ang = math.pi * step / 90.0
            r = _rot(pts, ang, ox, oy)
            xs_ = [q[0] for q in r]; ys_ = [q[1] for q in r]
            bw = max(xs_) - min(xs_); bh = max(ys_) - min(ys_)
            if bw <= 1e-6 or bh <= 1e-6:
                continue
            # log-distance from the target ratio, plus a light area term...
            sc = abs(math.log(bw / bh) - math.log(ASPECT_TARGET)) + 0.15 * math.log(bw * bh)
            # ...and a push towards the UPPER-LEFT TRIANGLE. The chromosome
            # panel is a row of bars on a shared baseline sorted short to tall,
            # so its silhouette is a staircase rising to the right - a lower-
            # right triangle. Cutting the square along that diagonal gives the
            # graph the complementary upper-left triangle, and the two shapes
            # interlock instead of each claiming a full rectangle. Scored as the
            # mean depth of the ink past the diagonal, so a layout only pays for
            # how far into the chromosomes' half it actually reaches.
            x0_, y0_ = min(xs_), min(ys_)
            over = 0.0
            for qx, qy in r:
                t = (qx - x0_) / bw + (qy - y0_) / bh - 1.0
                if t > 0:
                    over += t
            sc += 8.0 * (over / len(r))
            if best_sc is None or sc < best_sc:
                best, best_sc = ang, sc
        return best

    for g in comps.values():
        ang = _best_angle([poly[n] for n in g])
        if abs(ang) < 1e-9:
            continue
        pts_g = [p for n in g for p in poly[n]]
        ox = sum(p[0] for p in pts_g) / len(pts_g)
        oy = sum(p[1] for p in pts_g) / len(pts_g)
        for n in g:
            poly[n] = _rot(poly[n], ang, ox, oy)

    groups = sorted(comps.values(), key=lambda g: -sum(by_name[n].length for n in g))
    gap = segment_thickness() * 1.5

    # Isolated components go into the HOLES the main component leaves, biased
    # towards the top left.
    #
    # They used to be stacked in a column beside the main component, which put
    # them in the one corner the figure cannot spare: the bottom right, where
    # the chromosome row goes. It also grew the graph's bounding box sideways
    # for no reason, since a spring layout always leaves slack somewhere inside
    # its own box. An organelle ring or an unplaced fragment is small and has no
    # links, so it can sit anywhere that is empty - there is no reason for it to
    # claim new canvas.
    #
    # Greedy placement on a coarse grid: mark where the main component has ink,
    # then for each remaining component take the free position that scores best.
    # The score prefers positions towards the TOP LEFT - away from the diagonal
    # the chromosomes fill - and penalises any position that would enlarge the
    # overall bounding box, so filling a hole always beats growing the figure.
    if len(groups) > 1:
        cell = max(segment_thickness() * 0.9, 12.0)
        clear = segment_thickness() * 1.15

        def _pts(g):
            return [q for n in g for q in poly[n]]

        occupied: Set[Tuple[int, int]] = set()

        def _mark(points):
            r = int(math.ceil(clear / cell))
            for x, y in points:
                gx, gy = int(x // cell), int(y // cell)
                for a in range(-r, r + 1):
                    for b in range(-r, r + 1):
                        occupied.add((gx + a, gy + b))

        def _free(points, dx, dy):
            return all(
                (int((x + dx) // cell), int((y + dy) // cell)) not in occupied
                for x, y in points
            )

        _mark(_pts(groups[0]))
        mx0, my0, mx1, my1 = bbox(groups[0])
        anchor: Optional[Tuple[float, float]] = None
        for g in groups[1:]:
            gx0, gy0, gx1, gy1 = bbox(g)
            gw_, gh_ = gx1 - gx0, gy1 - gy0
            pts_g = _pts(g)
            span_x = max(mx1 - mx0, 1.0)
            span_y = max(my1 - my0, 1.0)
            best = None
            y_ = my0 - gap
            while y_ <= my1 + gap:
                x_ = mx0 - gap
                while x_ <= mx1 + gap:
                    dx, dy = x_ - gx0, y_ - gy0
                    if _free(pts_g, dx, dy):
                        cx = (x_ + gw_ / 2.0 - mx0) / span_x
                        cy = (y_ + gh_ / 2.0 - my0) / span_y
                        # Towards the BOTTOM LEFT. Not merely "away from the
                        # chromosomes" - the top left is further from them - but
                        # out of the reading path. The eye goes across the graph
                        # and then down to the chromosome row, and unlinked
                        # fragments dropped along the top of the graph sit right
                        # in the middle of that, reading as part of the structure
                        # when they are precisely the things that have none. The
                        # low left corner is the one place a reader passes
                        # through last.
                        sc = cx + (1.0 - cy)
                        # Keep them TOGETHER once one has been placed: a cluster
                        # of unlinked pieces reads as a category, the same pieces
                        # scattered around the margin read as four separate
                        # accidents.
                        if anchor is not None:
                            ax_, ay_ = anchor
                            sc += 0.6 * math.hypot(
                                (x_ + gw_ / 2.0 - ax_) / span_x,
                                (y_ + gh_ / 2.0 - ay_) / span_y,
                            )
                        # ...but never at the cost of a bigger figure
                        grow = (
                            max(0.0, mx0 - x_) + max(0.0, (x_ + gw_) - mx1)
                            + max(0.0, my0 - y_) + max(0.0, (y_ + gh_) - my1)
                        )
                        sc += 4.0 * grow / max(span_x, span_y)
                        if best is None or sc < best[0]:
                            best = (sc, dx, dy)
                    x_ += cell
                y_ += cell
            if best is None:
                # nowhere free inside the box: fall back to beside it
                best = (0.0, mx1 + gap - gx0, my0 - gy0)
                mx1 += gw_ + gap
            _, dx, dy = best
            for n in g:
                poly[n] = [(x + dx, y + dy) for x, y in poly[n]]
            anchor = ((gx0 + gx1) / 2.0 + dx, (gy0 + gy1) / 2.0 + dy)
            _mark(_pts(g))
            mx0 = min(mx0, gx0 + dx)
            my0 = min(my0, gy0 + dy)
            mx1 = max(mx1, gx1 + dx)
            my1 = max(my1, gy1 + dy)


    ang = _best_angle([poly[n] for n in names])
    if abs(ang) > 1e-9:
        allp2 = [p for n in names for p in poly[n]]
        ox = sum(p[0] for p in allp2) / len(allp2)
        oy = sum(p[1] for p in allp2) / len(allp2)
        poly = {n: _rot(v, ang, ox, oy) for n, v in poly.items()}

    xs = [p[0] for n in names for p in poly[n]]
    ys = [p[1] for n in names for p in poly[n]]
    # Everything drawn extends beyond the bead it hangs off: a ribbon by half its
    # width, a ring by its radius plus half its stroke, a coverage label by its
    # offset plus the text itself. Padding from bead positions alone shaved
    # whichever contig happened to sit at the edge of the layout.
    thick_ = segment_thickness()
    ring_reach = max(
        [thick_] + [
            max(segment_draw_length(by_name[n].length, args) / (2 * math.pi),
                thick_ * 0.85) + thick_ * 0.3
            for n in names
        ]
    )
    label_reach = thick_ / 2 + 18 + 12 + 18 * 3.0
    pad = 18.0 + max(thick_ * 0.6, ring_reach, label_reach)
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

    def outward(seg: str, end: str, back: float) -> Tuple[float, float]:
        """
        Unit vector pointing OUT of `seg` at `end`, measured over `back` pixels
        of arclength rather than over one bead. Beads sit about 0.4 of a ribbon
        width apart, so a single-bead tangent is dominated by layout noise;
        averaging over most of a width gives the direction the ribbon actually
        appears to be travelling as it arrives at the join.
        """
        pts = geom.get(seg) or []
        if len(pts) < 2:
            return (1.0, 0.0)
        chain = pts if end == "e" else list(reversed(pts))
        tip = chain[-1]
        ref = chain[0]
        acc = 0.0
        for k in range(len(chain) - 1, 0, -1):
            acc += math.hypot(chain[k][0] - chain[k - 1][0], chain[k][1] - chain[k - 1][1])
            if acc >= back:
                ref = chain[k - 1]
                break
        dx, dy = tip[0] - ref[0], tip[1] - ref[1]
        d = math.hypot(dx, dy) or 1.0
        return (dx / d, dy / d)

    def terminal(seg: str, end: str) -> Tuple[float, float]:
        # a tip is joined AT its live end, on the rim of the dot; a short contig
        # that really is a bridge keeps its centre so both sides meet in it
        if seg in dot_centre and len(ends_used.get(seg, set())) == 1:
            return raw_end(seg, end)
        if seg in dot_centre:
            return dot_centre[seg]
        return geom[seg][0] if end == "s" else geom[seg][-1]

    # ---- one attachment point per link ----
    # "Hub" is the wrong word for what this is, and using it invites the wrong
    # reading. Nothing routes THROUGH here. What the graph records is that five
    # contig ends abut ONE END of contig 9 - its e end - while its other end is
    # free and carries the telomere array. A path that entered contig 9 would
    # have to leave by the end it arrived at, which is not a path, so none of
    # those five can reach any of the others.
    #
    # Drawn as five lines converging on a single point, that is unreadable: it
    # looks like a junction, and a junction is exactly what a reader would then
    # route contig 3 through to contig 7. So each link gets its OWN attachment
    # point, spread across the width of contig 9's end face. Every line then
    # visibly terminates ON contig 9 rather than at a shared node - which is
    # precisely the claim the graph supports.
    #
    # Attachment points are assigned in the order the lines arrive, by angle, so
    # spreading them cannot introduce a crossing.
    incoming: Dict[Tuple[str, str], List[Tuple[str, str]]] = defaultdict(list)
    for l in links:
        if l.a not in geom or l.b not in geom or l.a == l.b:
            continue
        ae = "e" if l.a_orient == "+" else "s"
        be = "s" if l.b_orient == "+" else "e"
        incoming[(l.a, ae)].append((l.b, be))
        incoming[(l.b, be)].append((l.a, ae))

    port: Dict[Tuple[str, str, str, str], Tuple[float, float]] = {}
    for (seg, end), others in incoming.items():
        if len(others) < 2:
            continue
        hx, hy = terminal(seg, end)
        ux, uy = outward(seg, end, w * 0.75)
        # across the end face, not along the contig
        pxn, pyn = -uy, ux
        # Ordered by the ANGLE each line arrives at, measured in the end face's
        # own frame, not by how far along the face it happens to project. Two
        # lines coming from the same side at different distances project to the
        # same place and then cross on the way in - which is what put contig 1's
        # line over contig 3's. Angular order around the face is the standard
        # fix and is crossing-free for lines that converge on one point.
        def _arrival(o):
            tx_, ty_ = terminal(*o)
            vx_, vy_ = tx_ - hx, ty_ - hy
            return math.atan2(vx_ * pxn + vy_ * pyn, vx_ * ux + vy_ * uy)

        ranked = sorted(others, key=_arrival)
        k = len(ranked)
        # Kept inside the ribbon: half a width is its edge, so 0.34 either side
        # leaves the outermost line clearly on the contig rather than clipping
        # its corner.
        span = min(w * 0.68, w * 0.30 * (k - 1))
        step_ = span / (k - 1) if k > 1 else 0.0
        for i, o in enumerate(ranked):
            off = (i - (k - 1) / 2.0) * step_
            # Tucked BACK along the contig's own axis rather than left sitting
            # on the spine's end point. The ribbon is stroked with a round cap
            # and the polyline is trimmed half a width before stroking, so a
            # point offset sideways from the spine end falls outside the cap -
            # which is why the lines stopped short with a sliver of white
            # between them and contig 9. Connectors are drawn beneath the
            # contigs, so running them under the cap is invisible and closes
            # the gap at every offset.
            port[(seg, end, o[0], o[1])] = (
                hx + pxn * off - ux * w * 0.55,
                hy + pyn * off - uy * w * 0.55,
            )

    def attach(seg: str, end: str, other: Tuple[str, str]) -> Tuple[float, float]:
        return port.get((seg, end, other[0], other[1])) or terminal(seg, end)

    # junction connectors, behind the contigs
    P.append(
        f'<g id="layer-links" fill="none" stroke="{PALETTE["bar_edge"]}" '
        f'stroke-linecap="round">'
    )
    for l in links:
        if l.a not in geom or l.b not in geom or l.a == l.b:
            continue
        a_end = "e" if l.a_orient == "+" else "s"
        b_end = "s" if l.b_orient == "+" else "e"
        ax, ay = attach(l.a, a_end, (l.b, b_end))
        bx, by = attach(l.b, b_end, (l.a, a_end))
        # A cubic whose control points lie ONE RIBBON WIDTH beyond each contig
        # end, along that contig's own final tangent. The connector then leaves
        # each ribbon in the direction the ribbon was already going, so contig -
        # connector - contig is a single continuous stroke with no kink at the
        # join. This is how Bandage draws its edges, and it is the whole reason
        # its junctions look drawn rather than wired together. See Bandage's
        # GraphicsItemEdge::calculateAndSetPath():
        # https://github.com/rrwick/Bandage/blob/main/graph/graphicsitemedge.cpp
        #
        # What this replaces: a quadratic bowed away from the drawing's
        # CENTROID, whose direction depended on where the connector sat in the
        # picture rather than on which way the contigs pointed. Every sharp
        # angle at a join came from that, and it was decoration keyed to a
        # global coordinate - the same objection as the sinusoidal pre-bow.
        ux, uy = outward(l.a, a_end, w * 0.75)
        vx, vy = outward(l.b, b_end, w * 0.75)
        chord = math.hypot(bx - ax, by - ay)
        ext = min(w, chord / 2.0) if chord > 1e-6 else w
        c1x, c1y = ax + ux * ext, ay + uy * ext
        c2x, c2y = bx + vx * ext, by + vy * ext
        P.append(
            f'<path d="M {ax:.1f} {ay:.1f} C {c1x:.1f} {c1y:.1f} '
            f'{c2x:.1f} {c2y:.1f} {bx:.1f} {by:.1f}" '
            f'stroke-width="{max(w * 0.085, 3.0):.1f}" stroke-opacity="0.85"/>'
        )
    # Where several connectors land on the SAME end of a contig, mark the point.
    # Five lines arriving near a short contig read as five lines arriving
    # somewhere near it, and which end they share is exactly what the reader
    # needs - edge_9 has all five of its neighbours on one end. Drawn after the
    # lines so it caps them.
    ends_here: Dict[Tuple[float, float], int] = defaultdict(int)
    for l in links:
        if l.a not in geom or l.b not in geom or l.a == l.b:
            continue
        for seg_, end_ in ((l.a, "e" if l.a_orient == "+" else "s"),
                           (l.b, "s" if l.b_orient == "+" else "e")):
            px_, py_ = terminal(seg_, end_)
            ends_here[(round(px_, 1), round(py_, 1))] += 1
    # No junction dot. It was there to say "these lines meet here", but that is
    # the very reading the ports exist to remove - the lines do not meet, they
    # each end on the same contig. A dot re-drew them as one node.
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
    # Everything already on the page that a depth label must not land on: every
    # bead of every ribbon, plus the junction points where the connectors meet.
    # Sampled rather than exact - a label only needs to be clearly off the ink,
    # not provably disjoint from it.
    obstacles: List[Tuple[float, float]] = [
        p for v in geom.values() for p in v
    ]
    for v in geom.values():
        if len(v) >= 2:
            obstacles.append(v[0])
            obstacles.append(v[-1])
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
            f'<text x="{mx:.1f}" y="{my + FS_PRIMARY * 0.35:.1f}" text-anchor="middle" '
            f'font-size="{FS_PRIMARY}" font-weight="700" fill="{_text_on(colour)}">'
            f'{esc(_segment_number(name))}</text>'
        )
        # Coverage is shown for EVERY contig. It used to be suppressed on short
        # ones for want of room, which quietly hid it on exactly the segments
        # where depth decides the call - a repeat, a tip, an organelle.
        if c.depth is not None:
            # Coverage labels crowd badly around a hub, where several short
            # contigs meet. Each label tries a ring of candidate positions and
            # takes the one furthest from the labels already placed.
            off = w / 2 + FS_ANNOT + 12
            best_pt, best_score = None, -1.0
            for k in range(24):
                ang = 2 * math.pi * k / 24
                ox_, oy_ = math.cos(ang), math.sin(ang)
                # bias towards the ribbon's normal, so a label still reads as
                # belonging to its own contig
                bias = 1.0 + 0.6 * (ox_ * nx + oy_ * ny)
                cxp, cyp = mx + ox_ * off, my + oy_ * off
                near = min(
                    (math.hypot(cxp - px, cyp - py) for px, py in placed_labels),
                    default=1e6,
                )
                # ...and away from the DRAWING, not only from other labels. The
                # old rule dodged labels but not ribbons, so at a hub - where the
                # contigs converge and there is nothing else to collide with -
                # the depth label landed squarely on the junction. A hub is
                # exactly where depth matters most, so it is the worst place to
                # lose it.
                clear = min(
                    (math.hypot(cxp - px, cyp - py) for px, py in obstacles),
                    default=1e6,
                )
                score = min(min(near, clear * 1.4), 400.0) * bias
                if score > best_score:
                    best_pt, best_score = (cxp, cyp), score
            lx, ly = best_pt  # type: ignore
            placed_labels.append((lx, ly))
            # Muted, not full-strength. Depth is evidence about a contig, not the
            # contig's name, and setting it in the same ink as the identity
            # labels let it compete with them for the reader's first pass.
            P.append(
                f'<text x="{lx:.1f}" y="{ly:.1f}" text-anchor="middle" '
                f'font-size="{FS_PRIMARY}" font-weight="600" '
                f'fill="{PALETTE["muted"]}">{c.depth:.0f}x</text>'
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


def _essential_names(calls, hypothesis_joins):
    """Segments the drawing must not drop: the backbone, and the routes between it."""
    names = [c.name for c in calls if c.cls == "backbone"]
    for j in (hypothesis_joins or []):
        names.append(j.a)
        names.append(j.b)
        names.extend(j.via)
    return names


def prune_for_drawing(calls, links, args, log: Log, hypothesis_joins=None):
    """
    Keep the graph drawable by DROPPING segments, not by refusing to draw.

    A human or plant assembly graph has thousands of segments; the old rule
    printed "too many, figure skipped", which is the least useful thing it could
    do - the whole reason to look at a picture of a 5,000-segment graph is to
    find the handful of pieces that matter.

    What is kept is the longest `--max-graph-nodes` segments, because length is
    what chromosome structure is made of, plus any segment linked to one of
    those so a kept contig never appears to float free. What is dropped is
    reported by count and by span, so the reader knows the picture is partial
    and by how much. Nothing here touches the INFERENCE - classification and
    hypotheses have already run over the whole graph. This is a drawing filter.
    """
    limit = int(getattr(args, "max_graph_nodes", 300) or 300)
    if len(calls) <= limit:
        return calls, links, None
    ranked = sorted(calls, key=lambda c: -c.length)

    # Keep THE GRAPH THE INFERENCE USED first: the backbone segments and the
    # segments the drawn hypothesis's joins run through. On a real assembly,
    # "the 300 longest" is 300 unconnected dots that have nothing to do with
    # the chromosome panel beside them. Length fills whatever room is left.
    # This only bites once the graph is over the limit; below it, every segment
    # is drawn exactly as before.
    keep = set()
    for name in _essential_names(calls, hypothesis_joins):
        if len(keep) >= limit:
            break
        keep.add(name)
    for c in ranked:
        if len(keep) >= limit:
            break
        keep.add(c.name)
    # ...and anything directly attached to a kept segment, so no kept contig is
    # drawn with a link going nowhere.
    for l in links:
        if l.a in keep or l.b in keep:
            keep.add(l.a)
            keep.add(l.b)
    kept = [c for c in calls if c.name in keep]
    if len(kept) > limit * 2:      # attachment closure ran away; fall back
        kept = ranked[:limit]
        keep = {c.name for c in kept}
    dropped = [c for c in calls if c.name not in keep]
    span = sum(c.length for c in dropped)
    total = sum(c.length for c in calls) or 1
    note = (
        f"{len(dropped):,} of {len(calls):,} segments omitted from the figure "
        f"({human_bp(span)}, {100.0 * span / total:.1f}% of the assembly); "
        f"the {len(kept):,} drawn are the chromosome-sized pieces, the segments the "
        f"candidate joins run through, and then the longest of the rest"
    )
    log.warn(
        note + f". Raise --max-graph-nodes to draw more. Classification and the "
        f"chromosome hypotheses used ALL {len(calls):,} segments; only the picture is partial."
    )
    kept_names = {c.name for c in kept}
    kept_links = [l for l in links if l.a in kept_names and l.b in kept_names]
    return kept, kept_links, note


def render_graph_figure(
    calls: List[SegmentCall],
    links: List[GfaLink],
    colours: Dict[str, str],
    title: str,
    path: str,
    args,
    log: Log,
    hypothesis_joins=None,
) -> Optional[str]:
    """The Bandage graph redrawn with our own, fixed colours."""
    calls, links, _note = prune_for_drawing(calls, links, args, log, hypothesis_joins)
    with open(path, "w") as fh:
        fh.write(graph_svg_for_style(calls, links, title, colours, args, log))
    return path


def graph_svg_for_style(calls, links, title, colours, args, log) -> str:
    if args.graph_style == "bandage":
        return render_bandage_style_svg(calls, links, title, colours, args, log)
    return render_graph_svg(calls, links, title, colours)
