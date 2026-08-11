"""Enumerating and ranking chromosome hypotheses."""
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
from statistics import median
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
    human_bp,
)
from .graph import (
    OTHER_END,
)
from .calls import (
    SegmentCall,
    telomeric_segments,
)



# --------------------------------------------------------------------------
# chromosome hypotheses
# --------------------------------------------------------------------------
@dataclass
class Join:
    a: str
    b: str
    via: List[str]  # intermediate segment names, in order
    # which physical end of a and of b the route leaves from / arrives at. A
    # chain may consume each end at most once, which is what stops two different
    # joins from both hanging off the same side of a contig.
    a_end: str = "e"
    b_end: str = "s"
    # True when the route is NOT supported by a traversable path: the two
    # segments merely end in the same one-sided repeat. Kept as a declared
    # alternative rather than dropped, because it is often the biologically
    # right answer - it is just not something this graph establishes.
    speculative: bool = False

    @property
    def key(self) -> Tuple[str, str]:
        return tuple(sorted((self.a, self.b)))  # type: ignore

    @property
    def ends(self) -> Tuple[Tuple[str, str], Tuple[str, str]]:
        return ((self.a, self.a_end), (self.b, self.b_end))

    def describe(self) -> str:
        route = " - ".join([self.a] + self.via + [self.b])
        return route if self.via else f"{self.a} - {self.b} (direct link)"


def chain_end_status(
    chain: List[str],
    internal: Set[str],
    adj: Dict[str, Set[str]],
    telomeric: Dict[str, int],
    end_adj: Optional[Dict[Tuple[str, str], Set[str]]] = None,
) -> Tuple[int, int, List[str]]:
    """
    How many of a chain's two ends look finished.

    An end counts as capped when the terminal backbone segment has a telomere-
    bearing neighbour that is not already used inside this chain: a telomere in
    the middle of a chain is not capping anything. When end_adj is given, the
    two physical ends of each terminal segment are assessed independently, so
    a segment that abuts the same telomeric segment at BOTH its ends is
    credited with two capped ends, not one. Returns (capped, open, notes).
    """
    members = set(chain)
    notes: List[str] = []

    def telo_at(seg: str, side: str) -> List[str]:
        ext = {
            n
            for n in (end_adj or {}).get((seg, side), ())
            if n not in internal and n not in members
        }
        return sorted(n for n in ext if n in telomeric)

    if len(chain) == 1:
        seg = chain[0]
        if end_adj is not None:
            # Each physical end of the single segment is its own molecule end,
            # so the ends are assessed independently: one telomeric neighbour
            # sitting at both ends caps both of them.
            per_side = [telo_at(seg, side) for side in ("s", "e")]
            capped = sum(1 for t in per_side if t)
            if capped == 2 and per_side[0] == per_side[1]:
                notes.append(
                    f"{seg} abuts the telomeric segment {', '.join(per_side[0])} at both ends"
                )
            else:
                for t in per_side:
                    if t:
                        notes.append(f"{seg} abuts the telomeric segment {', '.join(t)}")
            return capped, 2 - capped, notes
        # No orientation information: bounded by distinct telomeric neighbours.
        ext = {n for n in adj.get(seg, ()) if n not in internal and n not in members}
        telo = sorted(n for n in ext if n in telomeric)
        capped = min(len(telo), 2)
        if telo:
            notes.append(f"{seg} abuts the telomeric segment {', '.join(telo)}")
        return capped, 2 - capped, notes

    capped = 0
    for end, inner in ((chain[0], chain[1]), (chain[-1], chain[-2])):
        if end_adj is not None:
            # The end of the terminal segment that faces into the chain
            # (towards the next backbone member, possibly through join
            # segments) cannot cap this molecule end; only the outward end can.
            inward = {
                side
                for side in ("s", "e")
                if any(n == inner or n in internal for n in end_adj.get((end, side), ()))
            }
            if len(inward) == 1:
                outer = "e" if inward == {"s"} else "s"
                telo = telo_at(end, outer)
                if telo:
                    capped += 1
                    notes.append(f"{end} abuts the telomeric segment {telo[0]}")
                continue
            # Both or neither end faces inward, so the orientation cannot be
            # resolved; fall through to the name-based check.
        ext = {n for n in adj.get(end, ()) if n not in internal and n not in members}
        telo = sorted(n for n in ext if n in telomeric)
        if telo:
            capped += 1
            notes.append(f"{end} abuts the telomeric segment {telo[0]}")
    return capped, 2 - capped, notes


@dataclass
class Hypothesis:
    rank: int
    chains: List[List[str]]
    joins: List[Join]
    score: float
    supporting: List[str]
    contradicting: List[str]
    resolve_with: List[str]
    capped_ends: int = 0
    open_ends: int = 0

    def chain_length(self, chain: List[str], lengths: Dict[str, int]) -> int:
        return sum(lengths.get(s, 0) for s in chain)


def find_joins(
    calls: List[SegmentCall],
    end_links: Dict[Tuple[str, str], Set[Tuple[str, str]]],
    max_hops: int,
) -> List[Join]:
    """
    Every route between two backbone segments that passes only through
    non-backbone segments, using at most max_hops intermediates. Low-depth
    segments are deliberately included: they look ignorable but they change the
    topology.

    The traversal is END-AWARE, and that is the whole point. A GFA link joins a
    specific end of one segment to a specific end of another, so a route that
    passes THROUGH an intermediate must arrive at one of its ends and leave by
    the opposite one. Walking a segment-level adjacency instead - which is what
    this function used to do - invents routes that enter and leave through the
    same end, which no assembly graph permits, and it lets a segment with links
    on one end only masquerade as a bridge when it is really a tip.
    """
    backbone = {c.name for c in calls if c.cls == "backbone"}
    joins: Dict[Tuple[str, str, str, str, Tuple[str, ...]], Join] = {}
    for start in backbone:
        for start_end in ("s", "e"):
            # each stack entry: the end we are about to leave from, and the
            # intermediates consumed so far
            stack: List[Tuple[Tuple[str, str], List[str]]] = [((start, start_end), [])]
            while stack:
                (node, exit_end), via = stack.pop()
                for nb, nb_end in sorted(end_links.get((node, exit_end), ())):
                    if nb == start or nb in via:
                        continue
                    if nb in backbone:
                        j = Join(
                            a=start, b=nb, via=list(via),
                            a_end=start_end, b_end=nb_end,
                        )
                        key = (start, start_end, nb, nb_end, tuple(via))
                        # the same physical route found from the other direction
                        rev = (nb, nb_end, start, start_end, tuple(reversed(via)))
                        if rev not in joins:
                            joins.setdefault(key, j)
                    elif len(via) < max_hops:
                        # enter nb at nb_end, so we may only leave by its far end
                        stack.append(((nb, OTHER_END[nb_end]), via + [nb]))
    return list(joins.values())


def enumerate_hypotheses(
    calls: List[SegmentCall],
    joins: List[Join],
    adj: Dict[str, Set[str]],
    args,
    log: Log,
    end_adj: Optional[Dict[Tuple[str, str], Set[str]]] = None,
) -> List[Hypothesis]:
    """
    Enumerate every way of chaining the backbone segments into linear paths
    using a subset of the available joins, then rank. Nothing is chosen
    silently: the full ranked list is returned.
    """
    lengths = {c.name: c.length for c in calls}
    cn = {c.name: c.copy_number for c in calls}
    cls = {c.name: c.cls for c in calls}
    backbone = [c.name for c in calls if c.cls == "backbone"]
    if not backbone:
        return []
    # ---- centromere-like bridges -------------------------------------
    # A regional centromere in many filamentous fungi (and some plants) is a
    # long, markedly AT-rich block. It assembles badly: coverage drops, the
    # assembler cannot read through it, and it lands in the graph as a
    # low-depth segment with links on ONE side only. That last property is
    # exactly what makes a join through it "speculative" - so for this class of
    # segment, the missing through-path is the EXPECTED observation rather than
    # evidence against the join. Recognising the pattern lets the tool propose
    # the join on its own, instead of needing to be told the chromosome count.
    #
    # This is a hypothesis-raiser, never a diagnosis. AT-rich also describes
    # organelle-derived sequence, repeat families and plain compositional bias,
    # and centromeres are not AT-rich in every lineage. Nothing here calls the
    # segment a centromere; it only stops the tool from dismissing the join.
    gc_vals = [c.gc for c in calls if c.cls == "backbone" and c.gc is not None]
    baseline_gc = median(gc_vals) if gc_vals else None
    centromere_like: Dict[str, float] = {}
    if baseline_gc is not None:
        for c in calls:
            if c.gc is None or c.length < args.centromere_min_length:
                continue
            deficit = baseline_gc - c.gc
            if deficit < args.at_rich_delta:
                continue
            k = c.copy_number
            if k is None or k >= 1.5:
                continue  # a multi-copy block could sit anywhere; no constraint
            centromere_like[c.name] = deficit
    bridging = {v for j in joins for v in j.via}
    centromere_like = {k: v for k, v in centromere_like.items() if k in bridging}
    if centromere_like:
        log.info(
            "AT-rich low-copy bridge candidate(s), treated as possible centromeric "
            "sequence when they join exactly two backbone ends: "
            + ", ".join(
                f"{k} (GC {baseline_gc - v:.0%} vs baseline {baseline_gc:.0%})"
                for k, v in sorted(centromere_like.items())
            )
        )

    telomeric = telomeric_segments(calls, args)
    if telomeric:
        log.info(
            "telomere-bearing segment(s): "
            + ", ".join(
                f"{k} (array of {v} repeat units)" for k, v in sorted(telomeric.items())
            )
        )
    else:
        log.warn(
            "no segment carries a convincing telomere repeat array, so chromosome ends "
            "cannot be recognised and the number of chromosomes is only weakly constrained. "
            "Check --telomere-motif matches your organism, and that the GFA stores sequence."
        )

    # Collapse to the best route per pair, but remember the alternatives. "Best"
    # is not simply the shortest: a route through a segment that sits beside
    # three backbone segments says less than one through a segment that sits
    # beside exactly two, even at the same hop count.
    backbone_set = set(backbone)

    def direct_backbone_degree(seg: str) -> int:
        return len({n for n in adj.get(seg, ()) if n in backbone_set})

    def route_rank(j: Join) -> Tuple:
        worst = max((direct_backbone_degree(v) for v in j.via), default=0)
        return (len(j.via), worst, -min((lengths.get(v, 0) for v in j.via), default=0))

    best_per_pair: Dict[Tuple[str, str], Join] = {}
    alt_routes: Dict[Tuple[str, str], List[Join]] = defaultdict(list)
    for j in joins:
        alt_routes[j.key].append(j)
        cur = best_per_pair.get(j.key)
        if cur is None or route_rank(j) < route_rank(cur):
            best_per_pair[j.key] = j
    alt_count = {k: len(v) for k, v in alt_routes.items()}
    edges = list(best_per_pair.values())

    if len(edges) > args.max_join_edges:
        edges.sort(key=lambda j: (len(j.via), -min(lengths.get(j.a, 0), lengths.get(j.b, 0))))
        log.warn(
            f"{len(edges)} candidate joins exceed --max-join-edges "
            f"({args.max_join_edges}); only the {args.max_join_edges} shortest routes are "
            f"enumerated, so the hypothesis list is not exhaustive"
        )
        edges = edges[: args.max_join_edges]

    # How many backbone segments does each connector touch directly? A segment
    # sitting next to three backbone segments cannot say which two belong
    # together, so it is weak evidence for any particular pairing.
    backbone_set = set(backbone)
    connector_reach: Dict[str, Set[str]] = {
        v: {n for n in adj.get(v, ()) if n in backbone_set}
        for j in edges
        for v in j.via
    }

    results: List[Hypothesis] = []
    n_edges = len(edges)
    for mask in range(1 << n_edges):
        chosen = [edges[i] for i in range(n_edges) if mask >> i & 1]
        chains = _linear_forest(backbone, chosen)
        if chains is None:
            continue  # not a valid set of disjoint linear paths
        score, sup, con, res, capped, opened = _score_hypothesis(
            chains, chosen, lengths, cn, cls, connector_reach, alt_count, alt_routes,
            adj, telomeric, args, end_adj, centromere_like
        )
        spec = [j for j in chosen if j.speculative]
        if spec:
            # not established by the graph, so it must not be allowed to win on
            # score alone; it stays in the list, clearly labelled
            for j in spec:
                # The penalty exists because "no traversable path" normally means
                # the graph does not support the join. For an AT-rich centromeric
                # block the assembler is EXPECTED to fail to read through, so the
                # one-sidedness is explained rather than damning. Discounted, not
                # waived - the pattern is suggestive, not diagnostic.
                mid = j.via[0] if j.via else None
                if mid in centromere_like:
                    score -= args.speculative_penalty * args.centromere_speculative_discount
                    con = list(con) + [
                        f"the join {j.a} - {j.b} runs through {mid}, which has links on one "
                        f"side only, so this graph does not establish it. That is the "
                        f"expected appearance of an AT-rich centromeric block, which is why "
                        f"it is still ranked, but it remains unproven here."
                    ]
                else:
                    score -= args.speculative_penalty
                    con = list(con) + [
                        f"the join {j.a} - {j.b} is NOT supported by a traversable path: both "
                        f"segments simply end in {j.via[0]}, which has links on one side only. "
                        f"Resolving it needs evidence from outside this graph."
                    ]
                res = list(res) + [
                    f"long reads spanning {j.via[0]}, or Hi-C contact between {j.a} and "
                    f"{j.b}, would settle whether they join"
                ]
        results.append(
            Hypothesis(0, chains, chosen, score, sup, con, res, capped, opened)
        )
    results.sort(key=lambda h: (-h.score, len(h.joins)))
    for i, h in enumerate(results[: args.max_hypotheses], 1):
        h.rank = i
    top = results[: args.max_hypotheses]

    # flag ties explicitly - this is the honest part
    tol = args.tie_threshold
    if len(top) > 1 and abs(top[0].score - top[1].score) <= tol:
        tied = [h for h in top if abs(h.score - top[0].score) <= tol]
        note = (
            f"{len(tied)} hypotheses (ranks {', '.join(str(h.rank) for h in tied)}) score within "
            f"{tol} of each other. The graph cannot separate them; treat the top-ranked one as "
            f"one option among several, not as an answer."
        )
        for h in tied:
            h.contradicting.insert(0, note)
    log.info(
        f"{len(results)} valid chromosome hypotheses; reporting the top "
        f"{min(len(results), args.max_hypotheses)}"
    )
    return top


def _linear_forest(vertices: List[str], edges: List[Join]) -> Optional[List[List[str]]]:
    """
    Return the chains formed by these edges, or None if they do not form a set
    of vertex-disjoint simple paths (degree <= 2 everywhere, no cycles).
    """
    deg: Dict[str, int] = defaultdict(int)
    nbr: Dict[str, List[str]] = defaultdict(list)
    parent = {v: v for v in vertices}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    used_ends: Set[Tuple[str, str]] = set()
    for e in edges:
        if e.a not in parent or e.b not in parent:
            return None
        # A contig has two ends. Two joins cannot both attach to the same one,
        # so an end is consumed the first time a join uses it. Degree <= 2 alone
        # does not catch this: both joins at a vertex could be on one side.
        for end in e.ends:
            if end in used_ends:
                return None
            used_ends.add(end)
        deg[e.a] += 1
        deg[e.b] += 1
        if deg[e.a] > 2 or deg[e.b] > 2:
            return None
        ra, rb = find(e.a), find(e.b)
        if ra == rb:
            return None  # cycle
        parent[ra] = rb
        nbr[e.a].append(e.b)
        nbr[e.b].append(e.a)

    chains: List[List[str]] = []
    visited: Set[str] = set()
    ends = [v for v in vertices if deg[v] <= 1]
    for v in ends:
        if v in visited:
            continue
        chain, cur, prev = [v], v, None
        visited.add(v)
        while True:
            nxt = [x for x in nbr[cur] if x != prev]
            if not nxt:
                break
            prev, cur = cur, nxt[0]
            if cur in visited:
                return None
            visited.add(cur)
            chain.append(cur)
        chains.append(chain)
    if len(visited) != len(vertices):
        return None
    return chains


def _score_hypothesis(
    chains, chosen, lengths, cn, cls, connector_reach, alt_count, alt_routes,
    adj, telomeric, args, end_adj=None, centromere_like=None
):
    centromere_like = centromere_like or {}
    score = 0.0
    sup: List[str] = []
    con: List[str] = []
    res: List[str] = []

    total = sum(lengths.get(s, 0) for c in chains for s in c)
    n = len(chains)

    # ---- chromosome ends, which is how the count gets inferred --------
    # A finished linear chromosome has a telomere array at each end. Counting
    # capped versus open ends lets the data decide how many chromosomes there
    # are, instead of the number being supplied. It also gets the direction of
    # the evidence right: joining two open ends is progress, joining two
    # telomere-capped ends would be destroying a finished chromosome.
    internal = {v for j in chosen for v in j.via}
    capped_total = open_total = 0
    cap_notes: List[str] = []
    for chain in chains:
        capped, opened, notes = chain_end_status(chain, internal, adj, telomeric, end_adj)
        capped_total += capped
        open_total += opened
        cap_notes.extend(notes)
    # Capped ends are rewarded. Open ends are NOT penalised: an open end means
    # the molecule is unfinished, which is not evidence for any particular join.
    # Penalising them would make the tool invent joins to tidy the picture, and
    # it would always prefer merging everything into one chromosome.
    score += args.telomere_bonus * capped_total - args.open_end_penalty * open_total

    if telomeric:
        complete = sum(
            1
            for chain in chains
            if chain_end_status(chain, internal, adj, telomeric, end_adj)[0] == 2
        )
        sup.append(
            f"{capped_total} of {2 * n} molecule ends are capped by a telomere repeat; "
            f"{complete} of {n} molecule(s) are capped at both ends"
        )
        if cap_notes:
            sup.append("; ".join(cap_notes[:6]))
        if open_total:
            con.append(
                f"{open_total} end(s) are open: no telomere array is adjacent, so those "
                f"molecules are unfinished and two of them could still be one chromosome"
            )
            res.append(
                "longer reads, or a telomere-to-telomere assembly, would close the open ends "
                "and fix the chromosome count outright"
            )

    if args.expected_chromosomes:
        diff = abs(n - args.expected_chromosomes)
        score -= 3.0 * diff
        if diff == 0:
            sup.append(f"{n} chains matches the expected chromosome count")
        else:
            con.append(
                f"{n} chains against {args.expected_chromosomes} expected "
                f"({'too many' if n > args.expected_chromosomes else 'too few'})"
            )
            res.append(
                "a karyotype, pulsed-field gel, or Hi-C contact map would settle the chromosome "
                "count directly"
            )
    if args.expected_genome_size:
        rel = abs(total - args.expected_genome_size) / float(args.expected_genome_size)
        score -= 6.0 * rel
        if rel <= 0.03:
            sup.append(
                f"backbone totals {human_bp(total)} against {human_bp(args.expected_genome_size)} "
                f"expected ({rel:.1%} difference), so the nuclear genome looks essentially complete"
            )
        else:
            con.append(
                f"backbone totals {human_bp(total)}, {rel:.0%} away from the expected "
                f"{human_bp(args.expected_genome_size)}"
            )

    for j in chosen:
        # A join is an assertion about chromosome structure, so it starts in
        # deficit and has to earn its way out on the evidence of the segment it
        # runs through.
        score -= args.join_cost
        score -= 0.25 * max(0, len(j.via) - 1)
        if not j.via:
            score += 0.3  # a direct link needs no intermediate to be believed
        detail = []
        ambiguous = False
        for v in j.via:
            reach = connector_reach.get(v, set())
            c = cn.get(v)
            unique_bridge = c is not None and c < 1.5 and len(reach) <= 2
            if unique_bridge:
                score += 0.7
                detail.append(
                    f"{v} is present in about {c:.2f} copies and touches only "
                    f"{len(reach)} backbone segment(s), so it is a unique bridge rather than a "
                    f"repeat that could sit anywhere"
                )
                if v in centromere_like and len(reach) == 2:
                    score += args.centromere_bonus
                    detail.append(
                        f"{v} is {human_bp(lengths.get(v, 0))} and markedly AT-rich "
                        f"({centromere_like[v]:.0%} below the assembly's GC baseline), and it "
                        f"bridges exactly two backbone ends. In lineages whose centromeres are "
                        f"long AT-rich regions - many filamentous fungi among them - that is "
                        f"the shape of a CENTROMERE sitting between the two arms of one "
                        f"chromosome, and it also explains the low depth. Treat it as a "
                        f"hypothesis to test, not a call: AT-rich equally describes "
                        f"organelle-derived sequence, a repeat family, or compositional bias"
                    )
            elif c is not None and 1.5 <= c <= 3.5 and len(reach) <= 2:
                score += 0.45
                detail.append(
                    f"{v} at {c:.1f} copies is consistent with joining exactly two loci"
                )
            elif c is not None and c > 3.5:
                score -= 0.5
                ambiguous = True
                detail.append(
                    f"{v} at {c:.1f} copies could sit at many loci, so it constrains little"
                )
            if len(reach) > 2:
                score -= 0.5
                ambiguous = True
                detail.append(
                    f"{v} sits directly beside {len(reach)} backbone segments "
                    f"({', '.join(sorted(reach))}), so this particular pairing is one of several"
                )
            if cls.get(v) == "low_coverage":
                detail.append(
                    f"{v} is below single-copy depth; it is easy to dismiss but it does change "
                    f"the topology"
                )
        if alt_count.get(j.key, 1) > 1:
            others = [
                " - ".join(o.via) or "direct"
                for o in alt_routes.get(j.key, [])
                if o.via != j.via
            ]
            detail.append(
                f"{alt_count[j.key]} distinct routes exist between {j.a} and {j.b}; the others "
                f"run through {', '.join(others[:3])}"
            )
        sup.append(f"join {j.describe()}" + (": " + "; ".join(detail) if detail else ""))
        if ambiguous:
            longest_via = max((lengths.get(v, 0) for v in j.via), default=0)
            res.append(
                f"a read or scaffold spanning {human_bp(longest_via)} of "
                f"{', '.join(j.via) or 'the junction'} would test the {j.a}-{j.b} join directly; "
                f"Hi-C or a reference alignment would do the same"
            )

    if not chosen:
        sup.append("no joins asserted: every backbone segment is treated as its own molecule")

    # deduplicate while preserving order
    def uniq(xs: List[str]) -> List[str]:
        seen, out = set(), []
        for x in xs:
            if x not in seen:
                seen.add(x)
                out.append(x)
        return out

    return score, uniq(sup), uniq(con), uniq(res), capped_total, open_total
