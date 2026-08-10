"""Assembly graph topology: adjacency, ends, components."""
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
from .records import (
    GfaLink,
    GfaSegment,
)



# ==========================================================================
# graph analysis
# ==========================================================================
def build_adjacency(links: List[GfaLink]) -> Dict[str, Set[str]]:
    adj: Dict[str, Set[str]] = defaultdict(set)
    for l in links:
        if l.a == l.b:
            continue
        adj[l.a].add(l.b)
        adj[l.b].add(l.a)
    return adj


def _link_end(orient: str, is_source: bool) -> str:
    """
    Which physical end of a segment a GFA link touches. A link joins the end
    of its source segment (as oriented) to the start of its target (as
    oriented): '+' as source means the forward end ('e') and '-' the forward
    start ('s'); the reverse for the target.
    """
    if is_source:
        return "e" if orient == "+" else "s"
    return "s" if orient == "+" else "e"


def build_end_adjacency(links: List[GfaLink]) -> Dict[Tuple[str, str], Set[str]]:
    """
    (segment, end) -> the NAMES of the segments attached at that end.

    A view over build_end_links for callers that only need to know which
    segments touch an end, not which of their ends does the touching. Derived
    rather than computed separately so the orientation rule lives in one place.
    """
    return {
        key: {name for name, _end in partners}
        for key, partners in build_end_links(links).items()
    }


OTHER_END = {"s": "e", "e": "s"}


def build_end_links(
    links: List[GfaLink],
) -> Dict[Tuple[str, str], Set[Tuple[str, str]]]:
    """
    The full bidirected adjacency: (segment, end) -> {(neighbour, neighbour_end)}.

    build_end_adjacency keeps which of OUR ends a link touches but throws away
    which of the NEIGHBOUR's ends it lands on, which makes it useless for
    traversal: to walk through a segment you must know the end you arrived at so
    you can leave by the opposite one. Losing that is what let the tool propose
    routes that enter and leave a segment through the same end.
    """
    adj: Dict[Tuple[str, str], Set[Tuple[str, str]]] = defaultdict(set)
    for l in links:
        if l.a == l.b:
            continue
        ae = (l.a, _link_end(l.a_orient, True))
        be = (l.b, _link_end(l.b_orient, False))
        adj[ae].add(be)
        adj[be].add(ae)
    return adj


def dead_end_repeats(
    end_links: Dict[Tuple[str, str], Set[Tuple[str, str]]], segs: Iterable[str]
) -> Dict[str, str]:
    """
    Segments with links on one end only. Such a segment CANNOT be traversed, so
    it is not a bridge between two contigs: it is a tip. Biologically it is the
    signature of a repeat whose far side the assembler never resolved - very
    often a telomeric or subtelomeric array sitting at the ends of several
    chromosomes at once.

    Returns {segment: the end that carries the links}.
    """
    out: Dict[str, str] = {}
    for name in segs:
        has_s = bool(end_links.get((name, "s")))
        has_e = bool(end_links.get((name, "e")))
        if has_s != has_e:
            out[name] = "s" if has_s else "e"
    return out


def find_circular(links: List[GfaLink]) -> Set[str]:
    """Segments with a self-link in a consistent orientation, i.e. a circle."""
    return {l.a for l in links if l.a == l.b and l.a_orient == l.b_orient}


def find_inverted(links: List[GfaLink]) -> Set[str]:
    """Self-link with flipped orientation: a hairpin / inverted repeat."""
    return {l.a for l in links if l.a == l.b and l.a_orient != l.b_orient}


def find_bubbles(adj: Dict[str, Set[str]], segs: Dict[str, GfaSegment], max_len: int) -> List[Tuple[str, str, List[str]]]:
    """
    Superbubble-lite: pairs of short segments that share exactly the same two
    neighbours. That is the classic diploid heterozygosity signature and is
    cheap to compute without a full superbubble algorithm.
    """
    by_neighbours: Dict[Tuple[str, ...], List[str]] = defaultdict(list)
    for name, nbrs in adj.items():
        if len(nbrs) == 2 and segs.get(name) and segs[name].length <= max_len:
            by_neighbours[tuple(sorted(nbrs))].append(name)
    out = []
    for pair, members in by_neighbours.items():
        if len(members) >= 2:
            out.append((pair[0], pair[1], sorted(members)))
    return out


def components(adj: Dict[str, Set[str]], names: Seq[str]) -> Dict[str, int]:
    """Connected components over the undirected graph, self-loops ignored."""
    comp: Dict[str, int] = {}
    cid = 0
    for n in names:
        if n in comp:
            continue
        stack, cid = [n], cid + 1
        while stack:
            cur = stack.pop()
            if cur in comp:
                continue
            comp[cur] = cid
            stack.extend(x for x in adj.get(cur, ()) if x not in comp)
    return comp
