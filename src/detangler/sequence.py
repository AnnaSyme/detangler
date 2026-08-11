"""Sequence-level measurements: GC, telomere arrays, periodicity."""
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
    human_bp,
    revcomp,
)


# Organelle size envelopes. Deliberately wide: plant mitogenomes reach several
# hundred kb, and some animal plastid-free lineages sit at the bottom edge.
MITO_RANGE = (11_000, 2_000_000)
MITO_TYPICAL = (14_000, 700_000)
PLASTID_RANGE = (70_000, 250_000)
PLASTID_TYPICAL = (110_000, 180_000)

# Canonical telomere repeat motifs. These are the published consensus motifs for
# the listed groups; add your own with --telomere-motif rather than assuming the
# organism is covered here.
TELOMERE_MOTIFS = {
    "TTAGGG": "vertebrates, most fungi incl. Ascomycota, many others",
    "TTTAGGG": "most land plants",
    "TTTTAGGG": "some green algae",
    "TTAGG": "many insects and other arthropods",
    "TTAGGC": "Caenorhabditis and other nematodes",
    "TTGGGG": "Tetrahymena and some other ciliates",
}


def gc_fraction(seq: str) -> Optional[float]:
    up = seq.upper()
    acgt = len(up) - up.count("N")
    if acgt <= 0:
        return None
    return (up.count("G") + up.count("C")) / acgt


def count_telomere_motifs(seq: str, motifs: Seq[str]) -> Dict[str, int]:
    """
    Non-overlapping counts of each motif and its reverse complement, ANYWHERE in
    the sequence.

    These are raw occurrences and are NOT evidence of a telomere. A 9 Mb contig
    contains roughly 2*L/4**k copies of any k-mer by chance - about 17,600
    TTAGGs - so this count says almost nothing on its own, and the shortest
    motif in the set always wins a naive argmax. Use find_telomere_arrays for
    anything that has to be true. Kept only to report background composition.
    """
    up = seq.upper()
    out: Dict[str, int] = {}
    for m in motifs:
        n = up.count(m.upper()) + up.count(revcomp(m))
        if n:
            out[m.upper()] = n
    return out


def find_telomere_arrays(
    seq: str, motifs: Seq[str], window: int, min_units: int
) -> Dict[str, Tuple[str, int, int]]:
    """
    Real telomere evidence: a TANDEM ARRAY of a motif, near an END of the
    sequence. Returns {end: (motif, units, offset)} for ends 's' (sequence
    start) and 'e' (sequence end), where offset is the distance in bp from that
    end to the array.

    Three things distinguish this from counting k-mers. The search is confined
    to a window at each end, because a telomere that is not at an end is not a
    telomere. Only consecutive, perfect repeats count, so a scattered handful of
    hits scores nothing. And candidates are ranked by the LENGTH IN BASES of the
    array rather than by the number of hits, which stops a short motif from
    beating a long one purely by being short.

    A run of even three units is ~4**-18 per position by chance, so min_units
    can be low without admitting noise; the returned unit count is what tells
    you how convincing a call really is.
    """
    up = seq.upper()
    n = len(up)
    if not n:
        return {}
    # The two windows must not overlap. On a segment SHORTER than twice the
    # window, an uncapped window makes both regions cover the whole sequence,
    # so one array near the start gets reported at BOTH ends - which is how a
    # one-sided telomeric tip came to look capped on both sides.
    w = min(max(window, 1), max(n // 2, 1))
    regions = {"s": (0, w), "e": (max(n - w, 0), n)}
    out: Dict[str, Tuple[str, int, int]] = {}
    for end, (lo, hi) in regions.items():
        sub = up[lo:hi]
        best: Optional[Tuple[int, str, int, int]] = None  # (bp, motif, units, offset)
        for m in motifs:
            mu = m.upper()
            for variant in {mu, revcomp(mu)}:
                if not variant:
                    continue
                for match in re.finditer(f"(?:{re.escape(variant)})+", sub):
                    units = len(match.group(0)) // len(variant)
                    if units < min_units:
                        continue
                    bp = units * len(variant)
                    start = lo + match.start()
                    offset = start if end == "s" else n - (start + len(match.group(0)))
                    cand = (bp, mu, units, offset)
                    if best is None or cand[0] > best[0]:
                        best = cand
        if best:
            out[end] = (best[1], best[2], best[3])
    return out


def estimate_period(seq: str, max_period: int = 3000, min_identity: float = 0.85) -> Optional[Tuple[int, float]]:
    """
    Smallest period p at which the sequence repeats itself, by direct comparison
    of s[i] with s[i+p] over sampled positions. Returns (period, identity).

    This is a screening heuristic, not a substitute for TRF or ULTRA. It will
    miss periods with indels, because it does no alignment.
    """
    up = seq.upper()
    n = len(up)
    if n < 40:
        return None
    limit = min(max_period, n // 2)
    step = max(1, (n - limit) // 4000)  # cap the work regardless of input size
    best: Optional[Tuple[int, float]] = None
    for p in range(2, limit + 1):
        span = n - p
        if span < 20:
            break
        same = total = 0
        for i in range(0, span, step):
            total += 1
            if up[i] == up[i + p]:
                same += 1
        if total and same / total >= min_identity:
            best = (p, same / total)
            break  # smallest period wins
    return best


def find_inverted_repeat_pair(
    seq: str,
    min_block: int = 10_000,
    k: int = 25,
    step: int = 200,
    max_scan: int = 400_000,
) -> Optional[Tuple[int, int, int]]:
    """
    Coarse scan for a large inverted-repeat pair: a block of at least min_block
    whose reverse complement also occurs elsewhere in the same sequence, the
    way the IRa/IRb pair does in most land-plant plastomes.

    Method: every k-mer of the reverse complement is indexed exactly, then the
    forward sequence is sampled every `step` bp and looked up. A genuine
    inverted pair puts its hits on one diagonal (constant offset between the
    forward position and the reverse-complement position), so hits are grouped
    into diagonal bands and the longest run wins. Returns
    (approximate_block_length, first_copy_start, second_copy_start) or None.

    Deliberately coarse and dependency-free: it answers "does this sequence
    contain a duplicated inverted block of plastome-IR scale" and nothing
    finer. Resolution is ~step bp, copies are assumed near-identical (plastome
    IR copies are), and sequences longer than max_scan are skipped (None)
    rather than scanned, keeping it O(n) in time and memory.
    """
    n = len(seq)
    if n < 2 * min_block + 1 or n > max_scan:
        return None
    fwd = seq.upper()
    rc = revcomp(fwd)
    index: Dict[str, List[int]] = {}
    for p in range(n - k + 1):
        index.setdefault(rc[p : p + k], []).append(p)

    # hits per diagonal band; band width = step so small indels stay together
    bands: Dict[int, List[int]] = defaultdict(list)
    for i in range(0, n - k + 1, step):
        hits = index.get(fwd[i : i + k])
        if not hits or len(hits) > 25:  # skip low-complexity k-mers
            continue
        for p in hits:
            partner = n - p - k  # start of the partner copy, forward coordinates
            if partner < i + k:  # partner must lie strictly downstream:
                continue  # a duplicated pair, not the centre of one palindrome
            bands[(p - i) // step].append(i)

    best: Optional[Tuple[int, int, int]] = None
    max_gap = 6 * step  # tolerate a few broken samples inside one block
    for band, starts in bands.items():
        starts.sort()
        run_start = prev = starts[0]
        for i in starts[1:] + [starts[-1] + max_gap + 1]:  # sentinel closes the last run
            if i - prev > max_gap:
                block = prev - run_start + k
                if block >= min_block and (best is None or block > best[0]):
                    d = band * step
                    second = n - (d + prev) - k  # approx start of the partner copy
                    best = (block, run_start, max(second, run_start + block))
                run_start = i
            prev = i
    return best


ORGANELLE_SUBTYPE_LABEL = {
    "mitochondrion": "mitochondrion-like",
    "plastid": "plastid-like",
    "unresolved": "type unresolved",
}


def classify_organelle_subtype(
    length: int,
    sequence: str,
    gc: Optional[float],
    baseline_gc: Optional[float],
) -> Tuple[str, str, Optional[Tuple[int, int, int]]]:
    """
    Say which organelle an organelle-candidate segment resembles, from
    graph-intrinsic and sequence-intrinsic evidence only. Returns
    (subtype, evidence_sentence, ir_block) with subtype one of
    "mitochondrion", "plastid", "unresolved".

    The one signal treated as strong is the plastome inverted-repeat pair:
    most land-plant plastomes are ~110-200 kb circles carrying two large
    inverted repeats (roughly 10-30 kb each) that separate the large and small
    single-copy regions. Size is never allowed to decide on its own, because
    mitogenome size varies enormously by lineage: ~15-20 kb in most animals,
    tens of kb in many fungi, and in plants frequently larger than a plastome.
    The figures above are broad documented ranges used as heuristics, not
    measurements, and the subtype stays a candidate either way - definitive
    organelle identification is delegated to dedicated external tools.
    """
    ir = find_inverted_repeat_pair(sequence) if sequence else None
    in_plastome_range = PLASTID_RANGE[0] <= length <= PLASTID_RANGE[1]
    plastome_typical = PLASTID_TYPICAL[0] <= length <= PLASTID_TYPICAL[1]

    if not sequence:
        sub = "unresolved"
        why = (
            "the GFA stores no sequence for this segment, so the plastome "
            "inverted-repeat test could not run, and size plus copy number alone "
            "cannot separate a mitogenome from a plastome"
        )
    elif ir and in_plastome_range:
        sub = "plastid"
        why = (
            f"a duplicated inverted block of ~{human_bp(ir[0])} was found (copies "
            f"near {human_bp(ir[1])} and {human_bp(ir[2])}) - matching the large "
            f"inverted-repeat pair, typically ~10-30 kb per copy, that divides most "
            f"land-plant plastomes into their single-copy regions - and "
            f"{human_bp(length)} sits inside the documented plastome size range"
        )
    elif ir:
        sub = "unresolved"
        why = (
            f"a duplicated inverted block of ~{human_bp(ir[0])} was found, but "
            f"{human_bp(length)} falls outside the documented plastome size range; "
            f"large inverted repeats also occur in some mitogenomes, so the two "
            f"signals disagree"
        )
    elif plastome_typical:
        sub = "unresolved"
        why = (
            f"no large (>= ~10 kb) inverted-repeat pair was found, yet "
            f"{human_bp(length)} is plastome-typical; this could be a plastome from "
            f"one of the lineages that lost an inverted-repeat copy, or a compact "
            f"plant mitogenome - the sequence alone cannot say"
        )
    else:
        sub = "mitochondrion"
        why = (
            f"no large (>= ~10 kb) inverted-repeat pair was found, so the hallmark "
            f"plastome architecture is absent (a minority of plastome lineages do "
            f"lack it), and {human_bp(length)} is compatible with a mitogenome; "
            f"mitogenome size varies too much by lineage (~15-20 kb in most "
            f"animals, tens of kb in many fungi, often plastome-sized or larger in "
            f"plants) for size to decide on its own"
        )
    if gc is not None and baseline_gc is not None:
        why += (
            f"; GC {gc:.0%} against an assembly baseline of {baseline_gc:.0%}, "
            f"recorded as context only - organellar GC varies too widely across "
            f"lineages to score"
        )
    return sub, f"organelle subtype: {ORGANELLE_SUBTYPE_LABEL[sub]} - {why}", ir
