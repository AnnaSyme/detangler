"""The drawable model, and reading or writing it as config."""
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
    VERSION,
    _natural_key,
    human_bp,
    median,
    nx_stat,
)
from .records import (
    Anchor,
    CoverageWindow,
    Evidence,
    GfaLink,
    GfaSegment,
    Placement,
    SeqRecord,
    Tangle,
)
from .graph import (
    build_adjacency,
    build_end_adjacency,
    find_bubbles,
    find_circular,
    find_inverted,
)
from .sequence import (
    MITO_RANGE,
    MITO_TYPICAL,
    ORGANELLE_SUBTYPE_LABEL,
    PLASTID_RANGE,
    PLASTID_TYPICAL,
)
from .palette import (
    CLASS_LABEL,
    assign_segment_colours,
)
from .calls import (
    SegmentCall,
    telomeric_segments,
)
from .hypotheses import (
    Hypothesis,
    chain_end_status,
)



# ==========================================================================
# classification
# ==========================================================================
NAME_PATTERNS = [
    # (regex, role, weight, human description)
    (r"^(chr)?(\d{1,2}|[IVXL]{1,5})$", "chromosome", 3.0, "name looks like a chromosome number"),
    (r"^chr[_\-]?\w{1,4}$", "chromosome", 2.0, "name has a chr prefix"),
    (r"(^|[_\W])(chr)?(mt|mito|mitochondri\w*|chrM|MT)([_\W]|$)", "mitochondrion", 6.0,
     "name contains a mitochondrial keyword"),
    (r"(^|[_\W])(chloroplast|plastid|chrC|cp|pltd|ptg_cp)([_\W]|$)", "plastid", 6.0,
     "name contains a plastid keyword"),
    (r"(^|[_\W])(chrX|chrY|chrZ|chrW|X|Y|Z|W)([_\W]|$)", "chromosome", 2.5,
     "name looks like a sex chromosome"),
    (r"(scaffold|scf|ctg|contig|utg|tig|ptg|jcf|unplaced|unloc|random|chrUn)", "unplaced", 2.5,
     "name looks like an unplaced contig/scaffold"),
    (r"(^|[_\W])(super[_\-]?scaffold|SUPER)", "chromosome", 1.0,
     "name uses a super-scaffold convention"),
]


def classify_sequences(
    seqs: List[SeqRecord],
    circular: Set[str],
    seg_depth_by_seq: Dict[str, float],
    args,
    log: Log,
) -> None:
    """
    Assign a role to every sequence with a recorded evidence trail.

    The rule that does most of the work is a length break: sequences at or above
    a threshold derived from the assembly's own size distribution are treated as
    chromosome candidates. The threshold is min(--min-chrom-length, a fraction
    of the longest sequence) so that it adapts to genomes of any size.
    """
    if not seqs:
        return

    lengths = [s.length for s in seqs]
    longest = max(lengths)
    total = sum(lengths)
    n50 = nx_stat(lengths, 50)

    if args.min_chrom_length is not None:
        chrom_cut = args.min_chrom_length
        cut_reason = f"--min-chrom-length {chrom_cut:,}"
    else:
        chrom_cut = max(int(longest * args.chrom_length_fraction), 100_000)
        cut_reason = (
            f"{args.chrom_length_fraction:.0%} of the longest sequence "
            f"({human_bp(longest)}) = {human_bp(chrom_cut)}"
        )
    log.info(f"chromosome length threshold: {cut_reason}")

    # nuclear depth baseline = depth-weighted median over long sequences
    long_depths = [
        seg_depth_by_seq[s.name]
        for s in seqs
        if s.name in seg_depth_by_seq and s.length >= chrom_cut
    ]
    nuclear_depth = median(long_depths) if long_depths else median(list(seg_depth_by_seq.values()))

    gc_values = [s.gc for s in seqs if s.gc is not None and s.length >= chrom_cut]
    nuclear_gc = median(gc_values) if gc_values else None

    for s in seqs:
        if s.manual:  # role already fixed by assembly report or config
            continue
        scores: Dict[str, float] = defaultdict(float)
        ev: List[Evidence] = []

        # --- name -------------------------------------------------------
        for pattern, role, weight, desc in NAME_PATTERNS:
            if re.search(pattern, s.name, re.IGNORECASE):
                scores[role] += weight
                ev.append(Evidence("name", f"{desc} ('{s.name}')", weight))

        # --- length -----------------------------------------------------
        if s.length >= chrom_cut:
            scores["chromosome"] += 4.0
            ev.append(
                Evidence("length", f"{human_bp(s.length)} >= threshold {human_bp(chrom_cut)}", 4.0)
            )
        else:
            scores["unplaced"] += 2.5
            ev.append(
                Evidence("length", f"{human_bp(s.length)} < threshold {human_bp(chrom_cut)}", 2.5)
            )

        # --- organelle size envelope -----------------------------------
        if MITO_RANGE[0] <= s.length <= MITO_RANGE[1]:
            w = 2.0 if MITO_TYPICAL[0] <= s.length <= MITO_TYPICAL[1] else 0.8
            scores["mitochondrion"] += w
            ev.append(Evidence("size", f"{human_bp(s.length)} is within the mitogenome range", w))
        if PLASTID_RANGE[0] <= s.length <= PLASTID_RANGE[1]:
            w = 2.0 if PLASTID_TYPICAL[0] <= s.length <= PLASTID_TYPICAL[1] else 0.8
            scores["plastid"] += w
            ev.append(Evidence("size", f"{human_bp(s.length)} is within the plastome range", w))

        # --- circularity from the graph --------------------------------
        if s.name in circular:
            s.circular = True
            scores["mitochondrion"] += 2.5
            scores["plastid"] += 2.0
            scores["chromosome"] -= 1.0
            ev.append(
                Evidence(
                    "graph",
                    "segment carries a self-link, i.e. the graph closes it into a circle",
                    2.5,
                )
            )

        # --- depth ------------------------------------------------------
        d = seg_depth_by_seq.get(s.name)
        if d is not None:
            s.depth = d
            if nuclear_depth and nuclear_depth > 0:
                ratio = d / nuclear_depth
                if ratio >= args.organelle_depth_ratio:
                    w = min(4.0, 1.5 * math.log10(max(ratio, 1.01)) + 1.5)
                    scores["mitochondrion"] += w
                    scores["plastid"] += w * 0.8
                    ev.append(
                        Evidence(
                            "depth",
                            f"{d:.1f}x is {ratio:.1f}x the nuclear baseline "
                            f"({nuclear_depth:.1f}x)",
                            w,
                        )
                    )
                elif 0.6 <= ratio <= 1.6:
                    scores["chromosome"] += 1.0
                    ev.append(
                        Evidence("depth", f"{d:.1f}x matches the nuclear baseline", 1.0)
                    )

        # --- GC ---------------------------------------------------------
        if s.gc is not None and nuclear_gc is not None and s.length < chrom_cut:
            delta = abs(s.gc - nuclear_gc)
            if delta >= 0.08:
                scores["mitochondrion"] += 1.0
                scores["plastid"] += 1.0
                ev.append(
                    Evidence(
                        "gc",
                        f"GC {s.gc:.1%} differs from the nuclear median {nuclear_gc:.1%} "
                        f"by {delta:.1%}",
                        1.0,
                    )
                )

        if not scores:
            scores["unplaced"] = 0.1
        role = max(scores.items(), key=lambda kv: kv[1])[0]

        # An organelle call needs a plausible size no matter what the name says.
        if role == "mitochondrion" and not (MITO_RANGE[0] <= s.length <= MITO_RANGE[1]):
            ev.append(
                Evidence(
                    "veto",
                    f"mitochondrion rejected: {human_bp(s.length)} is outside "
                    f"{human_bp(MITO_RANGE[0])}-{human_bp(MITO_RANGE[1])}",
                    -5.0,
                )
            )
            scores["mitochondrion"] = -5.0
            role = max(scores.items(), key=lambda kv: kv[1])[0]
        if role == "plastid" and not (PLASTID_RANGE[0] <= s.length <= PLASTID_RANGE[1]):
            ev.append(
                Evidence(
                    "veto",
                    f"plastid rejected: {human_bp(s.length)} is outside the plastome range",
                    -5.0,
                )
            )
            scores["plastid"] = -5.0
            role = max(scores.items(), key=lambda kv: kv[1])[0]

        s.role = role
        s.role_scores = dict(scores)
        s.evidence = ev

    # Only one mitochondrion and one plastid can be right; keep the best.
    for role in ("mitochondrion", "plastid"):
        cands = [s for s in seqs if s.role == role and not s.manual]
        if len(cands) > 1:
            cands.sort(key=lambda s: s.role_scores.get(role, 0), reverse=True)
            keep = cands[0]
            for s in cands[1:]:
                s.evidence.append(
                    Evidence(
                        "uniqueness",
                        f"demoted: '{keep.name}' scored higher for {role} "
                        f"({keep.role_scores.get(role, 0):.1f} vs "
                        f"{s.role_scores.get(role, 0):.1f}); could be a duplicate haplotype, "
                        f"a NUMT/NUPT, or a genuine multipartite organelle genome",
                        -3.0,
                    )
                )
                s.role = "chromosome" if s.length >= chrom_cut else "unplaced"
            log.warn(
                f"{len(cands)} sequences looked like a {role}; kept '{keep.name}'. "
                f"Check the report - multipartite organelle genomes and NUMTs both do this."
            )

    n_chrom = sum(1 for s in seqs if s.role == "chromosome")
    if n_chrom == 0:
        log.warn(
            "no sequence was called a chromosome. The assembly may be contig-level; "
            "consider --min-chrom-length or a manual config."
        )
    log.info(
        f"classified {len(seqs)} sequences: {n_chrom} chromosome(s), "
        f"{sum(1 for s in seqs if s.role == 'mitochondrion')} mitochondrion, "
        f"{sum(1 for s in seqs if s.role == 'plastid')} plastid, "
        f"{sum(1 for s in seqs if s.role == 'unplaced')} unplaced "
        f"(assembly {human_bp(total)}, N50 {human_bp(n50)})"
    )


def build_placements(
    segs: Dict[str, GfaSegment],
    seqs: List[SeqRecord],
    agp: List[Placement],
    paf: List[Placement],
    manual_map: List[Placement],
    log: Log,
) -> Dict[str, List[Placement]]:
    """
    Locate each graph segment on the assembled sequences. Precedence:
    explicit map > AGP > PAF > name identity. A segment may legitimately have
    many placements; that is what a repeat looks like.
    """
    seqnames = {s.name for s in seqs}
    placements: Dict[str, List[Placement]] = defaultdict(list)

    for label, source in (("map", manual_map), ("agp", agp), ("paf", paf)):
        if not source:
            continue
        kept = 0
        for p in source:
            if p.seqname not in seqnames:
                continue
            placements[p.segment].append(p)
            kept += 1
        if kept:
            log.info(f"placed {len(placements)} graph segments from {label} ({kept} records)")
            break  # highest-precedence source that produced anything wins outright
        log.warn(
            f"{label} source had {len(source)} records but none referenced a known "
            f"sequence name - check that names match between files"
        )

    if not placements:
        identical = [n for n in segs if n in seqnames]
        for n in identical:
            rec = next(s for s in seqs if s.name == n)
            placements[n].append(
                Placement(segment=n, seqname=n, start=0, end=rec.length, source="identity")
            )
        if identical:
            log.info(
                f"no AGP/PAF/map supplied; matched {len(identical)}/{len(segs)} graph segments "
                f"to sequences by name"
            )

    unplaced = [n for n in segs if n not in placements]
    if unplaced:
        log.warn(
            f"{len(unplaced)} of {len(segs)} graph segments could not be located on any "
            f"sequence. Tangles involving only those segments cannot be drawn. Supply "
            f"--agp or --paf (minimap2 the graph segments against the assembly) to fix this."
        )
    return placements


def analyse_graph(
    segs: Dict[str, GfaSegment],
    links: List[GfaLink],
    placements: Dict[str, List[Placement]],
    seqs: List[SeqRecord],
    args,
    log: Log,
) -> List[Tangle]:
    """
    Turn the assembly graph into a list of located, typed tangles.

    Four independent detectors run; a segment can raise more than one.
      1. multi-mapping segment          -> shared repeat between sequences
      2. link-based junction            -> neighbours sit on different sequences
      3. depth outlier + high degree    -> collapsed repeat
      4. self-links                     -> circular or inverted repeat
    """
    tangles: List[Tangle] = []
    adj = build_adjacency(links)
    role_of = {s.name: s.role for s in seqs}
    drawable = {s.name for s in seqs if s.role in ("chromosome", "mitochondrion", "plastid")}

    # depth baseline over unique-looking segments (degree <= 2, decent length)
    unique_depths = [
        s.depth
        for s in segs.values()
        if s.depth is not None and len(adj.get(s.name, ())) <= 2 and s.length >= args.min_segment_length
    ]
    base_depth = median(unique_depths) if unique_depths else None
    if base_depth:
        log.info(f"unique-segment depth baseline: {base_depth:.1f}x")

    def anchors_for(seg_name: str, restrict: Optional[Set[str]] = None) -> List[Anchor]:
        out = []
        for p in placements.get(seg_name, []):
            if restrict is not None and p.seqname not in restrict:
                continue
            out.append(Anchor(seqname=p.seqname, start=p.start, end=p.end, segment=seg_name))
        return out

    counter = 0

    def next_id(kind: str) -> str:
        nonlocal counter
        counter += 1
        return f"T{counter:03d}_{kind}"

    # --- 1. multi-mapping segments ------------------------------------
    for name, seg in segs.items():
        places = [p for p in placements.get(name, []) if p.seqname in drawable]
        if len(places) < 2 or seg.length < args.min_segment_length:
            continue
        distinct_seqs = {p.seqname for p in places}
        anchors = [Anchor(p.seqname, p.start, p.end, name) for p in places]
        depth_ratio = (seg.depth / base_depth) if (seg.depth and base_depth) else None
        if len(distinct_seqs) >= 2:
            ttype = "interchromosomal_junction"
            desc = (
                f"Segment {name} ({human_bp(seg.length)}) is placed on "
                f"{len(distinct_seqs)} sequences ({', '.join(sorted(distinct_seqs))}). "
                f"A single graph node shared between chromosomes is a repeat that the "
                f"assembler could not separate."
            )
        else:
            ttype = "intrachromosomal_repeat"
            only = next(iter(distinct_seqs))
            desc = (
                f"Segment {name} ({human_bp(seg.length)}) is placed {len(places)} times on "
                f"{only}, i.e. a repeat family within one chromosome."
            )
        evidence = [f"{len(places)} placements across {len(distinct_seqs)} sequence(s)"]
        if depth_ratio:
            evidence.append(f"depth {seg.depth:.1f}x = {depth_ratio:.1f}x the unique baseline")
        tangles.append(
            Tangle(
                id=next_id(ttype),
                type=ttype,
                segments=[name],
                anchors=anchors,
                multiplicity=round(depth_ratio) if depth_ratio else float(len(places)),
                depth_ratio=depth_ratio,
                length=seg.length,
                description=desc,
                evidence=evidence,
            )
        )

    # --- 2. link-based junctions --------------------------------------
    # Segments already reported above; a node whose neighbour is one of these is
    # only "joining" chromosomes because of that repeat, so reporting it again
    # would list the same event once per flanking unitig.
    repeat_segments = {
        s for t in tangles if t.type == "interchromosomal_junction" for s in t.segments
    }
    seen_pairs: Set[Tuple[str, ...]] = set()
    for name, nbrs in adj.items():
        if name not in segs or segs[name].length < args.min_segment_length:
            continue
        if name in repeat_segments or (nbrs & repeat_segments):
            continue
        # which sequences do this node's neighbours live on?
        nbr_seqs: Dict[str, List[Anchor]] = defaultdict(list)
        for nb in nbrs:
            for a in anchors_for(nb, drawable):
                nbr_seqs[a.seqname].append(a)
        own = {a.seqname for a in anchors_for(name, drawable)}
        touched = set(nbr_seqs) | own
        if len(touched) < 2:
            continue
        key = tuple(sorted(touched) + [name])
        if key in seen_pairs:
            continue
        seen_pairs.add(key)
        anchors = anchors_for(name, drawable)
        for seqname, al in nbr_seqs.items():
            anchors.extend(al[:2])
        if len({a.seqname for a in anchors}) < 2:
            continue
        # skip if the multi-mapping detector already reported this segment
        if any(name in t.segments and t.type == "interchromosomal_junction" for t in tangles):
            continue
        deg = len(nbrs)
        tangles.append(
            Tangle(
                id=next_id("interchromosomal_junction"),
                type="interchromosomal_junction",
                segments=[name] + sorted(nbrs),
                anchors=anchors,
                multiplicity=None,
                depth_ratio=(segs[name].depth / base_depth) if (segs[name].depth and base_depth) else None,
                length=segs[name].length,
                description=(
                    f"Graph node {name} (degree {deg}) links sequences "
                    f"{', '.join(sorted(touched))}. In Bandage this is the junction where "
                    f"those chromosomes appear to merge."
                ),
                evidence=[f"degree {deg}", f"neighbours placed on {len(touched)} sequences"],
            )
        )

    # --- 3. depth outliers --------------------------------------------
    if base_depth:
        for name, seg in segs.items():
            if seg.depth is None or seg.length < args.min_segment_length:
                continue
            ratio = seg.depth / base_depth
            if ratio < args.collapse_depth_ratio:
                continue
            if any(name in t.segments for t in tangles):
                continue
            anchors = anchors_for(name, drawable)
            if not anchors:
                continue
            # Organelles are legitimately present at high copy number per cell;
            # their depth says nothing about a collapsed nuclear repeat.
            if all(role_of.get(a.seqname) in ("mitochondrion", "plastid") for a in anchors):
                continue
            tangles.append(
                Tangle(
                    id=next_id("collapsed_repeat"),
                    type="collapsed_repeat",
                    segments=[name],
                    anchors=anchors,
                    multiplicity=round(ratio),
                    depth_ratio=ratio,
                    length=seg.length,
                    description=(
                        f"Segment {name} ({human_bp(seg.length)}) carries {seg.depth:.1f}x "
                        f"coverage, {ratio:.1f}x the unique baseline. Consistent with roughly "
                        f"{round(ratio)} copies collapsed into one node."
                    ),
                    evidence=[
                        f"depth {seg.depth:.1f}x vs baseline {base_depth:.1f}x",
                        f"degree {len(adj.get(name, ()))}",
                    ],
                )
            )

    # --- 4. self-links -------------------------------------------------
    for name in find_inverted(links):
        anchors = anchors_for(name, drawable)
        if not anchors:
            continue
        tangles.append(
            Tangle(
                id=next_id("inverted_repeat"),
                type="inverted_repeat",
                segments=[name],
                anchors=anchors,
                length=segs[name].length if name in segs else None,
                description=(
                    f"Segment {name} links to itself in the opposite orientation - an "
                    f"inverted repeat or hairpin. These are a common cause of assembly "
                    f"breaks and of the loops seen in Bandage."
                ),
                evidence=["self-link with opposite orientation"],
            )
        )

    for name in find_circular(links):
        anchors = anchors_for(name, drawable)
        if not anchors:
            continue
        tangles.append(
            Tangle(
                id=next_id("circular"),
                type="circular",
                segments=[name],
                anchors=anchors,
                length=segs[name].length if name in segs else None,
                description=f"Segment {name} closes on itself: the graph supports a circular molecule.",
                evidence=["self-link in consistent orientation"],
            )
        )

    # --- 5. bubbles -----------------------------------------------------
    for a, b, members in find_bubbles(adj, segs, args.max_bubble_length):
        anchors: List[Anchor] = []
        for m in members + [a, b]:
            anchors.extend(anchors_for(m, drawable)[:1])
        if not anchors:
            continue
        tangles.append(
            Tangle(
                id=next_id("bubble"),
                type="bubble",
                segments=members,
                anchors=anchors[:2],
                length=max((segs[m].length for m in members if m in segs), default=None),
                description=(
                    f"{len(members)} alternative short paths between {a} and {b} "
                    f"({', '.join(members)}). Typically residual heterozygosity in a "
                    f"collapsed diploid assembly, or a small structural variant."
                ),
                evidence=[f"{len(members)} parallel paths sharing the same two neighbours"],
            )
        )

    log.info(f"detected {len(tangles)} graph tangle(s)")
    return tangles


def summarise_coverage(
    windows: List[CoverageWindow], seqs: List[SeqRecord], args, log: Log
) -> Tuple[Dict[str, List[CoverageWindow]], Dict[str, float], List[Dict]]:
    """Bin coverage per sequence, find the per-sequence median, flag outliers."""
    by_seq: Dict[str, List[CoverageWindow]] = defaultdict(list)
    known = {s.name for s in seqs}
    for w in windows:
        if w.seqname in known:
            by_seq[w.seqname].append(w)
    for name in by_seq:
        by_seq[name].sort(key=lambda w: w.start)

    medians = {name: median([w.depth for w in ws]) for name, ws in by_seq.items()}
    chrom_medians = [
        medians[s.name] for s in seqs if s.role == "chromosome" and s.name in medians
    ]
    global_median = median(chrom_medians) if chrom_medians else median(list(medians.values()))

    anomalies: List[Dict] = []
    if global_median > 0:
        for name, ws in by_seq.items():
            run: Optional[Dict] = None
            for w in ws:
                r = w.depth / global_median
                kind = (
                    "high"
                    if r >= args.coverage_high_ratio
                    else ("low" if r <= args.coverage_low_ratio else None)
                )
                if kind and run and run["kind"] == kind and w.start <= run["end"] + 1:
                    run["end"] = w.end
                    run["peak"] = max(run["peak"], r) if kind == "high" else min(run["peak"], r)
                elif kind:
                    if run:
                        anomalies.append(run)
                    run = {
                        "seqname": name,
                        "start": w.start,
                        "end": w.end,
                        "kind": kind,
                        "peak": r,
                    }
                elif run:
                    anomalies.append(run)
                    run = None
            if run:
                anomalies.append(run)
    anomalies = [a for a in anomalies if a["end"] - a["start"] >= args.min_anomaly_length]
    if anomalies:
        log.info(
            f"coverage: {len(anomalies)} region(s) outside "
            f"{args.coverage_low_ratio}x-{args.coverage_high_ratio}x the genome median "
            f"({global_median:.1f}x)"
        )
    return by_seq, medians, anomalies


# ==========================================================================
# model assembly, config round trip
# ==========================================================================
class Model:
    """Everything the renderers and the report need."""

    def __init__(self) -> None:
        self.sequences: List[SeqRecord] = []
        self.tangles: List[Tangle] = []
        self.coverage: Dict[str, List[CoverageWindow]] = {}
        self.coverage_median: Dict[str, float] = {}
        self.coverage_anomalies: List[Dict] = []
        self.annotations: List[Dict] = []
        self.inputs: Dict[str, str] = {}
        self.settings: Dict[str, object] = {}
        self.warnings: List[str] = []
        self.title: str = "Assembly ideogram"
        # graph-first mode
        self.segment_calls: List["SegmentCall"] = []
        self.hypotheses: List["Hypothesis"] = []
        self.baseline_depth: Optional[float] = None
        self.baseline_basis: str = ""
        self.chosen_hypothesis: int = 0
        self.segment_colours: Dict[str, str] = {}
        self.foreign_gc_delta: float = 0.05
        self.candidate_reasons: Dict[str, List[str]] = {}

    # -- ordering ------------------------------------------------------
    def drawable(self) -> List[SeqRecord]:
        order = str(self.settings.get("order", "length"))
        chroms = [s for s in self.sequences if s.role == "chromosome"]
        if order == "name":
            chroms.sort(key=lambda s: _natural_key(s.display))
        elif order == "file":
            pass
        else:
            chroms.sort(key=lambda s: -s.length)
        organelles = [s for s in self.sequences if s.role in ("mitochondrion", "plastid")]
        organelles.sort(key=lambda s: s.role)
        return chroms + organelles

    def unplaced(self) -> List[SeqRecord]:
        return [s for s in self.sequences if s.role == "unplaced"]

    def count_range(self) -> Optional[Tuple[int, int, int]]:
        """
        (best, low, high) number of linear molecules. The range spans the
        hypotheses the graph cannot tell apart, so it collapses to a single
        number only when the evidence really does pin it down.

        Two kinds of "cannot tell apart" count. Hypotheses within the tie
        threshold of the best, and hypotheses that differ only by a speculative
        join - two contigs ending in the same one-sided repeat. The second kind
        scores lower on purpose, but excluding it from the range would state a
        chromosome count more confidently than the evidence allows.
        """
        if not self.hypotheses:
            return None
        best = self.hypotheses[0]
        near = [
            h
            for h in self.hypotheses
            if any("score within" in c for c in h.contradicting)
            or any(j.speculative for j in h.joins)
        ] or [best]
        counts = [len(h.chains) for h in near] + [len(best.chains)]
        return len(best.chains), min(counts), max(counts)

    def unassigned(self) -> List[SeqRecord]:
        """Sequences that fit no molecule and are shown in their own panel."""
        return sorted(
            (s for s in self.sequences if s.role == "unassigned"), key=lambda s: -s.length
        )

    def summary_sentence(self) -> str:
        chroms = [s for s in self.sequences if s.role == "chromosome"]
        parts = [f"{len(chroms)} chromosome-scale sequence{'s' if len(chroms) != 1 else ''}"]
        rng = self.count_range()
        if rng:
            best, low, high = rng
            complete = sum(1 for s in self.sequences if s.note == "telomere to telomere")
            span = "" if low == high else f", though the graph supports anywhere from {low} to {high}"
            parts[0] = (
                f"{best} linear molecule{'s' if best != 1 else ''} "
                f"({complete} capped by telomeres at both ends){span}"
            )
        ua = self.unassigned()
        if ua:
            parts.append(
                f"{len(ua)} sequence{'s' if len(ua) != 1 else ''} that fit no chromosome, "
                f"totalling {human_bp(sum(s.length for s in ua))}"
            )
        for role, word in (("mitochondrion", "a mitochondrial genome"), ("plastid", "a plastid genome")):
            if any(s.role == role for s in self.sequences):
                parts.append(word)
        up = self.unplaced()
        if up:
            parts.append(
                f"{len(up)} unplaced sequence{'s' if len(up) != 1 else ''} "
                f"totalling {human_bp(sum(s.length for s in up))}"
            )
        return "This assembly resolves into " + ", ".join(parts[:-1]) + (
            f" and {parts[-1]}." if len(parts) > 1 else parts[0] + "."
        )


def model_to_config(model: Model) -> Dict:
    return {
        "detangler_version": VERSION,
        "note": (
            "Edit this file and re-run with --config to override any call. "
            "'role' accepts chromosome, mitochondrion, plastid, unplaced or exclude. "
            "Set a tangle's 'show' to false to hide it, or edit its 'description' - "
            "your text is used verbatim in the report and the HTML."
        ),
        "inputs": model.inputs,
        "settings": model.settings,
        "sequences": [
            {
                "name": s.name,
                "label": s.display,
                "length": s.length,
                "role": s.role,
                "circular": s.circular,
                "depth": round(s.depth, 2) if s.depth is not None else None,
                "gc": round(s.gc, 4) if s.gc is not None else None,
                "call_confidence": _confidence(s),
                "why": [e.as_text() for e in s.evidence],
            }
            for s in model.sequences
        ],
        "tangles": [
            {
                "id": t.id,
                "type": t.type,
                "show": True,
                "segments": t.segments,
                "length": t.length,
                "multiplicity": t.multiplicity,
                "depth_ratio": round(t.depth_ratio, 2) if t.depth_ratio else None,
                "description": t.description,
                "why": t.evidence,
                "anchors": [
                    {"seqname": a.seqname, "start": a.start, "end": a.end, "segment": a.segment}
                    for a in t.anchors
                ],
            }
            for t in model.tangles
        ],
    }


def _confidence(s: SeqRecord) -> str:
    if s.manual:
        return "asserted"
    if not s.role_scores:
        return "unknown"
    ordered = sorted(s.role_scores.values(), reverse=True)
    top = ordered[0]
    second = ordered[1] if len(ordered) > 1 else 0.0
    gap = top - second
    if gap >= 4:
        return "high"
    if gap >= 2:
        return "medium"
    return "low"


def write_config(model: Model, path: str) -> str:
    cfg = model_to_config(model)
    if HAVE_YAML and path.endswith((".yaml", ".yml")):
        with open(path, "w") as fh:
            yaml.safe_dump(cfg, fh, sort_keys=False, width=100, allow_unicode=True)
    else:
        if path.endswith((".yaml", ".yml")):
            path = os.path.splitext(path)[0] + ".json"
        with open(path, "w") as fh:
            json.dump(cfg, fh, indent=2)
    return path


def load_config(path: str) -> Dict:
    with open(path) as fh:
        if path.endswith((".yaml", ".yml")):
            if not HAVE_YAML:
                raise SystemExit(
                    "This config is YAML but PyYAML is not installed. "
                    "Install PyYAML, or re-run without --config to emit a JSON config."
                )
            return yaml.safe_load(fh)
        return json.load(fh)


def apply_overrides(model: Model, cfg: Dict, log: Log) -> None:
    """Config wins over inference, everywhere."""
    by_name = {s.name: s for s in model.sequences}
    keep: List[SeqRecord] = []
    seen: Set[str] = set()
    for entry in cfg.get("sequences") or []:
        name = entry.get("name")
        if name is None:
            continue
        rec = by_name.get(name)
        if rec is None:
            rec = SeqRecord(name=name, length=int(entry.get("length") or 0))
            model.sequences.append(rec)
        seen.add(name)
        role = entry.get("role")
        if role == "exclude":
            continue
        if role and role != rec.role:
            rec.evidence.append(
                Evidence("override", f"role set to '{role}' in the config file", 99.0)
            )
            rec.role = role
            rec.manual = True
        if entry.get("label"):
            rec.label = entry["label"]
        if entry.get("circular") is not None:
            rec.circular = bool(entry["circular"])
        keep.append(rec)
    if seen:
        model.sequences = [s for s in model.sequences if s.name in seen and s in keep] or keep

    cfg_tangles = {t.get("id"): t for t in (cfg.get("tangles") or [])}
    if cfg_tangles:
        rebuilt: List[Tangle] = []
        existing = {t.id: t for t in model.tangles}
        for tid, entry in cfg_tangles.items():
            if entry.get("show") is False:
                continue
            t = existing.get(tid)
            if t is None:
                anchors = [
                    Anchor(a["seqname"], int(a["start"]), int(a["end"]), a.get("segment", ""))
                    for a in (entry.get("anchors") or [])
                ]
                if not anchors:
                    continue
                t = Tangle(
                    id=tid or f"cfg{len(rebuilt)}",
                    type=entry.get("type", "collapsed_repeat"),
                    segments=entry.get("segments") or [],
                    anchors=anchors,
                    multiplicity=entry.get("multiplicity"),
                    depth_ratio=entry.get("depth_ratio"),
                    length=entry.get("length"),
                    evidence=entry.get("why") or [],
                )
            if entry.get("description"):
                t.description = entry["description"]
            if entry.get("type"):
                t.type = entry["type"]
            rebuilt.append(t)
        model.tangles = rebuilt
    if cfg.get("settings"):
        model.settings.update(cfg["settings"])
    log.info(f"applied config overrides from {cfg.get('detangler_version', cfg.get('karyoglyph_version', '?'))}")



# --------------------------------------------------------------------------
# turn a hypothesis into the karyotype the ideogram draws
# --------------------------------------------------------------------------
def model_from_hypothesis(
    hyp: Hypothesis,
    calls: List[SegmentCall],
    links: List[GfaLink],
    args,
    log: Log,
) -> Model:
    by_name = {c.name: c for c in calls}
    adj = build_adjacency(links)
    end_adj = build_end_adjacency(links)
    colours = assign_segment_colours(calls)
    telomeric = telomeric_segments(calls, args)
    model = Model()
    model.settings = {"order": args.order, "coverage": False}
    model.segment_colours = colours

    used_in_chains: Set[str] = set()
    placements: Dict[str, List[Anchor]] = defaultdict(list)
    # For each anchor laid down inside a chain: which member precedes and which
    # follows it, so hanging features can be anchored to the correct end.
    anchor_context: Dict[Tuple[str, int, int, str], Tuple[Optional[str], Optional[str]]] = {}

    for idx, chain in enumerate(hyp.chains, 1):
        join_by_pair = {j.key: j for j in hyp.joins}
        members: List[str] = []
        for i, seg in enumerate(chain):
            members.append(seg)
            if i + 1 < len(chain):
                j = join_by_pair.get(tuple(sorted((seg, chain[i + 1]))))  # type: ignore
                if j:
                    members.extend(j.via)
        name = f"chain_{idx}"
        label = f"chain {idx}: " + " + ".join(chain)
        cursor = 0
        blocks: List[Tuple[int, int, str, str]] = []
        for i, seg in enumerate(members):
            L = by_name[seg].length if seg in by_name else 0
            placements[seg].append(Anchor(name, cursor, cursor + L, seg))
            anchor_context[(name, cursor, cursor + L, seg)] = (
                members[i - 1] if i > 0 else None,
                members[i + 1] if i + 1 < len(members) else None,
            )
            blocks.append((cursor, cursor + L, seg, colours.get(seg, "#cfcfcf")))
            used_in_chains.add(seg)
            cursor += L
        rec = SeqRecord(name=name, length=cursor, role="chromosome", label=label, manual=True)
        rec.blocks = blocks
        rec.blocks_tile = True

        # Repeats hanging off the FREE ends of this molecule. They are not part
        # of the chain - the graph does not let you walk through them - but they
        # are the most informative thing at a chromosome end (a telomeric array,
        # an rDNA block), and leaving them out of the figure made the right panel
        # look like the assembly had no repeats at all.
        used_ends = {e for j in hyp.joins for e in j.ends}
        backbone_names = {c.name for c in calls if c.cls == "backbone"}
        rec.caps = {}
        if len(chain) == 1:
            terminals = [(chain[0], "s", "top"), (chain[0], "e", "bottom")]
        else:
            terminals = [
                (chain[0], end, "top") for end in ("s", "e")
            ] + [(chain[-1], end, "bottom") for end in ("s", "e")]
        for seg, end, side in terminals:
            if (seg, end) in used_ends:
                continue
            nb = sorted(
                n for n in end_adj.get((seg, end), ())
                if n not in backbone_names and n not in members
            )
            if nb:
                rec.caps.setdefault(side, []).extend(
                    (n, colours.get(n, "#cfcfcf")) for n in nb
                )
        capped, opened, cap_notes = chain_end_status(
            chain, {v for j in hyp.joins for v in j.via}, adj, telomeric, end_adj
        )
        rec.evidence.append(
            Evidence(
                "hypothesis",
                f"hypothesis {hyp.rank}: {' + '.join(chain)}"
                + (f" joined via {', '.join(j.describe() for j in hyp.joins)}" if hyp.joins else ""),
                0.0,
            )
        )
        rec.evidence.append(
            Evidence(
                "ends",
                f"{capped} of 2 ends capped by a telomere repeat"
                + (f" ({'; '.join(cap_notes)})" if cap_notes else "")
                + (
                    f"; {opened} end(s) open, so this molecule may be part of a larger one"
                    if opened
                    else "; complete end to end"
                ),
                0.0,
            )
        )
        rec.note = (
            "telomere to telomere"
            if capped == 2
            else f"{opened} open end{'s' if opened != 1 else ''}"
        )
        model.sequences.append(rec)

    # organelle candidates get their own molecule
    for c in calls:
        if c.cls == "organelle_candidate":
            subtype = c.organelle_subtype or "unresolved"
            rec = SeqRecord(
                name=c.name,
                length=c.length,
                # the role slot only knows mitochondrion/plastid; an unresolved
                # candidate is filed under mitochondrion for layout, and the
                # label carries the honest call
                role="plastid" if subtype == "plastid" else "mitochondrion",
                label=f"{c.name} (organelle candidate, "
                f"{ORGANELLE_SUBTYPE_LABEL[subtype]})",
                circular=c.self_loop_same_orient,
                depth=c.depth,
                manual=True,
            )
            rec.evidence.append(Evidence("class", "; ".join(c.reasons), 0.0))
            model.sequences.append(rec)
            used_in_chains.add(c.name)

    # Anything left over sits in its own panel rather than being forced into a
    # chromosome or quietly dropped. The reason it did not fit is recorded, and
    # contamination is only ever offered as a candidate explanation.
    nuclear_component = None
    backbone_calls = [c for c in calls if c.cls == "backbone"]
    if backbone_calls:
        counts: Dict[int, int] = defaultdict(int)
        for c in backbone_calls:
            counts[c.component] += c.length
        nuclear_component = max(counts.items(), key=lambda kv: kv[1])[0]
    nuclear_gc = median([c.gc for c in backbone_calls if c.gc is not None]) or None

    for c in calls:
        if c.name in used_in_chains:
            continue
        why: List[str] = []
        isolated = c.degree == 0 or (
            nuclear_component is not None and c.component != nuclear_component
        )
        if isolated:
            why.append("does not connect to the nuclear part of the graph")
        gc_off = (
            nuclear_gc is not None
            and c.gc is not None
            and abs(c.gc - nuclear_gc) >= args.foreign_gc_delta
        )
        if gc_off:
            why.append(f"GC {c.gc:.0%} against {nuclear_gc:.0%} for the backbone")
        depth_off = c.copy_number is not None and (
            c.copy_number < args.low_coverage_max_copy or c.copy_number > 3.0
        )
        if depth_off:
            why.append(f"copy number {c.copy_number:.2f}, unlike the single-copy backbone")
        if isolated and (gc_off or depth_off):
            note = "candidate contaminant or foreign sequence: " + "; ".join(why)
        elif isolated:
            note = "unknown: " + "; ".join(why)
        elif c.length < args.backbone_min_length:
            note = f"too short to anchor a chromosome ({human_bp(c.length)})"
            why.append(note)
        else:
            note = "connected to the graph but not placed on any molecule"
            why.append(note)

        rec = SeqRecord(
            name=c.name, length=c.length, role="unassigned", depth=c.depth, manual=True
        )
        rec.note = note
        rec.evidence.append(Evidence("class", "; ".join(c.reasons) or CLASS_LABEL[c.cls], 0.0))
        rec.evidence.append(Evidence("unassigned", note, 0.0))
        model.sequences.append(rec)

    # tangles: repeats and arrays, positioned where the graph allows
    tangles: List[Tangle] = []
    counter = 0
    # Low-coverage segments are drawn too. They look ignorable on depth, but they
    # carry the topology, so hiding them would misrepresent the graph.
    class_to_tangle = {
        "repeat": "repeat_segment",
        "tandem_array": "tandem_array",
        "at_rich": "at_rich_region",
        "low_coverage": "low_coverage_region",
    }
    def facing_side(seg: str, other: Optional[str]) -> Optional[str]:
        """The physical end of seg that carries a link to other, if unambiguous."""
        if other is None:
            return None
        sides = [s for s in ("s", "e") if other in end_adj.get((seg, s), ())]
        return sides[0] if len(sides) == 1 else None

    for c in calls:
        if c.cls not in class_to_tangle:
            continue
        anchors = list(placements.get(c.name, []))
        approximate = False
        if not anchors:
            # Not inside a chain: hang it off each backbone neighbour instead,
            # at the end of that neighbour the link actually touches. Members
            # tile along the chain in order, so the end of a member that faces
            # the next member sits at anchor.end, the end facing the previous
            # member at anchor.start, and a single-member chain is drawn in
            # its forward orientation.
            for nb in sorted(adj.get(c.name, ())):
                touching = [s for s in ("s", "e") if c.name in end_adj.get((nb, s), ())]
                for a in placements.get(nb, []):
                    prev_seg, next_seg = anchor_context.get(
                        (a.seqname, a.start, a.end, a.segment), (None, None)
                    )
                    to_prev = facing_side(nb, prev_seg)
                    to_next = facing_side(nb, next_seg)
                    positions: List[int] = []
                    for side in touching:
                        if next_seg is not None and side == to_next:
                            positions.append(a.end)
                        elif prev_seg is not None and side == to_prev:
                            positions.append(a.start)
                        elif prev_seg is None and next_seg is not None and to_next is not None:
                            positions.append(a.start)  # outer end of the first member
                        elif next_seg is None and prev_seg is not None and to_prev is not None:
                            positions.append(a.end)  # outer end of the last member
                        elif prev_seg is None and next_seg is None:
                            positions.append(a.start if side == "s" else a.end)
                        else:
                            positions.append(a.end)  # orientation unresolved
                    if not positions:
                        positions.append(a.end)  # no orientation recorded for this link
                    for pos in sorted(set(positions)):
                        anchors.append(
                            Anchor(a.seqname, max(pos - 1, a.start), max(pos, a.start + 1), c.name)
                        )
                    approximate = True
        if not anchors:
            continue
        counter += 1
        ttype = class_to_tangle[c.cls]
        spans = sorted({a.seqname for a in anchors})
        desc = "; ".join(c.reasons)
        if len(spans) > 1:
            desc += (
                f". Present on {len(spans)} of the hypothesised molecules "
                f"({', '.join(spans)}), which is why they appear to join in the graph"
            )
        if c.at_rich and c.cls != "at_rich":
            desc += ". Also AT-rich"
        if approximate:
            desc += (
                ". Position is approximate: the graph places this segment beside a backbone "
                "segment but does not say where inside it."
            )
        tangles.append(
            Tangle(
                id=f"S{counter:03d}_{c.cls}",
                type=ttype,
                segments=[c.name],
                anchors=anchors,
                multiplicity=round(c.copy_number) if c.copy_number else None,
                depth_ratio=c.copy_number,
                length=c.length,
                description=f"{c.name} ({CLASS_LABEL[c.cls]}): {desc}",
                evidence=[f"copy number {c.copy_number:.1f}" if c.copy_number else "no depth"],
            )
        )
    model.tangles = tangles

    # A segment that is not its own molecule is still drawn, as a feature on the
    # chain it belongs to. Record that, so "unplaced" is not read as "ignored".
    drawn_on: Dict[str, Set[str]] = defaultdict(set)
    for t in tangles:
        for seg in t.segments:
            drawn_on[seg].update(a.seqname for a in t.anchors)
    for rec in model.sequences:
        if rec.role == "unassigned" and drawn_on.get(rec.name):
            rec.role = "unplaced"  # it is drawn on a molecule, just not as one
            rec.evidence.append(
                Evidence(
                    "placement",
                    f"not a molecule in its own right, but drawn as a feature on "
                    f"{', '.join(sorted(drawn_on[rec.name]))}",
                    0.0,
                )
            )

    model.title = args.title or "Chromosome hypothesis from the assembly graph"
    return model
