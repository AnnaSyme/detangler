"""Segment classification from graph-intrinsic evidence."""
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
    human_bp,
    median,
)
from .records import (
    ContigInfo,
    GfaLink,
    GfaSegment,
)
from .parsers import (
    resolve_path_element,
)
from .graph import (
    build_adjacency,
    components,
    find_circular,
    find_inverted,
)
from .sequence import (
    ORGANELLE_SUBTYPE_LABEL,
    TELOMERE_MOTIFS,
    classify_organelle_subtype,
    count_telomere_motifs,
    estimate_period,
    find_telomere_arrays,
    gc_fraction,
)
from .palette import (
    CLASS_LABEL,
)



@dataclass
class SegmentCall:
    """A segment with its observations kept separate from its inferences."""

    name: str
    # observations, straight from the files
    length: int
    depth: Optional[float]
    degree: int
    self_loop_same_orient: bool
    self_loop_flipped: bool
    component: int
    component_size: int
    gc: Optional[float] = None
    telomere_motifs: Dict[str, int] = field(default_factory=dict)
    # {end: (motif, units, offset_bp)} - the only telomere evidence worth acting on
    telomere_arrays: Dict[str, Tuple[str, int, int]] = field(default_factory=dict)
    path_terminal: int = 0
    path_interior: int = 0
    path_consecutive_repeats: int = 0
    contigs: List[str] = field(default_factory=list)
    # derived estimates
    copy_number: Optional[float] = None
    period: Optional[Tuple[int, float]] = None
    # inference
    cls: str = "unclassified"
    at_rich: bool = False
    organelle_subtype: Optional[str] = None  # mitochondrion | plastid | unresolved
    ir_block: Optional[Tuple[int, int, int]] = None  # (approx len, copy1 start, copy2 start)
    reasons: List[str] = field(default_factory=list)
    identity_hits: List[Dict] = field(default_factory=list)

    @property
    def one_line(self) -> str:
        bits = [f"{human_bp(self.length)}"]
        if self.depth is not None:
            bits.append(f"{self.depth:.0f}x")
        if self.copy_number is not None:
            bits.append(f"~{self.copy_number:.1f} copies")
        return f"{self.name}: " + ", ".join(bits)


def call_segments(
    segs: Dict[str, GfaSegment],
    links: List[GfaLink],
    contigs: List[ContigInfo],
    seq_by_segment: Dict[str, str],
    args,
    log: Log,
) -> Tuple[List[SegmentCall], float]:
    """
    Steps 1-5 of the design notes: baseline depth, copy number, classification,
    composition features, and position within each contig's graph path.
    """
    adj = build_adjacency(links)
    circ_same = find_circular(links)
    circ_flip = find_inverted(links)
    comp = components(adj, list(segs))
    comp_size: Dict[int, int] = defaultdict(int)
    for n in segs:
        comp_size[comp.get(n, 0)] += 1

    # ---- step 1: baseline single-copy depth ----
    long_depths = [
        s.depth for s in segs.values() if s.depth is not None and s.length >= args.baseline_min_length
    ]
    if long_depths:
        baseline = median(long_depths)
        basis = (
            f"median depth of the {len(long_depths)} segment(s) at least "
            f"{human_bp(args.baseline_min_length)} long"
        )
    else:
        ranked = sorted((s for s in segs.values() if s.depth is not None), key=lambda s: -s.length)
        top = ranked[: max(1, len(ranked) // 10)]
        baseline = median([s.depth for s in top]) if top else 0.0
        basis = (
            f"no segment reached {human_bp(args.baseline_min_length)}; fell back to the median "
            f"depth of the longest {len(top)} segment(s)"
        )
        log.warn(
            f"baseline depth derived from a fallback ({basis}). Copy numbers below are less "
            f"reliable than usual; set --baseline-min-length to something appropriate."
        )
    if not baseline:
        log.warn("no depth information in the GFA; copy number cannot be estimated")
    log.info(f"baseline single-copy depth: {baseline:.1f}x ({basis})")

    # ---- step 5: position within graph_path ----
    terminal: Dict[str, int] = defaultdict(int)
    interior: Dict[str, int] = defaultdict(int)
    consec: Dict[str, int] = defaultdict(int)
    seg_contigs: Dict[str, List[str]] = defaultdict(list)
    for c in contigs:
        resolved = [(el, resolve_path_element(el, segs)) for el in c.path]
        real = [(i, name) for i, (el, name) in enumerate(resolved) if name]
        if not real:
            continue
        first_i, last_i = real[0][0], real[-1][0]
        prev_name = None
        run = 1
        for i, name in real:
            if c.name not in seg_contigs[name]:
                seg_contigs[name].append(c.name)
            if i in (first_i, last_i):
                terminal[name] += 1
            else:
                interior[name] += 1
            if name == prev_name:
                run += 1
                consec[name] = max(consec[name], run)
            else:
                run = 1
            prev_name = name

    calls: List[SegmentCall] = []
    gc_all = []
    for name, seg in segs.items():
        sequence = seq_by_segment.get(name, "")
        gc = seg.gc if seg.gc is not None else (gc_fraction(sequence) if sequence else None)
        if gc is not None and seg.length >= args.baseline_min_length:
            gc_all.append(gc)
        calls.append(
            SegmentCall(
                name=name,
                length=seg.length,
                depth=seg.depth,
                degree=len(adj.get(name, ())),
                self_loop_same_orient=name in circ_same,
                self_loop_flipped=name in circ_flip,
                component=comp.get(name, 0),
                component_size=comp_size.get(comp.get(name, 0), 1),
                gc=gc,
                telomere_motifs=count_telomere_motifs(sequence, args.telomere_motif)
                if sequence
                else {},
                telomere_arrays=find_telomere_arrays(
                    sequence, args.telomere_motif,
                    args.telomere_window, args.min_telomere_units,
                )
                if sequence
                else {},
                path_terminal=terminal.get(name, 0),
                path_interior=interior.get(name, 0),
                path_consecutive_repeats=consec.get(name, 0),
                contigs=seg_contigs.get(name, []),
                copy_number=(seg.depth / baseline) if (seg.depth and baseline) else None,
            )
        )
    if gc_all:
        baseline_gc = median(gc_all)
    else:
        # No long segment carried sequence (many GFAs store '*' for the big ones).
        # Fall back to a length-weighted mean over whatever sequence there is, and
        # say so, because short segments are exactly the biased ones.
        weighted = [(c.gc, c.length) for c in calls if c.gc is not None]
        tot_w = sum(w for _, w in weighted)
        baseline_gc = (sum(g * w for g, w in weighted) / tot_w) if tot_w else None
        if baseline_gc is not None:
            log.warn(
                f"no segment at least {human_bp(args.baseline_min_length)} long carried "
                f"sequence, so the GC baseline ({baseline_gc:.0%}) is a length-weighted mean "
                f"over short segments only. AT-rich calls made against it are weaker than usual."
            )

    # ---- steps 3-4: classification ----
    for c in calls:
        cn = c.copy_number
        r = c.reasons
        isolated = c.component_size == 1 or c.degree == 0

        if (
            (c.self_loop_same_orient or isolated)
            and cn is not None
            and cn >= args.organelle_min_copy
            and args.organelle_min_length <= c.length <= args.organelle_max_length
        ):
            c.cls = "organelle_candidate"
            r.append(
                f"{'circular (self-link)' if c.self_loop_same_orient else 'no links'}, "
                f"{'isolated component' if isolated else f'component of {c.component_size}'}, "
                f"{human_bp(c.length)}, {cn:.1f}x relative copy number - within the size range "
                f"of an organellar genome and far above single copy"
            )
            c.organelle_subtype, subtype_reason, c.ir_block = classify_organelle_subtype(
                c.length, seq_by_segment.get(c.name, ""), c.gc, baseline_gc
            )
            r.append(subtype_reason)
            log.info(
                f"{c.name}: organelle candidate, "
                f"{ORGANELLE_SUBTYPE_LABEL[c.organelle_subtype]}"
            )
        elif c.self_loop_same_orient and cn is not None and cn >= args.tandem_min_copy:
            c.cls = "tandem_array"
            r.append(
                f"self-link plus {cn:.0f}x relative copy number: a tandem array of roughly "
                f"{cn:.0f} copies of a {human_bp(c.length)} unit"
            )
        elif c.path_consecutive_repeats >= 2 and cn is not None and cn >= args.tandem_min_copy:
            c.cls = "tandem_array"
            r.append(
                f"traversed {c.path_consecutive_repeats} times consecutively in a contig path "
                f"and carries {cn:.0f}x copy number"
            )
        elif cn is not None and cn >= args.repeat_min_copy:
            c.cls = "repeat"
            r.append(f"copy number {cn:.1f} (>= {args.repeat_min_copy}) relative to single copy")
        elif cn is not None and cn < args.low_coverage_max_copy:
            # Around HALF the baseline is not "foreign" - in a primary assembly
            # it is the expected depth of a haplotig the purge step left behind,
            # and that is the commonest thing in the file. Calling it
            # contamination was the tool's least defensible rule. What decides
            # the reading is not the organism's ploidy but what the ASSEMBLER
            # emitted, which is why --assembly-type is a declaration.
            atype = getattr(args, "assembly_type", "primary")
            half_band = args.haplotig_band
            looks_haplotig = (
                atype in ("primary", "phased")
                and half_band[0] <= cn <= half_band[1]
                and c.length >= args.haplotig_min_length
            )
            if looks_haplotig:
                c.cls = "haplotig"
                r.append(
                    f"copy number {cn:.2f}, about half the single-copy baseline, over "
                    f"{human_bp(c.length)}. In a {atype} assembly that is what an unpurged "
                    f"HAPLOTIG looks like - the second allele of a heterozygous region that "
                    f"the assembler kept but did not merge. It is not evidence of "
                    f"contamination, and it should not be counted as a chromosome of its own"
                )
            else:
                c.cls = "low_coverage"
                r.append(
                    f"copy number {cn:.2f} (< {args.low_coverage_max_copy}): below single copy, "
                    f"so not simply a low-confidence unique segment - candidates include a "
                    f"haplotype-specific region, contamination, or a real sub-stoichiometric "
                    f"molecule"
                )
        elif c.length >= args.backbone_min_length and (cn is None or args.backbone_copy_range[0] <= cn <= args.backbone_copy_range[1]):
            c.cls = "backbone"
            r.append(
                f"{human_bp(c.length)} at "
                + (f"copy number {cn:.2f}" if cn is not None else "unknown copy number")
                + f" (>= {human_bp(args.backbone_min_length)}, near single copy)"
            )
        else:
            c.cls = "short_single_copy"
            r.append(
                f"{human_bp(c.length)}, "
                + (f"copy number {cn:.2f}" if cn is not None else "no depth")
                + ": too short to anchor a chromosome"
            )

        # composition, added as extra evidence rather than replacing the class
        if c.gc is not None and baseline_gc is not None and c.gc <= baseline_gc - args.at_rich_delta:
            c.at_rich = True
            r.append(
                f"GC {c.gc:.0%} against a genome baseline of {baseline_gc:.0%}: markedly AT-rich. "
                f"What that means depends on the lineage - in many fungi it marks centromeric or "
                f"subtelomeric sequence, but it can equally be a repeat family, organelle-derived "
                f"sequence, or plain compositional bias - so treat it as a flag, not a diagnosis"
            )
            if c.cls in ("short_single_copy", "unclassified"):
                c.cls = "at_rich"
        if c.telomere_arrays:
            for end in ("s", "e"):
                if end not in c.telomere_arrays:
                    continue
                motif, units, off = c.telomere_arrays[end]
                where = "start" if end == "s" else "end"
                r.append(
                    f"telomere repeat array at the {where}: ({motif})x{units}, "
                    f"{units * len(motif)} bp, {off:,} bp from the {where} "
                    f"({TELOMERE_MOTIFS.get(motif, 'user-supplied motif')})"
                )
        elif c.telomere_motifs:
            tot = sum(c.telomere_motifs.values())
            r.append(
                f"no telomere repeat array at either end; the {tot} scattered motif "
                f"occurrence(s) are consistent with chance and are not evidence of a "
                f"chromosome end"
            )
        if c.path_terminal or c.path_interior:
            if c.path_terminal > c.path_interior:
                r.append(
                    f"appears at the end of a contig path {c.path_terminal} time(s) versus "
                    f"{c.path_interior} interior: positionally subtelomeric"
                )
            elif c.path_interior:
                r.append(
                    f"appears interior to a contig path {c.path_interior} time(s) versus "
                    f"{c.path_terminal} terminal: not subtelomeric"
                )
        if c.self_loop_flipped:
            r.append("self-link in the opposite orientation: inverted repeat or hairpin")

        # periodicity, only where it is cheap and meaningful
        if (
            c.cls in ("tandem_array", "repeat")
            and seq_by_segment.get(c.name)
            and c.length <= args.max_period_scan_length
        ):
            c.period = estimate_period(seq_by_segment[c.name], args.max_period)
            if c.period:
                r.append(
                    f"internal periodicity at {c.period[0]:,} bp ({c.period[1]:.0%} identity "
                    f"between offset copies), so the {human_bp(c.length)} segment is itself "
                    f"built of a shorter unit"
                )

    n_by_class: Dict[str, int] = defaultdict(int)
    for c in calls:
        n_by_class[c.cls] += 1
    log.info(
        "segment classes: "
        + ", ".join(f"{CLASS_LABEL[k]}={v}" for k, v in sorted(n_by_class.items()))
    )
    return calls, baseline


def telomeric_segments(calls: List[SegmentCall], args) -> Dict[str, int]:
    """
    Segments carrying a genuine telomere repeat ARRAY at one or both ends.
    The value is the number of repeat units in the best array, which is the
    honest measure of how strong the call is.

    This used to threshold on whole-sequence motif density, which counts chance
    k-mers: a 9 Mb contig contains ~17,600 TTAGGs by accident and sailed through.
    Only tandem arrays near a sequence end count now.
    """
    out: Dict[str, int] = {}
    for c in calls:
        if not c.telomere_arrays:
            continue
        out[c.name] = max(units for _m, units, _off in c.telomere_arrays.values())
    return out
