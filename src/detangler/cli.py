"""Argument parsing and the two top-level pipelines."""
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
    parse_size,
)
from .records import (
    GfaLink,
    GfaSegment,
    SeqRecord,
)
from .parsers import (
    parse_agp,
    parse_annotation,
    parse_assembly_report,
    parse_coverage_bed,
    parse_fai,
    parse_fasta,
    parse_flye_info,
    parse_gfa,
    parse_paf,
    parse_segment_map,
    read_segment_sequences,
)
from .graph import (
    build_adjacency,
    build_end_adjacency,
    build_end_links,
    dead_end_repeats,
    find_circular,
)
from .sequence import (
    TELOMERE_MOTIFS,
)
from .palette import (
    assign_segment_colours,
)
from .calls import (
    call_segments,
)
from .hypotheses import (
    Join,
    enumerate_hypotheses,
    find_joins,
)
from .model import (
    Model,
    analyse_graph,
    apply_overrides,
    build_placements,
    classify_sequences,
    load_config,
    model_from_hypothesis,
    summarise_coverage,
    write_config,
)
from .render_ideogram import (
    render_html,
    render_svg,
)
from .render_paired import (
    write_figures,
)
from .report import (
    render_report,
)
from .blast import (
    identify_repeats,
)
from .demo import (
    write_demo,
)



def build_model(args, log: Log) -> Tuple[Model, Dict[str, str]]:
    model = Model()
    model.settings = {
        "order": args.order,
        "coverage": not args.no_coverage_track,
    }

    # ---- sequences ---------------------------------------------------
    if args.fasta:
        seqs = parse_fasta(args.fasta)
        model.inputs["fasta"] = args.fasta
    elif args.fai:
        seqs = parse_fai(args.fai)
        model.inputs["fai"] = args.fai
    elif args.assembly_report:
        seqs = parse_assembly_report(args.assembly_report)
        model.inputs["assembly_report"] = args.assembly_report
    else:
        raise SystemExit(
            "Need sequence lengths: pass one of --fasta, --fai or --assembly-report "
            "(or use --demo)."
        )
    log.info(f"loaded {len(seqs)} sequences, {human_bp(sum(s.length for s in seqs))} total")

    # ---- graph -------------------------------------------------------
    segs: Dict[str, GfaSegment] = {}
    links: List[GfaLink] = []
    if args.gfa:
        segs, links = parse_gfa(args.gfa, log, args.assembler)
        model.inputs["gfa"] = args.gfa
        log.info(f"graph: {len(segs)} segments, {len(links)} links")

    agp = parse_agp(args.agp, log) if args.agp else []
    paf = parse_paf(args.paf, args.min_paf_block, log) if args.paf else []
    smap = parse_segment_map(args.segment_map, log) if args.segment_map else []
    for key, val in (("agp", args.agp), ("paf", args.paf), ("segment_map", args.segment_map)):
        if val:
            model.inputs[key] = val

    placements = build_placements(segs, seqs, agp, paf, smap, log) if segs else {}

    # depth and circularity projected from graph segments onto sequences
    circular_segments = find_circular(links)
    circular_seqs: Set[str] = set()
    depth_by_seq: Dict[str, List[float]] = defaultdict(list)
    for seg_name, places in placements.items():
        seg = segs.get(seg_name)
        for p in places:
            if seg_name in circular_segments:
                circular_seqs.add(p.seqname)
            if seg and seg.depth is not None:
                depth_by_seq[p.seqname].append(seg.depth)
    seq_depth = {k: median(v) for k, v in depth_by_seq.items() if v}

    classify_sequences(seqs, circular_seqs, seq_depth, args, log)
    model.sequences = seqs

    # ---- tangles -----------------------------------------------------
    if segs:
        model.tangles = analyse_graph(segs, links, placements, seqs, args, log)

    # ---- coverage ----------------------------------------------------
    if args.coverage:
        windows = parse_coverage_bed(args.coverage, log)
        model.inputs["coverage"] = args.coverage
        model.coverage, model.coverage_median, model.coverage_anomalies = summarise_coverage(
            windows, seqs, args, log
        )

    # ---- annotations -------------------------------------------------
    if args.annotation:
        model.annotations = parse_annotation(args.annotation, log)
        model.inputs["annotation"] = args.annotation
        log.info(f"loaded {len(model.annotations)} annotation features")

    # ---- graph segments: classify, colour, draw, and pick out what to identify
    extra: Dict[str, str] = {}
    if segs:
        base = os.path.join(args.out_dir, args.prefix)
        os.makedirs(args.out_dir, exist_ok=True)
        seq_by_segment = read_segment_sequences(args.gfa, args.segment_fasta, log)
        contigs = parse_flye_info(args.flye_info, log) if args.flye_info else []
        calls, baseline = call_segments(segs, links, contigs, seq_by_segment, args, log)
        model.segment_calls = calls
        model.baseline_depth = baseline
        model.baseline_basis = (
            f"median depth of segments at least {human_bp(args.baseline_min_length)} long"
        )
        colours = assign_segment_colours(calls, links)
        model.segment_colours = colours

        # the same colour a segment has in the graph figure, drawn where it lands
        # on the assembled chromosome
        by_seq: Dict[str, List[Tuple[int, int, str, str]]] = defaultdict(list)
        for seg_name, places in placements.items():
            for p in places:
                by_seq[p.seqname].append(
                    (p.start, p.end, seg_name, colours.get(seg_name, "#9aa0a6"))
                )
        for rec in seqs:
            if rec.role in ("chromosome", "mitochondrion", "plastid") and by_seq.get(rec.name):
                rec.blocks = sorted(by_seq[rec.name])
                rec.blocks_tile = False

        adj = build_adjacency(links)
        extra.update(
            identify_repeats(
                model, calls, model.tangles, adj, seq_by_segment, base, args, log
            )
        )
        model.title = args.title or f"{args.prefix} - assembly ideogram"
        extra.update(write_figures(model, calls, links, colours, base, args, log))

    model.title = args.title or f"{args.prefix} - assembly ideogram"
    return model, extra


def build_model_graph_first(args, log: Log) -> Tuple[Model, Dict[str, str]]:
    """
    The path taken when there is no chromosome-level assembly: work from the
    graph the assembler emitted, and produce ranked hypotheses.
    """
    segs, links = parse_gfa(args.gfa, log, args.assembler)
    log.info(f"graph: {len(segs)} segments, {len(links)} links")

    seq_by_segment = read_segment_sequences(args.gfa, args.segment_fasta, log)
    contigs = parse_flye_info(args.flye_info, log) if args.flye_info else []
    if not contigs and args.flye_info is None:
        log.info(
            "no assembly_info.txt supplied; path-position evidence (subtelomeric versus "
            "interior) is unavailable"
        )

    calls, baseline = call_segments(segs, links, contigs, seq_by_segment, args, log)
    adj = build_adjacency(links)
    end_links = build_end_links(links)
    joins = find_joins(calls, end_links, args.max_join_hops)
    log.info(f"{len(joins)} candidate route(s) between backbone segments")

    # Segments with links on one end only cannot be traversed, so they are tips,
    # not bridges. Report them as what they are: evidence about chromosome ENDS.
    backbone_names = {c.name for c in calls if c.cls == "backbone"}
    # only for the small stuff: a backbone contig with links on one end only is
    # simply a molecule with a free end, which is what a chromosome end looks
    # like. It is not a repeat and calling it a tip would be nonsense.
    tips = dead_end_repeats(
        end_links, [c.name for c in calls if c.name not in backbone_names]
    )
    for tip, live_end in sorted(tips.items()):
        touching = sorted({n for n, _e in end_links.get((tip, live_end), ())})
        if touching:
            log.info(
                f"{tip}: links on one end only ({len(touching)} neighbour(s): "
                f"{', '.join(touching)}); a tip, not a bridge, so it caps those ends "
                f"rather than joining them"
            )

    # Shared tips. Two backbone segments attached to the SAME END of the same
    # intermediate cannot both be joined through it - you would have to leave
    # that intermediate by an end you also arrived at. But this is exactly what a
    # repeat sitting at the end of two different chromosomes looks like, so each
    # pairing is offered as a declared, penalised alternative rather than dropped.
    # The test is per END, not per segment: edge_3 has a link on its far side and
    # so is not a dead end, yet edge_5 and edge_6 both hang off edge_3.R.
    seen_pairs: Set[Tuple[str, str, str]] = set()
    for (mid, mid_end), neighbours in sorted(end_links.items()):
        if mid in backbone_names:
            continue
        at_end = sorted({n for n, _e in neighbours if n in backbone_names})
        if len(at_end) < 2:
            continue
        log.warn(
            f"{', '.join(at_end)} all attach to the same end of {mid}. The graph does NOT "
            f"resolve whether any of them join through it, because a route would have to "
            f"leave {mid} by the end it arrived at. Treat such a join as an untested "
            f"alternative, not a result."
        )
        ends_at_mid = {n: e for n, e in neighbours if n in backbone_names}
        for i, n1 in enumerate(at_end):
            for n2 in at_end[i + 1:]:
                sig = (min(n1, n2), max(n1, n2), mid)
                if sig in seen_pairs:
                    continue
                seen_pairs.add(sig)
                joins.append(
                    Join(
                        a=n1, b=n2, via=[mid],
                        a_end=ends_at_mid.get(n1, "e"),
                        b_end=ends_at_mid.get(n2, "s"),
                        speculative=True,
                    )
                )

    hypotheses = enumerate_hypotheses(calls, joins, adj, args, log, build_end_adjacency(links))

    extra: Dict[str, str] = {}
    base = os.path.join(args.out_dir, args.prefix)
    os.makedirs(args.out_dir, exist_ok=True)

    if not hypotheses:
        log.warn(
            "no backbone segments were found, so no chromosome hypotheses could be built. "
            "Lower --backbone-min-length, or check that the GFA carries depth tags."
        )
        model = Model()
        model.sequences = [
            SeqRecord(name=c.name, length=c.length, role="unplaced", depth=c.depth, manual=True)
            for c in calls
        ]
        model.title = args.title or "Assembly graph (no backbone segments found)"
        pick = 1
    else:
        pick = min(max(args.hypothesis, 1), len(hypotheses))
        model = model_from_hypothesis(hypotheses[pick - 1], calls, links, args, log)
        model.chosen_hypothesis = pick

    model.segment_calls = calls
    model.hypotheses = hypotheses
    model.range_slack = args.speculative_penalty
    model.baseline_depth = baseline
    model.baseline_basis = (
        f"median depth of segments at least {human_bp(args.baseline_min_length)} long"
    )
    model.inputs["gfa"] = args.gfa
    if args.flye_info:
        model.inputs["flye_info"] = args.flye_info

    colours = model.segment_colours or assign_segment_colours(calls, links)
    model.segment_colours = colours
    extra.update(
        identify_repeats(model, calls, model.tangles, adj, seq_by_segment, base, args, log)
    )
    extra.update(write_figures(model, calls, links, colours, base, args, log))

    # Alternative hypotheses, each as its own figure. The point is honesty about
    # ambiguity: where the top two score within a whisker of each other, showing
    # one drawing and calling it the answer is the single easiest way for this
    # tool to mislead. Drawn from the SAME calls and the SAME colours, so a
    # contig is the same colour in every figure and the reader is comparing
    # topologies rather than re-learning the palette each time.
    n_draw = max(1, min(int(getattr(args, "draw_hypotheses", 1) or 1), 5))
    if n_draw > 1 and hypotheses:
        # N is the TOTAL number of figures, so N=3 means the top one plus two
        # alternatives - not the top one plus three.
        for rank in range(pick + 1, min(pick + n_draw - 1, len(hypotheses)) + 1):
            h = hypotheses[rank - 1]
            alt = model_from_hypothesis(h, calls, links, args, log)
            alt.chosen_hypothesis = rank
            alt.segment_calls = calls
            alt.hypotheses = hypotheses
            alt.range_slack = args.speculative_penalty
            alt.baseline_depth = baseline
            alt.segment_colours = colours
            alt.expected_chromosomes = getattr(model, "expected_chromosomes", None)
            alt_base = f"{base}_h{rank}"
            got = write_figures(alt, calls, links, colours, alt_base, args, log)
            for k, v in got.items():
                if "PAIRED" in k:
                    extra[f"hypothesis {rank} figure (score {h.score:.2f})"] = v
        log.info(
            f"drew {min(n_draw, len(hypotheses))} hypotheses; compare them before treating "
            f"any one as the answer"
        )
    return model, extra


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(
        prog="detangler",
        description=(
            "Infer a karyotype from assembly outputs and draw it as an annotated ideogram, "
            "with assembly-graph tangles mapped onto chromosome coordinates."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Typical run:\n"
            "  detangler.py --fai assembly.fa.fai --gfa assembly.gfa --paf segs.paf \\\n"
            "                --coverage cov.regions.bed.gz --out-dir out --prefix assembly\n\n"
            "Then edit out/assembly_karyotype.yaml and re-run with --config to lock in your calls."
        ),
    )
    g = p.add_argument_group("sequence input (one required)")
    g.add_argument("--fasta", help="assembly FASTA (also yields GC and N content)")
    g.add_argument("--fai", help="samtools faidx index of the assembly")
    g.add_argument("--assembly-report", help="NCBI *_assembly_report.txt")

    g = p.add_argument_group("assembly graph")
    g.add_argument("--gfa", help="GFA v1 assembly graph (the file Bandage opens)")
    g.add_argument("--flye-info", help="Flye assembly_info.txt (adds graph_path evidence)")
    g.add_argument("--segment-fasta", help="segment sequences, if the GFA stores '*'")
    g.add_argument("--agp", help="AGP placing contigs into scaffolds/chromosomes")
    g.add_argument("--paf", help="PAF of graph segments (query) vs assembly (target)")
    g.add_argument("--segment-map", help="TSV: segment, seqname, start, end[, orient]")

    g = p.add_argument_group("optional tracks")
    g.add_argument("--coverage", help="BED4 of per-window depth, e.g. mosdepth *.regions.bed.gz")
    g.add_argument("--annotation", help="BED or GFF3 of features to draw on the chromosomes")

    g = p.add_argument_group("output")
    g.add_argument("--out-dir", default=".", help="directory for outputs (default: .)")
    g.add_argument("--prefix", default="detangler", help="prefix for output filenames")
    g.add_argument("--title", help="figure title")
    g.add_argument("--config", help="karyotype config to apply (overrides all inference)")
    g.add_argument("--demo", metavar="DIR", help="write synthetic input data to DIR and use it")
    g.add_argument("--no-coverage-track", action="store_true", help="do not draw the coverage track")
    g.add_argument("--order", choices=["length", "name", "file"], default="length",
                   help="chromosome ordering (default: length)")
    g.add_argument("--quiet", action="store_true")

    g = p.add_argument_group("thresholds")
    g.add_argument("--min-chrom-length", type=int,
                   help="sequences at least this long are chromosome candidates "
                        "(default: derived from the assembly)")
    g.add_argument("--chrom-length-fraction", type=float, default=0.15,
                   help="fallback threshold as a fraction of the longest sequence (default 0.15)")
    g.add_argument("--organelle-depth-ratio", type=float, default=3.0,
                   help="depth relative to nuclear baseline that suggests an organelle (default 3)")
    g.add_argument("--collapse-depth-ratio", type=float, default=2.5,
                   help="segment depth ratio that flags a collapsed repeat (default 2.5)")
    g.add_argument("--min-segment-length", type=int, default=1000,
                   help="ignore graph segments shorter than this (default 1000)")
    g.add_argument("--max-bubble-length", type=int, default=100_000,
                   help="longest segment considered part of a bubble (default 100000)")
    g.add_argument("--min-paf-block", type=int, default=2000,
                   help="minimum PAF alignment block length to keep (default 2000)")
    g.add_argument("--coverage-high-ratio", type=float, default=1.8)
    g.add_argument("--coverage-low-ratio", type=float, default=0.5)
    g.add_argument("--min-anomaly-length", type=int, default=50_000)

    g = p.add_argument_group(
        "graph-first mode",
        "Used when there is no chromosome-level assembly: classify graph segments by copy "
        "number and rank hypotheses about which of them form chromosomes.",
    )
    g.add_argument("--expected-genome-size", type=parse_size,
                   help="OPTIONAL. Expected nuclear genome size, e.g. 36.5m. Used only as a "
                        "sanity check; leave it out and the tool works from the data alone.")
    g.add_argument("--expected-chromosomes", type=int,
                   help="OPTIONAL and not needed. The chromosome count is inferred from "
                        "telomere-capped ends. Supply this only to see how a hypothesis scores "
                        "against a karyotype you already trust.")
    g.add_argument("--telomere-bonus", type=float, default=1.2,
                   help="score added per telomere-capped molecule end (default 1.2)")
    g.add_argument("--open-end-penalty", type=float, default=0.0,
                   help="score subtracted per uncapped molecule end. Zero by default: an open "
                        "end means unfinished, not 'join something to it'.")
    g.add_argument("--join-cost", type=float, default=0.6,
                   help="score a join must earn back before it is asserted (default 0.6). "
                        "Raise it to make the tool more conservative about merging contigs.")
    g.add_argument("--baseline-min-length", type=parse_size, default=1_000_000,
                   help="segments at least this long set the single-copy depth baseline "
                        "(default 1m)")
    g.add_argument("--backbone-min-length", type=parse_size, default=500_000,
                   help="shortest segment that may anchor a chromosome (default 500k)")
    g.add_argument("--backbone-copy-range", type=float, nargs=2, default=(0.6, 1.7),
                   metavar=("LOW", "HIGH"), help="copy-number window for single copy")
    g.add_argument("--repeat-min-copy", type=float, default=1.7,
                   help="copy number at or above this is a repeat. The default sits below 2 "
                        "because depth-derived copy number is noisy: a true two-copy repeat "
                        "commonly estimates at 1.8-1.9 (default 1.7)")
    g.add_argument("--tandem-min-copy", type=float, default=20.0)
    g.add_argument("--organelle-min-copy", type=float, default=3.0)
    g.add_argument("--organelle-min-length", type=parse_size, default=15_000)
    g.add_argument("--organelle-max-length", type=parse_size, default=200_000)
    g.add_argument("--low-coverage-max-copy", type=float, default=0.6,
                   help="copy number below this is its own class, not a low-confidence single "
                        "copy (default 0.6)")
    g.add_argument("--foreign-gc-delta", type=float, default=0.05,
                   help="GC this far from the backbone marks an unassigned sequence as a "
                        "contaminant candidate (default 0.05)")
    g.add_argument("--at-rich-delta", type=float, default=0.08,
                   help="GC this far below the genome baseline counts as AT-rich (default 0.08)")
    g.add_argument("--telomere-motif", action="append", default=None, metavar="MOTIF",
                   help="telomere repeat motif to count; repeatable. Defaults to the canonical "
                        "motifs for vertebrates, plants, insects, nematodes and ciliates.")
    g.add_argument("--telomere-window", type=parse_size, default=20_000,
                   help="how far in from each end of a segment to look for a telomere repeat "
                        "array. A telomere that is not at an end is not a telomere "
                        "(default 20k)")
    g.add_argument("--min-telomere-units", type=int, default=3,
                   help="consecutive perfect repeats needed to call a telomere array. Three "
                        "units is already ~4**-18 by chance; the reported unit count is what "
                        "tells you how strong a call is (default 3)")
    g.add_argument("--max-period", type=int, default=3000,
                   help="largest tandem repeat unit to screen for (default 3000)")
    g.add_argument("--max-period-scan-length", type=parse_size, default=200_000)
    g.add_argument("--max-join-hops", type=int, default=2,
                   help="intermediate segments allowed between two backbone segments (default 2)")
    g.add_argument("--max-join-edges", type=int, default=18,
                   help="cap on candidate joins before enumeration (default 18)")
    g.add_argument("--max-hypotheses", type=int, default=10)
    g.add_argument("--tie-threshold", type=float, default=0.75,
                   help="hypotheses within this score of the best are reported as tied "
                        "(default 0.75)")
    g.add_argument("--assembler",
                   choices=["flye", "hifiasm", "verkko", "spades", "miniasm", "canu",
                            "unknown"],
                   default="unknown",
                   help="which assembler wrote the GFA. Assemblers disagree on things the "
                        "inference depends on: which tag carries depth and what it means "
                        "(hifiasm's rd:i:n is coverage n+1), and whether L-line CIGARs are 0M "
                        "or real overlaps that must be subtracted from a chain's length. "
                        "Overlaps are detected from the file either way; this flag mainly "
                        "fixes the depth reading and lets the report state what it assumed")
    g.add_argument("--assembly-type", choices=["primary", "phased", "collapsed"],
                   default="primary",
                   help="what the ASSEMBLER produced, which is what decides how depth should "
                        "be read - not the organism's ploidy. 'primary' (default): one "
                        "haplotype per locus, alternates written to a separate file, so a "
                        "segment at about half baseline is an unpurged haplotig. 'phased': "
                        "both haplotypes in this graph, so expect roughly twice the "
                        "chromosome count. 'collapsed': haplotypes merged, so the long "
                        "contigs sit at twice the haploid depth and every copy number below "
                        "is understated by a factor of two")
    g.add_argument("--haplotig-band", type=float, nargs=2, default=(0.35, 0.65),
                   metavar=("LOW", "HIGH"),
                   help="copy-number range that reads as a haplotig rather than as foreign "
                        "sequence (default 0.35 0.65)")
    g.add_argument("--haplotig-min-length", type=parse_size, default=20_000,
                   help="shorter than this, a half-depth segment is not called a haplotig - "
                        "there is not enough of it to tell (default 20k)")
    g.add_argument("--placement-tolerance", type=float, default=0.5,
                   help="how far a segment's depth-derived copy number may fall short "
                        "of the number of places it is drawn before that is reported as a "
                        "contradiction (default 0.5)")
    g.add_argument("--overcopy-factor", type=float, default=2.5,
                   help="how many times the number of places a segment is drawn its "
                        "depth-derived copy number may reach before that is reported as "
                        "a contradiction. This is the other direction of the same check: "
                        "a segment drawn once at 65 copies is a tandem array or a "
                        "collapsed repeat, not a single copy (default 2.5)")
    g.add_argument("--graph-triangle", dest="graph_triangle", action="store_true",
                   help="push the assembly graph out of the lower-right triangle of its panel. "
                        "Off by default: it existed so a chromosome row drawn BELOW the graph "
                        "could interlock with it, and the panels are now side by side, where "
                        "they cannot collide and the graph should use its whole rectangle")
    g.set_defaults(graph_triangle=False)
    g.add_argument("--direct-link-bonus", type=float, default=0.55,
                   help="score for a join made by a DIRECT link between two backbone segments, "
                        "over and above cancelling --join-cost. A direct link is the least "
                        "ambiguous evidence a graph carries - the assembler saw these two ends "
                        "adjoin, with no intermediate to misread - so it should be able to "
                        "carry a join on its own (default 0.55)")
    g.add_argument("--centromere-bonus", type=float, default=1.2,
                   help="score added when a join runs through a long, markedly AT-rich, "
                        "low-copy segment that bridges exactly two backbone ends. In lineages "
                        "with AT-rich regional centromeres (many filamentous fungi) that is "
                        "what a centromere between two chromosome arms looks like. Set to 0 to "
                        "switch the rule off (default 1.2)")
    g.add_argument("--centromere-min-length", type=parse_size, default=10_000,
                   help="shortest AT-rich bridge that can be read as centromeric; below this "
                        "there is not enough sequence to mean anything (default 10k)")
    g.add_argument("--centromere-speculative-discount", type=float, default=0.4,
                   help="fraction of --speculative-penalty applied to a join through an "
                        "AT-rich centromere candidate. An assembler is EXPECTED to fail to "
                        "read through such a block, so the missing through-path is explained "
                        "rather than damning - discounted, not waived (default 0.4)")
    g.add_argument("--speculative-penalty", type=float, default=1.5,
                   help="score penalty per join that is not supported by a traversable "
                        "path, i.e. where two segments merely end in the same one-sided "
                        "repeat. Such joins are reported, never silently used (default 1.5)")
    g.add_argument("--hypothesis", type=int, default=1,
                   help="which ranked hypothesis to draw as the ideogram (default 1)")
    g.add_argument("--draw-hypotheses", type=int, default=1, metavar="N",
                   help="draw the top N ranked hypotheses instead of only one (default 1, "
                        "max 5). Extra figures are written as <prefix>_paired_h2.svg and so "
                        "on, each headed with its rank and score. Worth using whenever the "
                        "report says the top hypotheses score within --tie-threshold of each "
                        "other: at that point the graph does not choose between them and one "
                        "picture on its own overstates the case")

    g = p.add_argument_group(
        "repeat identification",
        "Candidate sequences and ready-to-run commands are always written. BLAST itself only "
        "runs if you point it at a database, a subject FASTA, or NCBI.",
    )
    g.add_argument("--blast-db", help="local BLAST nucleotide database to search")
    g.add_argument("--blast-subject", help="FASTA to search against directly (no makeblastdb)")
    g.add_argument("--blast-remote", action="store_true",
                   help="search NCBI over the network (slow, rate-limited)")
    g.add_argument("--blast-remote-db", default="nt")
    g.add_argument("--blast-hits", help="existing outfmt-6 results to read instead of searching")
    g.add_argument("--blast-task", default="megablast",
                   choices=["megablast", "dc-megablast", "blastn", "blastn-short"])
    g.add_argument("--blast-evalue", default="1e-20")
    g.add_argument("--blast-max-target-seqs", type=int, default=5)
    g.add_argument("--blast-threads", type=int, default=4)
    g.add_argument("--blast-timeout", type=int, default=3600)
    g.add_argument("--blast-min-identity", type=float, default=85.0)
    g.add_argument("--blast-min-coverage", type=float, default=50.0)
    g.add_argument("--blast-report-hits", type=int, default=3)
    g.add_argument("--branch-degree", type=int, default=3,
                   help="a segment linking this many others counts as a branch point worth "
                        "identifying (default 3)")
    g.add_argument("--max-candidate-length", type=parse_size, default=1_000_000,
                   help="do not export segments longer than this: a multi-megabase contig is a "
                        "chromosome arm, not a repeat to look up (default 1m)")
    g.add_argument("--max-leader-lines", type=int, default=40,
                   help="leader lines drawn between the two panels of the paired figure "
                        "(default 40); colour still identifies segments beyond that")
    g.add_argument("--graph-style", choices=["bandage", "layered"], default="bandage",
                   help="'bandage' draws thick tapered segments in a force-directed layout, "
                        "like Bandage; 'layered' uses boxes in BFS layers (default: bandage)")
    g.add_argument("--graph-length-scale", type=float, default=0.10,
                   help="drawn pixels per sqrt(bp) for a segment (default 0.10)")
    g.add_argument("--graph-max-segment-px", type=float, default=300.0,
                   help="longest a segment may be drawn (default 300)")
    g.add_argument("--graph-label-limit", type=int, default=40,
                   help="label every segment up to this many; above it, only the notable ones")
    g.add_argument("--bandage-image", metavar="FILE",
                   help="a real Bandage export (PNG, JPEG or SVG) to use as the left panel of "
                        "the paired figure instead of our redraw. Load the emitted "
                        "*_bandage_colours.csv in Bandage first so the two panels agree on "
                        "colour.")
    g.add_argument("--bandage-max-width", type=float, default=1400.0,
                   help="cap on the width of an embedded Bandage image (default 1400)")
    g.add_argument("--rotate-graph", action=argparse.BooleanOptionalAction, default=False,
                   help="turn our graph redraw a quarter turn so the paired figure is narrower "
                        "and taller (default: on). Labels then read bottom-to-top.")
    g.add_argument("--max-graph-nodes", type=int, default=300,
                   help="skip the graph figure above this many segments (default 300)")
    g.add_argument("--blast-classes", nargs="+",
                   default=["repeat", "tandem_array", "organelle_candidate", "at_rich",
                            "low_coverage"],
                   help="segment classes to export and search")

    args = p.parse_args(argv)
    if args.telomere_motif is None:
        args.telomere_motif = list(TELOMERE_MOTIFS)
    log = Log(args.quiet)

    if args.demo:
        demo = write_demo(args.demo, log)
        args.fai = args.fai or demo["fai"]
        args.gfa = args.gfa or demo["gfa"]
        args.paf = args.paf or demo["paf"]
        args.coverage = args.coverage or demo["coverage"]
        args.annotation = args.annotation or demo["annotation"]
        args.title = args.title or "Demo assembly - synthetic data, not a real genome"

    cfg = None
    if args.config:
        cfg = load_config(args.config)
        for key in ("fasta", "fai", "assembly_report", "gfa", "agp", "paf",
                    "segment_map", "coverage", "annotation"):
            if not getattr(args, key, None) and (cfg.get("inputs") or {}).get(key):
                path = cfg["inputs"][key]
                if os.path.exists(path):
                    setattr(args, key, path)
                else:
                    log.warn(f"config references {key}={path}, which no longer exists")

    extra_files: Dict[str, str] = {}
    have_seqs = any([args.fasta, args.fai, args.assembly_report])

    if not have_seqs and args.gfa:
        # graph-first: no chromosome-level assembly, so hypothesise one
        log.info("no chromosome-level sequences supplied; running in graph-first mode")
        model, extra_files = build_model_graph_first(args, log)
    elif not have_seqs and cfg and cfg.get("sequences"):
        model = Model()
        model.settings = {"order": args.order, "coverage": not args.no_coverage_track}
        model.title = args.title or "Assembly ideogram (from config)"
        log.info("no input files available; rendering purely from the config file")
    elif not have_seqs:
        p.error(
            "no input. Pass --fasta/--fai/--assembly-report for an ideogram of an assembled "
            "genome, --gfa alone for graph-first hypotheses, or --config / --demo."
        )
    else:
        model, extra_files = build_model(args, log)

    if cfg:
        apply_overrides(model, cfg, log)

    if not model.sequences:
        raise SystemExit("no sequences to draw")

    os.makedirs(args.out_dir, exist_ok=True)
    base = os.path.join(args.out_dir, args.prefix)
    # keep the model's own warnings - the copy-number-versus-placement checks are
    # raised while the model is built, and assigning here used to wipe them
    model.warnings = list(log.warnings) + [
        w for w in model.warnings if w not in log.warnings
    ]

    svg_path = base + "_ideogram.svg"
    with open(svg_path, "w") as fh:
        fh.write(render_svg(model))

    html_path = base + "_ideogram.html"
    with open(html_path, "w") as fh:
        fh.write(render_html(model))

    cfg_path = write_config(model, base + ("_karyotype.yaml" if HAVE_YAML else "_karyotype.json"))

    files = {
        "publication SVG": svg_path,
        "interactive HTML": html_path,
        "editable karyotype config": cfg_path,
    }
    files.update(extra_files)
    report_path = base + "_report.md"
    with open(report_path, "w") as fh:
        fh.write(render_report(model, files))
    files["reasoning report"] = report_path

    if not args.quiet:
        print()
        print(model.summary_sentence())
        print()
        for k, v in files.items():
            print(f"  {k:28s} {v}")
        if not HAVE_YAML:
            print("\n  (PyYAML not found, so the config was written as JSON.)")
        print()
    return 0
if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, OSError) as exc:
        # Bad input files are expected; a stack trace helps nobody.
        print(f"[detangler] ERROR: {exc}", file=sys.stderr)
        sys.exit(2)
    except KeyboardInterrupt:
        sys.exit(130)
