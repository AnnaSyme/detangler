"""Optional similarity search over candidate segments."""
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
    TANGLE_LABEL,
    _maybe_float,
    smart_open,
)
from .records import (
    Tangle,
)
from .palette import (
    CLASS_LABEL,
)
from .calls import (
    SegmentCall,
)
from .model import (
    Model,
)



# --------------------------------------------------------------------------
# identifying what a repeat actually is
# --------------------------------------------------------------------------
BLAST_FMT = "6 qseqid sseqid pident length qcovhsp evalue bitscore stitle"
BLAST_COLS = ["qseqid", "sseqid", "pident", "length", "qcovhsp", "evalue", "bitscore", "stitle"]


def select_interesting_segments(
    calls: List[SegmentCall],
    tangles: List[Tangle],
    adj: Dict[str, Set[str]],
    args,
    log: Log,
) -> Dict[str, List[str]]:
    """
    Pick the segments worth identifying, and record why each was picked.

    "Interesting" means the segment is behaving like something other than plain
    single-copy sequence: it is present in more than one copy, it branches to
    several other segments, it closes on itself, it carries telomere repeat, or
    the tangle detector already flagged it. Very long segments are excluded -
    a nine-megabase contig is a chromosome arm, not a repeat to identify.
    """
    picked: Dict[str, List[str]] = {}

    def add(name: str, why: str) -> None:
        if why not in picked.setdefault(name, []):
            picked[name].append(why)

    by_name = {c.name: c for c in calls}
    for c in calls:
        if c.length < args.min_segment_length or c.length > args.max_candidate_length:
            continue
        deg = len(adj.get(c.name, ()))
        if c.copy_number is not None and c.copy_number >= args.repeat_min_copy:
            add(c.name, f"{c.copy_number:.1f}x copy number, so present in more than one place")
        if deg >= args.branch_degree:
            add(c.name, f"links {deg} other segments, so it sits at a branch in the graph")
        if c.self_loop_same_orient:
            add(c.name, "links to itself: circular molecule or tandem array")
        if c.self_loop_flipped:
            add(c.name, "links to itself inverted: hairpin or inverted repeat")
        if c.telomere_motifs:
            add(c.name, f"{sum(c.telomere_motifs.values())} telomere motif occurrences")
        if c.at_rich:
            add(c.name, "markedly AT-rich compared with the rest of the assembly")
        if c.cls in args.blast_classes:
            add(c.name, f"classified as {CLASS_LABEL.get(c.cls, c.cls)}")

    for t in tangles:
        for seg in t.segments:
            c = by_name.get(seg)
            if c and args.min_segment_length <= c.length <= args.max_candidate_length:
                add(seg, f"flagged as {TANGLE_LABEL.get(t.type, t.type).lower()}")

    if picked:
        log.info(f"{len(picked)} segment(s) selected for identification")
    return picked


def export_repeat_candidates(
    picked: Dict[str, List[str]],
    calls: List[SegmentCall],
    seq_by_segment: Dict[str, str],
    path: str,
    log: Log,
) -> Tuple[str, int]:
    """
    Write the selected segments to FASTA. The defline carries the observations
    and the reason it was picked, so the file explains itself if it is handed to
    a colleague or queued on a cluster.
    """
    by_name = {c.name: c for c in calls}
    n = missing = 0
    with open(path, "w") as fh:
        for name in sorted(picked, key=lambda k: -(by_name[k].length if k in by_name else 0)):
            c = by_name.get(name)
            seq = seq_by_segment.get(name)
            if c is None:
                continue
            if not seq:
                missing += 1
                continue
            cn = f"{c.copy_number:.1f}" if c.copy_number is not None else "NA"
            gc = f"{c.gc:.3f}" if c.gc is not None else "NA"
            why = "; ".join(picked[name])
            fh.write(
                f">{c.name} length={c.length} "
                f"depth={c.depth if c.depth is not None else 'NA'} copy_number={cn} gc={gc} "
                f"class={c.cls} why=\"{why}\"\n"
            )
            for i in range(0, len(seq), 60):
                fh.write(seq[i : i + 60] + "\n")
            n += 1
    if missing:
        log.warn(
            f"{missing} selected segment(s) have no sequence in the GFA (it stores '*'), so they "
            f"cannot be identified by similarity search. Re-export the graph with sequences, or "
            f"pass --segment-fasta."
        )
    if n == 0:
        log.warn("no candidate sequences could be written, so there is nothing to BLAST")
    else:
        log.info(f"wrote {n} candidate sequence(s) to {path}")
    return path, n


def attach_hits_to_tangles(calls: List[SegmentCall], tangles: List[Tangle], args) -> int:
    """Carry BLAST results through to the figure, so a repeat can be named on it."""
    best: Dict[str, Dict] = {}
    for c in calls:
        for h in c.identity_hits:
            if (h.get("pident") or 0) >= args.blast_min_identity and (
                h.get("qcovhsp") or 0
            ) >= args.blast_min_coverage:
                best[c.name] = h
                break
    n = 0
    for t in tangles:
        hits = [(s, best[s]) for s in t.segments if s in best]
        if not hits:
            continue
        seg, h = hits[0]
        title = (h.get("stitle") or h.get("sseqid") or "").strip()
        t.description += (
            f" Similarity search: {seg} matches \"{title[:110]}\" "
            f"({h['pident']:.1f}% identity over {h['qcovhsp']:.0f}% of the segment). "
            f"That is what it resembles, not a confirmed identity."
        )
        t.evidence.append(f"top BLAST hit for {seg}: {title[:70]}")
        n += 1
    return n


def identify_repeats(
    model: Model,
    calls: List[SegmentCall],
    tangles: List[Tangle],
    adj: Dict[str, Set[str]],
    seq_by_segment: Dict[str, str],
    base: str,
    args,
    log: Log,
) -> Dict[str, str]:
    """Selection, export, optional search, and feeding results back in."""
    extra: Dict[str, str] = {}
    picked = select_interesting_segments(calls, tangles, adj, args, log)
    model.candidate_reasons = picked
    if not picked:
        return extra

    fasta = base + "_repeat_candidates.fasta"
    _, n = export_repeat_candidates(picked, calls, seq_by_segment, fasta, log)
    if n:
        extra["candidate sequences to identify"] = fasta
        extra["ready-to-run BLAST commands"] = write_blast_commands(fasta, base, args)

    hits = args.blast_hits
    if not hits and n and (args.blast_db or args.blast_subject or args.blast_remote):
        hits = run_blast(fasta, base + "_repeat_hits.tsv", args, log)
        if hits:
            extra["similarity search hits"] = hits
    if hits and os.path.exists(hits):
        attach_blast_hits(calls, hits, args, log)
        named = attach_hits_to_tangles(calls, tangles, args)
        if named:
            log.info(f"{named} graph feature(s) now carry a similarity hit in the figure")
    return extra


def write_blast_commands(fasta: str, out_prefix: str, args) -> str:
    """
    Always emitted, whether or not BLAST is run here, so the search can be
    reproduced or moved to a cluster. These are the exact commands used.
    """
    path = out_prefix + "_blast_commands.sh"
    tsv = out_prefix + "_repeat_hits.tsv"
    with open(path, "w") as fh:
        fh.write("#!/bin/sh\n")
        fh.write("# Commands for identifying the repeat candidates exported by detangler.\n")
        fh.write("# Pick whichever line matches the resources you have. Requires BLAST+.\n\n")
        fh.write("# 1. Against a local nucleotide database (fastest, needs a formatted DB):\n")
        fh.write(
            f"blastn -task {args.blast_task} -query {fasta} -db nt "
            f"-outfmt '{BLAST_FMT}' -evalue {args.blast_evalue} "
            f"-max_target_seqs {args.blast_max_target_seqs} -num_threads {args.blast_threads} "
            f"> {tsv}\n\n"
        )
        fh.write("# 2. Against a repeat library or any FASTA, without formatting a database:\n")
        fh.write(
            f"blastn -task {args.blast_task} -query {fasta} -subject REPEAT_LIBRARY.fa "
            f"-outfmt '{BLAST_FMT}' -evalue {args.blast_evalue} > {tsv}\n\n"
        )
        fh.write("# 3. Against NCBI over the network (no local DB; slow and rate-limited,\n")
        fh.write("#    and -remote cannot be combined with -num_threads):\n")
        fh.write(
            f"blastn -task {args.blast_task} -query {fasta} -db nt -remote "
            f"-outfmt '{BLAST_FMT}' -evalue {args.blast_evalue} "
            f"-max_target_seqs {args.blast_max_target_seqs} > {tsv}\n\n"
        )
        fh.write("# Then feed the results back in:\n")
        fh.write(f"#   detangler.py ... --blast-hits {tsv}\n")
    os.chmod(path, 0o755)
    return path


def run_blast(fasta: str, out_tsv: str, args, log: Log) -> Optional[str]:
    """Run blastn if the user asked for it. Returns the TSV path, or None."""
    import shutil
    import subprocess

    exe = shutil.which("blastn")
    if not exe:
        log.warn(
            "blastn was requested but is not on PATH. Install BLAST+ (conda install -c bioconda "
            "blast) or run the commands in the emitted *_blast_commands.sh elsewhere."
        )
        return None

    cmd = [exe, "-task", args.blast_task, "-query", fasta, "-outfmt", BLAST_FMT,
           "-evalue", str(args.blast_evalue)]
    if args.blast_subject:
        cmd += ["-subject", args.blast_subject]
    elif args.blast_db:
        cmd += ["-db", args.blast_db, "-max_target_seqs", str(args.blast_max_target_seqs)]
    elif args.blast_remote:
        cmd += ["-db", args.blast_remote_db, "-remote",
                "-max_target_seqs", str(args.blast_max_target_seqs)]
    else:
        return None
    # -remote is incompatible with -num_threads; -subject ignores it
    if not args.blast_remote and not args.blast_subject:
        cmd += ["-num_threads", str(args.blast_threads)]

    log.info("running: " + " ".join(cmd))
    if args.blast_remote:
        log.info("remote BLAST against NCBI can take many minutes; be patient or use a local DB")
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=args.blast_timeout)
    except subprocess.TimeoutExpired:
        log.warn(f"blastn timed out after {args.blast_timeout}s; no hits recorded")
        return None
    if proc.returncode != 0:
        log.warn(f"blastn exited {proc.returncode}: {proc.stderr.strip()[:400]}")
        return None
    with open(out_tsv, "w") as fh:
        fh.write(proc.stdout)
    log.info(f"blastn wrote {len(proc.stdout.splitlines())} hit line(s) to {out_tsv}")
    return out_tsv


def attach_blast_hits(calls: List[SegmentCall], tsv: str, args, log: Log) -> None:
    """
    Read outfmt-6 results and attach them to the segments. Hits are recorded
    verbatim; the tool does not rename a segment on the strength of a hit, it
    only reports what matched and how well.
    """
    by_query: Dict[str, List[Dict]] = defaultdict(list)
    with smart_open(tsv) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < len(BLAST_COLS):
                f = f + [""] * (len(BLAST_COLS) - len(f))
            rec = dict(zip(BLAST_COLS, f))
            for k in ("pident", "qcovhsp", "bitscore"):
                rec[k] = _maybe_float(rec.get(k, ""))  # type: ignore
            rec["length"] = int(rec["length"]) if str(rec.get("length", "")).isdigit() else 0
            by_query[rec["qseqid"]].append(rec)

    n_annotated = 0
    for c in calls:
        hits = by_query.get(c.name) or []
        if not hits:
            continue
        hits.sort(key=lambda h: -(h.get("bitscore") or 0))
        c.identity_hits = hits[: args.blast_report_hits]
        n_annotated += 1
        good = [
            h
            for h in c.identity_hits
            if (h.get("pident") or 0) >= args.blast_min_identity
            and (h.get("qcovhsp") or 0) >= args.blast_min_coverage
        ]
        if good:
            top = good[0]
            c.reasons.append(
                f"top similarity hit: {top['stitle'][:110] or top['sseqid']} "
                f"({top['pident']:.1f}% identity over {top['qcovhsp']:.0f}% of the segment, "
                f"E={top['evalue']}) - this is what matched, not a taxonomic assignment"
            )
        else:
            best = c.identity_hits[0]
            c.reasons.append(
                f"best similarity hit fell below the reporting thresholds "
                f"({best.get('pident') or 0:.1f}% identity over "
                f"{best.get('qcovhsp') or 0:.0f}% coverage): "
                f"{best['stitle'][:90] or best['sseqid']}"
            )
    log.info(f"attached similarity hits to {n_annotated} segment(s)")
