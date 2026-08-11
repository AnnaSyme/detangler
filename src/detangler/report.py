"""The plain-language reasoning report."""
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
    TANGLE_LABEL,
    VERSION,
    human_bp,
    nx_stat,
)
from .records import (
    Tangle,
)
from .palette import (
    CLASS_LABEL,
)
from .model import (
    Model,
    _confidence,
)



# ==========================================================================
# Markdown report
# ==========================================================================
def _report_graph_sections(model: Model) -> List[str]:
    """
    Observations, derived estimates and hypotheses, kept in separate sections so
    a reader can see at a glance which is which.
    """
    calls = model.segment_calls
    L: List[str] = []

    L.append("## Segments: observations")
    L.append("")
    L.append(
        "Everything in this table is read directly from the GFA and, where supplied, "
        "assembly_info.txt. No inference has been applied."
    )
    L.append("")
    L.append("| Segment | Length | Depth | Degree | Self-link | Component | GC | Telomere motifs | Path position |")
    L.append("|---|---:|---:|---:|---|---:|---:|---:|---|")
    for c in sorted(calls, key=lambda c: -c.length):
        # Say what the file contains, not what it means: a same-orientation
        # self-link is compatible with a circle and with a tandem array, and
        # choosing between those is an inference made in the next section.
        loop = (
            "same orient" if c.self_loop_same_orient else ("flipped" if c.self_loop_flipped else "")
        )
        telo = sum(c.telomere_motifs.values()) or ""
        posn = (
            f"{c.path_terminal} terminal / {c.path_interior} interior"
            if (c.path_terminal or c.path_interior)
            else ""
        )
        L.append(
            f"| {c.name} | {c.length:,} | "
            f"{f'{c.depth:.0f}' if c.depth is not None else ''} | {c.degree} | {loop} | "
            f"{c.component} ({c.component_size}) | "
            f"{f'{c.gc:.0%}' if c.gc is not None else ''} | {telo} | {posn} |"
        )
    L.append("")

    L.append("## Segments: derived estimates and classification")
    L.append("")
    if model.baseline_depth:
        L.append(
            f"Baseline single-copy depth is **{model.baseline_depth:.1f}x**, taken as the "
            f"{model.baseline_basis}. Every copy number below is depth divided by that baseline, "
            f"so it is an estimate, not a measurement. It assumes uniform sequencing depth, which "
            f"GC bias, ploidy variation and sex chromosomes all break."
        )
    else:
        L.append(
            "No depth information was available, so copy number could not be estimated and the "
            "classification rests on topology and length alone."
        )
    L.append("")
    L.append("| Segment | Copy number (estimate) | Class (inference) | Unit period | Why |")
    L.append("|---|---:|---|---:|---|")
    for c in sorted(calls, key=lambda c: -c.length):
        per = f"{c.period[0]:,} bp" if c.period else ""
        L.append(
            f"| {c.name} | "
            f"{f'{c.copy_number:.2f}' if c.copy_number is not None else 'n/a'} | "
            f"{CLASS_LABEL[c.cls]} | {per} | {'; '.join(c.reasons)} |"
        )
    L.append("")

    if model.candidate_reasons:
        L.append("### Segments selected for identification")
        L.append("")
        L.append(
            "These were picked out as behaving like something other than plain single-copy "
            "sequence, written to the candidate FASTA, and are what the BLAST commands search."
        )
        L.append("")
        L.append("| Segment | Length | Depth | Why it was picked |")
        L.append("|---|---:|---:|---|")
        by_name = {c.name: c for c in calls}
        for name, reasons in sorted(
            model.candidate_reasons.items(),
            key=lambda kv: -(by_name[kv[0]].length if kv[0] in by_name else 0),
        ):
            c = by_name.get(name)
            if not c:
                continue
            L.append(
                f"| {name} | {human_bp(c.length)} | "
                f"{f'{c.depth:.0f}x' if c.depth is not None else ''} | {'; '.join(reasons)} |"
            )
        L.append("")

    hits = [c for c in calls if c.identity_hits]
    if hits:
        L.append("### Similarity search results")
        L.append("")
        L.append(
            "Reported verbatim from BLAST. A hit is evidence about what a sequence resembles; it "
            "is not a taxonomic assignment, and a repeat that hits many things equally is still "
            "unidentified."
        )
        L.append("")
        L.append("| Segment | Subject | % identity | % query covered | E-value | Description |")
        L.append("|---|---|---:|---:|---:|---|")
        for c in hits:
            for h in c.identity_hits:
                L.append(
                    f"| {c.name} | {h['sseqid']} | "
                    f"{h['pident'] if h['pident'] is not None else ''} | "
                    f"{h['qcovhsp'] if h['qcovhsp'] is not None else ''} | {h['evalue']} | "
                    f"{(h['stitle'] or '')[:120]} |"
                )
        L.append("")

    ua = model.unassigned()
    if ua:
        L.append("## Not assigned to any chromosome")
        L.append("")
        L.append(
            "These sequences were not forced into a chromosome. Contamination is offered as a "
            "candidate explanation only where the sequence is disconnected from the nuclear "
            "graph AND its GC or depth differs from the backbone; it is never asserted. "
            "Confirm by similarity search before discarding anything."
        )
        L.append("")
        L.append("| Sequence | Length | Depth | Why it is here |")
        L.append("|---|---:|---:|---|")
        for s in ua:
            L.append(
                f"| {s.display} | {human_bp(s.length)} | "
                f"{f'{s.depth:.0f}x' if s.depth is not None else ''} | {s.note} |"
            )
        L.append("")

    if model.hypotheses:
        top = model.hypotheses[0]
        best, low, high = model.count_range() or (len(top.chains),) * 3
        L.append("## How many chromosomes? (inferred, not supplied)")
        L.append("")
        L.append(
            f"Best estimate: **{best} linear molecule(s)**"
            + (
                f", but the graph cannot distinguish between **{low} and {high}**."
                if low != high
                else ", and the alternatives score clearly worse."
            )
        )
        L.append("")
        # This paragraph used to claim unconditionally that no karyotype was
        # used, while the scorer had already applied a penalty proportional to
        # the distance from --expected-chromosomes. That is a false statement in
        # the output, and it is false in the direction that flatters the tool.
        exp_n = getattr(model, "expected_chromosomes", None)
        if exp_n:
            L.append(
                f"**{exp_n} chromosomes were supplied** via `--expected-chromosomes`, and the "
                f"ranking used them: every hypothesis was penalised in proportion to how far "
                f"its molecule count sits from {exp_n}. So this is not an independent estimate "
                f"of the count - it is the graph's evidence weighted towards a number you "
                f"provided. Run without the flag to see what the graph says on its own."
            )
        else:
            L.append(
                "No expected karyotype was used to reach that. A finished linear chromosome "
                "carries a telomere repeat array at each end, so the count follows from how "
                "many ends are capped; a join between two contigs is only asserted when the "
                "segment bridging them is present in roughly one copy and touches nothing else."
            )
        L.append("")
        L.append(f"- {top.capped_ends} of {2 * len(top.chains)} ends are telomere-capped")
        L.append(f"- {top.open_ends} end(s) are open, so those molecules are unfinished")
        if low != high:
            L.append(
                f"- the range comes from joins the graph permits but does not require; each one "
                f"merges two molecules into one, which is why the count is a range and not a "
                f"number"
            )
            L.append(
                "- to close it: Hi-C contact data, a reference alignment, or reads long enough "
                "to span the bridging segments listed under each hypothesis"
            )
        if top.open_ends == 0 and low == high:
            L.append("- every molecule is closed end to end, so the count is well supported")
        organelles = [s for s in model.sequences if s.role in ("mitochondrion", "plastid")]
        for o in organelles:
            L.append(
                f"- plus {o.display}, {human_bp(o.length)}"
                + (", circular" if o.circular else "")
                + ", counted separately from the linear set"
            )
        L.append("")

        L.append("## Chromosome hypotheses (ranked)")
        L.append("")
        L.append(
            "These are hypotheses, not assignments. Topology alone is ambiguous wherever a repeat "
            "joins more than two backbone segments, so the alternatives are listed rather than "
            "one being chosen silently. Which chromosome is which is not addressed at all: size "
            "ordering is a hint, not evidence."
        )
        L.append("")
        lengths = {c.name: c.length for c in calls}
        for h in model.hypotheses:
            marker = "  <- drawn in the ideogram" if h.rank == model.chosen_hypothesis else ""
            L.append(f"### Hypothesis {h.rank} (score {h.score:.2f}){marker}")
            L.append("")
            L.append(f"{len(h.chains)} molecule(s) from the backbone segments:")
            L.append("")
            for i, chain in enumerate(h.chains, 1):
                L.append(
                    f"- chain {i}: {' + '.join(chain)} = "
                    f"{human_bp(sum(lengths.get(s, 0) for s in chain))}"
                )
            L.append("")
            if h.supporting:
                L.append("Supporting:")
                L.append("")
                for s in h.supporting:
                    L.append(f"- {s}")
                L.append("")
            if h.contradicting:
                L.append("Contradicting or unresolved:")
                L.append("")
                for s in h.contradicting:
                    L.append(f"- {s}")
                L.append("")
            if h.resolve_with:
                L.append("What would resolve this:")
                L.append("")
                for s in h.resolve_with:
                    L.append(f"- {s}")
                L.append("")
    return L


CAVEATS = """\
- Every call below is a heuristic over sequence length, name, graph topology and read
  depth. None of it is a substitute for an organelle-aware annotation tool or for
  aligning to a reference.
- A high-depth circular sequence of mitogenome size can also be a NUMT-rich contig, a
  plasmid, or a bacterial contaminant. Confirm organelle calls by annotating the
  expected gene set.
- A repeat shared between chromosomes in the graph may equally be an assembly artefact
  (chimeric join) or a real shared repeat family. The graph alone cannot distinguish
  these; orthogonal evidence such as Hi-C contact maps, long-read spanning, or an
  optical map is required.
- Depth-derived copy number assumes uniform sequencing depth. GC bias, ploidy variation
  and sex chromosomes all break that assumption.
- Segment placements taken from PAF inherit the aligner's mapping decisions; a repeat
  can be placed wherever the aligner chose to put it.
"""


def render_report(model: Model, files: Dict[str, str]) -> str:
    L: List[str] = []
    L.append(f"# {model.title}")
    L.append("")
    L.append(f"_Generated by detangler v{VERSION}._")
    L.append("")
    L.append("## Summary")
    L.append("")
    L.append(model.summary_sentence())
    L.append("")
    total = sum(s.length for s in model.sequences)
    lengths = [s.length for s in model.sequences]
    L.append(f"- Total assembly span: **{human_bp(total)}** across {len(model.sequences)} sequences")
    L.append(f"- N50: **{human_bp(nx_stat(lengths, 50))}**, longest {human_bp(max(lengths or [0]))}")
    placed = sum(s.length for s in model.drawable())
    L.append(
        f"- Anchored to drawn molecules: **{human_bp(placed)}** "
        f"({100.0 * placed / max(total, 1):.1f}% of the assembly)"
    )
    L.append("")

    if model.inputs:
        L.append("## Inputs")
        L.append("")
        for k, v in model.inputs.items():
            L.append(f"- `{k}`: {v}")
        L.append("")

    L.append("## Karyotype calls")
    L.append("")
    L.append("| Sequence | Length | Call | Confidence | Circular | Depth |")
    L.append("|---|---:|---|---|---|---:|")
    for s in model.drawable() + model.unplaced():
        L.append(
            f"| {s.display} | {human_bp(s.length)} | {s.role} | {_confidence(s)} | "
            f"{'yes' if s.circular else ''} | "
            f"{f'{s.depth:.1f}x' if s.depth is not None else ''} |"
        )
    L.append("")

    L.append("### Evidence behind each call")
    L.append("")
    for s in model.drawable():
        L.append(f"**{s.display}** -> `{s.role}`")
        L.append("")
        if s.evidence:
            for e in s.evidence:
                L.append(f"- {e.as_text()}")
        else:
            L.append("- no evidence recorded (role asserted)")
        L.append("")

    if model.segment_calls:
        L.extend(_report_graph_sections(model))

    L.append("## Assembly graph features")
    L.append("")
    if not model.tangles:
        L.append("No tangles were detected. Either no GFA was supplied, or the graph is linear "
                 "and unbranched at the thresholds used.")
        L.append("")
    else:
        by_type: Dict[str, List[Tangle]] = defaultdict(list)
        for t in model.tangles:
            by_type[t.type].append(t)
        for tt, group in by_type.items():
            L.append(f"### {TANGLE_LABEL.get(tt, tt)} ({len(group)})")
            L.append("")
            for t in group:
                where = ", ".join(
                    f"{a.seqname}:{a.start:,}-{a.end:,}" for a in t.anchors[:6]
                ) or "unplaced"
                L.append(f"- **{t.id}** - {t.description}")
                L.append(f"  - located at: {where}")
                if t.evidence:
                    L.append(f"  - evidence: {'; '.join(t.evidence)}")
            L.append("")

    if model.coverage_anomalies:
        L.append("## Coverage outliers")
        L.append("")
        L.append("| Region | Length | Direction | Depth vs genome median |")
        L.append("|---|---:|---|---:|")
        for a in sorted(model.coverage_anomalies, key=lambda a: -abs(a["peak"] - 1))[:50]:
            L.append(
                f"| {a['seqname']}:{a['start']:,}-{a['end']:,} | "
                f"{human_bp(a['end'] - a['start'])} | {a['kind']} | {a['peak']:.1f}x |"
            )
        L.append("")

    if model.warnings:
        L.append("## Warnings raised during this run")
        L.append("")
        for w in model.warnings:
            L.append(f"- {w}")
        L.append("")

    L.append("## How to read this, and what it cannot tell you")
    L.append("")
    L.append(CAVEATS)
    L.append("")
    L.append("## Files written")
    L.append("")
    for k, v in files.items():
        L.append(f"- {k}: `{v}`")
    L.append("")
    L.append(
        "To change any call, edit the karyotype config and re-run with `--config`. "
        "Config values always win over inference."
    )
    L.append("")
    return "\n".join(L)
