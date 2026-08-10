#!/usr/bin/env python3
"""
detangler - turn assembler output into an annotated ideogram.

Two modes, chosen by what you give it.

ASSEMBLED MODE (--fasta / --fai / --assembly-report, plus --gfa)
  You already have chromosome-scale sequences. The tool infers the karyotype -
  how many chromosomes, which sequence is the mitochondrion or plastid, what is
  unplaced - and locates assembly-graph tangles (collapsed repeats,
  inter-chromosomal junctions, inverted repeats, bubbles) onto chromosome
  coordinates, so you can see which repeat is making two chromosomes appear to
  join in Bandage.

GRAPH-FIRST MODE (--gfa alone, optionally --flye-info)
  You only have contigs and a graph. The tool sets a single-copy depth baseline,
  estimates a copy number for every segment, classifies segments (backbone,
  repeat, tandem array, organelle candidate, low-coverage/foreign), reads
  position within each contig's graph_path, and then enumerates and RANKS
  hypotheses about which segments form which chromosome. It does not pick one
  silently, and it never claims to know which chromosome is which.

Both modes write: a publication SVG, an interactive HTML view, an editable
karyotype config, a Markdown report that keeps observations, derived estimates
and hypotheses in separate sections, a redraw of the assembly graph coloured by
inferred class, a FASTA of the segments worth identifying, and ready-to-run
BLAST commands. The graph and the chromosome figure share one colour map, so a
node in the graph and a block on a chromosome match by construction.

Prior art worth knowing before you use this: telomere motif work is better done
with tidk, which discovers the motif rather than assuming it; contamination
calling belongs to BlobToolKit; MetagenomeScope is a far better interactive
graph browser. What is not covered elsewhere, as far as we can tell, is the
ranked contig-to-chromosome hypotheses with their evidence, and inferring the
chromosome count from telomere-capped ends.

Dependencies: Python 3.8+ standard library. PyYAML is used for the config file
if it is importable; otherwise the config is written as JSON. BLAST+ is only
needed if you ask for a similarity search.

Usage
-----
  # assembled genome
  detangler.py --fai asm.fa.fai --gfa asm.gfa --paf segs_to_asm.paf \\
                --coverage cov.regions.bed.gz --out-dir results --prefix myasm

  # graph only, with expectations to constrain the hypotheses
  detangler.py --gfa assembly_graph.gfa --flye-info assembly_info.txt \\
                --expected-genome-size 36.5m --expected-chromosomes 4 \\
                --out-dir results --prefix myasm

  # identify the repeats it found
  detangler.py --gfa assembly_graph.gfa --blast-subject repeats.fa ...
  # ...or run the emitted results/myasm_blast_commands.sh on a cluster and
  #    feed the table back with --blast-hits

  # apply your own calls: edit results/myasm_karyotype.yaml, then
  detangler.py --config results/myasm_karyotype.yaml --out-dir results --prefix myasm

  # synthetic data to kick the tyres
  detangler.py --demo demo_data --out-dir results --prefix demo

Every call carries its evidence, and the report ends with what the method cannot
tell you. Read that before putting any of it in a manuscript.
"""

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

VERSION = "1.0"

# --------------------------------------------------------------------------
# palette (Okabe-Ito derived, colour-vision-deficiency safe)
# --------------------------------------------------------------------------
PALETTE = {
    "chromosome": "#4c72b0",
    "mitochondrion": "#d55e00",
    "plastid": "#009e73",
    "unplaced": "#9aa0a6",
    "unassigned": "#8d8d8d",
    "bar_edge": "#2b2b2b",
    "bg": "#ffffff",
    "text": "#1a1a1a",
    "muted": "#6b6b6b",
    "grid": "#e3e3e3",
}

TANGLE_STYLE = {
    # type                         colour     dash
    "collapsed_repeat": ("#cc79a7", ""),
    "interchromosomal_junction": ("#d55e00", ""),
    "intrachromosomal_repeat": ("#0072b2", "6 3"),
    "inverted_repeat": ("#7b3294", "2 3"),
    "bubble": ("#56b4e9", "1 3"),
    "circular": ("#009e73", ""),
    # graph-first mode: same colours as the assembly graph figure, so the two
    # figures can be read side by side
    "repeat_segment": ("#c2703d", ""),
    "tandem_array": ("#b5487f", ""),
    "at_rich_region": ("#b0a04a", ""),
    "low_coverage_region": ("#a8a8a8", "4 3"),
}

TANGLE_LABEL = {
    "collapsed_repeat": "Collapsed repeat (multi-copy node)",
    "interchromosomal_junction": "Inter-chromosomal join (shared repeat)",
    "intrachromosomal_repeat": "Intra-chromosomal repeat",
    "inverted_repeat": "Inverted / hairpin repeat",
    "bubble": "Bubble (heterozygosity or small variant)",
    "circular": "Circular segment (self-link)",
    "repeat_segment": "Repeat segment (multi-copy)",
    "tandem_array": "Tandem array",
    "at_rich_region": "AT-rich segment",
    "low_coverage_region": "Low-coverage / foreign segment",
}


# ==========================================================================
# data model
# ==========================================================================
@dataclass
class Evidence:
    """One reason contributing to a classification call."""

    signal: str
    detail: str
    weight: float

    def as_text(self) -> str:
        sign = "+" if self.weight >= 0 else ""
        return f"{self.signal}: {self.detail} ({sign}{self.weight:.1f})"


@dataclass
class SeqRecord:
    name: str
    length: int
    gc: Optional[float] = None
    n_frac: Optional[float] = None
    role: str = "unknown"  # chromosome | mitochondrion | plastid | unplaced
    role_scores: Dict[str, float] = field(default_factory=dict)
    evidence: List[Evidence] = field(default_factory=list)
    circular: bool = False
    depth: Optional[float] = None
    label: Optional[str] = None
    manual: bool = False  # role came from a config override
    note: str = ""  # short reason, shown beside unassigned sequences
    # (start, end, segment_name, colour) blocks making up this molecule, so a
    # segment can be traced from the assembly graph figure to this one
    blocks: List[Tuple[int, int, str, str]] = field(default_factory=list)
    # True when the blocks tile the whole molecule (a chain built from segments),
    # False when they are placements dotted along an already-assembled sequence
    blocks_tile: bool = False
    # {'top'|'bottom': [(segment, colour)]} - repeats the graph attaches to a FREE
    # end of this molecule. They are not part of the chain and are not on the Mb
    # scale, so they are drawn hanging off the end of the bar rather than in it.
    caps: Dict[str, List[Tuple[str, str]]] = field(default_factory=dict)

    @property
    def display(self) -> str:
        return self.label or self.name


@dataclass
class GfaSegment:
    name: str
    length: int
    depth: Optional[float] = None
    gc: Optional[float] = None


@dataclass
class GfaLink:
    a: str
    a_orient: str
    b: str
    b_orient: str
    overlap: str = "*"


@dataclass
class Placement:
    """A graph segment located on an assembled sequence."""

    segment: str
    seqname: str
    start: int
    end: int
    orient: str = "+"
    identity: Optional[float] = None
    source: str = "identity"  # agp | paf | map | identity


@dataclass
class Anchor:
    """Where a tangle touches the ideogram."""

    seqname: str
    start: int
    end: int
    segment: str


@dataclass
class Tangle:
    id: str
    type: str
    segments: List[str]
    anchors: List[Anchor]
    multiplicity: Optional[float] = None
    depth_ratio: Optional[float] = None
    length: Optional[int] = None
    description: str = ""
    evidence: List[str] = field(default_factory=list)

    @property
    def sequences(self) -> List[str]:
        seen: List[str] = []
        for a in self.anchors:
            if a.seqname not in seen:
                seen.append(a.seqname)
        return seen


@dataclass
class CoverageWindow:
    seqname: str
    start: int
    end: int
    depth: float


# ==========================================================================
# small helpers
# ==========================================================================
def smart_open(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def human_bp(n: float) -> str:
    if n >= 1e9:
        return f"{n / 1e9:.2f} Gb"
    if n >= 1e6:
        return f"{n / 1e6:.2f} Mb"
    if n >= 1e3:
        return f"{n / 1e3:.1f} kb"
    return f"{int(n)} bp"


def median(values: Seq[float]) -> float:
    vals = sorted(v for v in values if v is not None)
    if not vals:
        return 0.0
    mid = len(vals) // 2
    if len(vals) % 2:
        return float(vals[mid])
    return (vals[mid - 1] + vals[mid]) / 2.0


def nx_stat(lengths: Seq[int], x: float = 50.0) -> int:
    """N50-style statistic. x is a percentage."""
    if not lengths:
        return 0
    total = sum(lengths)
    target = total * x / 100.0
    run = 0
    for L in sorted(lengths, reverse=True):
        run += L
        if run >= target:
            return L
    return min(lengths)


def esc(s: object) -> str:
    return (
        str(s)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


class Log:
    def __init__(self, quiet: bool = False):
        self.quiet = quiet
        self.warnings: List[str] = []

    def info(self, msg: str) -> None:
        if not self.quiet:
            print(f"[detangler] {msg}", file=sys.stderr)

    def warn(self, msg: str) -> None:
        self.warnings.append(msg)
        print(f"[detangler] WARNING: {msg}", file=sys.stderr)


# ==========================================================================
# parsers
# ==========================================================================
def parse_fai(path: str) -> List[SeqRecord]:
    """samtools faidx index: name, length, offset, linebases, linewidth."""
    recs: List[SeqRecord] = []
    with smart_open(path) as fh:
        for line_no, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 2:
                raise ValueError(f"{path}:{line_no}: expected at least 2 tab-separated fields")
            try:
                recs.append(SeqRecord(name=f[0], length=int(f[1])))
            except ValueError:
                raise ValueError(f"{path}:{line_no}: column 2 is not an integer length")
    if not recs:
        raise ValueError(f"{path}: no records parsed")
    return recs


def parse_fasta(path: str) -> List[SeqRecord]:
    """Streams a FASTA, recording length, GC and N fraction per sequence."""
    recs: List[SeqRecord] = []
    name: Optional[str] = None
    length = gc = n_count = 0

    def flush() -> None:
        nonlocal name, length, gc, n_count
        if name is None:
            return
        acgt = length - n_count
        recs.append(
            SeqRecord(
                name=name,
                length=length,
                gc=(gc / acgt) if acgt > 0 else None,
                n_frac=(n_count / length) if length else None,
            )
        )
        name, length, gc, n_count = None, 0, 0, 0

    with smart_open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                name = line[1:].strip().split()[0] if line[1:].strip() else "unnamed"
            elif name is not None:
                s = line.strip()
                length += len(s)
                up = s.upper()
                gc += up.count("G") + up.count("C")
                n_count += up.count("N")
    flush()
    if not recs:
        raise ValueError(f"{path}: no FASTA records found")
    return recs


def parse_assembly_report(path: str) -> List[SeqRecord]:
    """
    NCBI *_assembly_report.txt. Uses the column header comment line to locate
    fields rather than assuming positions, because the layout has changed
    between NCBI releases.
    """
    header: List[str] = []
    recs: List[SeqRecord] = []
    with smart_open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#"):
                stripped = line.lstrip("#").strip()
                if "Sequence-Name" in stripped and "Sequence-Length" in stripped:
                    header = [c.strip() for c in stripped.split("\t")]
                continue
            if not line.strip():
                continue
            f = line.split("\t")
            if header and len(f) >= len(header):
                row = dict(zip(header, f))
                name = row.get("RefSeq-Accn") or row.get("GenBank-Accn") or row.get("Sequence-Name")
                if name in (None, "na", ""):
                    name = row.get("Sequence-Name")
                try:
                    length = int(row.get("Sequence-Length", "0"))
                except ValueError:
                    continue
                rec = SeqRecord(name=str(name), length=length, label=row.get("Sequence-Name"))
                role = _role_from_assembly_report(row)
                if role:
                    rec.role = role
                    rec.manual = True
                    rec.evidence.append(
                        Evidence(
                            "assembly-report",
                            f"Sequence-Role={row.get('Sequence-Role')}, "
                            f"Assigned-Molecule-Location/Type="
                            f"{row.get('Assigned-Molecule-Location/Type')}",
                            10.0,
                        )
                    )
                recs.append(rec)
    if not recs:
        raise ValueError(f"{path}: could not parse an assembly report (no data rows matched header)")
    return recs


def _role_from_assembly_report(row: Dict[str, str]) -> Optional[str]:
    loc = (row.get("Assigned-Molecule-Location/Type") or "").lower()
    srole = (row.get("Sequence-Role") or "").lower()
    if "mitochond" in loc:
        return "mitochondrion"
    if "chloroplast" in loc or "plastid" in loc:
        return "plastid"
    if srole == "assembled-molecule" and "chromosome" in loc:
        return "chromosome"
    if srole in ("unplaced-scaffold", "unlocalized-scaffold"):
        return "unplaced"
    return None


_TAG_RE = re.compile(r"^([A-Za-z][A-Za-z0-9]):([AifZJHB]):(.*)$")


def _gfa_tags(fields: Seq[str]) -> Dict[str, str]:
    tags: Dict[str, str] = {}
    for f in fields:
        m = _TAG_RE.match(f)
        if m:
            tags[m.group(1)] = m.group(3)
    return tags


_SPADES_COV = re.compile(r"_cov_([0-9.]+)")


def parse_gfa(path: str, log: Log) -> Tuple[Dict[str, GfaSegment], List[GfaLink]]:
    """
    GFA v1 (the flavour Bandage reads). S and L lines only; P/W/C lines are
    ignored. Depth is taken from dp/DP if present, else derived from KC/RC
    divided by segment length, else from a SPAdes-style name.
    """
    segs: Dict[str, GfaSegment] = {}
    links: List[GfaLink] = []
    saw_gfa2 = False

    with smart_open(path) as fh:
        for line_no, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line:
                continue
            f = line.split("\t")
            rec = f[0]
            if rec == "H":
                tags = _gfa_tags(f[1:])
                if tags.get("VN", "").startswith("2"):
                    saw_gfa2 = True
            elif rec == "S":
                if len(f) < 3:
                    log.warn(f"{path}:{line_no}: malformed S line, skipped")
                    continue
                name, sequence = f[1], f[2]
                tags = _gfa_tags(f[3:])
                if "LN" in tags:
                    try:
                        length = int(tags["LN"])
                    except ValueError:
                        length = 0
                elif sequence != "*":
                    length = len(sequence)
                else:
                    log.warn(f"{path}:{line_no}: segment {name} has no sequence and no LN tag")
                    length = 0
                gc = None
                if sequence != "*" and sequence:
                    up = sequence.upper()
                    acgt = len(up) - up.count("N")
                    if acgt > 0:
                        gc = (up.count("G") + up.count("C")) / acgt
                depth = _segment_depth(name, length, tags)
                segs[name] = GfaSegment(name=name, length=length, depth=depth, gc=gc)
            elif rec == "L":
                if len(f) < 5:
                    log.warn(f"{path}:{line_no}: malformed L line, skipped")
                    continue
                links.append(
                    GfaLink(
                        a=f[1],
                        a_orient=f[2],
                        b=f[3],
                        b_orient=f[4],
                        overlap=f[5] if len(f) > 5 else "*",
                    )
                )

    if saw_gfa2:
        log.warn(
            "GFA header declares version 2. Only GFA1 S/L semantics are parsed; "
            "results may be incomplete."
        )
    if not segs:
        raise ValueError(f"{path}: no S (segment) lines found - is this a GFA?")

    # links referencing unknown segments would corrupt the graph analysis
    unknown = {l.a for l in links if l.a not in segs} | {l.b for l in links if l.b not in segs}
    if unknown:
        log.warn(
            f"{len(unknown)} segment name(s) referenced by L lines are not defined by S lines; "
            "those links are dropped"
        )
        links = [l for l in links if l.a in segs and l.b in segs]
    return segs, links


def _segment_depth(name: str, length: int, tags: Dict[str, str]) -> Optional[float]:
    for key in ("dp", "DP"):
        if key in tags:
            try:
                return float(tags[key])
            except ValueError:
                pass
    for key in ("KC", "RC", "FC"):
        if key in tags and length > 0:
            try:
                return float(tags[key]) / float(length)
            except ValueError:
                pass
    m = _SPADES_COV.search(name)
    if m:
        try:
            return float(m.group(1))
        except ValueError:
            pass
    return None


def parse_agp(path: str, log: Log) -> List[Placement]:
    """
    AGP 2.x. Only component lines (component_type not in {N, U}) are used, which
    is exactly the contig-into-scaffold placement we want.
    Columns: object, object_beg, object_end, part_number, component_type,
             component_id, component_beg, component_end, orientation
    """
    out: List[Placement] = []
    with smart_open(path) as fh:
        for line_no, line in enumerate(fh, 1):
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                log.warn(f"{path}:{line_no}: fewer than 9 columns, skipped")
                continue
            if f[4].upper() in ("N", "U"):
                continue
            try:
                out.append(
                    Placement(
                        segment=f[5],
                        seqname=f[0],
                        start=int(f[1]) - 1,  # AGP is 1-based inclusive
                        end=int(f[2]),
                        orient=f[8] if f[8] in ("+", "-") else "+",
                        source="agp",
                    )
                )
            except ValueError:
                log.warn(f"{path}:{line_no}: non-integer coordinates, skipped")
    return out


def parse_paf(path: str, min_block: int, log: Log) -> List[Placement]:
    """
    PAF from aligning graph segments (query) against the assembly (target).
    All alignments at or above min_block are kept - multi-mapping is the signal
    we are looking for, so we must not collapse to a best hit.
    """
    out: List[Placement] = []
    with smart_open(path) as fh:
        for line_no, line in enumerate(fh, 1):
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                log.warn(f"{path}:{line_no}: fewer than 12 PAF columns, skipped")
                continue
            try:
                strand, tname = f[4], f[5]
                tstart, tend = int(f[7]), int(f[8])
                matches, block = int(f[9]), int(f[10])
            except ValueError:
                log.warn(f"{path}:{line_no}: non-numeric PAF fields, skipped")
                continue
            if block < min_block:
                continue
            out.append(
                Placement(
                    segment=f[0],
                    seqname=tname,
                    start=tstart,
                    end=tend,
                    orient=strand if strand in ("+", "-") else "+",
                    identity=(matches / block) if block else None,
                    source="paf",
                )
            )
    return out


def parse_segment_map(path: str, log: Log) -> List[Placement]:
    """Manual TSV: segment, seqname, start, end[, orient]. 0-based half-open."""
    out: List[Placement] = []
    with smart_open(path) as fh:
        for line_no, line in enumerate(fh, 1):
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                log.warn(f"{path}:{line_no}: need at least 4 columns, skipped")
                continue
            try:
                out.append(
                    Placement(
                        segment=f[0],
                        seqname=f[1],
                        start=int(f[2]),
                        end=int(f[3]),
                        orient=f[4] if len(f) > 4 and f[4] in ("+", "-") else "+",
                        source="map",
                    )
                )
            except ValueError:
                log.warn(f"{path}:{line_no}: non-integer coordinates, skipped")
    return out


def parse_coverage_bed(path: str, log: Log) -> List[CoverageWindow]:
    """BED4: chrom, start, end, depth. mosdepth *.regions.bed.gz fits directly."""
    out: List[CoverageWindow] = []
    with smart_open(path) as fh:
        for line_no, line in enumerate(fh, 1):
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            try:
                out.append(
                    CoverageWindow(seqname=f[0], start=int(f[1]), end=int(f[2]), depth=float(f[3]))
                )
            except ValueError:
                continue
    if not out:
        log.warn(f"{path}: no coverage windows parsed (expected BED4 chrom/start/end/depth)")
    return out


def parse_annotation(path: str, log: Log) -> List[Dict]:
    """BED3+ or GFF3. Returns dicts with seqname, start, end, name, kind."""
    out: List[Dict] = []
    is_gff = path.lower().endswith((".gff", ".gff3", ".gtf", ".gff.gz", ".gff3.gz", ".gtf.gz"))
    with smart_open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            try:
                if is_gff:
                    if len(f) < 5:
                        continue
                    attrs = f[8] if len(f) > 8 else ""
                    m = re.search(r"(?:^|;)\s*(?:Name|ID)=([^;]+)", attrs)
                    out.append(
                        {
                            "seqname": f[0],
                            "start": int(f[3]) - 1,
                            "end": int(f[4]),
                            "name": m.group(1) if m else f[2],
                            "kind": f[2],
                        }
                    )
                else:
                    if len(f) < 3:
                        continue
                    out.append(
                        {
                            "seqname": f[0],
                            "start": int(f[1]),
                            "end": int(f[2]),
                            "name": f[3] if len(f) > 3 else "",
                            "kind": f[6] if len(f) > 6 else (f[3] if len(f) > 3 else "feature"),
                        }
                    )
            except ValueError:
                continue
    if not out:
        log.warn(f"{path}: no annotation features parsed")
    return out


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

# Organelle size envelopes. Deliberately wide: plant mitogenomes reach several
# hundred kb, and some animal plastid-free lineages sit at the bottom edge.
MITO_RANGE = (11_000, 2_000_000)
MITO_TYPICAL = (14_000, 700_000)
PLASTID_RANGE = (70_000, 250_000)
PLASTID_TYPICAL = (110_000, 180_000)


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
    Like build_adjacency, but keyed by (segment, end) with end 's' (forward
    start) or 'e' (forward end), so the two ends of a segment can be assessed
    independently.
    """
    adj: Dict[Tuple[str, str], Set[str]] = defaultdict(set)
    for l in links:
        if l.a == l.b:
            continue
        adj[(l.a, _link_end(l.a_orient, True))].add(l.b)
        adj[(l.b, _link_end(l.b_orient, False))].add(l.a)
    return adj


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


def _natural_key(s: str):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", s)]


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


# ==========================================================================
# layout + SVG rendering
# ==========================================================================
MARGIN_L = 86
MARGIN_R = 40
MARGIN_T = 96
BAR_W = 48
COV_W = 20
GAP = 96
MAX_BAR_H = 620
MIN_ORG_H = 52
LEGEND_H = 150
PANEL_W = 300


def wrap_text(text: str, max_chars: int) -> List[str]:
    """Greedy word wrap. max_chars is derived from an average glyph width."""
    words, lines, cur = text.split(), [], ""
    for w in words:
        if cur and len(cur) + 1 + len(w) > max_chars:
            lines.append(cur)
            cur = w
        else:
            cur = f"{cur} {w}".strip()
    if cur:
        lines.append(cur)
    return lines or [""]


class Layout:
    def __init__(self, model: Model, show_coverage: bool, header_h: float = MARGIN_T):
        self.header_h = header_h
        self.order = model.drawable()
        self.x: Dict[str, float] = {}
        self.top: Dict[str, float] = {}
        self.height: Dict[str, float] = {}
        self.length: Dict[str, int] = {}
        self.not_to_scale: Set[str] = set()
        self.show_coverage = show_coverage

        chrom_lengths = [s.length for s in self.order if s.role == "chromosome"] or [
            s.length for s in self.order
        ] or [1]
        self.max_len = max(chrom_lengths)
        self.scale = MAX_BAR_H / float(self.max_len)

        step = BAR_W + GAP + (COV_W + 6 if show_coverage else 0)
        self.step = step
        # Chain labels can be long ("chain 2: edge_7 + edge_2"), so wrap them to
        # the column width rather than letting neighbours collide.
        self.label_lines: Dict[str, List[str]] = {
            s.name: wrap_text(s.display, max(int(step / 6.2), 9))[:3] for s in self.order
        }
        self.max_label_lines = max((len(v) for v in self.label_lines.values()), default=1)
        for i, s in enumerate(self.order):
            self.x[s.name] = MARGIN_L + i * step
            h = s.length * self.scale
            if h < MIN_ORG_H:
                h = MIN_ORG_H
                self.not_to_scale.add(s.name)
            self.height[s.name] = h
            self.top[s.name] = header_h
            self.length[s.name] = s.length

        n = max(len(self.order), 1)
        self.panel_x = MARGIN_L + n * step + 16
        self.panel = bool(model.unassigned())
        panel_w = PANEL_W if self.panel else 0
        self.width = max(MARGIN_L + n * step - GAP + MARGIN_R + panel_w, 940)
        self.height_total = header_h + MAX_BAR_H + LEGEND_H  # refined by render_svg

    @property
    def text_cols(self) -> int:
        """Characters that fit on one line at ~12px in the drawing area."""
        return max(int((self.width - MARGIN_L - MARGIN_R) / 6.4), 40)

    @property
    def panel_cols(self) -> int:
        return max(int((PANEL_W - 20) / 5.6), 24)

    def y(self, seqname: str, pos: float) -> float:
        L = max(self.length.get(seqname, 1), 1)
        frac = min(max(pos / float(L), 0.0), 1.0)
        return self.top[seqname] + frac * self.height[seqname]

    def cx(self, seqname: str) -> float:
        return self.x[seqname] + BAR_W / 2.0


def _bar_path(x: float, y: float, w: float, h: float, rt: float, rb: float) -> str:
    """
    Rectangle with independently rounded top and bottom ends.

    Used instead of a clipPath because clipping is unevenly supported outside
    browsers - ImageMagick, for one, silently drops the clipped group, which
    would lose the segment blocks in any PNG conversion.
    """
    rt = max(0.0, min(rt, w / 2, h / 2))
    rb = max(0.0, min(rb, w / 2, h / 2))
    p = [f"M {x:.1f} {y + rt:.1f}"]
    if rt:
        p.append(f"A {rt:.1f} {rt:.1f} 0 0 1 {x + rt:.1f} {y:.1f}")
    else:
        p.append(f"L {x:.1f} {y:.1f}")
    p.append(f"L {x + w - rt:.1f} {y:.1f}")
    if rt:
        p.append(f"A {rt:.1f} {rt:.1f} 0 0 1 {x + w:.1f} {y + rt:.1f}")
    p.append(f"L {x + w:.1f} {y + h - rb:.1f}")
    if rb:
        p.append(f"A {rb:.1f} {rb:.1f} 0 0 1 {x + w - rb:.1f} {y + h:.1f}")
    else:
        p.append(f"L {x + w:.1f} {y + h:.1f}")
    p.append(f"L {x + rb:.1f} {y + h:.1f}")
    if rb:
        p.append(f"A {rb:.1f} {rb:.1f} 0 0 1 {x:.1f} {y + h - rb:.1f}")
    p.append("Z")
    return " ".join(p)


def _annotation_colour(kind: str) -> str:
    palette = ["#b07aa1", "#59a14f", "#edc948", "#e15759", "#76b7b2", "#ff9da7", "#9c755f"]
    return palette[abs(hash(kind)) % len(palette)]


def _arc_path(x1: float, y1: float, x2: float, y2: float) -> str:
    """Cubic Bezier between two ideogram points; bows outward if same column."""
    if abs(x1 - x2) < 1e-6:
        bow = 46 + min(abs(y2 - y1) * 0.35, 90)
        return f"M {x1:.1f} {y1:.1f} C {x1 + bow:.1f} {y1:.1f}, {x2 + bow:.1f} {y2:.1f}, {x2:.1f} {y2:.1f}"
    dx = (x2 - x1) * 0.45
    return (
        f"M {x1:.1f} {y1:.1f} C {x1 + dx:.1f} {y1:.1f}, "
        f"{x2 - dx:.1f} {y2:.1f}, {x2:.1f} {y2:.1f}"
    )


def representative_anchors(t: Tangle, drawn: Set[str]) -> List[Anchor]:
    """One anchor per sequence (the first), restricted to sequences we draw."""
    out: List[Anchor] = []
    seen: Set[str] = set()
    for a in t.anchors:
        if a.seqname in drawn and a.seqname not in seen:
            seen.add(a.seqname)
            out.append(a)
    if len(out) <= 1:
        # intra-sequence repeat: keep up to 4 distinct positions on one sequence
        pos: List[Anchor] = []
        for a in t.anchors:
            if a.seqname in drawn and not any(abs(a.start - p.start) < 1000 for p in pos):
                pos.append(a)
        if len(pos) > 1:
            return pos[:4]
    return out[:4]


def ideogram_geometry(model: Model) -> Tuple[Layout, List[str], bool, float]:
    """
    The ideogram's layout, header text and header height. A thin view over
    _ideogram_frame so there is only ever one layout calculation: if the
    renderer and the paired figure disagreed by a pixel, every leader line
    would point to the wrong place.
    """
    f = _ideogram_frame(model)
    return f["lay"], f["head_lines"], f["show_cov"], f["header_h"]  # type: ignore


def ideogram_block_anchors(model: Model) -> Dict[str, Tuple[float, float]]:
    """Left edge and vertical centre of each segment block, in ideogram coordinates."""
    lay, _, _, _ = ideogram_geometry(model)
    out: Dict[str, Tuple[float, float]] = {}
    for s in lay.order:
        for b_start, b_end, seg, _colour in s.blocks:
            y = (lay.y(s.name, b_start) + lay.y(s.name, b_end)) / 2.0
            out.setdefault(seg, (lay.x[s.name], y))
    return out


def _ideogram_frame(model: Model) -> Dict[str, object]:
    """
    Geometry for the chromosome figure, computed once and shared by the renderer
    and by anything that needs to point at a block from outside - the paired
    figure draws leader lines to these exact coordinates.
    """
    show_cov = bool(model.coverage) and model.settings.get("coverage", True)
    probe = Layout(model, show_cov)
    drawn = {s.name for s in probe.order}

    # v9: no summary paragraph under the panel title. The reasoning belongs in the
    # report; the figure carries only what points at something it draws.
    head_lines: List[str] = []
    header_h = 54 + 26 + 15 * (probe.max_label_lines - 1)

    lay = Layout(model, show_cov, header_h)
    legend_svg, legend_bottom = _legend_svg(model, lay)
    total_h = max(legend_bottom + 26, header_h + MAX_BAR_H + 90)
    return {
        "lay": lay,
        "head_lines": head_lines,
        "drawn": drawn,
        "show_cov": show_cov,
        "header_h": header_h,
        "legend_svg": legend_svg,
        "total_h": total_h,
    }


def render_svg(model: Model, interactive: bool = False) -> str:
    frame = _ideogram_frame(model)
    lay, head_lines, drawn = frame["lay"], frame["head_lines"], frame["drawn"]
    show_cov, header_h = frame["show_cov"], frame["header_h"]
    legend_svg, total_h = frame["legend_svg"], frame["total_h"]
    P: List[str] = []
    add = P.append

    add(
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{lay.width:.0f}" '
        f'height="{total_h:.0f}" viewBox="0 0 {lay.width:.0f} {total_h:.0f}" '
        f'font-family="Helvetica, Arial, sans-serif">'
    )
    add(
        '<defs>'
        '<pattern id="nts" width="6" height="6" patternUnits="userSpaceOnUse" '
        'patternTransform="rotate(45)">'
        '<rect width="6" height="6" fill="#ffffff" fill-opacity="0"/>'
        '<line x1="0" y1="0" x2="0" y2="6" stroke="#ffffff" stroke-width="2" stroke-opacity="0.55"/>'
        '</pattern>'
        '</defs>'
    )
    add(f'<rect width="100%" height="100%" fill="{PALETTE["bg"]}"/>')

    # ---- title ----
    add(
        f'<text x="{MARGIN_L}" y="40" font-size="{FS_HEADING}" font-weight="600" '
        f'fill="{PALETTE["text"]}">{esc(model.title)}</text>'
    )

    # ---- scale ruler ----
    add(f'<g id="ruler" stroke="{PALETTE["grid"]}">')
    tick = _nice_tick(lay.max_len)
    v = 0
    while v <= lay.max_len:
        y = header_h + v * lay.scale
        add(f'<line x1="{MARGIN_L - 46}" y1="{y:.1f}" x2="{MARGIN_L - 12}" y2="{y:.1f}"/>')
        add(
            f'<text x="{MARGIN_L - 50}" y="{y + 4:.1f}" font-size="10.5" text-anchor="end" '
            f'fill="{PALETTE["muted"]}" stroke="none">{_tick_label(v)}</text>'
        )
        v += tick
    add("</g>")

    # ---- tangle arcs (behind the bars) ----
    add('<g id="layer-tangles">')
    # v9 drops the arcs joining one chromosome to another. They were the biggest
    # source of visual noise and duplicated what the shared segment colours
    # already say: a repeat linking two chains is drawn as a copy on each. Only
    # single-point features survive, as a marker beside the bar.
    for t in model.tangles:
        anchors = representative_anchors(t, drawn)
        if len(anchors) >= 2:
            continue
        colour, _dash = TANGLE_STYLE.get(t.type, ("#888888", ""))
        attrs = ""
        if interactive:
            attrs = (
                f' class="tangle" data-id="{esc(t.id)}" data-type="{esc(t.type)}"'
                f' data-desc="{esc(t.description)}"'
            )
        for a in anchors:
            y = lay.y(a.seqname, (a.start + a.end) / 2.0)
            x = lay.x[a.seqname] + BAR_W + 5
            add(
                f'<path d="M {x:.1f} {y:.1f} l 9 -5 l 0 10 z" fill="{colour}" '
                f'fill-opacity="0.9"{attrs}/>'
            )
    add("</g>")

    # ---- chromosome bars ----
    add('<g id="layer-chromosomes">')
    for s in lay.order:
        x, top, h = lay.x[s.name], lay.top[s.name], lay.height[s.name]
        fill = PALETTE.get(s.role, PALETTE["chromosome"])
        rx = BAR_W / 2.0
        battrs = ""
        if interactive:
            battrs = (
                f' class="chrom" data-name="{esc(s.name)}" data-role="{esc(s.role)}"'
                f' data-length="{s.length}" data-depth="{s.depth if s.depth is not None else ""}"'
            )

        # v9: a circular molecule is drawn as a ring, not as a bar with a little
        # circle underneath it. Nothing circular is to scale against the nuclear
        # chromosomes anyway, so a ring is both truer and less misleading.
        if s.circular:
            # the segment's OWN colour, the one it has in the graph panel. Falling
            # back to the role colour here broke the figure's one promise: edge_11
            # came out cyan on the left and orange on the right.
            seg_colour = (
                s.blocks[0][3] if s.blocks
                else model.segment_colours.get(s.name, fill)
            )
            # an organelle is not on the nuclear scale, so it is drawn thinner
            # than a chromosome bar as well as round: nothing about it should
            # invite being read off the Mb axis
            ring_w = BAR_W * 0.45
            r = BAR_W * 0.95
            ccx, ccy = x + BAR_W / 2.0, top + r + 6
            add(
                f'<circle cx="{ccx:.1f}" cy="{ccy:.1f}" r="{r:.1f}" fill="none" '
                f'stroke="{PALETTE["bar_edge"]}" stroke-width="{ring_w + 2.0:.1f}"{battrs}/>'
            )
            add(
                f'<circle cx="{ccx:.1f}" cy="{ccy:.1f}" r="{r:.1f}" fill="none" '
                f'stroke="{seg_colour}" stroke-width="{ring_w:.1f}"/>'
            )
            add(
                f'<text x="{ccx:.1f}" y="{top - 10:.1f}" font-size="{FS_SUB + 2}" '
                f'text-anchor="middle" fill="{PALETTE["text"]}" font-weight="600">'
                f'{esc(s.role)}</text>'
            )
            add(
                f'<text x="{ccx:.1f}" y="{ccy + r * 2 + FS_SUB + 6:.1f}" '
                f'font-size="{FS_SUB}" text-anchor="middle" fill="{PALETTE["muted"]}">'
                f'{human_bp(s.length)}</text>'
            )
            continue

        add(
            f'<path d="{_bar_path(x, top, BAR_W, h, 0.0 if s.caps.get("top") else rx, 0.0 if s.caps.get("bottom") else rx)}" '
            f'fill="{fill}" fill-opacity="0.82" stroke="{PALETTE["bar_edge"]}" '
            f'stroke-width="0.9"{battrs}/>'
        )

        # Segment blocks: the same colour the segment has in the graph figure, so
        # a node over there can be found on a chromosome over here.
        if s.blocks:
            add('<g class="blocks">')
            last = len(s.blocks) - 1
            inset = 0.0 if s.blocks_tile else 3.5
            bw = BAR_W - 2 * inset
            for bi, (b_start, b_end, seg, colour) in enumerate(s.blocks):
                y1 = lay.y(s.name, b_start)
                y2 = max(lay.y(s.name, b_end), y1 + 1.6)
                bl = ""
                if interactive:
                    bl = (
                        f' class="block" data-desc="{esc(seg)}, {human_bp(b_end - b_start)}, '
                        f'on {esc(s.display)}"'
                    )
                # Tiled blocks ARE the bar, so the outer ones keep its rounded
                # ends - UNLESS a cap sits against that end, in which case the
                # molecule continues and the corner must be square. Only the
                # outermost piece of the whole molecule is rounded.
                rt = rx if (s.blocks_tile and bi == 0 and not s.caps.get("top")) else 0
                rb = rx if (s.blocks_tile and bi == last and not s.caps.get("bottom")) else 0
                add(
                    f'<path d="{_bar_path(x + inset, y1, bw, y2 - y1, rt, rb)}" '
                    f'fill="{colour}" fill-opacity="{0.95 if s.blocks_tile else 0.92}" '
                    f'stroke="#ffffff" stroke-width="0.6"{bl}/>'
                )
                # v9: the segment NUMBER, set inside the block, inked white or
                # dark for contrast against that block's own colour
                if y2 - y1 >= FS_LABEL + 4:
                    add(
                        f'<text x="{x + BAR_W / 2:.1f}" y="{(y1 + y2) / 2 + FS_LABEL * 0.35:.1f}" '
                        f'font-size="{FS_LABEL}" text-anchor="middle" fill="{_text_on(colour)}" '
                        f'font-weight="700">{esc(_segment_number(seg))}</text>'
                    )
            add("</g>")
        if s.name in lay.not_to_scale:
            add(
                f'<rect x="{x:.1f}" y="{top:.1f}" width="{BAR_W}" height="{h:.1f}" rx="{rx:.1f}" '
                f'ry="{rx:.1f}" fill="url(#nts)" stroke="none"/>'
            )

        # annotation bands
        for feat in model.annotations:
            if feat["seqname"] != s.name:
                continue
            y1 = lay.y(s.name, feat["start"])
            y2 = max(lay.y(s.name, feat["end"]), y1 + 1.2)
            c = _annotation_colour(str(feat.get("kind", "feature")))
            fattrs = (
                f' class="annot" data-desc="{esc(feat.get("name") or feat.get("kind"))} '
                f'{esc(s.name)}:{feat["start"]:,}-{feat["end"]:,}"'
                if interactive
                else ""
            )
            add(
                f'<rect x="{x:.1f}" y="{y1:.1f}" width="{BAR_W}" height="{y2 - y1:.1f}" '
                f'fill="{c}" fill-opacity="0.85" stroke="none"{fattrs}/>'
            )

        # coverage anomaly stripes
        for an in model.coverage_anomalies:
            if an["seqname"] != s.name:
                continue
            y1 = lay.y(s.name, an["start"])
            y2 = max(lay.y(s.name, an["end"]), y1 + 1.5)
            c = "#d62728" if an["kind"] == "high" else "#1f77b4"
            add(
                f'<rect x="{x - 5:.1f}" y="{y1:.1f}" width="4" height="{y2 - y1:.1f}" '
                f'fill="{c}" fill-opacity="0.9"/>'
            )

        # re-stroke the outline so bands do not spill past the rounded ends
        add(
            f'<path d="{_bar_path(x, top, BAR_W, h, 0.0 if s.caps.get("top") else rx, 0.0 if s.caps.get("bottom") else rx)}" '
            f'fill="none" stroke="{PALETTE["bar_edge"]}" stroke-width="1.1"/>'
        )

        # Repeats attached to a free end, drawn hanging OFF the bar rather than
        # inside it: they are not part of the molecule and not on the Mb scale,
        # but they are what tells you this end is a telomere or an rDNA block.
        # Flush against the bar, not floating beside it, so a molecule reads as
        # one object: cap, backbone, cap. Only the OUTER corner of the outermost
        # cap is rounded; every join between blocks is square so they abut.
        cap_h = BAR_W * 0.62
        for side, entries in sorted(s.caps.items()):
            n_side = len(entries)
            for ci, (seg, colour) in enumerate(entries):
                if side == "top":
                    cy = top - cap_h * (ci + 1)
                    rt = rx if ci == n_side - 1 else 0.0
                    rb = 0.0
                else:
                    cy = top + h + cap_h * ci
                    rt = 0.0
                    rb = rx if ci == n_side - 1 else 0.0
                add(
                    f'<path d="{_bar_path(x, cy, BAR_W, cap_h, rt, rb)}" fill="{colour}" '
                    f'fill-opacity="0.95" stroke="{PALETTE["bar_edge"]}" stroke-width="1.1"/>'
                )
                add(
                    f'<text x="{x + BAR_W / 2:.1f}" y="{cy + cap_h / 2 + FS_SUB * 0.36:.1f}" '
                    f'font-size="{FS_SUB + 1}" text-anchor="middle" fill="{_text_on(colour)}" '
                    f'font-weight="700">{esc(_segment_number(seg))}</text>'
                )

        # size only. Chain headings are gone: which contigs belong together is
        # shown by the numbered blocks in the bar, not by a caption above it.
        n_top = len(s.caps.get("top", []))
        add(
            f'<text x="{x + BAR_W / 2:.1f}" y="{top + h + 20 + cap_h * len(s.caps.get("bottom", [])):.1f}" '
            f'font-size="{FS_SUB}" text-anchor="middle" fill="{PALETTE["muted"]}">'
            f'{human_bp(s.length)}</text>'
        )
    add("</g>")

    # ---- coverage track ----
    if show_cov:
        add('<g id="layer-coverage">')
        gm = median([v for v in model.coverage_median.values() if v]) or 1.0
        cap = gm * 2.0
        for s in lay.order:
            ws = model.coverage.get(s.name) or []
            if not ws:
                continue
            x0 = lay.x[s.name] + BAR_W + 16
            add(
                f'<line x1="{x0:.1f}" y1="{lay.top[s.name]:.1f}" x2="{x0:.1f}" '
                f'y2="{lay.top[s.name] + lay.height[s.name]:.1f}" stroke="{PALETTE["grid"]}"/>'
            )
            pts = []
            for w in ws:
                y = lay.y(s.name, (w.start + w.end) / 2.0)
                xx = x0 + min(w.depth / cap, 1.0) * COV_W
                pts.append(f"{xx:.1f},{y:.1f}")
            if pts:
                add(
                    f'<polyline points="{" ".join(pts)}" fill="none" stroke="#333333" '
                    f'stroke-width="1.3" stroke-opacity="0.9" stroke-linejoin="round"/>'
                )
            mid = x0 + min(gm / cap, 1.0) * COV_W
            add(
                f'<line x1="{mid:.1f}" y1="{lay.top[s.name]:.1f}" x2="{mid:.1f}" '
                f'y2="{lay.top[s.name] + lay.height[s.name]:.1f}" stroke="#999999" '
                f'stroke-dasharray="2 3"/>'
            )
        add(
            f'<text x="{MARGIN_L}" y="{header_h + MAX_BAR_H + 52:.1f}" font-size="10" '
            f'fill="{PALETTE["muted"]}">Coverage track (right of each bar): 0 to 2x the genome '
            f'median ({gm:.0f}x); dashed line = median. Red/blue ticks left of a bar mark '
            f'depth outliers.</text>'
        )
        add("</g>")

    # ---- unassigned panel ----
    if lay.panel:
        add(_unassigned_panel_svg(model, lay, interactive))

    # ---- legend ----
    add(legend_svg)
    add("</svg>")
    return "\n".join(P)


def _unassigned_panel_svg(model: Model, lay: Layout, interactive: bool) -> str:
    """
    Sequences that fit no chromosome, kept visibly separate rather than being
    forced into the karyotype or dropped from the figure.

    Drawn as upright bars like the chromosomes, but deliberately narrower and on
    their own side of a divider, so they read as the same kind of object without
    implying they belong to the karyotype. Labels only - no sentences.
    """
    items = model.unassigned()
    x, y = lay.panel_x, lay.header_h
    out = ['<g id="layer-unassigned" font-family="Helvetica, Arial, sans-serif">']
    out.append(
        f'<line x1="{x - 14:.1f}" y1="{y - 26:.1f}" x2="{x - 14:.1f}" '
        f'y2="{y + MAX_BAR_H:.1f}" stroke="{PALETTE["grid"]}" stroke-dasharray="3 4"/>'
    )
    out.append(
        f'<text x="{x:.1f}" y="{y - 24:.1f}" font-size="{FS_SUB + 2}" font-weight="600" '
        f'fill="{PALETTE["text"]}">Not assigned</text>'
    )

    col_w = BAR_W * 2.4
    bw = BAR_W * 0.62
    max_len = max((s.length for s in items), default=1)
    top = y + 34.0
    max_h = 210.0
    shown = items[: max(int((lay.width - x) / col_w), 1)]
    for i, s in enumerate(shown):
        cx = x + i * col_w
        # log height: these span kb to tens of kb and would otherwise vanish
        h = max(
            18.0,
            max_h * math.log10(max(s.length, 10)) / math.log10(max(max_len, 100)),
        )
        colour = model.segment_colours.get(s.name, PALETTE["unassigned"])
        attrs = ""
        if interactive:
            attrs = (
                f' class="chrom" data-name="{esc(s.name)}" data-role="unassigned"'
                f' data-length="{s.length}"'
                f' data-depth="{s.depth if s.depth is not None else ""}"'
            )
        out.append(
            f'<rect x="{cx:.1f}" y="{top:.1f}" width="{bw:.1f}" height="{h:.1f}" '
            f'rx="{bw / 2:.1f}" ry="{bw / 2:.1f}" fill="{colour}" fill-opacity="0.9" '
            f'stroke="{PALETTE["bar_edge"]}" stroke-width="1.1"{attrs}/>'
        )
        if h >= FS_LABEL + 4:
            out.append(
                f'<text x="{cx + bw / 2:.1f}" y="{top + h / 2 + FS_LABEL * 0.35:.1f}" '
                f'font-size="{FS_LABEL}" text-anchor="middle" fill="{_text_on(colour)}" '
                f'font-weight="700">{esc(_segment_number(s.name))}</text>'
            )
        out.append(
            f'<text x="{cx + bw / 2:.1f}" y="{top + h + FS_SUB + 4:.1f}" '
            f'font-size="{FS_SUB}" text-anchor="middle" fill="{PALETTE["muted"]}">'
            f'{human_bp(s.length)}</text>'
        )
    out.append("</g>")
    return "\n".join(out)


def _nice_tick(max_len: int) -> int:
    raw = max_len / 8.0
    mag = 10 ** int(math.floor(math.log10(max(raw, 1))))
    for m in (1, 2, 5, 10):
        if raw <= m * mag:
            return int(m * mag)
    return int(10 * mag)


def _tick_label(v: int) -> str:
    if v == 0:
        return "0"
    if v >= 1e6:
        return f"{v / 1e6:g} Mb"
    if v >= 1e3:
        return f"{v / 1e3:g} kb"
    return str(v)


def _legend_svg(model: Model, lay: Layout) -> Tuple[str, float]:
    """
    v9: no floating key, no footnote block. Everything the reader needs is a
    label attached to the thing it describes, so this now draws nothing. The
    function survives because the layout asks it where the figure ends.
    """
    return "", lay.header_h + MAX_BAR_H + 40


def _legend_svg_unused(model: Model, lay: Layout) -> Tuple[str, float]:
    y0 = lay.header_h + MAX_BAR_H + 76
    out = [f'<g id="legend" font-size="11" fill="{PALETTE["text"]}">']
    x = MARGIN_L
    if any(s.blocks for s in model.sequences):
        out.append(
            f'<text x="{x}" y="{y0 + 1}" fill="{PALETTE["muted"]}">'
            f'Blocks within each bar are graph segments, coloured as in the assembly graph '
            f'figure so they can be traced between the two.</text>'
        )
        y0 += 20
    roles = [r for r in ("chromosome", "mitochondrion", "plastid") if any(s.role == r for s in model.sequences)]
    for r in roles:
        out.append(
            f'<rect x="{x}" y="{y0 - 9}" width="13" height="13" rx="6" fill="{PALETTE[r]}" '
            f'fill-opacity="0.82" stroke="{PALETTE["bar_edge"]}" stroke-width="0.8"/>'
        )
        out.append(f'<text x="{x + 19}" y="{y0 + 1}">{r}</text>')
        x += 24 + 7.2 * len(r)

    types = []
    for t in model.tangles:
        if t.type not in types:
            types.append(t.type)
    y = y0 + 24
    x = MARGIN_L
    for tt in types:
        colour, dash = TANGLE_STYLE.get(tt, ("#888888", ""))
        label = TANGLE_LABEL.get(tt, tt)
        dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
        out.append(
            f'<line x1="{x}" y1="{y - 4}" x2="{x + 26}" y2="{y - 4}" stroke="{colour}" '
            f'stroke-width="2.6"{dash_attr}/>'
        )
        out.append(f'<text x="{x + 32}" y="{y}">{esc(label)}</text>')
        x += 46 + 6.6 * len(label)
        if x > lay.width - 260:
            x = MARGIN_L
            y += 20

    notes: List[str] = []
    if lay.not_to_scale:
        notes.append(
            "* drawn at a minimum height so it stays visible: hatched bars are not to scale."
        )
    up = model.unplaced()
    if up:
        total = max(sum(s.length for s in model.sequences), 1)
        notes.append(
            f"{len(up)} unplaced sequence(s) are not drawn, totalling "
            f"{human_bp(sum(s.length for s in up))} "
            f"({100.0 * sum(s.length for s in up) / total:.1f}% of the assembly); "
            f"longest {human_bp(max(s.length for s in up))}."
        )
    for note in notes:
        for line in wrap_text(note, lay.text_cols):
            y += 18
            out.append(f'<text x="{MARGIN_L}" y="{y}" fill="{PALETTE["muted"]}">{esc(line)}</text>')
    out.append("</g>")
    return "\n".join(out), y


# ==========================================================================
# interactive HTML
# ==========================================================================
HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<title>__TITLE__</title>
<style>
  :root { --fg:#1a1a1a; --muted:#6b6b6b; --line:#e2e2e2; --panel:#fafafa; }
  * { box-sizing:border-box; }
  body { margin:0; font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif;
         color:var(--fg); background:#fff; }
  header { padding:16px 22px; border-bottom:1px solid var(--line); }
  h1 { font-size:18px; margin:0 0 4px; }
  .sub { color:var(--muted); font-size:13px; }
  .wrap { display:flex; align-items:flex-start; gap:0; }
  .canvas { flex:1 1 auto; overflow:auto; padding:10px 0 40px 0; }
  aside { width:360px; flex:0 0 360px; border-left:1px solid var(--line); height:calc(100vh - 78px);
          overflow:auto; background:var(--panel); padding:16px 18px; }
  .controls { padding:10px 22px; border-bottom:1px solid var(--line); display:flex; gap:18px;
              flex-wrap:wrap; align-items:center; font-size:13px; }
  label.chk { display:inline-flex; gap:6px; align-items:center; cursor:pointer; user-select:none; }
  h2 { font-size:13px; text-transform:uppercase; letter-spacing:.05em; color:var(--muted);
       margin:20px 0 8px; }
  h2:first-child { margin-top:0; }
  table { width:100%; border-collapse:collapse; font-size:12.5px; }
  th,td { text-align:left; padding:5px 6px; border-bottom:1px solid var(--line); vertical-align:top; }
  th { color:var(--muted); font-weight:600; }
  .card { border:1px solid var(--line); border-radius:7px; padding:9px 11px; margin-bottom:8px;
          background:#fff; cursor:pointer; }
  .card:hover, .card.on { border-color:#888; }
  .card .t { font-weight:600; font-size:12.5px; display:flex; align-items:center; gap:7px; }
  .swatch { width:11px; height:11px; border-radius:2px; flex:0 0 auto; }
  .card .d { color:var(--muted); font-size:12px; margin-top:4px; }
  .why { color:var(--muted); font-size:11.5px; margin-top:4px; font-style:italic; }
  .pill { display:inline-block; font-size:11px; padding:1px 6px; border-radius:9px;
          background:#eee; color:#444; margin-left:4px; }
  #tip { position:fixed; pointer-events:none; background:#111; color:#fff; padding:6px 9px;
         border-radius:5px; font-size:12px; max-width:340px; opacity:0; transition:opacity .1s;
         z-index:20; }
  .dim { opacity:.12 !important; }
  details { margin:6px 0; } summary { cursor:pointer; font-size:12.5px; }
  .warn { background:#fff5e6; border:1px solid #f0c987; border-radius:6px; padding:9px 11px;
          font-size:12.5px; margin-bottom:8px; }
  svg { display:block; margin:0 auto; }
</style></head><body>
<header><h1>__TITLE__</h1><div class="sub">__SUMMARY__</div></header>
<div class="controls">
  <label class="chk"><input type="checkbox" id="c-tangles" checked> Graph features</label>
  <label class="chk"><input type="checkbox" id="c-coverage" checked> Coverage</label>
  <label class="chk"><input type="checkbox" id="c-annot" checked> Annotations</label>
  <label class="chk"><input type="checkbox" id="c-legend" checked> Legend</label>
  <span style="margin-left:auto;color:var(--muted);font-size:12px">
    zoom <input type="range" id="zoom" min="50" max="220" value="100" style="vertical-align:middle">
    <span id="zv">100%</span></span>
</div>
<div class="wrap">
  <div class="canvas"><div id="svgbox">__SVG__</div></div>
  <aside>__SIDE__</aside>
</div>
<div id="tip"></div>
<script>
const tip = document.getElementById('tip');
function showTip(e, html){ tip.innerHTML = html; tip.style.opacity = 1;
  const x = Math.min(e.clientX + 14, window.innerWidth - 360);
  tip.style.left = x + 'px'; tip.style.top = (e.clientY + 16) + 'px'; }
function hideTip(){ tip.style.opacity = 0; }

document.querySelectorAll('.tangle').forEach(el => {
  el.style.cursor = 'pointer';
  el.addEventListener('mousemove', e => showTip(e,
    '<b>' + el.dataset.type.replace(/_/g,' ') + '</b><br>' + el.dataset.desc));
  el.addEventListener('mouseleave', hideTip);
  el.addEventListener('click', () => select(el.dataset.id));
});
document.querySelectorAll('.chrom').forEach(el => {
  el.addEventListener('mousemove', e => showTip(e, '<b>' + el.dataset.name + '</b><br>' +
    el.dataset.role + ', ' + Number(el.dataset.length).toLocaleString() + ' bp' +
    (el.dataset.depth ? '<br>depth ' + el.dataset.depth + 'x' : '')));
  el.addEventListener('mouseleave', hideTip);
});
document.querySelectorAll('.annot').forEach(el => {
  el.addEventListener('mousemove', e => showTip(e, el.dataset.desc));
  el.addEventListener('mouseleave', hideTip);
});

let current = null;
function select(id){
  current = (current === id) ? null : id;
  document.querySelectorAll('.tangle').forEach(el => {
    el.classList.toggle('dim', current !== null && el.dataset.id !== current); });
  document.querySelectorAll('.card').forEach(c => {
    c.classList.toggle('on', c.dataset.id === current); });
  if (current){ const c = document.querySelector('.card[data-id="'+current+'"]');
    if (c) c.scrollIntoView({block:'nearest', behavior:'smooth'}); }
}
document.querySelectorAll('.card').forEach(c =>
  c.addEventListener('click', () => select(c.dataset.id)));

function toggle(id, sel){ document.getElementById(id).addEventListener('change', e => {
  document.querySelectorAll(sel).forEach(el => el.style.display = e.target.checked ? '' : 'none');
}); }
toggle('c-tangles', '#layer-tangles');
toggle('c-coverage', '#layer-coverage');
toggle('c-annot', '.annot');
toggle('c-legend', '#legend');

const svg = document.querySelector('#svgbox svg');
const baseW = svg ? svg.getAttribute('width') : 0;
document.getElementById('zoom').addEventListener('input', e => {
  const z = e.target.value; document.getElementById('zv').textContent = z + '%';
  if (svg){ svg.style.width = (baseW * z / 100) + 'px'; svg.style.height = 'auto'; }
});
</script></body></html>
"""


def render_html(model: Model) -> str:
    svg = render_svg(model, interactive=True)
    side: List[str] = []

    if model.warnings:
        side.append("<h2>Warnings</h2>")
        for w in model.warnings:
            side.append(f'<div class="warn">{esc(w)}</div>')

    side.append("<h2>Karyotype calls</h2><table>")
    side.append("<tr><th>Sequence</th><th>Length</th><th>Call</th><th>Confidence</th></tr>")
    for s in model.drawable() + model.unplaced()[:15]:
        side.append(
            f"<tr><td>{esc(s.display)}</td><td>{human_bp(s.length)}</td>"
            f"<td>{esc(s.role)}</td><td>{esc(_confidence(s))}</td></tr>"
        )
    side.append("</table>")

    for s in model.drawable():
        if not s.evidence:
            continue
        side.append(
            f"<details><summary>Why {esc(s.display)} was called {esc(s.role)}</summary><ul>"
            + "".join(f"<li>{esc(e.as_text())}</li>" for e in s.evidence)
            + "</ul></details>"
        )

    side.append(f"<h2>Graph features ({len(model.tangles)})</h2>")
    if not model.tangles:
        side.append('<div class="sub">No tangles detected, or no assembly graph supplied.</div>')
    for t in model.tangles:
        colour = TANGLE_STYLE.get(t.type, ("#888", ""))[0]
        mult = (
            f'<span class="pill">~{t.multiplicity:g} copies</span>'
            if t.multiplicity
            else ""
        )
        side.append(
            f'<div class="card" data-id="{esc(t.id)}">'
            f'<div class="t"><span class="swatch" style="background:{colour}"></span>'
            f"{esc(TANGLE_LABEL.get(t.type, t.type))}{mult}</div>"
            f'<div class="d">{esc(t.description)}</div>'
            f'<div class="why">on {esc(", ".join(t.sequences)) or "unplaced"}'
            + (f" &middot; {esc('; '.join(t.evidence))}" if t.evidence else "")
            + "</div></div>"
        )

    if model.coverage_anomalies:
        side.append(f"<h2>Coverage outliers ({len(model.coverage_anomalies)})</h2><table>")
        side.append("<tr><th>Region</th><th>Type</th><th>vs median</th></tr>")
        for a in sorted(model.coverage_anomalies, key=lambda a: -abs(a["peak"] - 1))[:40]:
            side.append(
                f"<tr><td>{esc(a['seqname'])}:{a['start']:,}-{a['end']:,}</td>"
                f"<td>{esc(a['kind'])}</td><td>{a['peak']:.1f}x</td></tr>"
            )
        side.append("</table>")

    return (
        HTML_TEMPLATE.replace("__TITLE__", esc(model.title))
        .replace("__SUMMARY__", esc(model.summary_sentence()))
        .replace("__SVG__", svg)
        .replace("__SIDE__", "\n".join(side))
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
        L.append(
            "No expected karyotype was used to reach that. A finished linear chromosome carries "
            "a telomere repeat array at each end, so the count follows from how many ends are "
            "capped; a join between two contigs is only asserted when the segment bridging them "
            "is present in roughly one copy and touches nothing else."
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


# ==========================================================================
# GRAPH-FIRST MODE
#
# Everything below works from the assembler's own output only: a GFA and,
# optionally, Flye's assembly_info.txt. No chromosome-level assembly is
# assumed. It produces ranked hypotheses about which segments form which
# chromosome, and is deliberate about the difference between an observation,
# a derived estimate, and a hypothesis.
# ==========================================================================

# Palette carried over from the manual worked example.
CLASS_COLOUR = {
    "backbone": "#4a7ba7",
    "repeat": "#c2703d",
    "at_rich": "#b0a04a",
    "tandem_array": "#b5487f",
    "organelle_candidate": "#3f7d5c",
    "low_coverage": "#a8a8a8",
    "short_single_copy": "#8fa8bd",
    "unclassified": "#cfcfcf",
}
BACKBONE_COLOURS = ["#4a7ba7", "#5f9e6e", "#a87f4a", "#7d6ba7", "#3f8f9e"]

# The figure palette (v9 design). One maximally distinct colour per segment,
# reused in both panels, so a node in the graph can be found on a chromosome by
# colour alone. Class colours are NOT used in the figures any more; CLASS_COLOUR
# survives only for the text report.
SEGMENT_COLOURS = [
    "#e6194b",  # red
    "#4363d8",  # blue
    "#3cb44b",  # green
    "#f58231",  # orange
    "#911eb4",  # purple
    "#42d4f4",  # cyan
    "#f032e6",  # magenta
    "#9a6324",  # brown
    "#7ba428",  # olive
    "#ffe119",  # yellow
]
UNPLACED_GREY = "#b8b8b8"

# Type scale for the figures (v9): title, panel heading, label, sub-label.
FS_TITLE, FS_HEADING, FS_LABEL, FS_SUB = 40, 32, 27, 18

CLASS_LABEL = {
    "backbone": "single-copy backbone",
    "repeat": "repeat",
    "at_rich": "AT-rich",
    "tandem_array": "tandem array",
    "organelle_candidate": "organelle candidate",
    "low_coverage": "low coverage / foreign",
    "short_single_copy": "short, single copy",
    "unclassified": "unclassified",
}

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


def assign_segment_colours(calls: List["SegmentCall"]) -> Dict[str, str]:
    """
    One colour per segment, used by BOTH figures so a segment can be traced from
    the assembly graph to the chromosome it ends up in.

    EVERY segment gets its own maximally distinct colour, not a colour for its
    inferred class: the figure's job is correspondence between the two panels,
    and a shared class colour destroys that. Assignment is by descending length,
    so it is stable between runs. Once the palette is exhausted it is cycled and
    lightened, which keeps later segments distinguishable from earlier ones.
    """
    colours: Dict[str, str] = {}
    ordered = sorted(calls, key=lambda c: (-c.length, c.name))
    for i, c in enumerate(ordered):
        base = SEGMENT_COLOURS[i % len(SEGMENT_COLOURS)]
        cycle = i // len(SEGMENT_COLOURS)
        colours[c.name] = base if cycle == 0 else _lighten(base, 0.26 * cycle)
    return colours


def _text_on(hex_colour: str) -> str:
    """Black or white, whichever stays legible on the given fill."""
    try:
        r, g, b = (int(hex_colour[i : i + 2], 16) for i in (1, 3, 5))
    except (ValueError, IndexError):
        return "#1a1a1a"
    # sRGB relative luminance, good enough for picking ink
    lum = (0.2126 * r + 0.7152 * g + 0.0722 * b) / 255.0
    return "#1a1a1a" if lum > 0.55 else "#ffffff"


def _segment_number(name: str) -> str:
    """The digits in a segment name: 'edge_11' -> '11'. Falls back to the name."""
    m = re.findall(r"\d+", name)
    return m[-1] if m else name


def _lighten(hex_colour: str, amount: float) -> str:
    amount = min(max(amount, 0.0), 0.8)
    r, g, b = (int(hex_colour[i : i + 2], 16) for i in (1, 3, 5))
    r, g, b = (int(v + (255 - v) * amount) for v in (r, g, b))
    return f"#{r:02x}{g:02x}{b:02x}"


def revcomp(s: str) -> str:
    return s.upper()[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


def parse_size(text: str) -> int:
    """Accepts 36563796, 36.5m, 36.5Mb, 1.2g, 900k."""
    t = str(text).strip().lower().replace(",", "").replace("b", "")
    mult = 1
    if t.endswith("k"):
        mult, t = 1_000, t[:-1]
    elif t.endswith("m"):
        mult, t = 1_000_000, t[:-1]
    elif t.endswith("g"):
        mult, t = 1_000_000_000, t[:-1]
    return int(float(t) * mult)


@dataclass
class ContigInfo:
    """One row of Flye's assembly_info.txt."""

    name: str
    length: int
    cov: Optional[float]
    circular: bool
    repeat: bool
    mult: Optional[float]
    alt_group: str
    path: List[str]  # ordered element ids; "*" marks a dead end


def parse_flye_info(path: str, log: Log) -> List[ContigInfo]:
    """
    Flye assembly_info.txt. Columns are located from the header comment because
    the layout has changed between Flye versions (the 'telomere' column comes
    and goes). graph_path may repeat an id within one path, and that repetition
    is meaningful, so it is preserved.
    """
    rows: List[ContigInfo] = []
    header: List[str] = []
    with smart_open(path) as fh:
        for line_no, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line.strip():
                continue
            if line.startswith("#"):
                header = [c.strip().lower().rstrip(".") for c in line.lstrip("#").split("\t")]
                continue
            f = line.split("\t")
            if not header:
                log.warn(f"{path}:{line_no}: no header line seen; cannot parse assembly_info")
                break
            row = dict(zip(header, f))

            def get(*names, default=""):
                for n in names:
                    if row.get(n):
                        return row[n]
                return default

            try:
                length = int(get("length", default="0"))
            except ValueError:
                continue
            cov_raw = get("cov", "coverage")
            mult_raw = get("mult", "multiplicity")
            rows.append(
                ContigInfo(
                    name=get("seq_name", "#seq_name", "name"),
                    length=length,
                    cov=_maybe_float(cov_raw),
                    circular=get("circ", "circular").upper().startswith(("Y", "+", "T")),
                    repeat=get("repeat").upper().startswith(("Y", "+", "T")),
                    mult=_maybe_float(mult_raw),
                    alt_group=get("alt_group", default="*"),
                    path=[p.strip() for p in get("graph_path", "path").split(",") if p.strip()],
                )
            )
    if not rows:
        raise ValueError(f"{path}: no contig rows parsed - is this a Flye assembly_info.txt?")
    log.info(f"assembly_info: {len(rows)} contigs")
    return rows


def _maybe_float(text: str) -> Optional[float]:
    try:
        return float(text)
    except (TypeError, ValueError):
        return None


def resolve_path_element(element: str, segs: Dict[str, GfaSegment]) -> Optional[str]:
    """
    Map one graph_path element onto a GFA segment name. Flye writes signed
    integers; the GFA calls the same thing edge_N. The sign is orientation, not
    identity. '*' and '??' are dead-end markers, not segments.
    """
    e = element.strip()
    if e in ("*", "??", "", "?"):
        return None
    e = e.lstrip("+-")
    for candidate in (f"edge_{e}", e, f"utg{e}", f"contig_{e}"):
        if candidate in segs:
            return candidate
    return None


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
    w = min(max(window, 1), n)
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
            c.cls = "low_coverage"
            r.append(
                f"copy number {cn:.2f} (< {args.low_coverage_max_copy}): below single copy, so "
                f"not simply a low-confidence unique segment - candidates include a haplotype-"
                f"specific region, contamination, or a real sub-stoichiometric molecule"
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
            adj, telomeric, args, end_adj
        )
        spec = [j for j in chosen if j.speculative]
        if spec:
            # not established by the graph, so it must not be allowed to win on
            # score alone; it stays in the list, clearly labelled
            score -= args.speculative_penalty * len(spec)
            for j in spec:
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
    adj, telomeric, args, end_adj=None
):
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
    floor = segment_thickness() * 1.6
    px_per_bp = getattr(args, "graph_px_per_bp", None)
    if px_per_bp:
        # Same bases-per-pixel as the chromosome panel, so a contig is the same
        # size in both. Anything too short to draw at that scale is clamped to
        # the floor and is therefore drawn LARGER than true - unavoidable if a
        # 2 kb segment is to be visible beside a 9 Mb one, but it only ever
        # overstates the small ones.
        return max(length * px_per_bp, floor)
    return min(max(10.0 + args.graph_length_scale * math.sqrt(max(length, 1)), floor),
               args.graph_max_segment_px)


def segment_thickness(depth: Optional[float] = None) -> float:
    """
    Uniform, and equal to the chromosome bar width (v9 design). Thickness used to
    track read depth, but that put a second variable into the ribbon width and
    made the two panels hard to match up; depth is now carried by the label only.
    """
    return float(BAR_W)


def bandage_layout(
    calls: List[SegmentCall], links: List[GfaLink], args, log: Log
) -> Tuple[Dict[str, Tuple[float, float, float, float]], float, float]:
    """
    Force-directed layout over segment ENDS, not segment centres, which is what
    gives the Bandage look: each segment is a stiff spring of its own drawn
    length, and each link is a short spring tying one segment's end to another's.

    Returns {segment: (x1, y1, x2, y2)} plus the canvas size.
    """
    by_name = {c.name: c for c in calls}
    names = sorted(by_name)
    if not names:
        return {}, 100.0, 100.0

    # two point masses per segment: its start (+) and its end (-)
    pts: List[str] = []
    for n in names:
        pts += [n + "\x00s", n + "\x00e"]
        idx = {p: i for i, p in enumerate(pts)}

    n_pts = len(pts)
    rnd = _Rand(20260809)
    # deterministic ring start, jittered, so components unfold rather than
    # starting on top of one another
    radius = 40.0 + 9.0 * math.sqrt(n_pts)
    pos = []
    for i in range(n_pts):
        a = 2 * math.pi * i / n_pts
        pos.append([
            radius * math.cos(a) + rnd.uniform(-8, 8),
            radius * math.sin(a) + rnd.uniform(-8, 8),
        ])

    springs: List[Tuple[int, int, float, float]] = []  # a, b, rest, strength
    for n in names:
        springs.append((idx[n + "\x00s"], idx[n + "\x00e"],
                        segment_draw_length(by_name[n].length, args), 1.0))
    for l in links:
        if l.a not in by_name or l.b not in by_name:
            continue
        # a link leaves the end of a + oriented segment and enters the start of
        # the next; a - orientation flips which terminal is involved
        a_pt = l.a + ("\x00e" if l.a_orient == "+" else "\x00s")
        b_pt = l.b + ("\x00s" if l.b_orient == "+" else "\x00e")
        if a_pt == b_pt:
            continue
        # A junction is given a RADIUS rather than being a single point. Pulling
        # every linked end onto one coordinate made four contigs meeting at the
        # same hub (edge_9 has five ends on one side) collapse into a pile you
        # could not read. Held about a ribbon-width apart, they spread into an
        # arc and each join shows as its own short connector.
        springs.append((idx[a_pt], idx[b_pt], segment_thickness() * 1.15, 4.0))

    k = max(radius / max(math.sqrt(n_pts), 1.0), 135.0)
    iters = int(min(500, max(80, 9000 / max(n_pts, 1))))
    if n_pts > 400:
        log.info(f"graph layout: {n_pts} endpoints, {iters} iterations (this can take a moment)")
    temp = radius * 0.35

    # Endpoint pairs joined by a link must be free to touch. Left in the
    # all-pairs repulsion they settle at k^2/d against the spring, which for
    # k=46 parks them 40-135 px apart and turns every junction into a long bar.
    linked_pairs = {
        (min(a, b), max(a, b)) for a, b, rest, _s in springs
        if rest <= segment_thickness() * 1.2
    }

    for step in range(iters):
        disp = [[0.0, 0.0] for _ in range(n_pts)]
        # repulsion, all pairs except those a link is trying to hold together
        for i in range(n_pts):
            xi, yi = pos[i]
            for j in range(i + 1, n_pts):
                if (i, j) in linked_pairs:
                    continue
                dx = xi - pos[j][0]
                dy = yi - pos[j][1]
                d2 = dx * dx + dy * dy
                if d2 < 1e-6:
                    dx, dy, d2 = rnd.uniform(-1, 1), rnd.uniform(-1, 1), 1.0
                d = math.sqrt(d2)
                f = (k * k) / d
                ux, uy = dx / d, dy / d
                disp[i][0] += ux * f
                disp[i][1] += uy * f
                disp[j][0] -= ux * f
                disp[j][1] -= uy * f
        # springs
        for a, b, rest, strength in springs:
            dx = pos[a][0] - pos[b][0]
            dy = pos[a][1] - pos[b][1]
            d = math.hypot(dx, dy) or 1e-6
            f = strength * (d - rest) * 0.9
            ux, uy = dx / d, dy / d
            disp[a][0] -= ux * f
            disp[a][1] -= uy * f
            disp[b][0] += ux * f
            disp[b][1] += uy * f
        # pull to the centre so detached components (an unplaced contig, an
        # organelle) stay near the main mass instead of stranding themselves in a
        # far corner and stretching the canvas around a lot of white space
        for i in range(n_pts):
            disp[i][0] -= pos[i][0] * 0.038
            disp[i][1] -= pos[i][1] * 0.038
        # move, capped by the cooling temperature
        for i in range(n_pts):
            dx, dy = disp[i]
            d = math.hypot(dx, dy) or 1e-6
            lim = min(d, temp)
            pos[i][0] += dx / d * lim
            pos[i][1] += dy / d * lim
        temp = max(temp * 0.965, 0.6)

    raw = {
        n: (
            pos[idx[n + "\x00s"]][0], pos[idx[n + "\x00s"]][1],
            pos[idx[n + "\x00e"]][0], pos[idx[n + "\x00e"]][1],
        )
        for n in names
    }

    # ---- pack the connected components ----
    # A spring model left to itself flings a detached contig or an organelle into
    # a far corner, and the canvas then has to stretch around all that white
    # space. The components are laid out independently and then packed: the
    # largest keeps its position, the rest are set apart in a row beneath it.
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
        cx = [v for n in g for v in (raw[n][0], raw[n][2])]
        cy = [v for n in g for v in (raw[n][1], raw[n][3])]
        return min(cx), min(cy), max(cx), max(cy)

    thick = segment_thickness()
    groups = sorted(comps.values(), key=lambda g: -sum(by_name[n].length for n in g))
    placed: Dict[str, Tuple[float, float, float, float]] = {n: raw[n] for n in groups[0]}
    mx0, _my0, _mx1, my1 = bbox(groups[0])
    gap = thick * 3.0
    cur_x, row_y = mx0, my1 + gap
    for g in groups[1:]:
        gx0, gy0, gx1, _gy1 = bbox(g)
        dx, dy = cur_x - gx0, row_y - gy0
        for n in g:
            x1, y1, x2, y2 = raw[n]
            placed[n] = (x1 + dx, y1 + dy, x2 + dx, y2 + dy)
        cur_x += (gx1 - gx0) + gap

    # Orient the graph so its long axis is horizontal. A spring layout comes out
    # in an arbitrary rotation, and a tall one forces a tall figure: the
    # chromosome panel beside it is much shorter, so most of the canvas ends up
    # empty and everything has to be shrunk to fit a preview.
    xs0 = [v for t in placed.values() for v in (t[0], t[2])]
    ys0 = [v for t in placed.values() for v in (t[1], t[3])]
    if xs0 and (max(ys0) - min(ys0)) > (max(xs0) - min(xs0)):
        placed = {
            n: (y1, -x1, y2, -x2) for n, (x1, y1, x2, y2) in placed.items()
        }

    # normalise into a padded canvas
    xs = [v for t in placed.values() for v in (t[0], t[2])]
    ys = [v for t in placed.values() for v in (t[1], t[3])]
    pad = 18.0 + thick * 0.5
    minx, maxx, miny, maxy = min(xs), max(xs), min(ys), max(ys)
    width = (maxx - minx) + 2 * pad
    height = (maxy - miny) + 2 * pad
    out: Dict[str, Tuple[float, float, float, float]] = {}
    for n in names:
        x1, y1, x2, y2 = placed[n]
        out[n] = (
            x1 - minx + pad, y1 - miny + pad,
            x2 - minx + pad, y2 - miny + pad,
        )
    return out, width, height


def render_bandage_style_svg(
    calls: List[SegmentCall],
    links: List[GfaLink],
    title: str,
    colours: Dict[str, str],
    args,
    log: Log,
) -> str:
    """The graph drawn Bandage-fashion: thick tapered segments, force-directed."""
    by_name = {c.name: c for c in calls}
    geom, width, height = bandage_layout(calls, links, args, log)
    colours = colours or assign_segment_colours(calls)

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

    # Parallel segments - two contigs running between the same pair of junctions -
    # land on almost the same chord and the second one disappears underneath the
    # first. Bucket by the endpoints they share and fan the bow out, alternating
    # sign, so each is visible. Without this, edge_7 hides entirely behind edge_2.
    # Detected from the GRAPH, not from the drawn coordinates: two segments are
    # parallel when they have the same set of neighbours. Bucketing by pixel
    # position looked simpler but is far too brittle - edge_2 and edge_7 land in
    # adjacent buckets and stack anyway.
    neighbours_of: Dict[str, Set[str]] = defaultdict(set)
    for l in links:
        if l.a == l.b:
            continue
        neighbours_of[l.a].add(l.b)
        neighbours_of[l.b].add(l.a)

    parallel: Dict[Tuple[str, ...], List[str]] = defaultdict(list)
    for name in sorted(geom):
        nb = neighbours_of.get(name, set())
        if len(nb) >= 2:
            parallel[tuple(sorted(nb))].append(name)
    bow_slot: Dict[str, int] = {}
    label_t: Dict[str, float] = {}
    for group in parallel.values():
        if len(group) < 2:
            continue
        # 0, +1, -1, +2, -2 ... so the bundle spreads either side of the chord
        for i, nm in enumerate(sorted(group)):
            bow_slot[nm] = ((i + 1) // 2) * (1 if i % 2 else -1)
            # and stagger the labels ALONG the ribbons. Fanning separates the
            # middles but the ends still converge, so labels placed at the same
            # fraction of two bundled segments collide however wide the fan.
            label_t[nm] = 0.30 + 0.40 * (i / max(len(group) - 1, 1))

    # Segments carrying a self-link are circular molecules and are drawn as rings
    # rather than as ribbons with a loop hanging off one end.
    circular = {l.a for l in links if l.a == l.b}
    w = segment_thickness()

    # Junction stubs, behind the segments. Linked ends already abut after layout,
    # so a link is a short dark connector rather than a long thin line.
    P.append(
        f'<g id="layer-links" fill="none" stroke="{PALETTE["bar_edge"]}" '
        f'stroke-linecap="round">'
    )
    for l in links:
        if l.a not in geom or l.b not in geom or l.a == l.b:
            continue
        ax, ay = (geom[l.a][2], geom[l.a][3]) if l.a_orient == "+" else (geom[l.a][0], geom[l.a][1])
        bx, by = (geom[l.b][0], geom[l.b][1]) if l.b_orient == "+" else (geom[l.b][2], geom[l.b][3])
        # Thin. The ribbons are drawn with ROUND ends, so two of them meeting at
        # an angle no longer leave a white wedge at the corner and the connector
        # does not have to be wide enough to cover one. A junction should read as
        # a join, not as another piece of sequence.
        P.append(
            f'<line x1="{ax:.1f}" y1="{ay:.1f}" x2="{bx:.1f}" y2="{by:.1f}" '
            f'stroke-width="{w * 0.26:.1f}"/>'
        )
    P.append("</g>")

    # segments as thick ribbons of uniform width
    P.append('<g id="layer-segments" fill="none">')
    rings: Dict[str, Tuple[float, float, float]] = {}
    for name, (x1, y1, x2, y2) in sorted(geom.items(), key=lambda kv: -by_name[kv[0]].length):
        c = by_name[name]
        colour = colours.get(name, "#cfcfcf")
        if name in circular:
            # a ring whose circumference matches the drawn length of the segment
            seg_len = max(segment_draw_length(c.length, args), 30.0)
            r = max(seg_len / (2 * math.pi), 13.0)
            cx, cy = (x1 + x2) / 2, (y1 + y2) / 2
            rings[name] = (cx, cy, r)
            P.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r:.1f}" '
                     f'stroke="{PALETTE["bar_edge"]}" stroke-width="{w + 2.0:.1f}"/>')
            P.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r:.1f}" '
                     f'stroke="{colour}" stroke-width="{w:.1f}"/>')
            continue
        # a gentle bow, so nothing looks like a ruler
        mx, my = (x1 + x2) / 2, (y1 + y2) / 2
        nx, ny = -(y2 - y1), (x2 - x1)
        nlen = math.hypot(nx, ny) or 1.0
        bow = min(26.0, math.hypot(x2 - x1, y2 - y1) * 0.14)
        slot = bow_slot.get(name, 0)
        if slot:
            # fan a bundle of parallel segments apart rather than stacking them
            bow = slot * max(
                abs(bow),
                min(math.hypot(x2 - x1, y2 - y1) * 0.42, segment_thickness() * 6.0),
            )
        cx, cy = mx + nx / nlen * bow, my + ny / nlen * bow
        d = f"M {x1:.1f} {y1:.1f} Q {cx:.1f} {cy:.1f} {x2:.1f} {y2:.1f}"
        # dark casing then the colour, flat butt caps so segments abut cleanly
        P.append(f'<path d="{d}" stroke="{PALETTE["bar_edge"]}" stroke-width="{w + 2.0:.1f}" '
                 f'stroke-linecap="round"/>')
        P.append(f'<path d="{d}" stroke="{colour}" stroke-width="{w:.1f}" '
                 f'stroke-linecap="round"/>')
    P.append("</g>")

    # labels: the segment number sits INSIDE the ribbon, coverage beside it
    label_all = len(calls) <= args.graph_label_limit
    P.append(f'<g id="layer-labels" fill="{PALETTE["text"]}">')
    for name, (x1, y1, x2, y2) in geom.items():
        c = by_name[name]
        if not label_all and c.cls == "backbone" and c.length < args.backbone_min_length:
            continue
        colour = colours.get(name, "#cfcfcf")
        if name in rings:
            cx, cy, r = rings[name]
            nx, ny, mx, my = 0.0, -1.0, cx, cy - r
        else:
            cxm, cym = (x1 + x2) / 2, (y1 + y2) / 2
            nx, ny = -(y2 - y1), (x2 - x1)
            nlen = math.hypot(nx, ny) or 1.0
            nx, ny = nx / nlen, ny / nlen
            bow = min(26.0, math.hypot(x2 - x1, y2 - y1) * 0.14)
            slot = bow_slot.get(name, 0)
            if slot:
                bow = slot * max(
                    abs(bow),
                    min(math.hypot(x2 - x1, y2 - y1) * 0.42, segment_thickness() * 6.0),
                )
            # Evaluate the drawn quadratic Bezier at this segment's own t, and
            # take the normal from the tangent there. Bundled segments get
            # different t values so their labels never stack.
            qx, qy = cxm + nx * bow, cym + ny * bow  # the control point
            t = label_t.get(name, 0.5)
            u = 1.0 - t
            mx = u * u * x1 + 2 * u * t * qx + t * t * x2
            my = u * u * y1 + 2 * u * t * qy + t * t * y2
            tx = 2 * u * (qx - x1) + 2 * t * (x2 - qx)
            ty = 2 * u * (qy - y1) + 2 * t * (y2 - qy)
            tlen = math.hypot(tx, ty) or 1.0
            nx, ny = -ty / tlen, tx / tlen
        # number inside the ribbon, inked for contrast against its own colour
        P.append(
            f'<text x="{mx:.1f}" y="{my + FS_LABEL * 0.35:.1f}" text-anchor="middle" '
            f'font-size="{FS_LABEL}" font-weight="700" fill="{_text_on(colour)}">'
            f'{esc(_segment_number(name))}</text>'
        )
        # coverage only, set clear of the ribbon. Skipped on very short segments,
        # where there is no room for it to sit anywhere it would not collide.
        drawn_len = math.hypot(x2 - x1, y2 - y1)
        if c.depth is not None and (name in rings or drawn_len >= w * 2.2):
            off = w / 2 + FS_SUB + 8
            lx, ly = mx + nx * off, my + ny * off
            P.append(
                f'<text x="{lx:.1f}" y="{ly:.1f}" text-anchor="middle" '
                f'font-size="{FS_SUB}" fill="{PALETTE["muted"]}">{c.depth:.0f}x</text>'
            )
    P.append("</g>")

    P.append("</svg>")
    return "\n".join(P)


def _graph_legend_svg(calls: List[SegmentCall], height: float) -> str:
    out = [f'<g id="legend" font-size="10" fill="{PALETTE["text"]}">']
    x, y = 40.0, height - 24.0
    n_bb = sum(1 for c in calls if c.cls == "backbone")
    out.append(
        f'<text x="{x}" y="{y - 16:.1f}" fill="{PALETTE["muted"]}" font-size="9.5">'
        f'{n_bb} backbone segment(s), each its own colour and reused in the chromosome figure. '
        f'Other colours are by inferred class:</text>'
    )
    for cls_name in [c for c in CLASS_COLOUR if c != "backbone" and any(x2.cls == c for x2 in calls)]:
        out.append(f'<rect x="{x:.1f}" y="{y - 8:.1f}" width="11" height="11" rx="2" '
                   f'fill="{CLASS_COLOUR[cls_name]}"/>')
        out.append(f'<text x="{x + 16:.1f}" y="{y + 1:.1f}">{esc(CLASS_LABEL[cls_name])}</text>')
        x += 30 + 6.0 * len(CLASS_LABEL[cls_name])
    out.append("</g>")
    return "\n".join(out)


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
    colours = colours or assign_segment_colours(calls)
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


# ==========================================================================
# demo data
# ==========================================================================
def write_demo(out_dir: str, log: Log) -> Dict[str, str]:
    """
    Synthesises a small but structurally realistic dataset so the tool can be
    exercised without real inputs. Ground truth: 4 chromosomes, 1 circular
    high-depth mitogenome, 9 unplaced contigs, an rDNA array shared by chr1 and
    chr3, a transposon family on chr2/chr4, an inverted repeat on chr1, and a
    het bubble on chr2.
    """
    os.makedirs(out_dir, exist_ok=True)
    rnd = _Rand(20260806)

    chroms = [("chr1", 41_200_000), ("chr2", 33_800_000), ("chr3", 28_400_000), ("chr4", 19_600_000)]
    mito = ("MT", 16_540)
    unplaced = [(f"scaffold_{i:02d}", rnd.randint(18_000, 640_000)) for i in range(1, 10)]
    all_seqs = chroms + [mito] + unplaced

    fai = os.path.join(out_dir, "demo_assembly.fa.fai")
    with open(fai, "w") as fh:
        offset = 0
        for name, L in all_seqs:
            fh.write(f"{name}\t{L}\t{offset}\t60\t61\n")
            offset += L + len(name) + 2

    # ---- graph ----
    segs: List[Tuple[str, int, float]] = []
    links: List[Tuple[str, str, str, str]] = []
    placements: List[Tuple[str, str, int, int, str]] = []  # seg, chrom, start, end, strand

    for name, L in chroms:
        n_unique = 5
        step = L // n_unique
        prev = None
        for i in range(n_unique):
            sid = f"utg{name[-1]}{i:02d}"
            start, end = i * step, min((i + 1) * step, L)
            segs.append((sid, end - start, rnd.uniform(27.0, 33.0)))
            placements.append((sid, name, start, end, "+"))
            if prev:
                links.append((prev, "+", sid, "+"))
            prev = sid

    # rDNA array shared by chr1 and chr3 -> the classic inter-chromosomal join
    segs.append(("rDNA_array", 44_000, 168.0))
    placements.append(("rDNA_array", "chr1", 17_900_000, 17_944_000, "+"))
    placements.append(("rDNA_array", "chr3", 9_100_000, 9_144_000, "-"))
    links += [("utg102", "+", "rDNA_array", "+"), ("rDNA_array", "+", "utg103", "+"),
              ("utg301", "+", "rDNA_array", "-"), ("rDNA_array", "-", "utg302", "+")]

    # transposon family on chr2 and chr4
    segs.append(("TE_famA", 12_500, 96.0))
    placements.append(("TE_famA", "chr2", 24_300_000, 24_312_500, "+"))
    placements.append(("TE_famA", "chr4", 6_050_000, 6_062_500, "+"))
    links += [("utg201", "+", "TE_famA", "+"), ("TE_famA", "+", "utg202", "+"),
              ("utg400", "+", "TE_famA", "+"), ("TE_famA", "+", "utg401", "+")]

    # inverted repeat on chr1
    segs.append(("IR_chr1", 8_800, 58.0))
    placements.append(("IR_chr1", "chr1", 33_400_000, 33_408_800, "+"))
    links += [("utg104", "+", "IR_chr1", "+"), ("IR_chr1", "+", "IR_chr1", "-")]

    # heterozygous bubble on chr2
    for alt in ("het_a", "het_b"):
        segs.append((alt, 3_100, 15.5))
        placements.append((alt, "chr2", 12_000_000, 12_003_100, "+"))
        links += [("utg202", "+", alt, "+"), (alt, "+", "utg203", "+")]

    # circular mitogenome at very high depth
    segs.append(("mito_circ", mito[1], 1180.0))
    placements.append(("mito_circ", "MT", 0, mito[1], "+"))
    links.append(("mito_circ", "+", "mito_circ", "+"))

    # Short segments carry sequence so the composition screens and the BLAST
    # export can actually run; the multi-Mb ones use LN only, as real GFAs often
    # do. The bases are synthetic and mean nothing biologically.
    def fake_seq(n: int, gc: float, seed: int, unit: Optional[int] = None) -> str:
        r = _Rand(seed)
        def block(k: int) -> str:
            return "".join(
                ("GC" if r.uniform(0, 1) < gc else "AT")[r.randint(0, 1)] for _ in range(k)
            )
        if unit:
            core = block(unit)
            return (core * (n // unit + 1))[:n]
        return block(n)

    seq_recipe = {
        "rDNA_array": (0.55, 11, 4_400),
        "TE_famA": (0.42, 12, None),
        "IR_chr1": (0.44, 13, None),
        "het_a": (0.41, 14, None),
        "het_b": (0.41, 15, None),
        "mito_circ": (0.22, 16, None),
    }

    gfa = os.path.join(out_dir, "demo_assembly.gfa")
    with open(gfa, "w") as fh:
        fh.write("H\tVN:Z:1.0\n")
        for sid, L, dp in segs:
            if sid in seq_recipe:
                gc, seed, unit = seq_recipe[sid]
                seq = fake_seq(L, gc, seed, unit)
                if sid == "rDNA_array":  # give it something telomere-like to find
                    seq = "TTAGGG" * 60 + seq[360:]
                fh.write(f"S\t{sid}\t{seq}\tLN:i:{L}\tdp:f:{dp:.1f}\n")
            else:
                fh.write(f"S\t{sid}\t*\tLN:i:{L}\tdp:f:{dp:.1f}\n")
        for a, ao, b, bo in links:
            fh.write(f"L\t{a}\t{ao}\t{b}\t{bo}\t0M\n")

    paf = os.path.join(out_dir, "demo_segments_to_assembly.paf")
    seg_len = {s[0]: s[1] for s in segs}
    with open(paf, "w") as fh:
        tlen = dict(all_seqs)
        for sid, chrom, start, end, strand in placements:
            qlen = seg_len[sid]
            block = end - start
            fh.write(
                f"{sid}\t{qlen}\t0\t{qlen}\t{strand}\t{chrom}\t{tlen[chrom]}\t{start}\t{end}\t"
                f"{int(block * 0.985)}\t{block}\t60\n"
            )

    cov = os.path.join(out_dir, "demo_coverage.regions.bed")
    win = 200_000
    with open(cov, "w") as fh:
        for name, L in chroms:
            for start in range(0, L, win):
                end = min(start + win, L)
                d = rnd.uniform(27, 33)
                if name == "chr1" and 17_800_000 <= start <= 18_000_000:
                    d *= 5.4  # rDNA pile-up
                if name == "chr3" and 9_000_000 <= start <= 9_200_000:
                    d *= 5.1
                if name == "chr4" and 14_000_000 <= start <= 14_600_000:
                    d *= 0.36  # putative deletion / haploid region
                fh.write(f"{name}\t{start}\t{end}\t{d:.2f}\n")
        fh.write(f"MT\t0\t{mito[1]}\t1180.00\n")

    bed = os.path.join(out_dir, "demo_features.bed")
    with open(bed, "w") as fh:
        fh.write("chr1\t17900000\t17944000\trDNA_array\t0\t.\trDNA\n")
        fh.write("chr3\t9100000\t9144000\trDNA_array\t0\t.\trDNA\n")
        fh.write("chr1\t20100000\t20400000\tcentromere\t0\t.\tcentromere\n")
        fh.write("chr2\t16200000\t16500000\tcentromere\t0\t.\tcentromere\n")
        fh.write("chr3\t13800000\t14050000\tcentromere\t0\t.\tcentromere\n")
        fh.write("chr4\t9200000\t9450000\tcentromere\t0\t.\tcentromere\n")
        for name, L in chroms:
            fh.write(f"{name}\t0\t18000\ttelomere\t0\t.\ttelomere\n")
            fh.write(f"{name}\t{L - 18000}\t{L}\ttelomere\t0\t.\ttelomere\n")

    log.info(f"demo data written to {out_dir}")
    return {"fai": fai, "gfa": gfa, "paf": paf, "coverage": cov, "annotation": bed}


class _Rand:
    """Tiny deterministic LCG so demo output is byte-stable across Pythons."""

    def __init__(self, seed: int):
        self.s = seed & 0xFFFFFFFF

    def _next(self) -> int:
        self.s = (1103515245 * self.s + 12345) & 0x7FFFFFFF
        return self.s

    def uniform(self, a: float, b: float) -> float:
        return a + (b - a) * (self._next() / 0x7FFFFFFF)

    def randint(self, a: int, b: int) -> int:
        return a + self._next() % (b - a + 1)


# ==========================================================================
# pipeline
# ==========================================================================
PAIR_GUTTER = 96.0


def render_paired_svg(
    model: Model,
    calls: List[SegmentCall],
    links: List[GfaLink],
    colours: Dict[str, str],
    args,
    log: Log,
) -> str:
    """
    One figure: the assembly graph on the left, the chromosomes it resolves into
    on the right, and faint leader lines joining a graph node to the block it
    became. Both panels are nested <svg> elements, which keeps each renderer
    independent and avoids the two fighting over coordinates.
    """
    adj = build_adjacency(links)
    pos, gw, gh, _ = _graph_layout(calls, adj)
    # the combined figure carries the title, so the panels must not repeat it
    real_title = model.title
    model.title = "Hypothesis of chromosome structure"
    try:
        ideo_svg = render_svg(model)
        lay, _, _, _ = ideogram_geometry(model)
        anchors = ideogram_block_anchors(model)
    finally:
        model.title = real_title

    # Build the chromosome panel FIRST so the graph panel can borrow its scale:
    # a contig should be the same size in both halves of the figure.
    args.graph_px_per_bp = lay.scale
    graph_svg = graph_svg_for_style(calls, links, "", colours, args, log)
    if args.graph_style == "bandage":
        gw, gh = _svg_width(graph_svg), _svg_height(graph_svg)

    iw, ih = lay.width, _svg_height(ideo_svg)

    # ---- left panel: a real Bandage export if given, otherwise our redraw ----
    external = args.bandage_image and os.path.exists(args.bandage_image)
    if args.bandage_image and not external:
        log.warn(f"--bandage-image {args.bandage_image} not found; using our own redraw instead")
    if external:
        src = image_size(args.bandage_image)
        if src is None:
            log.warn(
                f"could not read the dimensions of {args.bandage_image}; expected PNG, JPEG or "
                f"SVG. Using our own redraw instead."
            )
            external = False
    if external:
        # scale the Bandage export to the height of the chromosome panel
        panel_h = ih - 70
        panel_w = panel_h * (src[0] / src[1])
        if panel_w > args.bandage_max_width:
            panel_w = args.bandage_max_width
            panel_h = panel_w * (src[1] / src[0])
        gw_eff, gh_eff = panel_w, panel_h
        left_label = "Assembly graph, as drawn by Bandage"
    elif args.rotate_graph:
        # rotating our redraw a quarter turn trades a very wide figure for a
        # taller, narrower one that fits a page or a slide
        gw_eff, gh_eff = gh, gw
        left_label = "Assembly graph"
    else:
        gw_eff, gh_eff = gw, gh
        left_label = "Assembly graph"

    ox = gw_eff + PAIR_GUTTER  # x offset of the ideogram panel
    width = ox + iw
    top = 84.0
    height = max(gh_eff, ih) + top

    P = [
        f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
        f'width="{width:.0f}" height="{height:.0f}" '
        f'viewBox="0 0 {width:.0f} {height:.0f}" font-family="Helvetica, Arial, sans-serif">',
        f'<rect width="100%" height="100%" fill="{PALETTE["bg"]}"/>',
        f'<text x="60" y="44" font-size="{FS_TITLE}" font-weight="600" '
        f'fill="{PALETTE["text"]}">{esc(model.title)}</text>',
        f'<text x="60" y="{top - 10:.0f}" font-size="{FS_HEADING}" font-weight="600" '
        f'fill="{PALETTE["text"]}">{esc(left_label)}</text>',
    ]

    if external:
        P.append(embed_image(args.bandage_image, 0, top, gw_eff, gh_eff))
    elif args.rotate_graph:
        P.append(_place_svg(graph_svg, 0, top, rotate=-90))
    else:
        P.append(_place_svg(graph_svg, 0, top))
    P.append(_place_svg(ideo_svg, ox, top))
    P.append("</svg>")
    return "\n".join(P)


def _svg_height(svg: str) -> float:
    m = re.search(r'<svg[^>]*\bheight="([0-9.]+)"', svg)
    return float(m.group(1)) if m else 800.0


def _svg_width(svg: str) -> float:
    m = re.search(r'<svg[^>]*\bwidth="([0-9.]+)"', svg)
    return float(m.group(1)) if m else 800.0


def image_size(path: str) -> Optional[Tuple[float, float]]:
    """Pixel dimensions of a PNG, JPEG or SVG, without needing an image library."""
    import struct

    with open(path, "rb") as fh:
        head = fh.read(4096)
    if head[:8] == b"\x89PNG\r\n\x1a\n":
        w, h = struct.unpack(">II", head[16:24])
        return float(w), float(h)
    if head[:2] == b"\xff\xd8":  # JPEG: walk the segment markers for SOFn
        with open(path, "rb") as fh:
            fh.read(2)
            while True:
                b = fh.read(1)
                while b and b != b"\xff":
                    b = fh.read(1)
                marker = fh.read(1)
                while marker == b"\xff":
                    marker = fh.read(1)
                if not marker:
                    return None
                if marker[0] in range(0xC0, 0xCF) and marker[0] not in (0xC4, 0xC8, 0xCC):
                    fh.read(3)
                    h, w = struct.unpack(">HH", fh.read(4))
                    return float(w), float(h)
                size = struct.unpack(">H", fh.read(2))[0]
                fh.read(size - 2)
    text = head.decode("utf-8", "replace")
    if "<svg" in text:
        w = re.search(r'<svg[^>]*\bwidth="([0-9.]+)', text)
        h = re.search(r'<svg[^>]*\bheight="([0-9.]+)', text)
        if w and h:
            return float(w.group(1)), float(h.group(1))
        vb = re.search(r'viewBox="[\s0-9.\-]*?([0-9.]+)\s+([0-9.]+)"', text)
        if vb:
            return float(vb.group(1)), float(vb.group(2))
    return None


def embed_image(path: str, x: float, y: float, w: float, h: float) -> str:
    """
    Place an external image as a panel, inlined as base64 so the figure stays a
    single portable file. An SVG export keeps its vectors; a PNG or JPEG is
    embedded as-is.
    """
    import base64
    import mimetypes

    if path.lower().endswith(".svg"):
        with open(path) as fh:
            svg = fh.read()
        body = svg.split(">", 1)[1].rsplit("</svg>", 1)[0]
        body = _BG_RECT_RE.sub("", body, count=1)
        src = image_size(path) or (w, h)
        scale = min(w / src[0], h / src[1])
        return (
            f'<g transform="translate({x:.1f},{y:.1f}) scale({scale:.4f})">{body}</g>'
        )
    mime = mimetypes.guess_type(path)[0] or "image/png"
    with open(path, "rb") as fh:
        data = base64.b64encode(fh.read()).decode("ascii")
    return (
        f'<image x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" '
        f'preserveAspectRatio="xMidYMin meet" xlink:href="data:{mime};base64,{data}"/>'
    )


def write_bandage_colour_csv(
    calls: List[SegmentCall], colours: Dict[str, str], path: str, log: Log
) -> str:
    """
    A CSV Bandage can load to colour the graph the way we do, so a real Bandage
    render and our chromosome figure agree.

    Bandage keys rows by segment name in the first column and recognises a
    colour column; the remaining columns show up as node labels. If the colours
    do not take, check the column name against Bandage's own colour-schemes
    documentation for your version rather than assuming this header is right.
    """
    import csv

    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Name", "Colour", "Color", "Class", "CopyNumber", "Depth", "Length"])
        for c in sorted(calls, key=lambda c: -c.length):
            colour = colours.get(c.name, "#cfcfcf")
            w.writerow([
                c.name, colour, colour, CLASS_LABEL.get(c.cls, c.cls),
                f"{c.copy_number:.2f}" if c.copy_number is not None else "",
                f"{c.depth:.1f}" if c.depth is not None else "",
                c.length,
            ])
    log.info(f"wrote Bandage colour CSV for {len(calls)} segments to {path}")
    return path


_BG_RECT_RE = re.compile(r'<rect width="100%" height="100%"[^/]*/>')


def _place_svg(svg: str, x: float, y: float, rotate: float = 0) -> str:
    """
    Re-position a complete SVG document as a panel of a larger figure.

    A translated <g> rather than a nested <svg>: nested-SVG positioning is
    another thing renderers disagree about - ImageMagick ignores the x/y and
    stacks both panels at the origin. The panel's own background rect is dropped
    so it cannot paint over its neighbour or over the leader lines.
    """
    body = svg.split(">", 1)[1].rsplit("</svg>", 1)[0]
    body = _BG_RECT_RE.sub("", body, count=1)
    if rotate == -90:
        # a quarter turn anticlockwise: the panel's own width becomes the height
        # of its footprint, so shift down by that width to land back in view
        w = _svg_width(svg)
        return f'<g transform="translate({x:.1f},{y + w:.1f}) rotate(-90)">{body}</g>'
    return f'<g transform="translate({x:.1f},{y:.1f})">{body}</g>'


def render_graph_figure(
    calls: List[SegmentCall],
    links: List[GfaLink],
    colours: Dict[str, str],
    title: str,
    path: str,
    args,
    log: Log,
) -> Optional[str]:
    """The Bandage graph redrawn with our own, fixed colours."""
    if len(calls) > args.max_graph_nodes:
        log.warn(
            f"the graph has {len(calls)} segments, more than --max-graph-nodes "
            f"({args.max_graph_nodes}), so the graph figure was skipped. It would be unreadable "
            f"at that size; raise the limit if you want it anyway."
        )
        return None
    with open(path, "w") as fh:
        fh.write(graph_svg_for_style(calls, links, title, colours, args, log))
    return path


def graph_svg_for_style(calls, links, title, colours, args, log) -> str:
    if args.graph_style == "bandage":
        return render_bandage_style_svg(calls, links, title, colours, args, log)
    return render_graph_svg(calls, links, title, colours)


def write_figures(
    model: Model,
    calls: List[SegmentCall],
    links: List[GfaLink],
    colours: Dict[str, str],
    base: str,
    args,
    log: Log,
) -> Dict[str, str]:
    """The graph panel, and the paired graph-plus-chromosomes figure."""
    out: Dict[str, str] = {}
    graph_path = render_graph_figure(
        calls, links, colours, model.title, base + "_graph.svg", args, log
    )
    if not graph_path:
        return out
    out["assembly graph figure"] = graph_path
    out["Bandage colour CSV (load this in Bandage)"] = write_bandage_colour_csv(
        calls, colours, base + "_bandage_colours.csv", log
    )
    paired = base + "_paired.svg"
    with open(paired, "w") as fh:
        fh.write(render_paired_svg(model, calls, links, colours, args, log))
    out["PAIRED FIGURE: graph and chromosomes"] = paired
    return out


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
        segs, links = parse_gfa(args.gfa, log)
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
        colours = assign_segment_colours(calls)
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
    segs, links = parse_gfa(args.gfa, log)
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
    else:
        pick = min(max(args.hypothesis, 1), len(hypotheses))
        model = model_from_hypothesis(hypotheses[pick - 1], calls, links, args, log)
        model.chosen_hypothesis = pick

    model.segment_calls = calls
    model.hypotheses = hypotheses
    model.baseline_depth = baseline
    model.baseline_basis = (
        f"median depth of segments at least {human_bp(args.baseline_min_length)} long"
    )
    model.inputs["gfa"] = args.gfa
    if args.flye_info:
        model.inputs["flye_info"] = args.flye_info

    colours = model.segment_colours or assign_segment_colours(calls)
    model.segment_colours = colours
    extra.update(
        identify_repeats(model, calls, model.tangles, adj, seq_by_segment, base, args, log)
    )
    extra.update(write_figures(model, calls, links, colours, base, args, log))
    return model, extra


def read_segment_sequences(
    gfa_path: str, extra_fasta: Optional[str], log: Log
) -> Dict[str, str]:
    """Segment sequences, from the GFA S lines or a supplementary FASTA."""
    out: Dict[str, str] = {}
    with smart_open(gfa_path) as fh:
        for line in fh:
            if not line.startswith("S\t"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3 and f[2] != "*":
                out[f[1]] = f[2]
    if extra_fasta:
        name, buf = None, []
        with smart_open(extra_fasta) as fh:
            for line in fh:
                if line.startswith(">"):
                    if name:
                        out[name] = "".join(buf)
                    name, buf = line[1:].strip().split()[0], []
                elif name:
                    buf.append(line.strip())
            if name:
                out[name] = "".join(buf)
    if not out:
        log.info("no segment sequences available (GFA stores '*'); composition screens are off")
    return out


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
            "  detangler.py --fai asm.fa.fai --gfa asm.gfa --paf segs.paf \\\n"
            "                --coverage cov.regions.bed.gz --out-dir out --prefix asm\n\n"
            "Then edit out/asm_karyotype.yaml and re-run with --config to lock in your calls."
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
    g.add_argument("--min-telomere-motifs", type=int, default=25,
                   help=argparse.SUPPRESS)  # retired: counted chance k-mers, see
                                            # --min-telomere-units
    g.add_argument("--min-telomere-fraction", type=float, default=0.02,
                   help=argparse.SUPPRESS)  # retired, as above
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
    g.add_argument("--speculative-penalty", type=float, default=1.5,
                   help="score penalty per join that is not supported by a traversable "
                        "path, i.e. where two segments merely end in the same one-sided "
                        "repeat. Such joins are reported, never silently used (default 1.5)")
    g.add_argument("--hypothesis", type=int, default=1,
                   help="which ranked hypothesis to draw as the ideogram (default 1)")

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
    model.warnings = list(log.warnings)

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
