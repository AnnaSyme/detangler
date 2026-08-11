"""Readers for every input format the tool accepts."""
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
    _maybe_float,
    smart_open,
)
from .records import (
    ContigInfo,
    CoverageWindow,
    Evidence,
    GfaLink,
    GfaSegment,
    Placement,
    SeqRecord,
)



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


def parse_gfa(
    path: str, log: Log, assembler: str = "unknown"
) -> Tuple[Dict[str, GfaSegment], List[GfaLink]]:
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
                depth = _segment_depth(name, length, tags, assembler)
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


def _segment_depth(
    name: str, length: int, tags: Dict[str, str], assembler: str = "unknown"
) -> Optional[float]:
    """
    Read depth off an S line. Assemblers do not agree on how to write it.

    `dp`/`DP` is a depth directly. `KC`/`RC`/`FC` are COUNTS - k-mers, reads,
    fragments - so they are divided by length; note that RC/length is reads per
    base, not coverage, and is out by roughly the mean read length. Ratios
    survive, which is what copy number needs, but the printed "depth" is then
    not a depth. hifiasm writes `rd`, and its documentation is explicit that
    `rd:i:n` means coverage n+1, so a segment supported by one read carries
    rd:i:0 - without the correction those segments look like zero coverage.
    https://hifiasm.readthedocs.io/en/latest/interpreting-output.html
    """
    if "rd" in tags:
        try:
            v = float(tags["rd"])
            return v + 1.0 if assembler == "hifiasm" else v
        except ValueError:
            pass
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
