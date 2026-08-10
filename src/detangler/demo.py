"""Synthetic demo data."""
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
    _Rand,
)



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
