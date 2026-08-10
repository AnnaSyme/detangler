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

from .cli import main
from .common import VERSION

__all__ = ["main", "VERSION"]
