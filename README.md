# detangler

Turn an assembler's output into a hypothesis of chromosome structure.

detangler reads a genome assembly graph (GFA), classifies every segment from
graph-intrinsic evidence (length, depth-derived copy number, links,
circularity, GC, telomere motif arrays), and proposes ranked hypotheses of
how the contigs resolve into linear chromosomes. It draws the result as a
paired figure: the assembly graph on the left, the inferred chromosomes on
the right, one colour per segment in both panels. No expected karyotype is
needed — the chromosome count is inferred from the structure of the graph.

## What it looks like

The paired figure, shown here on a real Flye assembly of *Fusarium
graminearum* — 11 graph segments resolving into 4 chromosomes plus a
mitochondrion (design mockup; the drawing code is being updated to match):

![Paired figure: assembly graph resolved into a chromosome
hypothesis](detangler_figure-mockup_v1.svg)

## How to run

The simplest run needs only the assembly graph (if the GFA carries
sequence, nothing else is required):

```
python detangler.py --gfa assembly_graph.gfa --out-dir results --prefix myasm
```

With more evidence (index, segment mapping, read coverage):

```
python detangler.py --fai asm.fa.fai --gfa asm.gfa --paf segs_to_asm.paf \
    --coverage cov.regions.bed.gz --out-dir results --prefix myasm
```

No data handy? Generate synthetic demo data and run on it:

```
python detangler.py --demo demo_data --out-dir results --prefix demo
```

To override any of the tool's calls, edit the emitted
`results/myasm_karyotype.yaml` and re-run from it:

```
python detangler.py --config results/myasm_karyotype.yaml --out-dir results --prefix myasm
```

## Required inputs

- `--gfa` the GFA v1 assembly graph (the file Bandage opens). If it carries
  sequence — Flye's does — this is the only input needed.

If the GFA lacks sequence, add one sequence source:

- `--fasta` assembly FASTA (also yields GC and N content), or
- `--fai` samtools faidx index, or
- `--assembly-report` NCBI assembly report

## Optional flags (the useful ones)

- `--flye-info assembly_info.txt` — adds Flye graph-path evidence
- `--coverage`, `--annotation`, `--paf`, `--agp` — extra evidence tracks
- `--telomere-motif MOTIF` — add a telomere motif beyond the built-in set
  (vertebrates, land plants, insects, nematodes, ciliates, fungi)
- `--expected-chromosomes N`, `--expected-genome-size N` — optional
  cross-checks only; never required
- `--hypothesis N` — draw a specific ranked hypothesis instead of the top one
- `--graph-style {bandage,layered}` — graph panel layout (default bandage)
- `--bandage-image FILE` — embed a real Bandage export as the left panel
- `--blast-db` / `--blast-subject` / `--blast-remote` — run BLAST on exported
  candidate segments and fold hits back in with `--blast-hits`
- `--out-dir`, `--prefix`, `--title` — output naming
- `--demo DIR` — generate demo data and run on it

Run `--help` for the full list, including all classification thresholds.

## Outputs

Written to `--out-dir` with `--prefix`:

- `<prefix>_report.md` — the evidence and ranked hypotheses, in plain language
- `<prefix>_paired.svg` / `.png` — assembly graph beside the chromosome
  hypothesis, shared colours
- `<prefix>_graph.svg` — the graph panel alone
- `<prefix>_ideogram.svg` / `.html` — the chromosome ideogram (HTML is
  interactive)
- `<prefix>_karyotype.yaml` — machine-readable result; reusable via `--config`
- `<prefix>_bandage_colours.csv` — colour map to load into Bandage
- `<prefix>_blast_commands.sh` and `<prefix>_repeat_candidates.fasta` —
  ready-to-run identification commands for the segments worth BLASTing
