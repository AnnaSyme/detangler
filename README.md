# detangler

Turn a genome assembly graph into a chromosome hypothesis.

detangler reads a genome assembly graph (GFA), classifies every segment, and
proposes ranked hypotheses of how the contigs resolve into linear chromosomes.
The graph on its own is enough: length, depth-derived copy number, links,
circularity, GC and telomere motif arrays all come out of the GFA. Outside
evidence is optional and folded in when given — read coverage, Flye's
`assembly_info.txt`, alignments, annotation, BLAST hits — as is
`--assembly-type`, which tells the tool what the assembler produced so it can
read depth correctly. It draws the result as a
paired figure: the assembly graph on the left, the inferred chromosomes on the
right, one colour per segment in both panels. No expected karyotype is
needed — the chromosome count is inferred from the structure of the graph.

## What it looks like

Real output, on a Flye assembly of *Fusarium graminearum* (ONT, ~22x). Eleven
graph segments on the left; on the right, the molecules they resolve into,
smallest to largest, each with the repeats the graph attaches to its ends. A
segment keeps the same colour in both panels.

![Paired figure: assembly graph resolved into a chromosome
hypothesis](detangler_figure_v3.svg)

A demo on a modest assembly, not a benchmark. What the tool infers depends on
how well the graph is resolved in the first place.

### Checking against a karyotype you already believe

If you know — or suspect — how many chromosomes there should be, pass
`--expected-chromosomes N`. The figure then shows the shortfall rather than
leaving you to count bars: every chromosome the graph did not produce gets a
dashed empty slot labelled *not found*, and a subheading names the number the
figure was given.

```
python detangler.py --gfa assembly_graph.gfa --expected-chromosomes 6 \
    --out-dir results --prefix myassembly
```

![The same assembly, told to expect six
chromosomes](detangler_figure_v3_expected6.svg)

Two things are worth understanding before you read this as a verdict.

**The flag is not a neutral overlay — it changes the answer.** Compare the two
figures above. Left to itself the tool asserts a join through the AT-rich
segment 8 and reports four chromosomes, with chr 1 built from contigs 2 + 8 + 7.
Told to expect six, it reports five separate molecules: leaving contigs apart
now scores better than joining them, so segment 8 falls back to being a cap on
the end of chr 5 rather than a centromere inside chr 1. Expecting more
chromosomes makes the tool more conservative about merging, which is the
intended behaviour but does mean you cannot use the flag to check a count
without also perturbing the inference that produced it.

**An empty slot means "the graph did not produce this", not "this is missing
from the genome".** The commonest reason a chromosome fails to appear is that
the assembler left it in pieces the graph gives no way to join — the tool says
so in the report — not that any sequence is absent. Reading a dashed slot as a
lost chromosome is the mistake this drawing most invites.

Nothing about the count is required. Without the flag the chromosome count is
inferred from the graph alone, which is the mode the tool is built around.

## How to run

The simplest run needs only the assembly graph (if the GFA carries
sequence, nothing else is required):

```
python detangler.py --gfa assembly_graph.gfa --out-dir results --prefix myassembly
```

With more evidence (index, segment mapping, read coverage):

```
python detangler.py --fai assembly.fa.fai --gfa assembly.gfa --paf segs_to_assembly.paf \
    --coverage cov.regions.bed.gz --out-dir results --prefix myassembly
```

No data handy? Generate synthetic demo data and run on it:

```
python detangler.py --demo demo_data --out-dir results --prefix demo
```

To override any of the tool's calls, edit the emitted
`results/myassembly_karyotype.yaml` and re-run from it:

```
python detangler.py --config results/myassembly_karyotype.yaml --out-dir results --prefix myassembly
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
- `<prefix>_paired.svg` — assembly graph beside the chromosome hypothesis,
  shared colours
- `<prefix>_graph.svg` — the graph panel alone
- `<prefix>_ideogram.svg` / `.html` — the chromosome ideogram (HTML is
  interactive)
- `<prefix>_karyotype.yaml` — machine-readable result; reusable via `--config`
- `<prefix>_bandage_colours.csv` — colour map to load into Bandage
- `<prefix>_blast_commands.sh` and `<prefix>_repeat_candidates.fasta` —
  ready-to-run identification commands for the segments worth BLASTing

## Credits

The assembly graph panel follows the approach used by
[Bandage](https://github.com/rrwick/Bandage) (Wick et al. 2015): each contig is
laid out as a chain of many small nodes rather than as a single rigid edge, so
that a contig's drawn length tracks its sequence length and long contigs can
curve around their neighbours. Bandage's edge construction is followed too —
each control point one contig-width beyond the end, along that contig's own
tangent, so a junction is continuous with the ribbons it joins. Both ideas are
reimplemented here in Python; no Bandage code is used, and detangler is not
derived from it.

Also drawn on for the graph layout:

- Initial placement comes from [Graphviz](https://graphviz.org) `sfdp` — Yifan
  Hu's multilevel force-directed algorithm
  ([paper](http://yifanhu.net/PUB/graph_draw_small.pdf)). Optional; there is a
  fallback.
- Angles at a junction are spread to a minimum gap by isotonic regression, using
  the pool-adjacent-violators algorithm of
  [Ayer et al. 1955](https://projecteuclid.org/euclid.aoms/1177728423).
- Contig ends are held at a fixed direction as well as a fixed position, so a
  contig settles as a minimum-bending-energy curve rather than a circular arc.
  Framing from [Levien & Séquin, *Interpolating splines: which is the fairest of
  them all?*](https://people.eecs.berkeley.edu/~sequin/PAPERS/2009_CAD_Levien_Sequin.pdf)
  (CAD & Applications 6, 2009).
- Each component is rotated to fill the panel, as Bandage does with
  `stepsForRotatingComponents`.

Segment colours are ColorBrewer's *Paired* qualitative scheme (Cynthia Brewer,
[colorbrewer2.org](https://colorbrewer2.org/#type=qualitative&scheme=Paired&n=10)),
reordered per figure so that segments drawn touching are as far apart in colour
as possible.

Co-built with [Claude](https://claude.ai).

## Licence

MIT. See [LICENSE](LICENSE).
