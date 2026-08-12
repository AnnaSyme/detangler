# Where to get assembly graphs to test on

Notes from a day of hunting, 11 Aug 2026. Recorded because the hunting was
harder than expected and the reason why is itself worth knowing.

## The headline problem

**Papers deposit the finished FASTA and throw the graph away.** Checked and came
up empty on all of these:

- **Zenodo 7962930**, Heng Li's own deposit of the plant assemblies from the
  hifiasm-UL paper - a paper *about* assembly graphs. Every file is `.fa.gz`.
  Includes `Arabidopsis.hifiasm.r587.*` and `Arabidopsis.verkko.1.3.1.*`, no GFA.
- **Zenodo 6947782**, the Galaxy Training assembly-QC dataset: yeast hifiasm and
  Flye assemblies. FASTAs and reads, no graph.
- Col-XJTU, Col-CEN and the Arabidopsis pan-genome papers: FASTAs.
- **Zenodo 15053552** (oatk) does hold 195 plant GFAs, but they are *organelle*
  genomes - single molecules, not a chromosome set.

No public Arabidopsis assembly graph was found. That is a mild argument for the
tool: the graph carries information the FASTA has already discarded, and almost
nobody keeps it.

## Where graphs do survive

**[GenomeArk](https://www.genomeark.org/) is the best source.** The VGP's public
S3 bucket, hundreds of vertebrate species, laid out as
`species/{Genus_species}/{ToLID}/assembly_{pipeline}_{ver}/`. Those
`assembly_hifiasm_*` and `assembly_verkko_*` directories hold the assembler's
raw output, which is where a GFA survives. Browsable at
`https://genomeark.s3.amazonaws.com/index.html?prefix=species/`

**[HPRC](https://github.com/human-pangenomics/HPP_Year1_Assemblies)** - raw
hifiasm output including GFAs is under each sample's
`assemblies/hifiasm_v0.14_raw/` in `s3://human-pangenomics/working/HPRC/`.

**[marbl/HG002](https://github.com/marbl/HG002)** publishes two Verkko graphs
directly, but both are `noseq` and homopolymer-compressed - see the caveats
below, they are close to useless here.

**Assemble your own.** The Galaxy Training QC dataset above includes
`SRR13577847`, real PacBio reads for yeast at 30x, ~237 MB. Yeast is only 12 Mb
but has **16 chromosomes**, all telomere-capped, which is a harder karyotype test
than a 4-chromosome fungus. For Arabidopsis, `SRR14728885` is the Col-XJTU HiFi
read set.

## Two caveats that decide whether a graph is usable

**`noseq` is disqualifying.** No sequence means no GC and no telomere arrays, and
telomere-capped ends are how the chromosome count is inferred. Without them
nothing caps anything, no join is asserted, and every backbone contig becomes its
own molecule. Verkko's published graphs are all `noseq`.

**Homopolymer-compressed is tolerable.** It shortens everything roughly
proportionally, so contig ratios, depth-derived copy number and topology all
survive; only the absolute Mb labels are wrong. Annoying, not fatal.

**And the one that surprised us: use the UNITIG graph, not the contig graph.**
See the README section. A `p_ctg.gfa` is post-cleanup and has almost no links
left. The bird graph tested here had 34 L-lines for 683 segments.

## What has actually been run

| data | segments | result |
|---|---|---|
| *Fusarium graminearum*, Flye `assembly_graph.gfa` | 11 | the worked example; 4 or 5 molecules |
| synthetic, human-like scale | 311 | ran in 1m43s, 109 molecules from 23 chromosomes, no sequence so no telomere evidence |
| *Baeolophus bicolor* `bBaeBic1` hap1 `p_ctg`, GenomeArk | 683 | 50 s, 242 molecules, zero joins - contig graph, nothing to detangle |
| synthetic 25-contig chain of direct links | 25 | found the direct-link scoring bug |
| synthetic bacterial chromosome + plasmid | 2 | found the "plasmid called a mitochondrion" bug |
| *Eucalyptus*, Flye `assembly_graph.gfa`, 645 Mb | 11,194 | found two bugs that made the tool hang or die silently; see below |

## The Eucalyptus graph, 11 Aug 2026

A real Flye repeat graph from Anna's own data: **11,194 segments, 7,557 links,
645.3 Mb, longest segment 2.20 Mb.** Eucalyptus has 11 chromosomes of roughly
40-60 Mb, so the longest contig in this graph is about 4% of a chromosome arm.
That is the whole story: there is no path from this graph's scale to chromosome
scale, and no inference can invent one.

**As a biological test it is uninformative.** Filtered to segments >= 500 kb
plus their bridging neighbours (716 segments, 225 Mb) it reports **218 molecules**
against a true count of 11 - one per backbone contig, exactly the artefact the
human-scale test predicted. Six telomere arrays were found in that slice (16 in
the >= 100 kb slice, 12-93 units each), and the baseline depth read correctly at
10x, so the sequence-derived evidence works; there is simply nothing to join.

**As a stress test it was very informative.** It found two bugs that no smaller
graph could:

1. **Silent death during hypothesis enumeration.** `enumerate_hypotheses` kept
   every valid hypothesis in a list, and each one holds a full linear forest
   over the backbone. At 2,991 segments that exhausted 3 GB and the process was
   killed with no message at all - the tool just vanished. Nothing downstream
   ever looks past `--max-hypotheses`, so it now keeps a bounded heap instead.
   Memory went from climbing past 1 GB and dying to flat at ~508 MB.
2. **An apparent hang in the colour assignment.** `assign_segment_colours` runs
   an exhaustive pairwise swap, O(4 x n^2) trials, each one copying the whole
   colour map and re-measuring every linked pair. At 716 segments that is
   roughly 1.4 billion operations and the run sat there long after the biology
   had finished. Found with `faulthandler.dump_traceback_later`, not by
   guessing. Graphs of 60 segments or fewer keep the exhaustive search
   unchanged, so the README figures are byte-identical; above that it repairs
   the single worst-contrast pair at a time.

Runtimes after both fixes: **716 segments 41 s, 2,991 segments 2m25s**, where
before neither finished. The full 11,194-segment file still needs more memory
than a 3 GB machine has, at the parse step.

There is now a scope warning at `SCOPE_SEGMENT_WARN = 400` segments, so a user
hears this before the run rather than discovering it from a 218-bar figure.

**One call to check.** The graph yields four organelle candidates, three of them
15-18 kb circles at 30-110x labelled "mitochondrion-like". Plant mitogenomes are
200-700 kb, so those are more likely rDNA or other high-copy circles. The
summary line also says "a mitochondrial genome", singular, while the table lists
four. Not yet fixed.

**In progress.** An Arabidopsis Col-0 assembly is running on Galaxy
(`gtntesting`, history `a2190c6d75c35794`) from `SRR14728885`, the Col-XJTU HiFi
read set. Submitted 11 Aug 2026:

- `fasterq-dump` gave one single-end file, 18.5 GB - about 130x for a 135 Mb
  genome, well past what hifiasm needs and enough to risk a memory or walltime
  failure hours in.
- `seqtk_sample` 1.5+galaxy0, fraction 0.3, seed 42 -> about 39x (hid 6).
- `Hifiasm` 0.25.0+galaxy3, standard mode, all option groups left at defaults
  (hid 7-12), chained so it queues behind the subsample.

The outputs to want are **hid 8, the raw unitig graph (`r_utg`)** and hid 9, the
processed unitig graph (`p_utg`). The wrapper writes the `noseq` versions to a
separate collection, so hid 8 should carry sequence - which is the thing that
disqualified every published Verkko graph.

Why this genome: 5 chromosomes, telomere motif `TTTAGGG` already in the built-in
set, and CEN180 centromere arrays that are the same situation as `edge_8` on the
Fusarium graph but larger and with a published answer to check against.

**Still not tested: a genuinely tangled graph from anything but a fungus.** That
is the test that would stress the inference rather than the drawing. A hifiasm
`p_utg.gfa` from a vertebrate is the obvious candidate.

## Why the Fusarium example flatters the tool

Flye's `assembly_graph.gfa` is a **repeat graph**: it keeps the tangles that a
hifiasm contig graph has already resolved away. So the worked example is close
to a best case. Worth saying so rather than implying the tool generalises to any
file with a `.gfa` extension.
