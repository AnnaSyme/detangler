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

**In progress.** An Arabidopsis Col-0 assembly is running on Galaxy
(`gtntesting`, history `a2190c6d75c35794`) from `SRR14728885`, the Col-XJTU HiFi
read set, via fasterq-dump then Hifiasm. The aim is a `p_ctg`/`p_utg` GFA *with*
sequence for a 5-chromosome genome whose telomere motif (`TTTAGGG`) is already
in the built-in set and whose CEN180 centromere arrays are the same situation as
`edge_8` but bigger and with a published answer.

**Still not tested: a genuinely tangled graph from anything but a fungus.** That
is the test that would stress the inference rather than the drawing. A hifiasm
`p_utg.gfa` from a vertebrate is the obvious candidate.

## Why the Fusarium example flatters the tool

Flye's `assembly_graph.gfa` is a **repeat graph**: it keeps the tangles that a
hifiasm contig graph has already resolved away. So the worked example is close
to a best case. Worth saying so rather than implying the tool generalises to any
file with a `.gfa` extension.
