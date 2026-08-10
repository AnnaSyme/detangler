# detangler — enhancement list

Things worth doing, kept out of the code so the code stays honest about what it
does today. BLAST-specific proposals live in `detangler_blast-enhancements_v1.md`.

## 1. Benchmark against published finished assemblies

**The idea.** Collect genomes that are finished to chromosome level AND have a
publication that reports an assembly graph (or from which the graph can be
regenerated). Run detangler on the graph and compare its inferred chromosome
complement with what the paper concluded.

**What it would tell us**, in rough order of value:

- Does the tool reach the published chromosome count? On how many genomes?
- Where it disagrees, *why*: is the graph genuinely ambiguous and the paper used
  outside evidence (Hi-C, optical maps, a related reference), or did the tool
  miss something the graph does support?
- Does the reported range bracket the truth even when the point estimate misses?
  That matters more than the point estimate — the tool's claim is about honest
  uncertainty, not about being right.
- Does it degrade gracefully at low coverage, high heterozygosity, or with
  organelles and contamination present?

**Why it matters.** Right now the only validation is one *Fusarium graminearum*
assembly, where the tool says 5 (range 3–5) against a published 4. One data point
cannot distinguish "the tool is conservative" from "the tool is wrong". A dozen
genomes across fungi, plants and animals would.

**Design notes / traps.**

- Chromosome-count semantics must be pinned down first: cytogenetic haploid *n*,
  NCBI `totalNumberOfChromosomes` (which includes organelles), and detangler
  backbone molecules are three different numbers. Comparing them naively will
  manufacture disagreements. This is the same reconciliation problem flagged in
  the BLAST enhancements doc.
- Prefer assemblies where the *graph* is archived, not just the final FASTA. A
  graph regenerated with a different assembler or version is a different object
  and is not a fair test of the same inference.
- Record the assembler and version per genome — the shape of the graph, and
  therefore what detangler can infer, depends heavily on it.
- Expect the honest answer on many genomes to be "the graph does not resolve
  this", and treat that as a result rather than a failure.

**Status:** not started. Raised 10 Aug 2026.

## 2. Review the literature on untangling assembly graphs

**The idea.** Read across the work on getting from an assembly graph to
chromosomes, and see which ideas are worth borrowing. detangler's inference is
currently hand-built: classify segments from depth and composition, enumerate
legal joins, score them. That is a reasonable first pass, but it is not informed
by what anyone else has tried.

**Seed reference: arXiv:2206.00668** (Vrček, Bresson, Schmitz, Šikić). Trains a
GatedGCN over an assembly graph to predict, per edge, whether it leads to an
optimal reconstruction, then greedy-decodes those probabilities into paths and
converts them to contigs. Trained on graphs from simulated chromosome 19 HiFi
reads, evaluated on real CHM13 data, and compared against Raven's Layout
heuristics on the *same* input graph — which is the right way to isolate the
untangling method from the graph construction. Reported fewer contigs at
comparable or better genome fraction and NG50/NGA50.

**Read it for the framing, not the model.** Two things transfer directly even
if we never train anything:

- *Compare on the same input graph.* They deliberately hold the graph fixed and
  vary only the untangling, so differences cannot be blamed on the assembler.
  That is exactly the design the benchmarking programme above needs.
- *Label what "correct" means.* They generate ground-truth edge labels by
  simulating reads from a finished genome and keeping positional information.
  An equivalent for detangler would give real training or evaluation data:
  simulate from a finished assembly, build the graph, and you know which joins
  are true.

**The granularity differs and this matters.** Their nodes are READS and they
work in the Layout phase of OLC, before contigs exist. detangler starts from an
assembler's finished GFA, where nodes are already contigs. So their method is
not a drop-in replacement — it is upstream of us. Worth being precise about that
in any writing, because "untangling assembly graphs" describes both.

**Others worth finding.** Not yet surveyed: the graph-native scaffolders
(Rukki and similar), the organelle path-finders (oatk `pathfinder`), Tapestry,
and whatever exists for karyotype inference from graph topology. See
[[genomeviz-novelty-position]] in the assistant's notes for what the earlier
adversarial sweeps already established, so this does not repeat them.

**Status:** not started. Raised 10 Aug 2026.
