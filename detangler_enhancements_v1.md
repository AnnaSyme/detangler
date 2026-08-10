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
