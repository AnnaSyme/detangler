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

## 3. Ploidy: let the user declare it

**The problem.** Everything downstream of the depth baseline assumes one thing:
that a segment at baseline depth is present once. The baseline is the median
depth of segments >= 1 Mb, so that assumption is only safe for a haploid
assembly, or a diploid one that has been fully purged. It is wrong for:

- a **collapsed diploid**, where the long segments are two haplotypes merged and
  sit at 2x - so every copy number is halved and real repeats look single-copy;
- an unpurged diploid, where **haplotigs sit at ~0.5x** and are currently
  classified `low coverage / foreign`, i.e. as contamination. This is the single
  most attackable rule in the tool: ~0.5x is exactly the heterozygous-haplotig
  expectation, and calling the commonest phenomenon in a diploid assembly
  "foreign" will be the first thing a reviewer picks on;
- **polyploids**, where the whole notion of one baseline breaks down.

**The proposal.** An optional `--ploidy N` (or `--ploidy haploid|diploid|N`).
Default: unset, and the tool says so rather than assuming 1. When given, it
would:

1. Scale the expectation: a single-copy segment is expected at
   `baseline`, a haplotig at `baseline / N`, a two-copy repeat at `2 * baseline`.
2. Reclassify the ~0.5x class. With `--ploidy 2` a segment at half baseline,
   of decent length, and sitting in a bubble or beside a same-length segment is
   a HAPLOTIG - not contamination. Say so, and offer to exclude it from the
   molecule count rather than counting it as a chromosome.
3. Adjust the expected molecule count. A diploid's chromosome complement is 2n;
   whether the tool should report n molecules or 2n depends on whether the
   assembly is phased, which the user knows and the graph does not.
4. Feed the figure: draw haplotig pairs side by side rather than as separate
   chromosomes.

**Keep it a declaration, not an inference.** Ploidy could be estimated from the
depth histogram - a second mode at half the main one is diagnostic - but that is
real computational work and it is the assembler's and purge_dups' job. Asking
the user, who almost always knows, is one flag and no computation. The tool
should still SAY what it assumed, so the assumption travels with the result.

**Also worth a flag: expected chromosome sizes.** The karyogram panel draws each
molecule as an outline filled by the contigs that compose it. With only the
graph, the outline is the contigs' own extent, so it always reads as full.
Given expected sizes - from a related assembly, a cytogenetic estimate, or a
previous version - the outline becomes the target and the shortfall is visible
as unfilled space. That is what would let the figure answer "are any expected
bits missing", which it cannot honestly answer today.

**Status:** DONE 11 Aug 2026, as `--assembly-type primary|phased|collapsed`
(default `primary`) rather than `--ploidy N` — what governs how depth should be
read is what the ASSEMBLER emitted, not the organism's ploidy. A segment in the
`--haplotig-band` (default 0.35-0.65) and at least `--haplotig-min-length`
(20 kb) is now classed `haplotig`, not `low coverage / foreign`. Items 3 and 4
above (rescaling for collapsed assemblies, drawing haplotig pairs side by side)
are NOT yet done. The `--expected-chromosome-sizes` idea at the end is also
still open.

## 4. Centromere-shaped bridges as positive evidence

**Status:** DONE 11 Aug 2026. A long (>= `--centromere-min-length`, 10 kb),
markedly AT-rich, low-copy segment that bridges exactly two backbone ends now
earns `--centromere-bonus` (1.2), and a join through it pays only
`--centromere-speculative-discount` (0.4) of the usual speculative penalty. The
reasoning: an assembler is EXPECTED to fail to read through an AT-rich
centromere, so the missing through-path is explained rather than damning.

**Not validated - tuned.** On the *Fusarium graminearum* test GFA it lifts the
4-chromosome answer to rank 1 with no `--expected-chromosomes` supplied, but the
two constants were chosen to make that happen: set either to its neutral value
and the answer reverts to 5. One free parameter fitted to one known answer is
not validation, and calling it that here was wrong. The rule stays silent on the
demo data, which is the only thing about it that has been tested independently.
See `detangler_review-actions_v1.md`. It is a hypothesis-raiser, not a diagnosis, and
the report says so in those words.

**Still open:** a join's bridging segment is not DRAWN inside the chromosome
bar, so chr 1 renders as edge_2 + edge_7 with the 51.5 kb centromere invisible.
Also unbuilt: a `centromere_candidate` segment class, and testing the rule on a
lineage whose centromeres are not AT-rich.
