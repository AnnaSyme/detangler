# Review actions, 11 Aug 2026

An external review of the whole repo, plus a scale test at human genome size.
This is the to-do list that came out of it, in the order I would do it.

## Fix before showing anyone

**1. Stop presenting a tuned result as an inference.** `--centromere-bonus`
(1.2) and `--centromere-speculative-discount` (0.4) together contribute +2.1 to
hypothesis 1. Hypothesis 1 scores 5.50, hypothesis 2 scores 4.80. Set either
constant to neutral and the answer reverts to 5 molecules. The tie threshold is
0.75 and the gap is 0.70, so the report already says the graph cannot separate
them - and then the first line of the same report says "This assembly resolves
into 4 linear molecules".

PARTLY DONE 11 Aug: the summary now leads with the range whenever the top two
score within `--tie-threshold`, then names the drawn reading and by how little
it won. STILL TO DO: Print each
hypothesis's score with and without the centromere terms so the rule's
contribution is visible. Delete the word "validated" from
`detangler_enhancements_v1.md` - one free parameter fitted to one known answer
is not validation.

**2. DONE 11 Aug. Stop quoting 0.3% as evidence for the centromere join.** The
shipped code and README never asserted it; the claim lived in the working notes.
`detangler_method-rationale_v1.md` now states the arithmetic both ways and says
the fit is better WITHOUT edge_8, and no longer calls edge_8 a centromere as
though established. Verified against
RefSeq: edge_2+edge_8+edge_7 = 11,735,737 bp against chr1 = 11,697,295, so
+0.33%. But edge_2+edge_7 WITHOUT edge_8 = 11,684,239, which is -0.11%. The
length agreement fits BETTER without the centromere. It is not corroboration.

What is defensible: edge_2 and edge_7 each carry a telomeric neighbour at one
end and edge_8 at the other; joining them yields a telomere-to-telomere molecule
matching the largest published chromosome to ~0.1-0.3%, and all four molecules
then match all four published chromosomes to ~1%. That is a good result. Lead
with it. Do not claim the graph establishes the join - it does not, both
attachments are to edge_8's same end - and do not call edge_8 a centromere. All
that is known is that a 51 kb 86%-AT block sits at the junction.

**3. DONE 11 Aug. The report lied when `--expected-chromosomes` was supplied.**
`report.py:190` prints "No expected karyotype was used to reach that"
unconditionally, while `hypotheses.py:512` has already applied
`score -= 3.0 * abs(n - expected)`. The sentence is now conditional and names the penalty when it applied.

## Correctness

**4. DONE 11 Aug. A direct GFA link scored net NEGATIVE.** `hypotheses.py:544` charges
`--join-cost` 0.6 and line 546 refunds only 0.3, so an unambiguous direct link
nets -0.3 while a speculative one-sided AT-rich bridge earns +1.3. A 25-contig
chain joined by 24 plain links reports 25 separate molecules. New `--direct-link-bonus` (0.55) on top of cancelling the join cost. The
25-contig chain now reports 7, the residue being the `--max-join-edges` cap.

**5. DONE 11 Aug. The one-sided contradiction text was wrong for half its
cases.**
`hypotheses.py:364` hardcodes "which has links on one side only". edge_3 has
links on both ends; the real reason its join fails is that both partners attach
to the SAME end. The message now distinguishes the two failures.

**6. DONE 11 Aug. Telomere credit was unbounded by copy number.** `hypotheses.py:488` credits
edge_9 - one segment at copy number 1.95 - with capping four molecule ends, 4.8
of hypothesis 1's 5.50 total. Five backbone ends attach to its single 'e' end;
if all were real you would expect ~95x depth, and the observed depth is 37x.
`chain_end_status` now returns which segment earned each cap; the score counts
at most `round(copy number)` of them and the rest become open ends, with the
shortfall stated. Hypothesis 1 fell 5.50 to 3.10, hypothesis 2 4.80 to 2.40.

**7. Circular replicons are outside the model and the tool does not say so.** On
a bacterial GFA - 5 Mb circular chromosome plus a 60 kb plasmid - the summary
prints "1 linear molecule and a mitochondrial genome". The plasmid satisfies
every organelle rule and the chromosome is drawn linear because `_linear_forest`
rejects cycles. Detect the prokaryote case and refuse rather than mislabel.

**8. Overlaps.** DONE 11 Aug: L-line CIGARs are now subtracted from chain
lengths, and non-zero overlaps are reported. Flye writes 0M so this was
invisible; on miniasm, canu or any unitig graph every molecule length was
inflated by the sum of its junctions.

**9. hifiasm depth.** DONE 11 Aug: `rd` added to the tag list, with hifiasm's
`rd:i:n = coverage n+1` correction, behind `--assembler`. Before this, hifiasm
GFAs gave `depth = None`, and `calls.py:295` then let everything over 500 kb
pass as backbone - no repeat, organelle, haplotig or centromere ever called,
silently.

**10. Drawing size guard.** DONE 11 Aug: `--max-graph-nodes` now PRUNES to the
longest N segments and states what was dropped, instead of refusing to draw. The
paired figure had no guard at all before this.

## Presentation

**11. The figure asserts what the report withdraws.** Bars are labelled
`chr 1..N` by size rank with no caveat. On this dataset that mislabels half the
karyotype - the tool's chr 3 is nearest published chr4 and its chr 4 is nearest
published chr3 - and the numbering changes when `--expected-chromosomes` is
supplied. A speculative join is drawn identically to an asserted one. Label the
bars `chain N`, or `chr N?`, and draw a speculative join as a break.

**12. `low coverage / foreign` as a class name.** The same document calls
edge_8 "foreign" in the segment table and a candidate centromere in the
hypothesis section. "sub-baseline depth" carries the same information without
the accusation.

**13. DONE 11 Aug. Stale AND duplicated docs.** The two byte-identical
`chromviz_*` leftovers are gone from the working tree and from `git ls-files`.
`detangler_method-rationale_v1.md` now carries a dated superseded banner and
three inline CORRECTED notes: centromere-shaped bridges DO score now, edge_8 has
no telomere array (edge_9 is the only telomeric segment), and the
`low coverage / foreign` class has been renamed.

**14. DONE 11 Aug. Raised to 10.** `--min-telomere-units 3`: The chance argument in `cli.py:446` is for one
specified 18-mer at one position, but the scan is 6 motifs x 2 strands x 2
windows x every segment, and the motifs are AT-skewed. Three units is also not
biologically a telomere. Raise to ~10; keep reporting the unit count.

## Scope - what the human-scale test showed

A 311-segment, 3.03 Gb graph with human-like structure ran in 1m43s and
produced a figure, but the figure is not usable and the inference is not either:

- **109 molecules from 23 chromosomes.** With no sequence in the GFA there is no
  telomere evidence, so nothing caps anything, no joins are asserted, and every
  backbone contig becomes its own molecule. The count is an artefact of the
  input format, not a finding.
- **The palette cycles.** Ten Paired colours across 300 segments means colour no
  longer ties the two panels together, which is the figure's whole premise.
- **262,144 hypotheses enumerated** - the 2^18 cap - to no purpose.
- **The chromosome panel is a 109-bar staircase**, unreadable at any zoom.

The honest conclusion is that this tool is built for SMALL graphs: fungal,
bacterial, organellar, a resolved plant chromosome set. That is a reasonable
scope and it should be stated in the README rather than discovered. A vertebrate
assembly graph needs a different figure - probably per-component, drawn one
chromosome at a time - and a different enumeration strategy.

Ceilings measured: 40 segments 3.0 s, 100 -> 18.3 s, 200 -> 74.2 s, and the
layout is O(n^2) and runs twice. Pruning keeps the drawing tractable; it does
not make the inference meaningful at that scale.

## Testing

**15. DONE 11 Aug. The regression suite does not cover the headline claim.**
`real_graph_checks` in the test script now runs the real GFA and pins three
things: the tool reaches 4 molecules by default, setting `--centromere-bonus 0`
changes the answer to 5, and the top two hypotheses stay inside the tie
threshold so the headline must not assert a single number. It skips cleanly if
`real_data/flye_assembly_graph.gfa` is absent. 23/23.

Original finding: Its fixture sets
edge_8 to GC 47% against a 48% baseline, so the centromere rule cannot fire, and
the suite reports 5 molecules while the shipped tool reports 4. It also puts
edge_9's telomere array at both ends where the real one is at the start only.
Add a test that runs the real GFA and asserts the score gap between hypotheses 1
and 2, with and without `--centromere-bonus`.

## Scale test 2: a real vertebrate graph, 11 Aug 2026

`bBaeBic1_hap1_contig_graph.gfa` from GenomeArk (VGP bird, hifiasm hap1). 683
segments, 1.09 Gb, sequence present, `rd:i:` depth tags. Ran in 50 s with
`--assembler hifiasm --max-graph-nodes 60`.

Result: **242 molecules, 441 unplaced, one hypothesis, zero joins asserted,
score 0.00.** Baseline depth 23x read correctly from the hifiasm tags; 65
telomere arrays found, so the sequence-derived evidence worked.

The reason is the input, not the tool: the file has **34 L-lines for 683
segments**. A `p_ctg` graph is post-cleanup - the assembler has already resolved
and popped every tangle it could. There was nothing to detangle.

Two things follow:

1. **Document which GFA to use.** Done, as a new README section: unitig graph
   (`r_utg` / `p_utg`), not contig graph. This is the single most consequential
   piece of user guidance in the tool and it was missing.
2. **The Fusarium graph is unusually good for this tool** - Flye's
   `assembly_graph.gfa` is a repeat graph, so it keeps the tangles a hifiasm
   contig graph has already discarded. Worth saying so rather than implying the
   tool generalises to any file named `.gfa`.

Still untested: a genuinely tangled vertebrate unitig graph. That is the test
that would actually stress the inference.
