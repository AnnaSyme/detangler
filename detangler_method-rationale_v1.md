# detangler — method rationale: assembly graph to linear chromosome hypothesis

> **Dated 9 Aug 2026. Partly superseded - read the corrections in place below.**
> Three things changed after this was written: centromere-shaped bridges now
> score (they did not on 9 Aug), the `low coverage / foreign` class was renamed,
> and `--min-telomere-units` went from 3 to 10. The corrections are marked
> CORRECTED inline. For the current behaviour see the README and
> `detangler_review-actions_v1.md`.

How detangler converts an assembler's GFA into a ranked hypothesis of linear
chromosomes. Captured 9 Aug 2026 from the working session on the real Flye
*Fusarium graminearum* assembly (`real_data/flye_assembly_graph.gfa`).

## The logic, step by step

- **Copy number from depth.** Each segment's read depth (`dp:i:` tag) is
  normalised against the long segments, so ~1.0x means single-copy, ~2x means
  present twice, 67x means 67 copies. This one number drives most decisions.

- **Classify every segment first.**
  - Long segments near 1.0x → **single-copy backbone** — the only material
    chromosomes are built from.
  - ≥ ~1.7x → **repeat** (edge_9).
  - Self-link plus very high copy number → **tandem array** (edge_10: ~67
    copies of a 7.7 kb unit).
  - Circular, isolated, organelle-typical size, elevated copy number,
    divergent GC → **organelle candidate** (edge_11).
  - < ~0.6x → **below single-copy depth** — haplotype-specific,
    sub-stoichiometric, or contamination (edge_3, edge_4, edge_8). CORRECTED:
    this class was called `low coverage / foreign` when this was written. The
    word "foreign" was an accusation the depth does not support, and the same
    document then called edge_8 a candidate centromere.

- **Joins are asserted conservatively.** Two backbones are chained into one
  chromosome only when the segment bridging them is present in roughly one
  copy and touches nothing else (edge_5 + edge_6 joined through the edge_3
  bridge — the only join asserted on this data). Multi-way hubs like edge_9,
  which touches five segment ends, are never traversed: a repeat at a junction
  is compatible with many different walks.

- **Telomeres decide the molecule count.** Sequence ends are scanned for
  canonical telomere repeat motif arrays (TTAGG-family here). A finished
  linear chromosome carries an array at each end, so the number of capped ends
  bounds how many molecules the graph supports — this is where "4, but the
  graph supports 4–5" comes from. Chromosome ends abutting a telomeric segment
  count as capped. CORRECTED: this said "the telomeric segments edge_8/edge_9".
  **edge_8 has no telomere array.** The only telomeric segment on this graph is
  edge_9, a perfect (TTAGGG)x17 array 4 bp from its free `s` end. CORRECTED
  also: telomere credit is now capped by copy number - edge_9 is present in
  ~1.95 copies and so cannot cap four ends, which is what it was being paid for.

- **Centromeres are not called.** CORRECTED: no longer true as written. Since
  11 Aug, a long, markedly AT-rich, low-copy segment that bridges exactly two
  backbone ends earns `--centromere-bonus`, and a join through it pays a reduced
  speculative penalty — the reasoning being that an assembler is expected to
  fail to read through an AT-rich centromere, so the missing through-path is
  explained rather than damning. It is still a hypothesis-raiser and not a
  diagnosis: the tool does not assert that such a segment IS a centromere, and
  the two constants involved were tuned on one known answer, not validated.
  Proper telomere/centromere identification is still delegated to dedicated
  tools (tidk etc.).

  On edge_8 specifically, all that is established is that a 51 kb, 86%-AT block
  sits at the junction between edge_2 and edge_7. The length agreement with
  published chr1 is **not** evidence for including it: with edge_8 the chain is
  +0.33% against the published length, without it −0.11%. The fit is better
  without. What is defensible is that edge_2 and edge_7 each carry a telomeric
  neighbour at their outer ends, so joining them yields a telomere-to-telomere
  molecule, and all four molecules then match the four published chromosomes to
  about 1%.

- **Nothing is deleted; it is set apart.**
  - Isolated low-coverage sequence with divergent GC that touches nothing
    (edge_4) is excluded from chromosomes and shown as a candidate
    contaminant.
  - The organelle candidate is separated from the nuclear genome.
  - Repeats are not given one definitive location: edge_9 is drawn at each
    chromosome end the graph links it to, with its depth (~2 copies) as the
    check on how many placements can be real.

## Generality

The method must not be fungi-specific — it is intended to work for animals and
plants too. The telomere motif table already covers vertebrates, land plants,
insects, nematodes and ciliates (plus `--telomere-motif` for anything else).
The AT-rich flag is now worded lineage-conditionally in the code (in many
fungi it marks centromeric/subtelomeric sequence, but it can be a repeat
family, organelle-derived, or plain compositional bias). Organelle candidates
are now subtyped mitochondrion-like vs plastid-like (chloroplast) using
size, copy number, circularity, GC and a large inverted-repeat pair — the
plastome hallmark — with "type unresolved" when evidence is mixed; size never
decides alone because plant mitogenomes can exceed plastome size.

## Known issues (found 9 Aug 2026, subagent review against the raw GFA)

Status: items 1, 2 and 4 (report footer) were FIXED in code on 9 Aug 2026
(per-end adjacency `build_end_adjacency`, `facing_side` feature anchoring,
footer rebrand; 20/20 tests pass). Other "karyoglyph" occurrences (argparse
prog, log prefix, config key) remain and need a deliberate rebrand decision.
Item 3 was a mockup-only drawing error, corrected in the mockup.

1. **Telomere-cap counting.** The GFA links edge_9 to BOTH ends of edge_1
   (`edge_1+ – edge_9-` and `edge_1- – edge_9-`), but the report counts only
   one, giving 5 of 8 ends capped where the graph supports 6, and 2 molecules
   capped at both ends where the graph supports 3.
2. **edge_10 placement.** The report places the tandem array at the edge_3
   junction (`chain_3:7,962,990–7,962,991`); the GFA links it to edge_5's
   opposite, terminal end.
3. **Graph-figure topology.** edge_9 and edge_8 are dead-end spurs (all links
   on one end) but were drawn as pass-through bridges, implying a cycle that
   does not exist.
4. **Branding.** Report footer says "karyoglyph v1.0"; everything else says
   detangler.
