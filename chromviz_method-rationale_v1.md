# chromviz — method rationale: assembly graph to linear chromosome hypothesis

How chromviz converts an assembler's GFA into a ranked hypothesis of linear
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
  - < ~0.6x → **low coverage / foreign** — haplotype-specific,
    sub-stoichiometric, or contamination (edge_3, edge_4, edge_8).

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
  graph supports 4–5" comes from. Chromosome ends abutting the telomeric
  segments edge_8/edge_9 count as capped.

- **Centromeres are not called.** The tool only flags markedly AT-rich
  segments (GC well below the genome baseline, 48% here), because in many
  fungi AT-rich blocks mark centromeric or subtelomeric sequence. That is a
  hint, not a call — proper telomere/centromere identification is parked for
  delegation to dedicated tools (tidk etc.).

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
   chromviz.
