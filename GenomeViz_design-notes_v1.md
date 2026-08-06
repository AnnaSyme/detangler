# Chromosome structure hypothesis from assembly graph output — design notes

Written 5 Aug 2026, from a worked manual example on *Fusarium graminearum* Nanopore data.
Everything below was derived by hand from real Flye output; the numbers are a usable test case.

_Moved into the `genome viz` project 6 Aug 2026 from the originating project. Content unchanged._

## What the tool is for

Take the files an assembler already produces, and return a ranked set of hypotheses about how the
contigs correspond to chromosomes, where the repeats are, and what is organellar or foreign.

The point is that assemblers already emit the evidence needed to do this. Most users never look at
it, because the graph is a 36 MB text file and the summary table has an undocumented column format.

## Inputs

`assembly_graph.gfa` — the primary source.

- `S` lines: `S<TAB>edge_N<TAB><sequence><TAB>dp:i:<depth>`. Sequence length and depth per segment.
- `L` lines: `L<TAB>edge_A<TAB><orient><TAB>edge_B<TAB><orient><TAB>0M`. Adjacency with orientation.
  A link from a segment to itself is a self-loop, which means a tandem array or a circular molecule.

`assembly_info.txt` — contig-level summary. Columns: `seq_name`, `length`, `cov.`, `circ.`,
`repeat`, `mult.`, `alt_group`, `graph_path`. `graph_path` is the ordered list of signed edge IDs a
contig traverses. `*` marks a dead end. Example: `*,6,-3,-9,*` means the contig runs edge_6, then
edge_3 reversed, then edge_9 reversed, with uncertain ends. `*,10,10,10,10,5,-3` means edge_10
traversed four times, i.e. a tandem array.

Optional but valuable: expected genome size, expected chromosome count, a reference assembly.

Note that contigs and edges are different objects with overlapping numbering. `contig_6` is not
necessarily `edge_6`. Do not assume the mapping; read it from `graph_path`.

## Inference steps

1. **Establish baseline single-copy depth.** Take the median depth of segments longer than some
   threshold (1 Mb worked here). In the example this was 19–20. Everything downstream is relative to
   this. Do not assume the assembler's own coverage estimate is on the same scale.

2. **Estimate copy number per segment** as depth / baseline. This is the single most informative
   derived value, and it is what separates the classes below.

3. **Classify segments.**

   | class | test |
   |---|---|
   | single-copy backbone | copy number ≈ 1, length above threshold |
   | repeat | copy number ≥ 2 |
   | tandem array | copy number high (>20), self-loop present |
   | organelle candidate | circular, copy number >> 1, length 15–200 kb, isolated component |
   | low-coverage / foreign | copy number < 0.5 |

4. **Add composition features.** Cheap and diagnostic:

   - GC content. Fungal centromeric and subtelomeric sequence is markedly AT-rich.
   - Telomere motif counts (TTAGGG / CCCTAA and organism-appropriate variants).
   - Periodicity for tandem arrays, to recover the repeat unit length.

5. **Use position within `graph_path`, not just presence.** A repeat appearing at path termini
   across several contigs is subtelomeric. A repeat appearing interior to a contig is something
   else. This distinction did most of the work in the example and is easy to compute.

6. **Build the graph, find connected components, enumerate traversals.** Then constrain: total
   backbone length should approximate expected genome size, and the number of chromosome-scale paths
   should approximate expected chromosome count.

7. **Return ranked hypotheses, not one answer.** See the honesty section below.

## Worked example — use as a test case

Flye 2.8.2, ONT R9 reads subsampled to ~22x, *Fusarium graminearum*. Reference for validation:
4 chromosomes, 36,563,796 bp nuclear; mitochondrion 95,638 bp.

Segments:

| segment | length | depth | copy no. | note |
|---|---:|---:|---:|---|
| edge_6 | 9,008,043 | 19 | 1 | backbone |
| edge_7 | 8,942,195 | 20 | 1 | backbone |
| edge_5 | 7,962,991 | 20 | 1 | backbone |
| edge_1 | 7,793,581 | 19 | 1 | backbone |
| edge_2 | 2,742,039 | 19 | 1 | backbone |
| edge_9 | 15,807 | 37 | ~2 | repeat, telomere motifs, 38% GC |
| edge_3 | 2,358 | 10 | <1 | 25% GC, strongly AT-rich |
| edge_10 | 7,685 | 1269 | ~67 | self-loop, tandem array, likely rDNA |
| edge_11 | 98,177 | 203 | ~11 | self-loop, circular, isolated — organelle |
| edge_8 | 51,498 | 4 | ~0.2 | low coverage, bridges two backbone segments |
| edge_4 | 13,947 | 5 | ~0.25 | low coverage, no links |

Links: 1–9 (both orientations), 2–8, 2–9, 3–5, 3–6, 3–9, 5–10, 7–8, 7–9, 10–10 self, 11–11 self.

**Correct conclusions the tool should reach:**

- Backbone totals 36.45 Mb against an expected 36.56 Mb, so the nuclear genome is essentially
  complete.
- `edge_11` is organellar: circular, isolated, 98 kb, ~11 copies.
- `edge_10` is a tandem array on the chromosome containing `edge_5`, ~67 copies of a 7.7 kb unit.
- `edge_9` is a repeat at contig ends, present in ~2 copies, carrying telomere motifs.
- Five backbone segments against four expected chromosomes means one chromosome is split.

**Conclusions the tool should not assert:**

- Which specific chromosome is which. Size ordering is a hint, not evidence.
- That `edge_6` + `edge_2` is one chromosome. It fits the size of the largest chromosome, but
  `edge_6` + `edge_5` through `edge_3`, and `edge_7` + `edge_2` through `edge_8`, are equally valid
  traversals.

## Honesty requirements

This is the part that will determine whether the tool is trusted.

- **Topology alone is ambiguous.** Where a repeat joins several backbone segments, the number of
  valid traversals grows quickly. The tool must enumerate them and rank, not pick one silently.
- **Separate observation from inference in the output.** Lengths, depths and links are observations.
  Copy number is a derived estimate. Chromosome assignment is a hypothesis. Label them differently
  and let the user see which is which.
- **Attach a reason to every call.** "edge_11 is organellar" is unhelpful. "edge_11: circular,
  isolated, 98 kb, 11x relative copy number, consistent with an organellar genome" is checkable.
- **Say what would resolve it.** For an ambiguous chromosome assignment, that is alignment to a
  reference or Hi-C data. Naming the missing evidence is more useful than a confident guess.

## Pitfalls found the hard way

- Bandage colours are random per render. Never key anything on colour, and never ask a user to match
  colours between two Bandage images.
- Default Bandage output overlaps its own labels when node lengths vary by three orders of
  magnitude. Consider emitting a layout rather than relying on the Bandage image.
- A single segment can appear at both ends of one contig (`edge_9` in `*,9,1,-9,*`). Path parsing
  must handle repeated IDs within a path.
- Low-depth segments can bridge two chromosome-scale contigs (`edge_8` here, at 4x). These look
  ignorable on depth but change the topology, so they must not be filtered out before traversal.
- A depth below baseline (`edge_3` at 10 against 19) is real and needs a class of its own. It is not
  simply a low-confidence single-copy segment.

## Suggested outputs

- A redrawn graph with segments coloured by class rather than at random, labels placed without
  overlap, and depth shown.
- A chromosome-level diagram with repeats marked in position, using the same colours.
- A table of segments with observation columns and inference columns visually separated.
- A ranked list of chromosome hypotheses, each with its supporting and contradicting evidence.

Colours used in the worked example, if useful as a starting palette: backbone segments in `#4a7ba7`,
`#5f9e6e`, `#a87f4a`, `#7d6ba7`, `#3f8f9e`; repeat `#c2703d`; AT-rich `#b0a04a`; tandem array
`#b5487f`; organelle `#3f7d5c`; low coverage `#a8a8a8`.

## Possible extensions

- Accept output from other assemblers. hifiasm and Verkko both emit GFA, though the depth tag and
  the auxiliary summary differ.
- Telomere motif detection at contig ends, to identify which chromosome ends are complete.
- Read-depth in windows along each contig, to spot collapsed repeats within a backbone segment,
  which the segment-level depth will not reveal.
- Optional reference alignment to promote hypotheses to assignments.
