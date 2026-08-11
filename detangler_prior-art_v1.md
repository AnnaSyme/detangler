# Prior art: does anything already do this?

Assessment 11 Aug 2026, prompted by finding Verkko-Fillet. Everything below was
read from the tools' own docs or source. Where something could not be verified
it says so.

## Verdict

Verkko-Fillet **dents but does not sink** the position. It does not produce a
chromosome complement from an assembly graph. It is an interactive curation
environment for editing Verkko's *path* files, and every chromosome-level
identity it produces comes from aligning contigs to a **reference genome** with
mashmap.

The claim that survives, and it is sharper than what the README says now:

> detangler infers a chromosome complement from a bare GFA, from any assembler,
> non-interactively, and reports it as a ranked list of alternatives.

The claim to **stop making**: that nobody else emits ranked alternatives.
Verkko-Fillet's `pp.searchNodes` returns a ranked table of candidate paths
scored by ONT read support. It is ranked *from reads*, for *one* junction the
user names, and a human still chooses - but "ranked, not a single answer" is
not by itself novel.

## Verkko-Fillet, in detail

**What it outputs.** `pp.writeFixedPaths` writes a corrected path/GAF file;
Verkko's own consensus module then builds the FASTA. Chromosome identity comes
from `tl.chrAssign`, whose docstring is "Run the script to align the assembly to
the given reference using mashmap and obtain the chromosome assignment results"
- `ref` is a required argument. The paper: "Each contig was renamed according to
its chromosomal assignment and haplotype, based on alignments to a matched
reference genome when available." With no reference, the fallback is to "sort
the chromosomes by length and orient them based on user-provided centromere
information". The giraffe chromosome count is cited from prior cytogenetics, not
derived.

**What it needs.** `pp.read_Verkko` hard-requires `assembly.paths.tsv`,
`assembly.scfmap`, `assembly.homopolymer-compressed.noseq.gfa` and
`assembly.colors.csv`; `checkFiles` demands the whole Verkko run tree. Note the
default graph is `noseq` - no sequence, so no GC and no telomere arrays from the
graph at all. The giraffe result also used parental Illumina, Hi-C, ONT
alignments, HiFi coverage, a prior giraffe reference AND a cattle reference for
orientation.

**Interactive, and the decisions are the user's.** `findGaps` only extracts the
`[...]` gap tokens Verkko already wrote. `fillGaps(obj, gapId, final_path, ...)`
takes the path as a literal string from the user and merely checks the flanks
match. `connectContigs` performs no evidence check at all. Its `cat` values
include `"random_assign"` and `"maynot_correct"`, which exist to flag joins
asserted without evidence.

**Verkko-specific, and the paper says so:** "it is currently limited to the
Verkko framework. Support for other assemblers could be a potential area of
future development."

**Two overlaps worth naming before a reviewer does.** `pp.calNodeDepth`
computes a normalised-depth proxy, and `pp.impute_depth` fits a regression
against user-supplied "real" values - user-calibrated, and the code never says
"copy number". And `tl.getT2T`, `pp.find_intra_telo`, `tl.detect_internal_telomere`
all detect telomere repeats, but for QC and trimming, not to cap a hypothesis.
It also writes a Bandage colour CSV, same idea as ours.

## Worth borrowing: a `--gaf` flag

Our own code already names the missing evidence. `hypotheses.py` prints "a read
or scaffold spanning ... would test the join directly". `searchNodes` is exactly
the machine that produces that number.

What it would need:

1. **A GAF parser.** `parsers.py` has `parse_paf` but nothing for GAF. Column 6
   is an oriented walk (`>utig4-12<utig4-88>utig4-7`); we need the ordered
   (segment, orientation) list plus query/path coordinates and mapq. Segment
   names must match the GFA's `S` names, which holds only if the GAF was made
   against the same graph - the only case to support.
2. **A spanning test on the Join we already build.** `find_joins` yields
   `Join(a, b, via, a_end, b_end)` and is end-aware. Count reads whose walk
   contains `a -> via -> b` contiguously, in either direction. Two anchors are
   essential: the alignment must extend at least k bp past the junction into
   BOTH `a` and `b` (a read that merely touches the repeat proves nothing), and
   must cover every `via` segment fully. Add `--gaf-min-mapq` and
   `--gaf-min-spanning`.
3. **Two scoring hooks.** A bounded bonus is the easy one. The one that matters
   is the tie-break: where `len(reach) > 2` we currently set `ambiguous` and
   take the same penalty for every competing pairing. With spanning counts we
   could rank the competitors instead - assert the supported one, penalise the
   ones with zero.
4. **A hard invariant.** Read support must never assert a join `find_joins` did
   not enumerate. The graph stays the legality gate; GAF only ranks within it.
   That keeps "a join is asserted only when the graph supports a traversable
   path" true.
5. **Two traps.** Verkko graphs are homopolymer-compressed, so GAF coordinates
   are in HPC space and a bp anchor threshold means something else - detect it
   or refuse. And GraphAligner multimaps, so without a best-alignment filter a
   repeat looks supported from every direction at once.

**Honest caveat.** On our own Fusarium case, `edge_8` is a 51 kb AT-rich block.
If it is longer than the reads, nothing spans it and `--gaf` changes nothing
there. It would help most on shorter ambiguous bridges. Say that rather than
overselling it.

## Everything else checked

| tool | overlap |
|---|---|
| **Rukki** | GFA plus a REQUIRED parental trio k-mer table; walks out one haplotype path set. Cannot run on topology alone, and gives one answer. Closest in spirit, furthest in evidence. |
| **gfastats** | graph metrics and scripted JOIN/SPLIT. Reports; does not infer which joins to make. |
| **Bandage / Bandage-NG** | visualiser plus BLAST-in-graph. No path inference. |
| **GraphAligner** | reads + GFA to GAF. Supplies the evidence above; not a competitor. |
| **GFAse** | phases a graph into haplotypes, but needs Hi-C/Pore-C BAM or trio k-mers. |
| **ntJoin / SALSA / YaHS** | scaffolders needing a reference or Hi-C. YaHS does use telomere motifs as an ancillary signal. |
| **gfatools, gaftools, MBG, Merqury** | parsing, conversion, construction, QV. No overlap. |

**The real conceptual neighbourhood is the plasmid/organelle family**, and it
should be cited:

- **plasmidSPAdes** separates chromosomal from plasmidic edges using coverage
  plus topology, no external database - but outputs a binary partition, not a
  ranked complement.
- **Recycler** extracts circular sequences from a FASTG using topology,
  coverage and length. Closest of all, but needs a paired-end BAM.
- **PlasBin-flow** uses GFA + flow-based depth + GC + a plasmid score in a MILP,
  structurally the nearest feature set to ours - but the plasmid score comes
  from BLASTing an external gene database, and the output is an optimised bin
  set, not ranked alternatives.
- **GetOrganelle** takes a GFA and can emit MULTIPLE path sequences, one per
  genome structure - genuine enumeration of alternatives - but needs a reference
  seed database and does not rank them.

Nothing found reports a chromosome count or karyotype from graph-intrinsic
evidence alone. Two independent searches; neither claims to be exhaustive.

## Verification status

Read directly: the verkko-fillet GitHub README, the readthedocs `chrAssign`
page, `docs/usage-principles.md`, the API docs, and the sources `_fill_gaps.py`,
`_searchNodes.py`, `_find_gaps.py`, `_read_wirte.py`, `_estLoop.py`,
`_graphAlign.py`, `_chrNaming.py`. The bioRxiv v3 full text gave Introduction,
Results and Discussion but **not the Abstract or STAR Methods**, so no claim
above comes from a methods section.

Could not confirm: PMC12621783 and the PubMed record (reCAPTCHA); the Cell
Genomics 2026 citation resolves to a Cell Press PII consistent with 2026 but the
page returned no readable text, so do not cite volume/issue/pages. The Zenodo
DOI 10.5281/zenodo.19867931 is the **software archive** of verkko-fillet
v0.1.25, not a paper or data deposit - it will not serve as a data citation.
`pp.get_NodeChr` is exported but its behaviour is unconfirmed.
