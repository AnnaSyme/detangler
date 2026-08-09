# detangler: design notes for sequence-matching enhancements

Date: 2026-08-09
Author: research notes for Anna
Status: design notes only, nothing implemented

## Context

detangler reads an assembler's GFA, classifies segments from graph-intrinsic
evidence (coverage, topology, composition, telomere motifs), and proposes a set
of linear chromosomes with a paired graph/ideogram figure and a text report. Its
philosophy is graph-first: extract cheap intrinsic signal, and delegate deep
sequence analysis to established external tools rather than reimplementing them.

It already does one piece of external delegation. When run, it exports repeat and
ambiguous segments to a candidate FASTA and writes a ready-to-run helper script
(`real_data/detangler_fusarium_blast_commands.sh`) offering three `blastn`
invocations (local `nt`, a `-subject` FASTA, or `-remote`), all using the fixed
output format `6 qseqid sseqid pident length qcovhsp evalue bitscore stitle`. The
results are fed back with `--blast-hits`, filtered by `--blast-min-identity` and
`--blast-min-coverage`, and the top hit becomes an evidence line and a label on
the figure. The classes exported for search are controlled by `--blast-classes`
(default: `repeat`, `tandem_array`, `organelle_candidate`, `at_rich`,
`low_coverage`).

These notes describe how to extend that pattern. They do not duplicate the
existing BLAST helper; they build on it. Each enhancement below states the
question it answers about chromosome structure, what slice of the assembly you
would extract, which verified tool or database it would use, how the result flows
back into the figure/report, and an honest usefulness/effort judgement.

A recurring tension worth stating up front: detangler is currently offline and
self-contained. Every enhancement here adds either an external database download
(often large) or a network call. None should become a hard dependency; all
should be optional, gated behind a flag, and degrade gracefully to "not run" the
way the current BLAST step already does.

The validation graph (Fusarium graminearum, Flye, 11 segments: 5 backbones,
telomeric repeat edge_9, telomeric low-coverage edge_8, tandem array edge_10,
AT-rich contaminant candidate edge_4, circular mito edge_11) is used as the
running example.

---

## Enhancement 1: identify unplaced or ambiguous segments by homology (extend the existing BLAST step)

- **Question it answers:** what is this segment, really? For segments detangler
  cannot confidently assign (short single-copy, unclassified, contaminant
  candidates like edge_4), a homology hit says whether it is host sequence, a
  known repeat family, an organelle fragment, or foreign DNA — which decides
  whether it belongs on a chromosome at all.
- **What you extract:** the segment nucleotide sequence (already done for the
  classes in `--blast-classes`).
- **Verified tool/DB:**
  - `blastn` against NCBI `core_nt`. NCBI made `core_nt` the default web
    nucleotide BLAST database in August 2024 and released it for standalone
    BLAST download; it is `nt` with most eukaryotic chromosome sequences
    removed (less than half the size) while keeping transcript and gene-related
    sequences. Note the tradeoff: because `core_nt` deliberately drops most
    eukaryotic chromosome arms, it can be *worse* than `nt` for the specific
    task of "which chromosome does this piece come from." For detangler's use
    (naming repeats, spotting contaminants, ID'ing organelles) `core_nt` is a
    reasonable, smaller default; keep `nt` selectable.
  - DIAMOND (`diamond blastx` / `blastp`) for protein-level search when
    nucleotide identity is too low to register — much faster than BLASTX against
    large protein sets. Verified as a real tool (bbuchfink/diamond), used this
    way inside the BlobToolKit pipeline against UniProt reference proteomes.
- **How it flows back:** identical mechanism to the current `--blast-hits` path
  — top hit becomes an evidence line and an optional label. A protein-level hit
  would need a second column set (DIAMOND's default outfmt 6 differs from the
  blastn one already hard-coded), so this is a real, if small, code change.
- **Usefulness/effort:** High usefulness, low-to-moderate effort. This is the
  smallest, safest step because it generalises machinery that already exists.
  The main new work is (a) offering `core_nt` vs `nt` as a choice and (b)
  optionally adding a DIAMOND protein path.

## Enhancement 2: name rDNA arrays and telomere-adjacent repeats explicitly

- **Question it answers:** is this tandem array (edge_10) a ribosomal DNA (rDNA)
  cluster, and is the telomeric segment (edge_9) genuinely subtelomeric? rDNA
  arrays and telomere-adjacent repeats are the two repeat types most likely to
  break or bloat a graph, and knowing which is which explains why a segment sits
  off the backbone.
- **What you extract:** tandem-array and telomere-flagged segments (detangler
  already classifies `tandem_array` and detects telomere motifs).
- **Verified tool/DB:**
  - **rRNA genes:** barrnap (tseemann/barrnap) annotates rRNA in bacteria,
    archaea and fungi (`--kingdom fun`), takes a FASTA and emits GFF3 with
    `16S`/`23S`/`5S` (and eukaryotic equivalents via Rfam models); it uses
    Infernal/Rfam internally. This is the direct way to confirm a tandem array
    is rDNA. (It reports rRNA genes, not full array copy number.)
  - **Telomeric repeats:** tidk (tolkit/telomeric-identifier, "Telomere
    Identification toolKit"). `tidk explore`/`tidk find`/`tidk search` identify
    telomeric repeat motifs de novo or against a clade's known motif and count
    them in windows across the assembly. This corroborates detangler's own motif
    detection with an established, citable tool (Bioinformatics, 2025).
- **How it flows back:** "rDNA (barrnap: 18S/5.8S/28S)" or "telomeric repeat
  confirmed (tidk)" as evidence lines; a distinct glyph or annotation for rDNA
  arrays and telomere ends on the ideogram. Both tools emit coordinates, so
  positions could in principle be drawn along a segment.
- **Usefulness/effort:** High usefulness for exactly the segment types in the
  Fusarium example; moderate effort. Both are conda-installable, fast, and
  file-in/file-out — a good fit for the delegate-to-external-tools philosophy.

## Enhancement 3: organelle identification (the mito segment, edge_11)

- **Question it answers:** is the circular low/variant-coverage segment the
  mitochondrion (or, in plants, the plastid), and should it be labelled as an
  organelle rather than a nuclear chromosome? detangler already flags edge_11 as
  a circular mito candidate; this would confirm it by gene content.
- **What you extract:** segments classed `organelle_candidate` (or the circular
  subgraph).
- **Verified tool/DB (three real options, different mechanisms):**
  - **oatk** (c-zhou/oatk) with **OatkDB** HMM profile databases. `oatk`/
    `hmmannot` run `nhmmscan` (HMMER) against organelle gene profiles and
    `pathfinder` parses the organelle subgraph. It can take an existing
    assembly GFA via `-G` (not just reads), which matches detangler's input.
    Caveat from its own docs: it relies heavily on sequence/arc coverage tags
    and is tuned for HiFi/`syncasm`-style graphs; on a Flye graph you would need
    to set the coverage tags correctly and results may be imperfect.
  - **MitoHiFi** (marcelauliano/MitoHiFi) can start from assembled contigs
    (`-c`), finds mito contigs by BLASTing against a close reference it fetches
    with `findMitoReference.py`, separates NUMTs, circularises and annotates.
    Designed for PacBio HiFi; supports `-a {animal,plant,fungi}`.
  - **GetOrganelle** (Kinggerm/GetOrganelle) has a `get_organelle_from_assembly.py`
    mode that takes a FASTG/GFA and extracts the organelle using seed and label
    databases; it ships `fungus_mt` and `fungus_nr` modes among others.
- **How it flows back:** an "organelle confirmed (tool, genes found)" evidence
  line and an organelle label/colour on the figure; confirmed organelles could
  be drawn separately from the proposed nuclear chromosome set.
- **Usefulness/effort:** Useful as confirmation, but moderate-to-high effort and
  arguably heavier than the question needs. detangler already infers "circular +
  off-backbone + coverage anomaly = organelle" from the graph. Full organelle
  assemblers are built to *assemble* organelles from reads, not to label a
  single pre-existing segment. The lightest honest option is often just
  Enhancement 1 (BLAST the segment against `core_nt`) or barrnap/HMM gene
  detection, reserving oatk/MitoHiFi/GetOrganelle for cases where the organelle
  is fragmented across several segments and needs reconstructing.

## Enhancement 4: contamination screening (the AT-rich candidate, edge_4)

- **Question it answers:** is this segment foreign (a different organism, an
  adaptor/vector), and should it be excluded from the chromosome set? Directly
  relevant to edge_4, detangler's AT-rich contaminant candidate.
- **What you extract:** whole assembly, or the `at_rich`/contaminant-candidate
  segments.
- **Verified tool/DB:**
  - **NCBI FCS** (ncbi/fcs): FCS-GX detects cross-species contamination by
    aligning sequences to a large NCBI genome database and assigning a taxonomic
    division, flagging sequences whose assignment differs from the declared
    taxon; FCS-adaptor screens for adaptor/vector contamination. Runs via
    Docker/Singularity. Real caveat: FCS-GX's database is very large (the
    documented resource requirement is on the order of ~470+ GB, ideally in
    RAM), so this is a cluster/server step, never something to bundle offline.
  - **BlobToolKit** (genomehubs/blobtoolkit, BlobTools2): builds a BlobDir from
    coverage + BLAST/DIAMOND taxonomy + BUSCO and produces blob/snail plots and
    filtered datasets to separate cobionts/contaminants. Heavier, more
    interactive; better as a full-assembly QC companion than a per-segment call.
- **How it flows back:** a "flagged as possible contaminant (FCS-GX: <taxon>)"
  evidence line; contaminant segments visually marked and excluded (or held
  aside) from the proposed chromosome list. This strengthens exactly the kind of
  call detangler already makes tentatively about edge_4.
- **Usefulness/effort:** High usefulness (contamination is a real reason a
  segment should not be placed), but high infrastructure effort. Best treated as
  "detangler emits the FASTA and a suggested FCS command; you run it on a server
  and feed a verdict back," mirroring the existing offline BLAST helper. Do not
  make it a dependency.

## Enhancement 5: organism identification from an assembly slice

- **Question it answers:** what organism is this assembly (or a suspicious
  segment) actually from? This underpins Anna's follow-up idea — you must know
  the species before you can look up its expected chromosome number — and it is
  also how you tell "off-backbone segment is a cobiont" from "off-backbone
  segment is host repeat."
- **What you extract:** a representative slice (e.g. the longest backbone
  segment) for whole-assembly ID, or an individual ambiguous segment.
- **Verified tool/DB:**
  - **Mash** (marbl/mash): MinHash sketching; `mash screen`/`mash dist` against
    a sketched RefSeq database gives fast approximate identity. Mash's own paper
    sketched all of RefSeq release 70. Fast and small, but sketch-only distances
    are less reliable for incomplete/medium-quality sequence.
  - **skani** (bluenote-1577/skani): ANI + aligned fraction via sparse chaining;
    more accurate than Mash for incomplete/MAG sequence, can classify an
    assembly against a preprocessed database of tens of thousands of genomes in
    seconds. Good middle ground.
  - **sourmash**: provides prepared databases including "NCBI Eukaryotes
    (Jan 2025)" and GTDB bacterial/archaeal sets, usable for k-mer based
    taxonomic ID/containment.
  - **BLAST against `core_nt`/`nt`** (Enhancement 1) remains the highest-
    resolution but slowest option, and is the standard "what is this" answer.
- **How it flows back:** a report header line ("assembly most closely matches
  <species> (skani ANI x%, aligned fraction y%)"); per-segment, a taxonomic tag
  that separates host from cobiont. It also produces the species name that
  Enhancement 6 needs.
- **Usefulness/effort:** Useful, moderate effort. The sketch tools (Mash, skani,
  sourmash) are light and fast; the cost is the reference database download.
  Sensible as an optional "identify" step that prints a best guess and its
  confidence, never overriding the user.

## Enhancement 6: expected chromosome number as a sanity check (Anna's follow-up idea)

- **Question it answers:** does the number of chromosomes detangler proposed match
  what is known for this species? Concretely, it could help resolve the
  4-vs-5 ambiguity in the Fusarium example. (For reference: the published
  completed PH-1 reference of *Fusarium graminearum* is organised as four
  chromosomes — King et al. 2015, BMC Genomics — so an expected-count lookup for
  this species should return 4, which is a useful check against a 5-backbone
  inference that may be splitting one chromosome or counting an organelle/repeat
  as a chromosome.)
- **What you extract:** nothing new sequence-wise; this consumes the species name
  from Enhancement 5 (or a user-supplied name) and queries a metadata database.
- **Verified DB/tool:**
  - **GoaT (Genomes on a Tree, genomehubs)** at goat.genomehubs.org, with a
    documented **goat-cli** (bioconda `goat`) and a public API
    (goat.genomehubs.org/api-docs). GoaT holds direct or estimated values for
    dozens of taxon attributes across eukaryotes, and the GoaT publication and
    GenomeHubs docs explicitly list **cytologically determined chromosome
    number** as a taxon attribute. This is the most fit-for-purpose source: it
    is queryable programmatically and organised on the taxonomy tree, so it can
    fall back to a genus/family estimate when the exact species is absent.
    (Uncertain: the exact machine attribute name(s). GoaT documentation
    describes a "chromosome number" attribute; I could not verify from the
    sources fetched whether the queryable field is literally `chromosome_number`
    and whether a separate `haploid_number` field exists. The API/`goat-cli`
    should be checked directly before coding against a field name.)
  - **NCBI Datasets** `datasets summary genome` / `dataformat`: the genome
    assembly data report includes a verified field **`totalNumberOfChromosomes`**
    (CLI flag `total-number-of-chromosomes`), defined as the count of nuclear
    chromosomes, organelles and plasmids in a *submitted assembly*. Important
    distinction: this is the chromosome count of existing *assemblies* in NCBI,
    not a cytogenetic expectation — useful as "how many chromosomes do published
    assemblies of this species have," and note it includes organelles/plasmids
    in the count.
  - **CCDB (Chromosome Counts Database)**, ccdb.tau.ac.il: real, large, and the
    biggest source of plant chromosome counts (covers ~18% of vascular plants),
    giving mode count and n vs 2n. **Plants only** — not applicable to the
    Fusarium example or to fungi/animals.
  - **Animal Genome Size Database**, genomesize.com: real, but it stores
    **C-values (genome size in pg)**, not chromosome number. It does *not*
    answer the expected-chromosome-count question and should not be cited as if
    it did.
- **How it flows back:** a report line such as "proposed 5 chromosomes; GoaT
  lists chromosome number = 4 for this species — review whether two backbones
  are one chromosome or whether an organelle/array is being counted." Never an
  automatic override; a flagged discrepancy for the user to adjudicate.
- **Usefulness/effort:** Genuinely useful as a *sanity check*, low-to-moderate
  effort (a single API call given a species name). But the risks are real and
  must be stated in any implementation:
  - **Karyotype variation within a species:** chromosome number is not always
    fixed. Aneuploidy, B chromosomes, supernumerary/dispensable chromosomes
    (well known in *Fusarium* and related fungi), and polyploidy mean the
    "expected" number can legitimately differ from a correct assembly. A
    mismatch is a prompt to look, not proof of error.
  - **Database coverage gaps:** most species have no cytogenetic record; GoaT
    will often fall back to a higher-rank estimate or return nothing.
  - **Misidentification upstream:** the check is only as good as the species ID
    from Enhancement 5; a wrong organism call yields a wrong expectation.
  - **Counting conventions differ:** cytogenetic haploid number, NCBI's
    `totalNumberOfChromosomes` (which includes organelles/plasmids), and
    detangler's "backbones" are three different things. Comparing them naively
    will manufacture false discrepancies. Any implementation must reconcile
    what is being counted (e.g. exclude the mito segment before comparing to a
    nuclear haploid number).
  - **Offline philosophy:** this step is inherently online. Keep it optional and
    behind a flag, and let detangler run fully without it.

---

## How these fit detangler's existing design

- All six reuse the established pattern: detangler classifies from the graph,
  optionally exports a FASTA/segment set plus a suggested command, and accepts a
  results file back that becomes evidence lines and figure annotations. The
  current `--blast-hits` plumbing is the template.
- Ordered by best effort-to-value for the Fusarium case: (1) extend BLAST to
  `core_nt`/DIAMOND, (2) barrnap + tidk for rDNA/telomeres, (4) FCS as an
  emit-command-only contaminant step, (5) skani/Mash for organism ID, (6) GoaT
  expected-count sanity check, (3) full organelle assemblers only when an
  organelle is fragmented.
- Nothing here should become a hard dependency or break the offline default.

---

## Verified sources (URLs actually fetched or read for these notes)

- NCBI FCS (FCS-GX, FCS-adaptor): https://github.com/ncbi/fcs
- oatk / OatkDB: https://github.com/c-zhou/oatk
- MitoHiFi: https://github.com/marcelauliano/MitoHiFi
- GetOrganelle: https://github.com/Kinggerm/GetOrganelle
- BlobToolKit (BlobTools2): https://github.com/genomehubs/blobtoolkit
- barrnap: https://github.com/tseemann/barrnap
- tidk (Telomere Identification toolKit): https://github.com/tolkit/telomeric-identifier
- minimap2 (asm5 assembly-to-assembly preset confirmed): https://github.com/lh3/minimap2 (README)
- skani: https://github.com/bluenote-1577/skani
- Mash: https://github.com/marbl/mash (via Mash / Mash Screen Genome Biology papers)
- sourmash prepared databases (NCBI Eukaryotes Jan 2025, GTDB): https://sourmash.readthedocs.io/en/latest/databases.html
- GoaT CLI (goat-cli, links to GoaT API docs): https://github.com/genomehubs/goat-cli (README)
- GoaT publication (chromosome number listed among taxon attributes): https://pmc.ncbi.nlm.nih.gov/articles/PMC9971660/
- NCBI Datasets genome assembly data report (`totalNumberOfChromosomes` field confirmed): https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-reports/genome-assembly/
- NCBI core_nt announcement: https://ncbiinsights.ncbi.nlm.nih.gov/2024/07/18/new-blast-core-nucleotide-database/
- CCDB (Chromosome Counts Database, plants): http://ccdb.tau.ac.il/ (via New Phytologist 2015, PMID 25423910)
- Animal Genome Size Database (C-values, not chromosome number): https://www.genomesize.com/
- Fusarium graminearum PH-1 four-chromosome reference: King et al. 2015, BMC Genomics 16:544, https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-015-1756-1

## Unverified / to check before implementing

- **GoaT exact attribute field name(s).** GoaT documentation describes a
  "chromosome number" taxon attribute, but I did not confirm the literal
  queryable field name (e.g. whether it is `chromosome_number`) or whether a
  distinct `haploid_number` attribute exists. The live GoaT API
  (goat.genomehubs.org/api-docs) or `goat-cli` should be queried directly. My
  attempt to hit the GoaT search API programmatically returned no body in this
  environment, so field names here are unconfirmed.
- **DIAMOND default outfmt columns** differ from the blastn outfmt 6 string
  detangler hard-codes; the exact column set for a protein path was not verified
  here and must be pinned before wiring results back.
- **oatk on a Flye GFA.** oatk documents GFA input via `-G` and coverage-tag
  options, but I did not verify it produces sensible output on a non-syncasm
  (Flye) fungal graph specifically; treat as needs-testing.
- **FCS-GX database size / RAM.** The ~470+ GB figure is from general knowledge
  of FCS-GX, not confirmed against the wiki in this session; check the current
  FCS-GX quickstart for exact requirements before recommending hardware.
- **Fungus-specific rRNA/organelle model coverage** in barrnap (`--kingdom fun`)
  and GetOrganelle (`fungus_mt`/`fungus_nr`) exists per their docs, but suitability
  for *Fusarium* specifically was not tested.
- **Chromosome-count semantics reconciliation** (cytogenetic haploid n vs NCBI
  `totalNumberOfChromosomes` including organelles vs detangler backbones) is a
  design problem, not a fact to verify, but must be solved before any automated
  comparison is trustworthy.
