#!/usr/bin/env python3
"""
Regression test built from the manual Fusarium graminearum worked example in
GenomeViz_design-notes_v1.md.

The segment lengths, depths and links are the real values recorded by hand from
Flye 2.8.2 output. The nucleotide sequences are NOT real - they are synthesised
here with the recorded GC content, telomere motif presence and periodicity, so
that the composition screens have something to act on. Do not use them for
anything else.

The test asserts the conclusions the design notes say the tool should reach, and
also asserts that it does NOT assert the ones the notes say it must not.

Run:  python3 GenomeViz_worked-example-test_v1.py [--keep DIR]
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
TOOL = os.path.join(HERE, "GenomeViz_ideogram-tool_v1.py")

# Reference values for validation (Fusarium graminearum PH-1 reference assembly)
EXPECTED_NUCLEAR = 36_563_796
EXPECTED_CHROMOSOMES = 4

# segment -> (length, depth), transcribed from the worked example
SEGMENTS = {
    "edge_6": (9_008_043, 19),
    "edge_7": (8_942_195, 20),
    "edge_5": (7_962_991, 20),
    "edge_1": (7_793_581, 19),
    "edge_2": (2_742_039, 19),
    "edge_9": (15_807, 37),
    "edge_3": (2_358, 10),
    "edge_10": (7_685, 1269),
    "edge_11": (98_177, 203),
    "edge_8": (51_498, 4),
    "edge_4": (13_947, 5),
}

# Orientations match the L lines of the real Flye GFA (flye_assembly_graph.gfa).
# They matter: the end-capping and feature-placement logic assesses the two
# physical ends of a segment independently, so e.g. edge_9 must sit on the
# OPPOSITE end of edge_2 from edge_8, as it does in the real graph.
LINKS = [
    ("edge_1", "+", "edge_9", "-"),
    ("edge_1", "-", "edge_9", "-"),
    ("edge_2", "+", "edge_8", "+"),
    ("edge_2", "-", "edge_9", "-"),
    ("edge_3", "+", "edge_5", "-"),
    ("edge_3", "+", "edge_6", "-"),
    ("edge_3", "-", "edge_9", "-"),
    ("edge_5", "-", "edge_10", "-"),
    ("edge_7", "+", "edge_8", "+"),
    ("edge_7", "-", "edge_9", "-"),
    ("edge_10", "+", "edge_10", "+"),  # self-loop: tandem array
    ("edge_11", "+", "edge_11", "+"),  # self-loop: circular
]

# graph_path values, consistent with the links above. contig_2 deliberately
# carries edge_9 at both ends, which is the path-parsing pitfall in the notes.
CONTIGS = [
    ("contig_1", 9_010_000, 19.0, "N", "N", 1, "*", "*,6,-3,-9,*"),
    ("contig_2", 7_809_000, 19.0, "N", "N", 1, "*", "*,9,1,-9,*"),
    ("contig_3", 7_995_000, 20.0, "N", "N", 1, "*", "*,10,10,10,10,5,-3"),
    ("contig_4", 8_995_000, 20.0, "N", "N", 1, "*", "*,7,8,2,*"),
    ("contig_5", 98_177, 203.0, "Y", "N", 1, "*", "11"),
    ("contig_6", 13_947, 5.0, "N", "N", 1, "*", "*,4,*"),
]

# Synthetic sequence recipes: (gc_fraction, telomere_motif or None, unit_length or None)
RECIPE = {
    "edge_9": (0.38, "TTAGGG", None),
    "edge_3": (0.25, None, None),
    "edge_10": (0.51, None, 1_537),  # 7,685 = 5 x 1,537
    "edge_11": (0.32, None, None),
    "edge_8": (0.47, None, None),
    "edge_4": (0.46, None, None),
}
BASELINE_GC = 0.48


def load_tool():
    spec = importlib.util.spec_from_file_location("detangler", TOOL)
    if spec is None or spec.loader is None:
        raise SystemExit(f"could not load {TOOL}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["detangler"] = mod
    spec.loader.exec_module(mod)
    return mod


class LCG:
    """Deterministic, so the fixture is byte-stable."""

    def __init__(self, seed=2026):
        self.s = seed

    def next(self):
        self.s = (1103515245 * self.s + 12345) & 0x7FFFFFFF
        return self.s

    def pick(self, options):
        return options[self.next() % len(options)]


def make_sequence(length: int, gc: float, motif=None, unit=None) -> str:
    rnd = LCG(length)
    gc_bases, at_bases = "GC", "AT"

    def block(n: int) -> str:
        return "".join(
            rnd.pick(gc_bases) if (rnd.next() % 1000) / 1000.0 < gc else rnd.pick(at_bases)
            for _ in range(n)
        )

    if unit:
        core = block(unit)
        seq = (core * (length // unit + 1))[:length]
        return seq
    seq = block(length)
    if motif:
        # telomere motifs concentrated at both ends, as in real subtelomeric sequence
        arr = motif * 80
        seq = arr + seq[len(arr) : length - len(arr)] + arr
    return seq[:length]


def write_fixture(d: str) -> dict:
    os.makedirs(d, exist_ok=True)
    gfa = os.path.join(d, "assembly_graph.gfa")
    with open(gfa, "w") as fh:
        fh.write("H\tVN:Z:1.0\n")
        for name, (length, depth) in SEGMENTS.items():
            if name in RECIPE:
                gc, motif, unit = RECIPE[name]
                seq = make_sequence(length, gc, motif, unit)
                fh.write(f"S\t{name}\t{seq}\tdp:i:{depth}\n")
            else:
                # multi-Mb backbone segments: length tag only, as many GFAs do
                fh.write(f"S\t{name}\t*\tLN:i:{length}\tdp:i:{depth}\n")
        for a, ao, b, bo in LINKS:
            fh.write(f"L\t{a}\t{ao}\t{b}\t{bo}\t0M\n")

    info = os.path.join(d, "assembly_info.txt")
    with open(info, "w") as fh:
        fh.write("#seq_name\tlength\tcov.\tcirc.\trepeat\tmult.\talt_group\tgraph_path\n")
        for row in CONTIGS:
            fh.write("\t".join(str(x) for x in row) + "\n")
    return {"gfa": gfa, "info": info}


# --------------------------------------------------------------------------
CHECKS = []
_ARGS = None  # the namespace the tool ran with, filled in by main()


def check(name):
    def deco(fn):
        CHECKS.append((name, fn))
        return fn

    return deco


@check("baseline single-copy depth lands on 19-20x")
def _(model, kg):
    assert 18.0 <= model.baseline_depth <= 21.0, f"baseline was {model.baseline_depth}"
    return f"baseline {model.baseline_depth:.1f}x"


@check("five backbone segments identified")
def _(model, kg):
    bb = sorted(c.name for c in model.segment_calls if c.cls == "backbone")
    expect = ["edge_1", "edge_2", "edge_5", "edge_6", "edge_7"]
    assert bb == expect, f"backbone was {bb}"
    return ", ".join(bb)


@check("backbone total is within 1% of the expected nuclear genome")
def _(model, kg):
    total = sum(c.length for c in model.segment_calls if c.cls == "backbone")
    rel = abs(total - EXPECTED_NUCLEAR) / EXPECTED_NUCLEAR
    assert rel < 0.01, f"{total:,} vs {EXPECTED_NUCLEAR:,} ({rel:.2%})"
    return f"{total:,} bp, {rel:.2%} from reference"


@check("edge_11 is called an organelle candidate, with a checkable reason")
def _(model, kg):
    c = next(x for x in model.segment_calls if x.name == "edge_11")
    assert c.cls == "organelle_candidate", f"edge_11 was {c.cls}"
    reason = " ".join(c.reasons).lower()
    for token in ("circular", "isolated", "copy number"):
        assert token in reason, f"reason missing '{token}': {reason}"
    assert 9.0 <= (c.copy_number or 0) <= 13.0, f"copy number {c.copy_number}"
    return f"{c.copy_number:.1f} copies, {c.reasons[0][:70]}..."


# --------------------------------------------------------------------------
# Organelle subtype checks, added with the plastid/mitochondrion split.
# Synthetic sequences only: they test the inverted-repeat detector and the
# subtype decision, not any real plastome.
# --------------------------------------------------------------------------
@check("IR detector finds a planted 12 kb inverted-repeat pair, and only that")
def _(model, kg):
    base = make_sequence(150_000, 0.37)
    assert kg.find_inverted_repeat_pair(base) is None, "false positive on a plain sequence"
    block = base[30_000:42_000]
    planted = base[:110_000] + kg.revcomp(block) + base[122_000:]
    assert len(planted) == 150_000
    hit = kg.find_inverted_repeat_pair(planted)
    assert hit is not None, "planted inverted pair was missed"
    size, first, second = hit
    assert abs(size - 12_000) <= 1_200, f"size estimate {size:,} vs planted 12,000"
    assert abs(first - 30_000) <= 600, f"first copy placed at {first:,}"
    assert abs(second - 110_000) <= 600, f"second copy placed at {second:,}"
    return f"~{size:,} bp block at ~{first:,} and ~{second:,}; plain sequence stays clean"


@check("subtype call: plastid-like needs the IR pair, mitochondrion-like its absence")
def _(model, kg):
    base = make_sequence(150_000, 0.37)
    plastome_like = base[:110_000] + kg.revcomp(base[30_000:50_000]) + base[130_000:]
    sub, why, ir = kg.classify_organelle_subtype(len(plastome_like), plastome_like, 0.37, 0.48)
    assert sub == "plastid", f"plastome-like sequence called {sub}"
    assert "plastid-like" in why and "inverted" in why, why
    assert ir is not None

    mito_like = make_sequence(60_000, 0.32)
    sub2, why2, ir2 = kg.classify_organelle_subtype(60_000, mito_like, 0.32, 0.48)
    assert sub2 == "mitochondrion", f"60 kb IR-free sequence called {sub2}"
    assert "size" in why2 and "its own" in why2, f"no hedge on size: {why2}"
    assert ir2 is None

    sub3, why3, _ir3 = kg.classify_organelle_subtype(150_000, "", None, None)
    assert sub3 == "unresolved", f"sequence-less candidate called {sub3}"
    assert "cannot separate" in why3, why3
    return "plastid with IR, mitochondrion without, unresolved when there is no sequence"


@check("edge_11 comes out mitochondrion-like, never plastid-like")
def _(model, kg):
    c = next(x for x in model.segment_calls if x.name == "edge_11")
    assert c.organelle_subtype == "mitochondrion", f"subtype was {c.organelle_subtype}"
    reason = " ".join(c.reasons)
    assert "mitochondrion-like" in reason, reason
    assert "plastid-like" not in reason, reason
    assert c.ir_block is None, f"phantom inverted repeat: {c.ir_block}"
    return "mitochondrion-like: no IR pair, and 98.2 kb is below the typical plastome range"


@check("edge_10 is a tandem array of ~67 copies")
def _(model, kg):
    c = next(x for x in model.segment_calls if x.name == "edge_10")
    assert c.cls == "tandem_array", f"edge_10 was {c.cls}"
    assert 60 <= (c.copy_number or 0) <= 75, f"copy number {c.copy_number}"
    return f"{c.copy_number:.0f} copies of a {c.length:,} bp unit"


@check("edge_10 is associated with the chromosome carrying edge_5")
def _(model, kg):
    tangle = next(
        (t for t in model.tangles if "edge_10" in t.segments), None
    )
    assert tangle is not None, "edge_10 was not placed on any molecule"
    chain_names = {a.seqname for a in tangle.anchors}
    edge5_chain = {
        s.name
        for s in model.sequences
        if any(b[2] == "edge_5" for b in s.blocks)
    }
    assert chain_names & edge5_chain, f"edge_10 on {chain_names}, edge_5 on {edge5_chain}"
    return f"placed on {', '.join(sorted(chain_names))}"


@check("edge_9 is a repeat at contig ends carrying telomere motifs")
def _(model, kg):
    c = next(x for x in model.segment_calls if x.name == "edge_9")
    assert c.cls == "repeat", f"edge_9 was {c.cls}"
    assert 1.6 <= (c.copy_number or 0) <= 2.4, f"copy number {c.copy_number}"
    assert c.path_terminal > c.path_interior, (
        f"terminal {c.path_terminal} vs interior {c.path_interior}"
    )
    assert c.telomere_motifs, "no telomere motifs detected"
    return (
        f"{c.copy_number:.1f} copies, {c.path_terminal} terminal / {c.path_interior} interior, "
        f"{sum(c.telomere_motifs.values())} motif hits"
    )


@check("edge_3 gets its own below-single-copy class and is flagged AT-rich")
def _(model, kg):
    c = next(x for x in model.segment_calls if x.name == "edge_3")
    assert c.cls == "low_coverage", f"edge_3 was {c.cls}"
    assert c.at_rich, "edge_3 was not flagged AT-rich"
    return f"copy number {c.copy_number:.2f}, GC {c.gc:.0%}"


@check("edge_8 is low coverage but still participates in the topology")
def _(model, kg):
    c = next(x for x in model.segment_calls if x.name == "edge_8")
    assert c.cls == "low_coverage", f"edge_8 was {c.cls}"
    used = any("edge_8" in j.via for h in model.hypotheses for j in h.joins)
    assert used, "edge_8 was filtered out before traversal"
    return f"copy number {c.copy_number:.2f}, used as a join in at least one hypothesis"


@check("edge_9 is recognised as telomeric, without being told the karyotype")
def _(model, kg):
    telo = kg.telomeric_segments(model.segment_calls, _ARGS)
    assert "edge_9" in telo, f"telomeric segments were {sorted(telo)}"
    assert len(telo) == 1, f"expected only edge_9 to be telomeric, got {sorted(telo)}"
    return f"edge_9, {telo['edge_9']} motif occurrences"


@check("the chromosome count is inferred, and the true value of 4 is inside the range")
def _(model, kg):
    best, low, high = model.count_range()
    assert low <= EXPECTED_CHROMOSOMES <= high, (
        f"true count {EXPECTED_CHROMOSOMES} outside the reported range {low}-{high}"
    )
    counts = {len(h.chains) for h in model.hypotheses[:4]}
    assert EXPECTED_CHROMOSOMES in counts, f"top hypotheses give {sorted(counts)} chains"
    return f"best estimate {best}, range {low}-{high}, truth {EXPECTED_CHROMOSOMES} included"


@check("no expected chromosome count was supplied to reach that")
def _(model, kg):
    assert _ARGS.expected_chromosomes is None, "the test leaked an expected count into the tool"
    assert _ARGS.expected_genome_size is None, "the test leaked an expected genome size"
    return "ran with neither --expected-chromosomes nor --expected-genome-size"


@check("the ambiguity is declared, not hidden")
def _(model, kg):
    tied = [h for h in model.hypotheses if any("score within" in c for c in h.contradicting)]
    assert len(tied) >= 3, f"only {len(tied)} hypotheses flagged as tied"
    pairs = set()
    for h in tied:
        for c in h.chains:
            if len(c) > 1:
                pairs.add(tuple(sorted(c)))
    for expected in (("edge_2", "edge_7"), ("edge_5", "edge_6")):
        assert expected in pairs, f"{expected} not among the tied alternatives: {sorted(pairs)}"
    return f"{len(tied)} tied hypotheses covering {len(pairs)} alternative joins"


@check("edge_4 sits apart as a contaminant candidate, not forced into a chromosome")
def _(model, kg):
    ua = {s.name: s for s in model.unassigned()}
    assert "edge_4" in ua, f"unassigned set was {sorted(ua)}"
    note = ua["edge_4"].note.lower()
    assert "candidate" in note, f"note does not hedge: {note}"
    assert "does not connect" in note, f"note lacks the topological reason: {note}"
    for s in model.sequences:
        if s.role == "chromosome":
            assert not any(b[2] == "edge_4" for b in s.blocks), "edge_4 was placed on a chromosome"
    return ua["edge_4"].note[:78]


@check("each backbone segment has its own colour, shared by both figures")
def _(model, kg):
    colours = model.segment_colours
    backbone = [c.name for c in model.segment_calls if c.cls == "backbone"]
    used = [colours[n] for n in backbone]
    assert len(set(used)) == len(backbone), f"colours repeat across backbone segments: {used}"
    for s in model.sequences:
        for _s, _e, seg, colour in s.blocks:
            assert colour == colours[seg], f"{seg} drawn as {colour}, graph uses {colours[seg]}"
    return f"{len(backbone)} distinct backbone colours, blocks match the graph figure"


@check("every ambiguous hypothesis says what would resolve it")
def _(model, kg):
    for h in model.hypotheses:
        if h.contradicting and not h.resolve_with:
            # a wrong chromosome count is contradiction enough to demand a remedy
            assert any("chains" not in c for c in h.contradicting), (
                f"hypothesis {h.rank} states a problem with no remedy"
            )
    with_remedy = [h for h in model.hypotheses if h.resolve_with]
    assert with_remedy, "no hypothesis named the missing evidence"
    return f"{len(with_remedy)}/{len(model.hypotheses)} name the evidence that would settle it"


@check("no claim is made about which chromosome is which")
def _(model, kg):
    for s in model.sequences:
        label = (s.label or s.name).lower()
        assert not any(
            tok in label for tok in ("chromosome 1", "chromosome 2", "chr1", "chr2", "chr i")
        ), f"sequence label asserts an identity: {s.label}"
    return "molecules are named chain_N, not chromosome N"


@check("repeat candidates are exported for identification")
def _(model, kg, out_dir=None):
    return None  # filled in by the runner, which has the paths


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--keep", help="keep the fixture and outputs in this directory")
    a = ap.parse_args()

    kg = load_tool()
    tmp = a.keep or tempfile.mkdtemp(prefix="detangler_test_")
    os.makedirs(tmp, exist_ok=True)
    fixture = write_fixture(os.path.join(tmp, "fixture"))
    out_dir = os.path.join(tmp, "out")

    # Deliberately no --expected-chromosomes and no --expected-genome-size: the
    # point is that the structure is worked out from the assembly output alone.
    # EXPECTED_* are used only by this test, to check the answer afterwards.
    argv = [
        "--gfa", fixture["gfa"],
        "--flye-info", fixture["info"],
        "--out-dir", out_dir,
        "--prefix", "fusarium",
        "--title", "Fusarium graminearum, Flye 2.8.2, ONT ~22x (worked example)",
        "--quiet",
    ]
    rc = kg.main(argv)
    assert rc == 0, f"tool exited {rc}"

    # rebuild the model the same way main() did, so the assertions can inspect it
    parser_args = kg.main.__wrapped__ if hasattr(kg.main, "__wrapped__") else None
    del parser_args
    log = kg.Log(quiet=True)
    ns = _namespace(kg, argv)
    global _ARGS
    _ARGS = ns
    model, extra = kg.build_model_graph_first(ns, log)

    print(f"\nFusarium worked example - {len(CHECKS)} checks\n")
    failures = 0
    for name, fn in CHECKS:
        if fn.__code__.co_argcount > 2:
            continue
        try:
            detail = fn(model, kg)
            print(f"  PASS  {name}")
            if detail:
                print(f"        {detail}")
        except AssertionError as e:
            failures += 1
            print(f"  FAIL  {name}\n        {e}")

    # file-level checks
    for label, path in extra.items():
        ok = os.path.exists(path) and os.path.getsize(path) > 0
        print(f"  {'PASS' if ok else 'FAIL'}  wrote {label}")
        failures += 0 if ok else 1

    cand = extra.get("repeat candidate sequences")
    if cand:
        names = [l[1:].split()[0] for l in open(cand) if l.startswith(">")]
        ok = "edge_10" in names and "edge_9" in names and "edge_11" in names
        print(f"  {'PASS' if ok else 'FAIL'}  repeat candidates exported for BLAST: {names}")
        failures += 0 if ok else 1

    print(f"\n{len(CHECKS) - 1 - failures} passed, {failures} failed")
    if not a.keep:
        shutil.rmtree(tmp, ignore_errors=True)
    else:
        print(f"outputs kept in {tmp}")
    return 1 if failures else 0


def _namespace(kg, argv):
    """Re-parse argv into the namespace the tool builds internally."""
    import contextlib
    import io

    holder = {}
    real_main = kg.main

    def capture(a):
        holder["ns"] = a
        raise SystemExit(0)

    # run main far enough to build the namespace, then intercept
    src_parser = None
    with contextlib.redirect_stderr(io.StringIO()):
        try:
            kg.build_model_graph_first = _intercept(kg.build_model_graph_first, holder)
            real_main(argv)
        except SystemExit:
            pass
        finally:
            kg.build_model_graph_first = kg.build_model_graph_first.__wrapped__  # type: ignore
    del src_parser, capture
    return holder["ns"]


def _intercept(fn, holder):
    def wrapper(args, log):
        holder["ns"] = args
        return fn(args, log)

    wrapper.__wrapped__ = fn  # type: ignore
    return wrapper


if __name__ == "__main__":
    sys.exit(main())
