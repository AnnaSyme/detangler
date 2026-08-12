#!/usr/bin/env python3
"""
Tune detangler's graph layout against measured legibility, not taste.

The problem this solves: the layout is governed by a handful of constants, and
picking them by looking at one figure overfits to that dataset. So instead -
sweep the constants, render several graphs of different sizes with each setting,
measure every figure with layout_metrics.py, and rank by the score ACROSS all
graphs. A setting only wins if it is good everywhere.

It never touches your repo. Each trial gets its own copy of the source under
/tmp, edited there.

Two stages, because the big graphs are slow:

    stage 1   sweep everything on the fast graphs
    stage 2   take the best few and check them on the slow ones

A setting that wins stage 1 and loses stage 2 is exactly the overfit this is
built to catch.

Usage
    python3 layout_tune.py --repo ~/genomeviz \\
        --fast  ~/genomeviz/real_data/flye_assembly_graph.gfa:flye \\
        --slow  ~/genomeviz-data/arabidopsis/ath_p_ctg.gfa:hifiasm \\
        --out   ~/genomeviz-data/tuning

Each graph is given as PATH:ASSEMBLER.
"""

from __future__ import annotations

import argparse
import itertools
import os
import shutil
import subprocess
import sys
import time
from typing import Dict, List, Optional, Tuple

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import layout_metrics as LM            # noqa: E402


# ---------------------------------------------------------------------------
# the constants being tuned
# ---------------------------------------------------------------------------
# Each entry is (label, exact source line, template with {v}). They are edited
# as text in a throwaway copy of render_graph.py, which keeps this script
# independent of the tool's internals - nothing has to be refactored to make it
# tunable, and the repo is never modified.
#
# Why these three:
#   spacing  sets how many beads a contig is made of. Fewer, wider-spaced beads
#            make a stiffer ribbon that cannot fold back on itself - the direct
#            lever on WINDING, which is the metric that tracks human judgement.
#   w_sep    the hard minimum distance between ribbons, in ribbon widths. The
#            direct lever on CLEARANCE.
#   rep_cut  how far repulsion reaches, in bead widths. Sets whether distant
#            contigs shove each other and inflate the canvas.

KNOBS = {
    "spacing": ("spacing = max(segment_thickness() * 0.42, 8.0)",
                "spacing = max(segment_thickness() * {v}, 8.0)",
                [0.42, 0.60, 0.85]),
    "w_sep":   ("w_sep = segment_thickness() * 1.12",
                "w_sep = segment_thickness() * {v}",
                [1.12, 1.6, 2.1]),
    "rep_cut": ("if d > k * 14.0:",
                "if d > k * {v}:",
                [14.0, 8.0, 22.0]),
}


def patch(src_root: str, settings: Dict[str, float]) -> None:
    p = os.path.join(src_root, "src", "detangler", "render_graph.py")
    s = open(p).read()
    for name, value in settings.items():
        old, tmpl, _ = KNOBS[name]
        new = tmpl.format(v=value)
        if old not in s:
            raise SystemExit(
                f"could not find the line for '{name}' in render_graph.py:\n  {old}\n"
                f"The source has changed - update KNOBS in this script.")
        s = s.replace(old, new)
    open(p, "w").write(s)


def run_one(repo: str, settings: Dict[str, float], graphs: List[Tuple[str, str]],
            workdir: str, outdir: str, label: str) -> Optional[Dict[str, Dict]]:
    """Render every graph with one setting, and measure each figure."""
    trial = os.path.join(workdir, label)
    shutil.rmtree(trial, ignore_errors=True)
    os.makedirs(trial, exist_ok=True)
    shutil.copytree(os.path.join(repo, "src"), os.path.join(trial, "src"))
    shutil.copy(os.path.join(repo, "detangler.py"), trial)
    patch(trial, settings)

    out: Dict[str, Dict] = {}
    for gfa, assembler in graphs:
        name = os.path.basename(gfa).split(".")[0]
        dest = os.path.join(outdir, label, name)
        os.makedirs(dest, exist_ok=True)
        r = subprocess.run(
            [sys.executable, "detangler.py", "--gfa", gfa, "--assembler", assembler,
             "--out-dir", dest, "--quiet"],
            cwd=trial, capture_output=True, text=True, timeout=1800)
        svg = os.path.join(dest, "detangler_graph.svg")
        if r.returncode != 0 or not os.path.exists(svg):
            print(f"    {name}: FAILED ({r.stderr.strip().splitlines()[-1:] or ''})")
            return None
        out[name] = LM.measure(LM.read_figure(svg))
    return out


def combined(res: Dict[str, Dict], refs: Dict[str, float]) -> float:
    """
    One number for a setting: the WORST graph's score.

    Deliberately the worst rather than the mean. A setting that renders the
    fungal graph beautifully and the plant graph as a hairball is not a good
    setting, and averaging would hide that.
    """
    if not res:
        return -99.0
    return min(LM.score(m, refs.get(name)) for name, m in res.items())


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--fast", nargs="+", required=True, help="PATH:ASSEMBLER, swept fully")
    ap.add_argument("--slow", nargs="*", default=[], help="PATH:ASSEMBLER, used to check the winners")
    ap.add_argument("--out", required=True)
    ap.add_argument("--keep", type=int, default=3, help="how many winners go to stage 2")
    args = ap.parse_args()

    def parse(xs):
        out = []
        for x in xs:
            path, _, asm = x.rpartition(":")
            out.append((os.path.expanduser(path), asm or "unknown"))
        return out

    repo = os.path.expanduser(args.repo)
    fast, slow = parse(args.fast), parse(args.slow)
    outdir = os.path.expanduser(args.out)
    workdir = os.path.join(outdir, "_work")
    os.makedirs(workdir, exist_ok=True)

    grid = [dict(zip(KNOBS, combo))
            for combo in itertools.product(*(v[2] for v in KNOBS.values()))]
    print(f"{len(grid)} settings x {len(fast)} fast graph(s)\n")

    # The baseline is the current code, and every area is scored relative to it,
    # so "made the picture 3x bigger" is visible as a cost rather than free.
    base = {k: v[2][0] for k, v in KNOBS.items()}
    print("baseline (current settings)")
    b = run_one(repo, base, fast + slow, workdir, outdir, "baseline")
    if b is None:
        return 1
    refs = {name: m["area_u2"] for name, m in b.items()}
    for name, m in b.items():
        print(f"    {name:<22} wind {m['wind_mean']:5.0f}  clear {m['clearance']:5.2f}  "
              f"score {LM.score(m, refs.get(name)):5.2f}")

    results = []
    for i, settings in enumerate(grid, 1):
        label = "t%02d" % i
        desc = "  ".join(f"{k}={v}" for k, v in settings.items())
        t0 = time.time()
        res = run_one(repo, settings, fast, workdir, outdir, label)
        sc = combined(res, refs)
        results.append((sc, settings, res, label))
        print(f"[{i:2d}/{len(grid)}] {desc:<52} worst-score {sc:6.2f}  "
              f"({time.time()-t0:.0f}s)")

    results.sort(key=lambda r: -r[0])
    print("\n--- stage 1, best on the fast graphs ---")
    for sc, settings, res, label in results[:args.keep]:
        print(f"  {sc:6.2f}  " + "  ".join(f"{k}={v}" for k, v in settings.items()))
        for name, m in res.items():
            print(f"            {name:<20} wind {m['wind_mean']:5.0f} "
                  f"clear {m['clearance']:5.2f} xr {m['r_cross']}")

    if not slow:
        print("\nNo --slow graphs given, so nothing has been checked for overfitting.")
        return 0

    print(f"\n--- stage 2, checking the top {args.keep} on the slow graphs ---")
    final = []
    for sc, settings, _res, label in results[:args.keep]:
        res = run_one(repo, settings, slow, workdir, outdir, label + "_slow")
        s2 = combined(res, refs)
        final.append((min(sc, s2), sc, s2, settings))
        print(f"  fast {sc:6.2f}   slow {s2:6.2f}   "
              + "  ".join(f"{k}={v}" for k, v in settings.items()))
        if res:
            for name, m in res.items():
                print(f"            {name:<20} wind {m['wind_mean']:5.0f} "
                      f"clear {m['clearance']:5.2f} xr {m['r_cross']}")

    final.sort(key=lambda r: -r[0])
    best = final[0]
    print("\n=== winner (best worst-case across every graph) ===")
    for k, v in best[3].items():
        print(f"  {KNOBS[k][1].format(v=v)}")
    print(f"  fast {best[1]:.2f}   slow {best[2]:.2f}")
    print("\nNothing has been written to the repo. Apply by hand if you like it.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
