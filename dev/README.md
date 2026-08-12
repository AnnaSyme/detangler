# dev

Tools for working on detangler. Not needed to run it.

## layout_metrics.py

Measures how legible a graph figure is, so "this looks crowded" becomes a
number.

```
python3 dev/layout_metrics.py detangler_paired.svg
```

Everything is dimensionless — distances in ribbon widths, angles in degrees —
so a figure of eleven contigs and one of six hundred are scored on the same
scale.

| metric | meaning |
|---|---|
| `wind` | mean degrees a contig turns along its length. **The metric that matches human judgement.** A gentle arc is ~30°; a contig folded on itself is 200+ |
| `wmax` | the worst single contig |
| `clear` | smallest gap between two contigs' centrelines, in ribbon widths. **Below 1.0 the ribbons overlap** |
| `xr` | ribbon × ribbon crossings. Must be 0 |
| `xrc` | ribbon × connector crossings. Must be 0 |
| `area` | the graph's bounding box, in ribbon widths squared |

Two things worth knowing before you trust it:

- **Measure `detangler_paired.svg`, not `detangler_graph.svg`.** They are laid
  out differently — `render_paired.py` makes the graph borrow the ideogram's
  scale — so one does not predict the other. The README uses the paired figure.
- A connector attaching to the end of a ribbon can register as an `xrc`
  crossing. Check how far the crossing is from the nearest ribbon end before
  believing it.

`hub_gap` is printed but not scored. It ranked a figure people call clear
*below* one they call crowded, so it does not measure what its name suggests.

### Use it as a regression check

`detangler_graph-layout-research_v1.md` states the invariant that ribbon
crossings must stay zero. Nothing enforced it, and the code was in fact drawing
ribbons at a clearance of 0.99 — overlapping — until August 2026. Run this after
any change to `render_graph.py`.

## layout_tune.py

Sweeps the layout constants, renders several graphs with each setting, measures
every figure, and ranks them.

```
python3 dev/layout_tune.py --repo . \
    --fast real_data/flye_assembly_graph.gfa:flye \
    --out /tmp/tuning
```

It scores each setting by its **worst** graph rather than the average, so a
setting that flatters one dataset and wrecks another is rejected. Give it more
than one graph or it will happily overfit to whichever one you supply.

Each trial gets its own copy of the source under `/tmp` and is edited there —
the repo is never modified. The constants it knows about are listed in `KNOBS`
at the top; if `render_graph.py` changes those lines it will say so rather than
silently doing nothing.

## What tuning these constants taught us

**Winding and separation pull against each other.** Pushing ribbons apart makes
them curl to fit, so raising the separation alone made the figure worse — 1.12
to 2.2 took winding from 78° to 104°. Never tune one without watching the other.

**Bead count, not bead spacing, controls curl.** A contig is drawn as a chain of
beads, and how much it can fold is set by how many it has. With a fixed spacing
a 12 Mb contig gets a long floppy chain and knots, while a 500 kb one gets three
beads and cannot bend at all.

**These constants were fitted on one graph.** The structures are scale-free —
a cap on a count, distances in ribbon widths — but the values came from
Fusarium. On a 197-segment Aspergillus graph they barely help: at that density
no spacing constant works, and the answer is to draw fewer contigs. Run the
harness on a new genome before trusting them there.
