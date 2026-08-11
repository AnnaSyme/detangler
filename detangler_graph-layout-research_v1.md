# Making the graph panel look drawn, not simulated

Research pass, 11 Aug 2026. **R1-R4 implemented the same day** (see Status
below); R5-R9 are still open.

## Status

| | done | note |
|---|---|---|
| R1 tangent-continuous connectors | yes | replaced the centroid bow |
| R2 clamped chain end tangents | yes | damped rotation, 6 alternating passes |
| R3 hub fanning by PAVA | yes | min gap 1.15/2 rad ~ 33 deg |
| R4 graphviz K fix | yes | K = spacing/72; sfdp kept, crossings re-measured |
| R5 Schneider Bezier fitting | no | |
| R6 coarser beads | no | needs R5 first |
| R7 curvature ceiling | no | |
| R8 centripetal Catmull-Rom | no | superseded by R5 if R5 is done |
| R9 raise MIN_DRAWN_PX | no | shared with the chromosome panel |

Follow-ups after R1-R4, same day: the rigid end rotation in R2 was replaced by
a DECAYING profile applied as nested suffix rotations, because rotating the
first few beads rigidly made a straight stretch at each end - a hairpin read as
two straight legs and a bend rather than as one curve. Connector stroke reduced
to 0.17 w. Depth labels now repel from the drawing, not just from other labels.
The panel rotation searches 90 angles and scores the bounding box, targeting a
wide box AND pushing ink out of the lower-right triangle, which is where the
chromosome panel sits.

Measured after R1-R4 on the real Flye GFA: drawn px per bp varies **1.48%**
across the five contigs above the length floor, so the length invariant holds.
Ribbon x ribbon crossings **0**, ribbon x connector crossings **0**; the 13
connector x connector crossings are all at the 5-way edge_9 hub, where
connectors share an endpoint and must overlap.

## The measurement that reframes the problem

Every ribbon in `detangler_figure_v2.svg` was measured. Nine linear contigs:

| contig | beads | drawn length | length / ribbon width | total turn | sd of turn / mean |
|---|---|---|---|---|---|
| 6 | 22 | 577 px | 9.0 | -45 deg | 0.40 |
| 7 | 21 | 538 px | 8.4 | -157 deg | 0.30 |
| 5 | 19 | 503 px | 7.9 | -63 deg | 0.32 |
| 1 | 18 | 454 px | 7.1 | -188 deg | 0.80 |
| 2, 8, 9, 4, 3 | 3-5 | 77-125 px | 1.2-2.0 | ~0 | ~0 |

Two facts follow.

**Three of the four long contigs are circular arcs.** Turn-per-bead is constant
to within 30-40%. A circular arc has zero curvature variation, which is the most
machine-looking shape there is.

**The five "straight" contigs are pinned at the minimum-length floor and cannot
curve.** `MIN_DRAWN_PX = RIBBON_W * 2.2 = 140.8`, minus `_trim_polyline(w/2)` at
each end, is 76.8 px - exactly the measured 76.7-76.9. At 1.2 ribbon widths there
is no room for a bend. That complaint has no layout fix.

## Three real bugs

**1. `sfdp` silently ignores `len=`.** Tested on graphviz 2.43: with `sfdp`, a
chain with one `len=3.0` edge among `len=0.2` edges comes out uniform; with
`neato` the long edge is 3.24 against 0.19. The Graphviz docs agree - `len` is
marked *neato, fdp only*. Chain length still survives through bead COUNT, so
nothing is grossly mis-scaled, but the absolute scale is sfdp's, not ours, so
`-Gsep=+18` is being applied against a spacing we did not choose. `-GK=1.4` is
in inches = 100.8 px, about 3.7x the intended bead spacing.

**2. The smoothing loop converges to circular arcs by construction.** Running the
exact loop (18 rounds, 0.62/0.38, 3 length sub-passes) on a wiggly chain drives
sd/|mean| of the turn angle from 47 to 5.3, and to 2.7 at 200 rounds. Mean
curvature is preserved; curvature VARIATION is destroyed. Wobble and organic
character are the same quantity at different frequencies, so this dial cannot
separate them. More rounds makes it worse, not better.

**3. The connectors bow away from the global centroid**, so their direction is a
function of where they sit in the picture rather than of which way the contigs
point. That is the source of every kink - and it is fake in the same way the
rejected sinusoid was fake.

## Framing error

Bandage drawings look organic because they are CROWDED: `nodeLengthPerMegabase
= 1000`, `averageNodeWidth = 5`, so a 9 Mb contig is ~1800 widths long - thin
spaghetti, hundreds of strands shoving each other. Ours is 7-9 widths. Ten fat
sausages in an empty canvas. A force model in an empty room correctly outputs
straight lines and gentle arcs.

**So the layout algorithm is not the problem. The wins are all in the drawing.**

## What Bandage actually does (read from source)

- **Layout**: OGDF FMMM, `unitEdgeLength(1.0)` with real lengths carried in a
  per-edge `EdgeArray`, which OGDF honours. `stepsForRotatingComponents(50)`,
  `minDistCC`. [graphlayoutworker.cpp](https://raw.githubusercontent.com/rrwick/Bandage/main/program/graphlayoutworker.cpp)
- **Contig rendering**: `remakePath()` is `moveTo` then `lineTo` in a loop.
  **No smoothing at all.** It reads as smooth only because segments are 4 ribbon
  widths long and there are hundreds of them.
  [graphicsitemnode.cpp](https://raw.githubusercontent.com/rrwick/Bandage/main/graph/graphicsitemnode.cpp)
- **Connector rendering** - the important one. Each control point sits one
  ribbon width beyond the contig end, ALONG THE CONTIG'S OWN FINAL TANGENT:
  ```cpp
  extensionLength = min(edgeLength, edgeDistance / 2.0);
  controlPoint1 = extendLine(beforeStart, start, extensionLength);
  controlPoint2 = extendLine(afterEnd,    end,   extensionLength);
  path.cubicTo(controlPoint1, controlPoint2, end);
  ```
  G1-continuous with both ribbons. [graphicsitemedge.cpp](https://raw.githubusercontent.com/rrwick/Bandage/main/graph/graphicsitemedge.cpp)

Ratios worth carrying over, as multiples of ribbon width `w`:

| | Bandage | ours now | suggested |
|---|---|---|---|
| bead spacing | 4.0 w | 0.42 w | 0.9-1.2 w |
| connector control extension | 1.0 w | n/a (centroid bow) | 1.0 w |
| connector stroke | 0.30 w | 0.15 w | 0.25-0.30 w |
| connector colour | 70% black | opaque `#2b2b2b` | grey 60-70% |

## Recommendations, ranked by value / cost

**R1. Rebuild connectors as tangent-continuous cubic Beziers.** Bandage's exact
construction, with one deviation: measure the tangent over ~0.75 w of arclength
rather than one bead, because our beads are 0.42 w apart and a single-bead
tangent is noisy. Deletes every kink in the figure. ~25 lines, no deps, layout
untouched.

**R2. Clamp the chain end tangents with a ghost bead.** Currently the smoothing
holds end POSITIONS and leaves end TANGENTS free, which is why the ends
straighten and the middle becomes a uniform arc. Pin the direction too: after
each smoothing round, rigidly rotate the first ~3 beads about the pinned end so
the chord lies along the target direction, then re-run the length passes (a
rotation is an isometry, so drawn length is preserved exactly). A chain with
pinned endpoints AND pinned end tangents minimising bending energy is an Euler
elastica, whose curvature VARIES along the stroke - the confident drawn sweep,
derived from real geometry rather than decoration. ~35 lines, no deps.

**R3. Fan the hub by isotonic regression.** For k contig ends meeting at a hub
at natural angles theta_i, solve `min sum (theta'_i - theta_i)^2` subject to a
minimum angular gap, via pool-adjacent-violators. Minimum gap from geometry: two
ribbons of width w are clear at radial distance d when the gap exceeds
`1.15 w / d`; at d = 2w that is ~33 deg. Because PAVA moves angles as little as
possible, the fan is not mechanically symmetric. Feeds the directions to R1 and
R2, and lets the junction dots shrink or go. ~50 lines, no deps. This is
engineering, not a cited method - PAVA-with-minimum-gap is the standard tool for
the identical 1-D axis-label problem.

**R4. Fix the graphviz call.** Set `K = spacing/72` so the initial layout is at
the scale everything downstream assumes. If `len` is wanted, `neato` must be used
(its default `mode=major` is stress majorization and honours `len`). The code
comment claiming sfdp beat neato 0 crossings to 2 was measured while sfdp was
silently discarding `len`, so it should be redone. Note `neato -Gmode=sgd` is
NOT available on graphviz 2.43. ~5 lines.

**R5. Re-draw each contig as 1-3 fitted cubic Beziers** instead of plotting 22
beads, via Schneider's algorithm with prescribed end tangents from R2/R3,
tolerance ~w/8. Converts "plot of a simulation" into "drawn with a pen".
Watch the length invariant: a fitted cubic is shorter than its polyline (0.01-0.3%
at our turn angles, 2.3% at 30 deg/bead), so measure and re-fit. ~120 lines, no
deps. [Graphics Gems FitCurves.c](https://github.com/erich666/GraphicsGems/blob/master/gems/FitCurves.c)

**R6. Coarsen beads to ~1 ribbon width.** The repulsion length scale is currently
0.48 w, less than one ribbon width, so contigs cannot push each other at the
scale where they visibly overlap. Also quarters the O(n^2) loops. **Do not ship
without R5** - at 9 beads the polyline is visible. Do not go to Bandage's 4 w;
our contigs are 8 widths long, not 1800.

**R7. Curvature ceiling in the constraint loop.** Cap turn per bead at `s/w`
radians so a ribbon can never eat itself on the inside of a bend. Nothing pinches
today, but R2/R3/R6 could cause it. ~20 lines.

**R8. Centripetal Catmull-Rom** (alpha = 0.5) instead of the current smoothing -
only if R5 is skipped. Centripetal is provably the only parameterisation in the
family with no cusps and no self-intersection within a segment (Yuksel, Schaefer
& Keyser, CAD 43(7):747-755, 2011). At current bead density it is visually
indistinguishable from what we have; it matters only after R6.

**R9. Stop expecting the floored contigs to look organic.** Either accept them as
tokens (short square-capped bars, which is already what happens), or raise
`MIN_DRAWN_PX` from 2.2 w to ~3.5 w so they can carry one gentle bend. Note
`MIN_DRAWN_PX` is shared with the chromosome panel by design, so this changes
both.

## Explicitly rejected

- **FMMM via `ogdf-python`** - heavy dep (cppyy + compiled OGDF wheel). FMMM
  differs from sfdp in coarsening, multipole repulsion, fine-tuning and component
  packing, none of which matters on an 11-node graph. Its one real advantage,
  per-edge desired lengths, is available free from `neato`.
- **Stress majorization / SMACOF on the bead graph** - would make it worse.
  Stress wants beads at graph-distance d to sit at Euclidean distance d, which is
  achievable only if the chain is STRAIGHT. It optimises directly for the thing
  being complained about.
- **PG-SGD (ODGI-style)** - same objection; it is stress with path-nucleotide
  distances. ODGI looks tangled because thousands of paths pull in conflicting
  directions. With 11 contigs there is no conflict.
  [odgi layout](https://odgi.readthedocs.io/en/latest/rst/commands/odgi_layout.html),
  [Guarracino et al. 2022](https://doi.org/10.1093/bioinformatics/btac308)
- **Force-directed edge bundling** - we have 10 links. Also: none of the seven
  genome-graph viewers surveyed does any bundling.
- **Chaikin corner-cutting** - `_smooth_path` already emits the Chaikin limit
  curve (the quadratic B-spline of the bead polygon). Strictly zero gain.
- **More smoothing rounds** - provably counterproductive, see bug 2. If anything
  reduce to ~8 once R2/R5 are in.
- **Any re-introduction of a global-coordinate bow** - same species as the
  rejected sinusoid.

## Other viewers, for the record

Bandage-NG: still FMMM, adds parallel per-component layout and OGDF rectangle
packing. gfaestus and Waragraph: compute no layout, consume `odgi layout` TSV,
draw one quad per node, no edges. GfaViz: OGDF Stress Minimization by default,
FMMM via `--fmmm`. VRPG: reference nodes pinned on a coordinate line, bubbles by
d3-force, Straight/Curved toggle. panGraphViewer: vis-network under 200 nodes,
else cytoscape.js fcose. MoMI-G: no force layout, metro-map ordering plus a
hand-written d3 Circos ring.

## Suggested order

R1 -> R3 -> R2 -> R4 first: mutually reinforcing, ~110 lines, no new deps, and
they turn every join into a single continuous stroke while stopping the long
contigs being circular arcs. Then R5 -> R6 -> R7. R8 only if R5 is skipped. R9
is a five-minute decision that needs the chromosome panel checked.

## Sources credited in the code and README

Each of these is named at the point of use in `src/detangler/render_graph.py`
and listed in the README's Credits section.

- Bandage, `GraphicsItemEdge::calculateAndSetPath()` - the connector
  construction. Also `stepsForRotatingComponents` for the per-component rotation.
  https://github.com/rrwick/Bandage
- Yifan Hu, "Efficient and high quality force-directed graph drawing",
  Mathematica Journal 10(1), 2005 - the algorithm behind graphviz `sfdp`.
  http://yifanhu.net/PUB/graph_draw_small.pdf
- Ayer, Brunk, Ewing, Reid & Silverman (1955), Ann. Math. Statist. 26:641-647 -
  pool-adjacent-violators, used to fan hub angles.
  https://projecteuclid.org/euclid.aoms/1177728423
- Levien & Sequin, "Interpolating splines: which is the fairest of them all?",
  CAD & Applications 6 (2009) - framing for why a clamped-endpoint
  minimum-bending-energy curve is an elastica rather than an arc.
  https://people.eecs.berkeley.edu/~sequin/PAPERS/2009_CAD_Levien_Sequin.pdf

Not used, so not credited: Schneider's curve fitting (R5), centripetal
Catmull-Rom (R8), OGDF/FMMM, ODGI PG-SGD, edge bundling.

## The interlocking-triangle layout: what worked and what did not

The idea: cut the canvas on a diagonal, graph in the upper-left triangle,
chromosomes filling the lower-right. The chromosome row is sorted short to tall
on a shared baseline, so its silhouette IS a rising staircase - a lower-right
triangle - which is why the two shapes ought to interlock.

**Done.** Row order is now unplaced, then organelle, then chromosomes ascending,
with the unplaced panel moved from the right-hand end to the left. Chromosome
numbering is by SIZE RANK rather than draw order, so chr 1 is still the largest
even though it now sits at the right. The panel is right-aligned and scaled up
to the graph panel's width (capped at 2x), which roughly doubled the bars. The
graph's rotation search penalises ink past the diagonal: mean overshoot falls
from 0.071 to 0.014 at a weight of 8.

**Attempted and backed out.** A hard keep-out was tried: project every bead
past the diagonal back across it, then restore chain lengths, repeatedly. It
does clear the corner, but bolted on AFTER the solver it fights the length pass
and loses - contigs came out visibly compressed and overlapping each other and
the organelle ring. It survives behind `--graph-triangle`, off by default, as a
record of what does not work. Doing it properly means adding the keep-out to the
main constraint loop alongside repulsion and separation, so all three negotiate
in the same iteration.

**Where it landed.** The COMPOSITION is diagonal: one canvas, graph upper-left,
chromosomes lower-right, overlapping vertically rather than stacked. The graph
is only nudged out of the corner, by the rotation search, so the composition
checks what it actually left there - `_ink_bottom_in_column` - and grows the
canvas if the two would collide. When the graph vacates the corner the panels
interlock and the figure is short; when it does not, it degrades back towards a
stack instead of drawing them on top of each other.

**Not fully solved.** On this graph one long contig still hangs into the lower
right at every angle, so the panels are still mostly The chromosome panel is allowed to rise into whatever the graph leaves
empty above it - `_ink_bottom_in_column` measures how far the graph's ink
actually reaches in the column the chromosomes occupy - but on this graph the
answer is "all the way down", because one long contig hangs into the lower
right whatever angle is chosen. Rotation can only reorient the shape the spring
layout produced; it cannot reshape it.

Genuinely interlocking would mean giving the FORCE MODEL a keep-out region -
the chromosome triangle as an obstacle the beads are pushed out of, in the same
constraint-projection pass that already enforces chain lengths and separation.
That is a real change, not a scoring tweak, and it is not attempted here.

## Composition: two panels, side by side

Settled here after trying a stack, then a diagonal overlap. Graph left,
chromosomes right, one heading CENTRED across the top reading "Assembly graph ->
Possible chromosomes" so the figure states what it is FOR rather than captioning
two halves separately, and a dotted rule around the whole thing including the
heading - without it the two panels read as two images that happen to have been
saved together, and the heading looks like a caption for whichever one it sits
nearest.

**Why the clever versions were abandoned.** Stacking, and then overlapping the
panels on a diagonal, both saved the white space each panel's empty corner
costs. Both also required the two panels to know about each other's shape: a
per-bar clearance test, a keep-out inside the graph's own solver, a rotation
that had to aim somewhere. Side by side needs none of it - each panel is a
rectangle, they cannot collide, and the reader gets two pictures instead of one
picture with two halves. The complexity was buying compactness, which is not
what the figure is for.

What survived from that work and is still in use:

- Chromosome row sorted short to tall, unplaced column at the left end,
  numbered by size rank so chr 1 is still the largest.
- Bar height derived from the space available rather than a fixed constant, so
  the chromosome panel scales with the graph beside it.
- Link-less components packed into the holes the main component leaves, scored
  towards one corner and penalised for growing the bounding box, instead of
  being stacked in a column beside it.
- The chromosome column butts against the graph's INK, not its declared width,
  which carries worst-case padding for rings and labels.

`--graph-triangle` remains, off by default, and is now pointless in the default
composition. Its measurements are kept below because they say something real
about how the solver's constraints interact.

## Where the triangle idea landed

The composition is diagonal: one canvas, graph upper left, chromosome row lower
right, with a clear band of white between them. Row sorted short to tall with
the unplaced column at the left end, right-justified against the canvas edge,
heading beneath the bars. Clearance is checked PER BAR in each bar's own column,
because a staircase's outline is nothing like its bounding box - a single-column
test let a contig run through a chromosome - and the required gap is two ribbon
widths, since a gap narrower than the things it separates does not read as one.

**Confining the graph to its triangle works, but only in the FORCE phase.**
Three versions were tried on the real Flye GFA, measuring the closest distance
between two ribbon centre-lines against the 64 px stroke:

| where the keep-out is applied | min ribbon gap | ink past the diagonal |
|---|---|---|
| projection after the solver | contigs visibly squashed | - |
| force phase AND constraint phase | 50 px (overlapping) | 0.023 |
| **force phase only** | **89 px** | **0.045** |
| none | 91 px | 0.071 |

The middle row is the trap: interleaving the keep-out with the constraint
projection kept shoving beads together after the separation pass had pulled them
apart, so separation never won. Left to the force phase alone, the 120 pure
length-and-separation passes that follow have the last word, and the keep-out
still halves the ink past the diagonal. On by default; `--no-graph-triangle`
turns it off.

Final figure: 1870x1474, against 2456 tall for the plain stack this replaced.

## Packing the isolated components

Components with no links - an organelle ring, an unplaced fragment - used to be
stacked in a COLUMN BESIDE the main component. That put them in the one corner
the figure cannot spare, the bottom right where the chromosome row goes, and it
grew the graph's bounding box sideways for no reason.

They are now placed greedily into the HOLES the main component leaves. Mark the
main component's ink on a coarse grid (cell ~0.9 ribbon widths, clearance
1.15), then for each remaining component take the free position scoring best,
where the score prefers the top left - away from the diagonal the chromosomes
fill - and adds a heavy penalty for any position that would enlarge the overall
bounding box. Filling a hole therefore always beats growing the figure. Falls
back to the old beside-placement when nothing inside the box is free.

On the real Flye GFA this moved edge_11 (mito) and edge_4 (unplaced) out of the
chromosome corner and into the empty upper left, and the canvas fell from
1870x1474 to 1695x1441 with ribbon separation unchanged at 88.8 px.
