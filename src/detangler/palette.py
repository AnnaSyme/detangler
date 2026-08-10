"""Per-segment colour assignment shared by both figure panels."""
from __future__ import annotations

import argparse
import glob
import gzip
import json
import math
import os
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field, asdict
from typing import Dict, Iterable, List, Optional, Sequence as Seq, Set, Tuple

try:  # optional, only affects config file format
    import yaml  # type: ignore

    HAVE_YAML = True
except Exception:  # pragma: no cover
    yaml = None  # type: ignore
    HAVE_YAML = False



# ==========================================================================
# GRAPH-FIRST MODE
#
# Everything below works from the assembler's own output only: a GFA and,
# optionally, Flye's assembly_info.txt. No chromosome-level assembly is
# assumed. It produces ranked hypotheses about which segments form which
# chromosome, and is deliberate about the difference between an observation,
# a derived estimate, and a hypothesis.
# ==========================================================================

# Palette carried over from the manual worked example.
CLASS_COLOUR = {
    "backbone": "#4a7ba7",
    "repeat": "#c2703d",
    "at_rich": "#b0a04a",
    "tandem_array": "#b5487f",
    "organelle_candidate": "#3f7d5c",
    "low_coverage": "#a8a8a8",
    "short_single_copy": "#8fa8bd",
    "unclassified": "#cfcfcf",
}
BACKBONE_COLOURS = ["#4a7ba7", "#5f9e6e", "#a87f4a", "#7d6ba7", "#3f8f9e"]

# The figure palette (v9 design). One maximally distinct colour per segment,
# reused in both panels, so a node in the graph can be found on a chromosome by
# colour alone. Class colours are NOT used in the figures any more; CLASS_COLOUR
# survives only for the text report.
SEGMENT_COLOURS = [
    # Vibrant but harmonious: saturated enough to be easy to follow between the
    # two panels, without the clashing you get from a maximally-distinct set.
    "#1982C4",  # azure
    "#FF595E",  # coral
    "#8AC926",  # lime
    "#FFCA3A",  # gold
    "#6A4C93",  # violet
    "#00BBF9",  # sky
    "#F15BB5",  # pink
    "#2A9D8F",  # teal
    "#E36414",  # burnt orange
    "#9B5DE5",  # purple
    "#00C49A",  # mint
    "#C1121F",  # crimson
    "#4361EE",  # indigo
    "#B5179E",  # magenta
    "#57CC99",  # spring green
    "#F9844A",  # tangerine
]
UNPLACED_GREY = "#b8b8b8"

CLASS_LABEL = {
    "backbone": "single-copy backbone",
    "repeat": "repeat",
    "at_rich": "AT-rich",
    "tandem_array": "tandem array",
    "organelle_candidate": "organelle candidate",
    "low_coverage": "low coverage / foreign",
    "short_single_copy": "short, single copy",
    "unclassified": "unclassified",
}


def _rgb(hex_colour: str) -> Tuple[int, int, int]:
    return tuple(int(hex_colour[i : i + 2], 16) for i in (1, 3, 5))  # type: ignore


def colour_distance(a: str, b: str) -> float:
    """
    How far apart two fills look. Weighted to match the eye's sensitivity -
    plain RGB distance rates red against orange as further apart than it looks.
    """
    r1, g1, b1 = _rgb(a)
    r2, g2, b2 = _rgb(b)
    rm = (r1 + r2) / 2.0
    dr, dg, db = r1 - r2, g1 - g2, b1 - b2
    return math.sqrt(
        (2 + rm / 256.0) * dr * dr + 4 * dg * dg + (2 + (255 - rm) / 256.0) * db * db
    )


def assign_segment_colours(
    calls: List["SegmentCall"], links: Optional[List] = None
) -> Dict[str, str]:
    """
    One colour per segment, used by BOTH figures so a segment can be traced from
    the assembly graph to the chromosome it ends up in.

    EVERY segment gets its own maximally distinct colour, not a colour for its
    inferred class: the figure's job is correspondence between the two panels,
    and a shared class colour destroys that. Assignment is by descending length,
    so it is stable between runs.

    When `links` are given, colours are then swapped between segments to push
    LINKED pairs as far apart as possible. Linked segments are exactly the ones
    drawn touching - a repeat is drawn capping the contig it attaches to - so
    they are where a near-miss in colour actually costs you. edge_3 caps edge_6,
    and the two were coming out as a red and a slightly lighter red.
    """
    colours: Dict[str, str] = {}
    ordered = sorted(calls, key=lambda c: (-c.length, c.name))
    for i, c in enumerate(ordered):
        base = SEGMENT_COLOURS[i % len(SEGMENT_COLOURS)]
        cycle = i // len(SEGMENT_COLOURS)
        colours[c.name] = base if cycle == 0 else _lighten(base, 0.26 * cycle)

    if not links:
        return colours

    names = [c.name for c in ordered]
    adjacent = [
        (l.a, l.b) for l in links
        if l.a != l.b and l.a in colours and l.b in colours
    ]
    if not adjacent:
        return colours

    def worst(cmap: Dict[str, str]) -> float:
        return min(colour_distance(cmap[a], cmap[b]) for a, b in adjacent)

    best = worst(colours)
    for _pass in range(4):
        improved = False
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                trial = dict(colours)
                trial[names[i]], trial[names[j]] = trial[names[j]], trial[names[i]]
                score = worst(trial)
                if score > best + 1e-9:
                    colours, best = trial, score
                    improved = True
        if not improved:
            break
    return colours


def _text_on(hex_colour: str) -> str:
    """Black or white, whichever stays legible on the given fill."""
    try:
        r, g, b = (int(hex_colour[i : i + 2], 16) for i in (1, 3, 5))
    except (ValueError, IndexError):
        return "#1a1a1a"
    # sRGB relative luminance, good enough for picking ink
    lum = (0.2126 * r + 0.7152 * g + 0.0722 * b) / 255.0
    return "#1a1a1a" if lum > 0.55 else "#ffffff"


def _segment_number(name: str) -> str:
    """The digits in a segment name: 'edge_11' -> '11'. Falls back to the name."""
    m = re.findall(r"\d+", name)
    return m[-1] if m else name


def _lighten(hex_colour: str, amount: float) -> str:
    amount = min(max(amount, 0.0), 0.8)
    r, g, b = (int(hex_colour[i : i + 2], 16) for i in (1, 3, 5))
    r, g, b = (int(v + (255 - v) * amount) for v in (r, g, b))
    return f"#{r:02x}{g:02x}{b:02x}"
