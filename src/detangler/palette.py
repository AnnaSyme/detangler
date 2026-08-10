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
    "#e6194b",  # red
    "#4363d8",  # blue
    "#3cb44b",  # green
    "#f58231",  # orange
    "#911eb4",  # purple
    "#42d4f4",  # cyan
    "#f032e6",  # magenta
    "#9a6324",  # brown
    "#7ba428",  # olive
    "#ffe119",  # yellow
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


def assign_segment_colours(calls: List["SegmentCall"]) -> Dict[str, str]:
    """
    One colour per segment, used by BOTH figures so a segment can be traced from
    the assembly graph to the chromosome it ends up in.

    EVERY segment gets its own maximally distinct colour, not a colour for its
    inferred class: the figure's job is correspondence between the two panels,
    and a shared class colour destroys that. Assignment is by descending length,
    so it is stable between runs. Once the palette is exhausted it is cycled and
    lightened, which keeps later segments distinguishable from earlier ones.
    """
    colours: Dict[str, str] = {}
    ordered = sorted(calls, key=lambda c: (-c.length, c.name))
    for i, c in enumerate(ordered):
        base = SEGMENT_COLOURS[i % len(SEGMENT_COLOURS)]
        cycle = i // len(SEGMENT_COLOURS)
        colours[c.name] = base if cycle == 0 else _lighten(base, 0.26 * cycle)
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
