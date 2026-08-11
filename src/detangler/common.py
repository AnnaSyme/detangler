"""Shared primitives: palette, small helpers, logging."""
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

VERSION = "1.0"

# --------------------------------------------------------------------------
# palette (Okabe-Ito derived, colour-vision-deficiency safe)
# --------------------------------------------------------------------------
PALETTE = {
    "chromosome": "#4c72b0",
    "mitochondrion": "#d55e00",
    "plastid": "#009e73",
    "unplaced": "#9aa0a6",
    "unassigned": "#8d8d8d",
    "bar_edge": "#2b2b2b",
    "bg": "#ffffff",
    "text": "#1a1a1a",
    "muted": "#6b6b6b",
    "grid": "#e3e3e3",
}

TANGLE_STYLE = {
    # type                         colour     dash
    "collapsed_repeat": ("#cc79a7", ""),
    "interchromosomal_junction": ("#d55e00", ""),
    "intrachromosomal_repeat": ("#0072b2", "6 3"),
    "inverted_repeat": ("#7b3294", "2 3"),
    "bubble": ("#56b4e9", "1 3"),
    "circular": ("#009e73", ""),
    # graph-first mode: same colours as the assembly graph figure, so the two
    # figures can be read side by side
    "repeat_segment": ("#c2703d", ""),
    "tandem_array": ("#b5487f", ""),
    "at_rich_region": ("#b0a04a", ""),
    "low_coverage_region": ("#a8a8a8", "4 3"),
}

TANGLE_LABEL = {
    "collapsed_repeat": "Collapsed repeat (multi-copy node)",
    "interchromosomal_junction": "Inter-chromosomal join (shared repeat)",
    "intrachromosomal_repeat": "Intra-chromosomal repeat",
    "inverted_repeat": "Inverted / hairpin repeat",
    "bubble": "Bubble (heterozygosity or small variant)",
    "circular": "Circular segment (self-link)",
    "repeat_segment": "Repeat segment (multi-copy)",
    "tandem_array": "Tandem array",
    "at_rich_region": "AT-rich segment",
    "low_coverage_region": "Below single-copy depth",
}


# ==========================================================================
# small helpers
# ==========================================================================
def smart_open(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def human_bp(n: float) -> str:
    if n >= 1e9:
        return f"{n / 1e9:.2f} Gb"
    if n >= 1e6:
        return f"{n / 1e6:.2f} Mb"
    if n >= 1e3:
        return f"{n / 1e3:.1f} kb"
    return f"{int(n)} bp"


def figure_bp(n: float) -> str:
    """Size as it appears ON a figure: whole units, no decimals."""
    if n >= 1e9:
        return f"{n / 1e9:.0f} Gb"
    if n >= 1e6:
        return f"{n / 1e6:.0f} Mb"
    if n >= 1e3:
        return f"{n / 1e3:.0f} kb"
    return f"{int(n)} bp"


def median(values: Seq[float]) -> float:
    vals = sorted(v for v in values if v is not None)
    if not vals:
        return 0.0
    mid = len(vals) // 2
    if len(vals) % 2:
        return float(vals[mid])
    return (vals[mid - 1] + vals[mid]) / 2.0


def nx_stat(lengths: Seq[int], x: float = 50.0) -> int:
    """N50-style statistic. x is a percentage."""
    if not lengths:
        return 0
    total = sum(lengths)
    target = total * x / 100.0
    run = 0
    for L in sorted(lengths, reverse=True):
        run += L
        if run >= target:
            return L
    return min(lengths)


def esc(s: object) -> str:
    return (
        str(s)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


class Log:
    def __init__(self, quiet: bool = False):
        self.quiet = quiet
        self.warnings: List[str] = []

    def info(self, msg: str) -> None:
        if not self.quiet:
            print(f"[detangler] {msg}", file=sys.stderr)

    def warn(self, msg: str) -> None:
        self.warnings.append(msg)
        print(f"[detangler] WARNING: {msg}", file=sys.stderr)


def _natural_key(s: str):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", s)]


def wrap_text(text: str, max_chars: int) -> List[str]:
    """Greedy word wrap. max_chars is derived from an average glyph width."""
    words, lines, cur = text.split(), [], ""
    for w in words:
        if cur and len(cur) + 1 + len(w) > max_chars:
            lines.append(cur)
            cur = w
        else:
            cur = f"{cur} {w}".strip()
    if cur:
        lines.append(cur)
    return lines or [""]


def revcomp(s: str) -> str:
    return s.upper()[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


def parse_size(text: str) -> int:
    """Accepts 36563796, 36.5m, 36.5Mb, 1.2g, 900k."""
    t = str(text).strip().lower().replace(",", "").replace("b", "")
    mult = 1
    if t.endswith("k"):
        mult, t = 1_000, t[:-1]
    elif t.endswith("m"):
        mult, t = 1_000_000, t[:-1]
    elif t.endswith("g"):
        mult, t = 1_000_000_000, t[:-1]
    return int(float(t) * mult)


def _maybe_float(text: str) -> Optional[float]:
    try:
        return float(text)
    except (TypeError, ValueError):
        return None


class _Rand:
    """Tiny deterministic LCG so demo output is byte-stable across Pythons."""

    def __init__(self, seed: int):
        self.s = seed & 0xFFFFFFFF

    def _next(self) -> int:
        self.s = (1103515245 * self.s + 12345) & 0x7FFFFFFF
        return self.s

    def uniform(self, a: float, b: float) -> float:
        return a + (b - a) * (self._next() / 0x7FFFFFFF)

    def randint(self, a: int, b: int) -> int:
        return a + self._next() % (b - a + 1)
