"""Geometry and helpers shared by the renderers."""
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
from .common import (
    wrap_text,
)
from .records import (
    Anchor,
    Tangle,
)
from .model import (
    Model,
)



# ==========================================================================
# layout + SVG rendering
# ==========================================================================
MARGIN_L = 86
MARGIN_R = 40
MARGIN_T = 96
BAR_W = 48
COV_W = 20
GAP = 96
MAX_BAR_H = 620
MIN_ORG_H = 52
LEGEND_H = 150
PANEL_W = 300


class Layout:
    def __init__(self, model: Model, show_coverage: bool, header_h: float = MARGIN_T):
        self.header_h = header_h
        self.order = model.drawable()
        self.x: Dict[str, float] = {}
        self.top: Dict[str, float] = {}
        self.height: Dict[str, float] = {}
        self.length: Dict[str, int] = {}
        self.not_to_scale: Set[str] = set()
        self.show_coverage = show_coverage

        chrom_lengths = [s.length for s in self.order if s.role == "chromosome"] or [
            s.length for s in self.order
        ] or [1]
        self.max_len = max(chrom_lengths)
        self.scale = MAX_BAR_H / float(self.max_len)

        step = BAR_W + GAP + (COV_W + 6 if show_coverage else 0)
        self.step = step
        # Chain labels can be long ("chain 2: edge_7 + edge_2"), so wrap them to
        # the column width rather than letting neighbours collide.
        self.label_lines: Dict[str, List[str]] = {
            s.name: wrap_text(s.display, max(int(step / 6.2), 9))[:3] for s in self.order
        }
        self.max_label_lines = max((len(v) for v in self.label_lines.values()), default=1)
        for i, s in enumerate(self.order):
            self.x[s.name] = MARGIN_L + i * step
            h = s.length * self.scale
            if h < MIN_ORG_H:
                h = MIN_ORG_H
                self.not_to_scale.add(s.name)
            self.height[s.name] = h
            self.top[s.name] = header_h
            self.length[s.name] = s.length

        n = max(len(self.order), 1)
        self.panel_x = MARGIN_L + n * step + 16
        self.panel = bool(model.unassigned())
        panel_w = PANEL_W if self.panel else 0
        self.width = max(MARGIN_L + n * step - GAP + MARGIN_R + panel_w, 940)
        self.height_total = header_h + MAX_BAR_H + LEGEND_H  # refined by render_svg

    @property
    def text_cols(self) -> int:
        """Characters that fit on one line at ~12px in the drawing area."""
        return max(int((self.width - MARGIN_L - MARGIN_R) / 6.4), 40)

    @property
    def panel_cols(self) -> int:
        return max(int((PANEL_W - 20) / 5.6), 24)

    def y(self, seqname: str, pos: float) -> float:
        L = max(self.length.get(seqname, 1), 1)
        frac = min(max(pos / float(L), 0.0), 1.0)
        return self.top[seqname] + frac * self.height[seqname]

    def cx(self, seqname: str) -> float:
        return self.x[seqname] + BAR_W / 2.0


# The smallest a contig may be drawn, in either panel. Bandage uses the same
# idea (a minimumNodeLength floor under a length-proportional scale); the point
# of pinning it here is that the graph panel and the chromosome panel must agree,
# or a 15 kb telomeric repeat looks like two different sizes in one figure.
# A full ribbon width. A round-capped stroke is never shorter than its own
# width however short the path is, so anything below this cannot be drawn
# honestly in the graph panel - and the chromosome panel then uses the same
# floor so the two agree.
MIN_DRAWN_PX = BAR_W


def drawn_length_px(length: int, px_per_bp: float) -> float:
    """A contig's drawn length, identical wherever it is drawn."""
    return max(float(length) * px_per_bp, MIN_DRAWN_PX)


def _bar_path(x: float, y: float, w: float, h: float, rt: float, rb: float) -> str:
    """
    Rectangle with independently rounded top and bottom ends.

    Used instead of a clipPath because clipping is unevenly supported outside
    browsers - ImageMagick, for one, silently drops the clipped group, which
    would lose the segment blocks in any PNG conversion.
    """
    rt = max(0.0, min(rt, w / 2, h / 2))
    rb = max(0.0, min(rb, w / 2, h / 2))
    p = [f"M {x:.1f} {y + rt:.1f}"]
    if rt:
        p.append(f"A {rt:.1f} {rt:.1f} 0 0 1 {x + rt:.1f} {y:.1f}")
    else:
        p.append(f"L {x:.1f} {y:.1f}")
    p.append(f"L {x + w - rt:.1f} {y:.1f}")
    if rt:
        p.append(f"A {rt:.1f} {rt:.1f} 0 0 1 {x + w:.1f} {y + rt:.1f}")
    p.append(f"L {x + w:.1f} {y + h - rb:.1f}")
    if rb:
        p.append(f"A {rb:.1f} {rb:.1f} 0 0 1 {x + w - rb:.1f} {y + h:.1f}")
    else:
        p.append(f"L {x + w:.1f} {y + h:.1f}")
    p.append(f"L {x + rb:.1f} {y + h:.1f}")
    if rb:
        p.append(f"A {rb:.1f} {rb:.1f} 0 0 1 {x:.1f} {y + h - rb:.1f}")
    p.append("Z")
    return " ".join(p)


def _annotation_colour(kind: str) -> str:
    palette = ["#b07aa1", "#59a14f", "#edc948", "#e15759", "#76b7b2", "#ff9da7", "#9c755f"]
    return palette[abs(hash(kind)) % len(palette)]


def _arc_path(x1: float, y1: float, x2: float, y2: float) -> str:
    """Cubic Bezier between two ideogram points; bows outward if same column."""
    if abs(x1 - x2) < 1e-6:
        bow = 46 + min(abs(y2 - y1) * 0.35, 90)
        return f"M {x1:.1f} {y1:.1f} C {x1 + bow:.1f} {y1:.1f}, {x2 + bow:.1f} {y2:.1f}, {x2:.1f} {y2:.1f}"
    dx = (x2 - x1) * 0.45
    return (
        f"M {x1:.1f} {y1:.1f} C {x1 + dx:.1f} {y1:.1f}, "
        f"{x2 - dx:.1f} {y2:.1f}, {x2:.1f} {y2:.1f}"
    )


def representative_anchors(t: Tangle, drawn: Set[str]) -> List[Anchor]:
    """One anchor per sequence (the first), restricted to sequences we draw."""
    out: List[Anchor] = []
    seen: Set[str] = set()
    for a in t.anchors:
        if a.seqname in drawn and a.seqname not in seen:
            seen.add(a.seqname)
            out.append(a)
    if len(out) <= 1:
        # intra-sequence repeat: keep up to 4 distinct positions on one sequence
        pos: List[Anchor] = []
        for a in t.anchors:
            if a.seqname in drawn and not any(abs(a.start - p.start) < 1000 for p in pos):
                pos.append(a)
        if len(pos) > 1:
            return pos[:4]
    return out[:4]

# Type scale for the figures (v9): title, panel heading, label, sub-label.
FS_TITLE, FS_HEADING, FS_LABEL, FS_SUB = 40, 32, 27, 18


def _svg_height(svg: str) -> float:
    m = re.search(r'<svg[^>]*\bheight="([0-9.]+)"', svg)
    return float(m.group(1)) if m else 800.0


def _svg_width(svg: str) -> float:
    m = re.search(r'<svg[^>]*\bwidth="([0-9.]+)"', svg)
    return float(m.group(1)) if m else 800.0


def image_size(path: str) -> Optional[Tuple[float, float]]:
    """Pixel dimensions of a PNG, JPEG or SVG, without needing an image library."""
    import struct

    with open(path, "rb") as fh:
        head = fh.read(4096)
    if head[:8] == b"\x89PNG\r\n\x1a\n":
        w, h = struct.unpack(">II", head[16:24])
        return float(w), float(h)
    if head[:2] == b"\xff\xd8":  # JPEG: walk the segment markers for SOFn
        with open(path, "rb") as fh:
            fh.read(2)
            while True:
                b = fh.read(1)
                while b and b != b"\xff":
                    b = fh.read(1)
                marker = fh.read(1)
                while marker == b"\xff":
                    marker = fh.read(1)
                if not marker:
                    return None
                if marker[0] in range(0xC0, 0xCF) and marker[0] not in (0xC4, 0xC8, 0xCC):
                    fh.read(3)
                    h, w = struct.unpack(">HH", fh.read(4))
                    return float(w), float(h)
                size = struct.unpack(">H", fh.read(2))[0]
                fh.read(size - 2)
    text = head.decode("utf-8", "replace")
    if "<svg" in text:
        w = re.search(r'<svg[^>]*\bwidth="([0-9.]+)', text)
        h = re.search(r'<svg[^>]*\bheight="([0-9.]+)', text)
        if w and h:
            return float(w.group(1)), float(h.group(1))
        vb = re.search(r'viewBox="[\s0-9.\-]*?([0-9.]+)\s+([0-9.]+)"', text)
        if vb:
            return float(vb.group(1)), float(vb.group(2))
    return None


def embed_image(path: str, x: float, y: float, w: float, h: float) -> str:
    """
    Place an external image as a panel, inlined as base64 so the figure stays a
    single portable file. An SVG export keeps its vectors; a PNG or JPEG is
    embedded as-is.
    """
    import base64
    import mimetypes

    if path.lower().endswith(".svg"):
        with open(path) as fh:
            svg = fh.read()
        body = svg.split(">", 1)[1].rsplit("</svg>", 1)[0]
        body = _BG_RECT_RE.sub("", body, count=1)
        src = image_size(path) or (w, h)
        scale = min(w / src[0], h / src[1])
        return (
            f'<g transform="translate({x:.1f},{y:.1f}) scale({scale:.4f})">{body}</g>'
        )
    mime = mimetypes.guess_type(path)[0] or "image/png"
    with open(path, "rb") as fh:
        data = base64.b64encode(fh.read()).decode("ascii")
    return (
        f'<image x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" '
        f'preserveAspectRatio="xMidYMin meet" xlink:href="data:{mime};base64,{data}"/>'
    )


_BG_RECT_RE = re.compile(r'<rect width="100%" height="100%"[^/]*/>')


def _place_svg(svg: str, x: float, y: float, rotate: float = 0) -> str:
    """
    Re-position a complete SVG document as a panel of a larger figure.

    A translated <g> rather than a nested <svg>: nested-SVG positioning is
    another thing renderers disagree about - ImageMagick ignores the x/y and
    stacks both panels at the origin. The panel's own background rect is dropped
    so it cannot paint over its neighbour or over the leader lines.
    """
    body = svg.split(">", 1)[1].rsplit("</svg>", 1)[0]
    body = _BG_RECT_RE.sub("", body, count=1)
    if rotate == -90:
        # a quarter turn anticlockwise: the panel's own width becomes the height
        # of its footprint, so shift down by that width to land back in view
        w = _svg_width(svg)
        return f'<g transform="translate({x:.1f},{y + w:.1f}) rotate(-90)">{body}</g>'
    return f'<g transform="translate({x:.1f},{y:.1f})">{body}</g>'
