"""Plain data records passed between stages."""
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
# data model
# ==========================================================================
@dataclass
class Evidence:
    """One reason contributing to a classification call."""

    signal: str
    detail: str
    weight: float

    def as_text(self) -> str:
        sign = "+" if self.weight >= 0 else ""
        return f"{self.signal}: {self.detail} ({sign}{self.weight:.1f})"


@dataclass
class SeqRecord:
    name: str
    length: int
    gc: Optional[float] = None
    n_frac: Optional[float] = None
    role: str = "unknown"  # chromosome | mitochondrion | plastid | unplaced
    role_scores: Dict[str, float] = field(default_factory=dict)
    evidence: List[Evidence] = field(default_factory=list)
    circular: bool = False
    depth: Optional[float] = None
    label: Optional[str] = None
    manual: bool = False  # role came from a config override
    note: str = ""  # short reason, shown beside unassigned sequences
    # (start, end, segment_name, colour) blocks making up this molecule, so a
    # segment can be traced from the assembly graph figure to this one
    blocks: List[Tuple[int, int, str, str]] = field(default_factory=list)
    # True when the blocks tile the whole molecule (a chain built from segments),
    # False when they are placements dotted along an already-assembled sequence
    blocks_tile: bool = False
    # {'top'|'bottom': [(segment, colour)]} - repeats the graph attaches to a FREE
    # end of this molecule. They are not part of the chain and are not on the Mb
    # scale, so they are drawn hanging off the end of the bar rather than in it.
    caps: Dict[str, List[Tuple[str, str]]] = field(default_factory=dict)

    @property
    def display(self) -> str:
        return self.label or self.name


@dataclass
class GfaSegment:
    name: str
    length: int
    depth: Optional[float] = None
    gc: Optional[float] = None


@dataclass
class GfaLink:
    a: str
    a_orient: str
    b: str
    b_orient: str
    overlap: str = "*"


@dataclass
class Placement:
    """A graph segment located on an assembled sequence."""

    segment: str
    seqname: str
    start: int
    end: int
    orient: str = "+"
    identity: Optional[float] = None
    source: str = "identity"  # agp | paf | map | identity


@dataclass
class Anchor:
    """Where a tangle touches the ideogram."""

    seqname: str
    start: int
    end: int
    segment: str


@dataclass
class Tangle:
    id: str
    type: str
    segments: List[str]
    anchors: List[Anchor]
    multiplicity: Optional[float] = None
    depth_ratio: Optional[float] = None
    length: Optional[int] = None
    description: str = ""
    evidence: List[str] = field(default_factory=list)

    @property
    def sequences(self) -> List[str]:
        seen: List[str] = []
        for a in self.anchors:
            if a.seqname not in seen:
                seen.append(a.seqname)
        return seen


@dataclass
class CoverageWindow:
    seqname: str
    start: int
    end: int
    depth: float


@dataclass
class ContigInfo:
    """One row of Flye's assembly_info.txt."""

    name: str
    length: int
    cov: Optional[float]
    circular: bool
    repeat: bool
    mult: Optional[float]
    alt_group: str
    path: List[str]  # ordered element ids; "*" marks a dead end
