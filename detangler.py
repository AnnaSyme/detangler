#!/usr/bin/env python3
"""Entry point for detangler.

The implementation currently lives in GenomeViz_ideogram-tool_v1.py; this
wrapper exists so the documented command `python detangler.py ...` works.
"""
import pathlib
import runpy

runpy.run_path(
    str(pathlib.Path(__file__).with_name("GenomeViz_ideogram-tool_v1.py")),
    run_name="__main__",
)
