#!/usr/bin/env python3
"""Entry point for detangler.

The implementation lives in the `detangler` package under src/. This wrapper
exists so the documented command `python detangler.py ...` keeps working from a
checkout, without the package needing to be installed.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))

from detangler.cli import main  # noqa: E402

if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, OSError) as exc:
        # Bad input files are expected; a stack trace helps nobody.
        print(f"[detangler] ERROR: {exc}", file=sys.stderr)
        sys.exit(2)
    except KeyboardInterrupt:
        sys.exit(130)
