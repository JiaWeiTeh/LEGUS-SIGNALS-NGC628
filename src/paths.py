# -*- coding: utf-8 -*-
"""Central path constants for the LEGUS-SIGNALS-NGC628 repository.

Everything resolves relative to the repo root, so the code runs from any clone with
no machine-specific absolute paths. Constants are strings ending in a separator, to
match the codebase's `path + "filename"` concatenation style.
"""
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]   # src/paths.py -> repo root


def _d(*parts):
    """Absolute directory path (with trailing separator) under the repo root."""
    return str(ROOT.joinpath(*parts)) + os.sep


DAT          = _d("src", "dat")
FIGS         = _d("src", "figs")
LEGUS        = _d("lib", "LEGUS")
SIGNALS      = _d("lib", "SIGNALS")
SLUG_CLUSTER = _d("lib", "SLUG2", "cluster_slug")
SLUG_LIB     = _d("lib", "SLUG2", "cluster_lib")
