#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_figures.py — master runner for the paper figures.

Runs every figure script in ``src/plots/`` (the top-level ``fg*.py`` files), each in
its own subprocess so their global state stays isolated and the saved figures are
bit-identical to running the scripts by hand. The scripts write their output to
``src/figs/`` themselves (via ``src.tools.plot_tools.save``).

Usage (run from anywhere — paths resolve relative to this file):

    python make_figures.py            # run all figure scripts
    python make_figures.py fg5 fg7    # run only scripts whose name contains fg5 or fg7
    python make_figures.py --list     # list the scripts that would run, then exit

A failing script does not stop the others: a summary is printed at the end and the
process exits non-zero if any script failed.

The scripts read their pre-computed arrays from ``src/dat/`` (each has ``reRun = False``);
they do not regenerate the heavy Monte-Carlo data. See DATA.md for the ``src/dat/`` bundle.
"""
import re
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PLOTS = ROOT / "src" / "plots"


def fig_key(path):
    """Sort key: numbered figures in order (fg1, fg1b, fg2, ... fg10, fg10b, fg11),
    with non-numeric appendix figures (e.g. fgb1) last."""
    m = re.match(r"fg(\d+)([a-z]*)", path.stem)
    if m:
        return (0, int(m.group(1)), m.group(2))
    return (1, 0, path.stem)


def discover(filters):
    """Top-level ``fg*.py`` scripts in src/plots/, in figure order, optionally
    filtered to those whose filename contains any of ``filters``."""
    scripts = sorted(PLOTS.glob("fg*.py"), key=fig_key)
    if filters:
        scripts = [s for s in scripts if any(f in s.name for f in filters)]
    return scripts


def main(argv):
    filters = [a for a in argv if not a.startswith("-")]
    scripts = discover(filters)

    if not scripts:
        print("No matching figure scripts found in src/plots/.")
        return 1

    if "--list" in argv:
        for s in scripts:
            print(s.relative_to(ROOT))
        return 0

    failures = []
    for i, script in enumerate(scripts, 1):
        rel = script.relative_to(ROOT)
        print(f"\n[{i}/{len(scripts)}] {rel}", flush=True)
        t0 = time.time()
        result = subprocess.run([sys.executable, str(script)], cwd=ROOT)
        dt = time.time() - t0
        if result.returncode == 0:
            print(f"    ok ({dt:.1f}s)")
        else:
            print(f"    FAILED (exit {result.returncode}, {dt:.1f}s)")
            failures.append(rel)

    print("\n" + "=" * 60)
    print(f"Done: {len(scripts) - len(failures)}/{len(scripts)} succeeded.")
    if failures:
        print("Failed:")
        for f in failures:
            print(f"  - {f}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
