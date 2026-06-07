#!/usr/bin/env python3
"""Run all release experiments."""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
for script in [
    "near_boundary_path.py",
    "regime_selection_all_shifts.py",
    "componentwise_radii.py",
    "precision_scaling.py",
]:
    print(f"running {script}")
    subprocess.check_call([sys.executable, str(ROOT / "experiments" / script)])
