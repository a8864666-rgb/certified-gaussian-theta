#!/usr/bin/env python3
"""Generate manuscript-style regime-selection table for all supported shifts."""
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from certified_gaussian_theta import evaluate

M_VALUES = ["1e-6", "1e-4", "1e-2", "1e-1", "1", "10", "100"]
SHIFTS = ["0", "1/2", "1/3", "1/4", "1/6"]


def main() -> None:
    out = ROOT / "results" / "regime_selection_all_shifts.csv"
    out.parent.mkdir(exist_ok=True)
    rows = []
    for a in SHIFTS:
        for m in M_VALUES:
            rows.append(evaluate(m, a=a, D=1000, guard=50).as_dict())
    with out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
