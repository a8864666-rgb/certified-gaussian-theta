#!/usr/bin/env python3
"""Generate representative componentwise radius diagnostics."""
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from certified_gaussian_theta import evaluate

CASES = [("1e-6", "0"), ("1e-4", "0"), ("1", "0"), ("10", "0"), ("100", "1/2"), ("1", "1/3")]


def main() -> None:
    out = ROOT / "results" / "componentwise_radii.csv"
    out.parent.mkdir(exist_ok=True)
    rows = [evaluate(m, a=a, D=1000, guard=50).as_dict() for m, a in CASES]
    with out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
