# Release notes - v2.1-certified-near-boundary

## Summary

This release supports the revised manuscript emphasis on the positive-real near-boundary path

```math
q=e^{-m}\to 1^-.
```

The algorithmic core remains the certified recurrence framework with Poisson-transformed small-`m` evaluation. This release adds clearer reproducibility support for the near-boundary interpretation.

## Changes from v2.0-certified

- Added near-boundary-path experiment script:
  - `experiments/near_boundary_path.py`
- Added generated CSV:
  - `results/near_boundary_path.csv`
- Added componentwise diagnostic CSV:
  - `results/componentwise_radii.csv`
- Added all-shift regime-selection CSV:
  - `results/regime_selection_all_shifts.csv`
- Added lightweight precision-scaling depth diagnostic:
  - `results/precision_scaling.csv`
- Updated README to explain the positive-real nome path `q=e^{-m}->1-`.
- Included manuscript snapshot as `docs/manuscript.pdf`.

## Main reproducibility claim

At `m=1e-6`, the direct planned depth is large, while the transformed planned depth is one retained term. This illustrates the Poisson-dual geometry of the near-boundary path.

## Not a new algorithmic claim

This release should be described as an artifact and reproducibility update, not as a different algorithm from the manuscript. The mathematical core remains:

- recurrence-generated Gaussian terms,
- explicit Gaussian-tail bounds,
- finite-precision roundoff enclosures,
- Poisson-transformed small-`m` evaluation,
- regime-aware representation selection.
