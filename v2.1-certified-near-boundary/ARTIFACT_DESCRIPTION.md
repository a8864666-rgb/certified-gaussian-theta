# Artifact description

## Purpose

The artifact provides a reproducible reference implementation for certified evaluation of real Gaussian theta sums

```math
S(m,a)=\sum_{n\in\mathbb Z} e^{-m(n+a)^2}
```

for positive real `m` and supported shifts `a in {0, 1/2, 1/3, 1/4, 1/6}`.

## Design principle

The code is written to be inspectable rather than maximally optimized. Each run returns an explicit enclosure

```math
S(m,a)\in \widehat S(m,a)\pm R_{total},
```

with

```math
R_{total}=R_{tail}+R_{rec}+R_{sum}+R_{trans}+R_{coef}.
```

For the supported shifts, `R_coef=0` because the transformed cosine weights are binary-exact dyadic rationals.

## Experiments

- `near_boundary_path.py` tests the positive-real path `q=e^{-m}->1-`.
- `regime_selection_all_shifts.py` regenerates regime-selection diagnostics for all supported shifts.
- `componentwise_radii.py` exports representative radius components.
- `precision_scaling.py` exports planned-depth scaling diagnostics.

## Reproducibility command

```bash
python experiments/run_all.py
```

## Known limitations

- Complex theta arguments are outside the present artifact.
- Non-dyadic rational shifts require a nonzero coefficient-roundoff radius.
- Irrational real shifts require nonperiodic cosine-coefficient certification.
- Higher-dimensional nonseparable lattice sums are future work.
