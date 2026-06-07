# Certified Gaussian Theta Sums - v2.1 certified near-boundary release

This release accompanies the manuscript:

**A Certified Recurrence Framework for Gaussian Theta Sums: Tail-Roundoff Enclosures and Regime-Aware Transformation**  
Thiat-ui Lau, Independent Researcher (Taiwan)

The artifact implements a transparent certified evaluator for real Gaussian theta sums

```math
S(m,a)=\sum_{n\in\mathbb Z} e^{-m(n+a)^2}, \qquad m>0,
```

for

```math
a\in\{0,1/2,1/3,1/4,1/6\}.
```

The release emphasizes the positive-real nome path

```math
q=e^{-m}\to 1^- ,
```

where direct Gaussian summation becomes inefficient, while Poisson transformation maps the problem to a rapidly decaying dual Gaussian series.

## What is included

- `src/certified_gaussian_theta/`: reference Python implementation.
- `experiments/`: reproducibility scripts for the manuscript tables and near-boundary diagnostics.
- `results/`: generated CSV files.
- `docs/manuscript.pdf`: manuscript snapshot associated with this release.
- `RELEASE_NOTES.md`: release summary.
- `ARTIFACT_DESCRIPTION.md`: artifact scope, limitations, and reproducibility notes.
- `CITATION.cff`: citation metadata template.
- `LICENSE_TODO.md`: license selection note.

## Install

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

The code uses `mpmath` for arbitrary-precision arithmetic.

## Run all experiments

```bash
python experiments/run_all.py
```

This regenerates:

```text
results/near_boundary_path.csv
results/regime_selection_all_shifts.csv
results/componentwise_radii.csv
results/precision_scaling.csv
```

## Single evaluation example

```python
from certified_gaussian_theta import evaluate

result = evaluate("1e-6", a="0", D=1000, guard=50)
print(result)
```

Expected qualitative behavior for the near-boundary case `m=1e-6`:

```text
planned_direct_depth      49171
planned_transformed_depth 1
selected_regime           transformed
R_coef                    0
```

## Radius components

The evaluator reports the componentwise error budget:

```text
R_tail   truncation tail radius
R_rec    recurrence and initialization roundoff radius
R_sum    balanced summation roundoff radius
R_trans  transformation-prefactor radius
R_coef   coefficient-representation radius; zero for supported binary-exact shifts
R_total  total certified radius
```

For the supported shifts, the Poisson cosine weights belong to `{1, -1, 0, 1/2, -1/2}` and are exactly representable in binary floating-point arithmetic, so `R_coef = 0`.

## Scope

This artifact is intentionally specialized. It does not replace general complex theta-function algorithms or ball-arithmetic libraries. It provides an inspectable, white-box certified framework for positive-real Gaussian theta sums with exact-weight rational shifts.

## Important limitation

For general rational shifts with non-dyadic cosine weights, or irrational real shifts, an additional nonzero coefficient-roundoff component `R_coef` is required. That extension is planned as future work.
