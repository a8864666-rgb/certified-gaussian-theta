# Certified Gaussian Theta Framework

This repository contains the reference implementation and reproducibility scripts for the manuscript:

**A Certified Recurrence Framework for Gaussian Theta Sums: Tail--Roundoff Enclosures and Regime-Aware Transformation**

The code evaluates the two canonical real Gaussian theta sums

\[
S_0(m)=\sum_{n\in\mathbb Z} e^{-m n^2}
\]

and

\[
S_{1/2}(m)=\sum_{n\in\mathbb Z} e^{-m(n+1/2)^2},
\]

corresponding to \(\theta_3(e^{-m})\) and \(\theta_2(e^{-m})\).

## Contents

- `v2-certified/evaluator.py`: certified evaluator
- `v2-certified/experiments.py`: scripts for generating numerical experiments
- `v2-certified/tables.py`: LaTeX table generator
- `v2-certified/outputs/`: CSV files used in the manuscript tables
- `v2-certified/latex_tables/`: LaTeX table files
- `manuscript/main.tex`: LaTeX manuscript source for Overleaf
- `manuscript/certified_gaussian_theta_main_v004.pdf`: compiled manuscript preview

## Requirements

Install dependencies with:

```bash
pip install -r requirements.txt
```

The main dependency is:

```text
mpmath>=1.3.0
```

## Reproducing the tables

From the repository root:

```bash
cd v2-certified
python experiments.py
python tables.py
```

This regenerates the CSV files in `outputs/` and LaTeX table files in `latex_tables/`.

## Main evaluator

The main function is:

```python
eval_theta(m, a, D, Delta, force_regime=None)
```

where:

- `a=0` evaluates \(S_0(m)=\theta_3(e^{-m})\)
- `a=0.5` evaluates \(S_{1/2}(m)=\theta_2(e^{-m})\)
- `D` is the requested decimal precision
- `Delta` is the guard margin

The function returns:

```python
(S_hat, R_total, regime, N, diagnostics)
```

with certified enclosure:

\[
S(m,a)\in \widehat S(m,a)\pm R_{\mathrm{total}}.
\]

## Release

The manuscript version is intended to be archived as:

```text
v2.0-certified
```

## License

MIT License.
