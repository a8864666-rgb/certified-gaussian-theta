"""Convert CSV experiment outputs into LaTeX table drafts."""
from __future__ import annotations

import csv
import os
from typing import List, Dict

OUTPUT_DIR = "outputs"
TABLE_DIR = "latex_tables"


def ensure_table_dir() -> None:
    os.makedirs(TABLE_DIR, exist_ok=True)


def read_csv(filename: str) -> List[Dict[str, str]]:
    path = os.path.join(OUTPUT_DIR, filename)
    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing CSV file: {path}")
    with open(path, "r", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def fmt_value(x: str) -> str:
    if x is None:
        return ""
    x = str(x).strip()
    if x.lower() in {"inf", "+inf", "infinity"}:
        return r"\mathrm{n/a}"
    replacements = {
        "1e-6": r"10^{-6}",
        "1e-4": r"10^{-4}",
        "1e-2": r"10^{-2}",
        "1e-1": r"10^{-1}",
        "1/2": r"\frac12",
        "pi": r"\pi",
    }
    return replacements.get(x, x)


def fmt_float_digits(x: str, digits: int = 1) -> str:
    x = str(x).strip()
    if x.lower() in {"inf", "+inf", "infinity"}:
        return r"\mathrm{n/a}"
    try:
        return f"{float(x):.{digits}f}"
    except Exception:
        return fmt_value(x)


def latex_table(caption: str, label: str, headers: List[str], rows: List[List[str]], align: str) -> str:
    lines = [r"\begin{table}[t]", r"\centering", r"\small", rf"\begin{{tabular}}{{{align}}}", r"\hline"]
    lines.append(" & ".join(headers) + r" \\")
    lines.append(r"\hline")
    for row in rows:
        lines.append(" & ".join(row) + r" \\")
    lines.extend([r"\hline", r"\end{tabular}", rf"\caption{{{caption}}}", rf"\label{{{label}}}", r"\end{table}", ""])
    return "\n".join(lines)


def write_latex(filename: str, content: str) -> None:
    ensure_table_dir()
    path = os.path.join(TABLE_DIR, filename)
    with open(path, "w", encoding="utf-8") as f:
        f.write(content)
    print(f"Wrote {path}")


def make_regime_selection_table(csv_filename: str, tex_filename: str, caption: str, label: str) -> None:
    data = read_csv(csv_filename)
    headers = [r"$m$", r"$N_{\mathrm{dir}}^{(0)}$", r"$N_{\mathrm{tr}}^{(0)}$", "regime", r"$N$", r"$-\log_{10}R_{\mathrm{total}}$"]
    rows = []
    for row in data:
        rows.append([rf"${fmt_value(row['m_label'])}$", fmt_value(row["N_dir"]), fmt_value(row["N_tr"]), fmt_value(row["selected_regime"]), fmt_value(row["final_N"]), fmt_float_digits(row["minus_log10_R_total"], 1)])
    content = latex_table(caption, label, headers, rows, "c c c c c c")
    write_latex(tex_filename, content)


def make_componentwise_table(csv_filename: str, tex_filename: str, caption: str, label: str) -> None:
    data = read_csv(csv_filename)
    headers = [r"$m$", "regime", r"$N$", r"$-\log_{10}R_{\mathrm{tail}}$", r"$-\log_{10}R_{\mathrm{rec}}$", r"$-\log_{10}R_{\mathrm{sum}}$", r"$-\log_{10}R_{\mathrm{trans}}$", r"$-\log_{10}R_{\mathrm{total}}$"]
    rows = []
    for row in data:
        rows.append([rf"${fmt_value(row['m_label'])}$", fmt_value(row["selected_regime"]), fmt_value(row["final_N"]), fmt_float_digits(row["minus_log10_R_tail"], 1), fmt_float_digits(row["minus_log10_R_rec"], 1), fmt_float_digits(row["minus_log10_R_sum"], 1), fmt_float_digits(row["minus_log10_R_trans"], 1), fmt_float_digits(row["minus_log10_R_total"], 1)])
    content = latex_table(caption, label, headers, rows, "c c c c c c c c")
    write_latex(tex_filename, content)


def make_crossover_table() -> None:
    data = read_csv("tableX_crossover_plateau.csv")
    headers = [r"$m$", r"$N_{\mathrm{dir}}^{(0)}$", r"$N_{\mathrm{tr}}^{(0)}$", "selected regime", r"$N$"]
    rows = []
    for row in data:
        rows.append([rf"${fmt_value(row['m_label'])}$", fmt_value(row["N_dir"]), fmt_value(row["N_tr"]), fmt_value(row["selected_regime"]), fmt_value(row["final_N"])])
    content = latex_table(r"Integer crossover plateau near the leading-order scale $m_{\mathrm{crit}}=\pi$ for $D=1000$ and $\Delta=50$.", "tab:crossover-plateau", headers, rows, "c c c c c")
    write_latex("tableX_crossover_plateau.tex", content)


def make_precision_scaling_table(csv_filename: str, tex_filename: str, caption: str, label: str) -> None:
    data = read_csv(csv_filename)
    headers = [r"$m$", r"$D$", "regime", r"$N$", r"$N/\sqrt{D}$", r"$-\log_{10}R_{\mathrm{total}}$"]
    rows = []
    for row in data:
        rows.append([rf"${fmt_value(row['m_label'])}$", fmt_value(row["D"]), fmt_value(row["selected_regime"]), fmt_value(row["final_N"]), fmt_float_digits(row["N_over_sqrt_D"], 3), fmt_float_digits(row["minus_log10_R_total"], 1)])
    content = latex_table(caption, label, headers, rows, "c c c c c c")
    write_latex(tex_filename, content)


def make_rec_vs_exp_table() -> None:
    data = read_csv("table7_rec_vs_exp.csv")
    headers = [r"$m$", r"$a$", r"$D$", "expected", "selected", r"$N$", r"$-\log_{10}\Delta_{\mathrm{rec-exp}}$", r"$-\log_{10}R_{\mathrm{total}}$", "status"]
    rows = []
    for row in data:
        rows.append([rf"${fmt_value(row['m_label'])}$", rf"${fmt_value(row['a_label'])}$", fmt_value(row["D"]), fmt_value(row["expected_regime"]), fmt_value(row["selected_regime"]), fmt_value(row["N"]), fmt_float_digits(row["minus_log10_Delta_rec_exp"], 1), fmt_float_digits(row["minus_log10_R_total"], 1), fmt_value(row["status"])])
    content = latex_table(r"Consistency check comparing recurrence-generated values against independently evaluated direct-exponential sums.", "tab:rec-vs-exp", headers, rows, "c c c c c c c c c")
    write_latex("table7_rec_vs_exp.tex", content)


def make_guard_stability_table() -> None:
    data = read_csv("table8_guard_stability.csv")
    headers = [r"$m$", r"$a$", r"$D$", r"$\Delta$", r"$\Delta_{\mathrm{high}}$", "base", "high", r"$-\log_{10}\Delta S$", r"$-\log_{10}(R_\Delta+R_{\Delta_{\mathrm{high}}})$", "status"]
    rows = []
    for row in data:
        rows.append([rf"${fmt_value(row['m_label'])}$", rf"${fmt_value(row['a_label'])}$", fmt_value(row["D"]), fmt_value(row["Delta"]), fmt_value(row["Delta_high"]), fmt_value(row["base_regime"]), fmt_value(row["high_regime"]), fmt_float_digits(row["minus_log10_Delta_S"], 1), fmt_float_digits(row["minus_log10_bound"], 1), fmt_value(row["status"])])
    content = latex_table(r"Guard-precision stability test comparing computations with $\Delta=50$ and $\Delta_{\mathrm{high}}=150$.", "tab:guard-stability", headers, rows, "c c c c c c c c c c")
    write_latex("table8_guard_stability.tex", content)


def make_reference_comparison_table() -> None:
    data = read_csv("table9_reference_comparison.csv")
    headers = [r"$m$", r"$a$", r"$D$", "regime", r"$N$", r"$-\log_{10}\Delta_{\mathrm{ref}}$", r"$-\log_{10}R_{\mathrm{total}}$", "status"]
    rows = []
    for row in data:
        rows.append([rf"${fmt_value(row['m_label'])}$", rf"${fmt_value(row['a_label'])}$", fmt_value(row["D"]), fmt_value(row["selected_regime"]), fmt_value(row["N"]), fmt_float_digits(row["minus_log10_Delta_ref"], 1), fmt_float_digits(row["minus_log10_R_total"], 1), fmt_value(row["status"])])
    content = latex_table(r"Numerical comparison with high-precision mpmath reference values. The comparison is used as an external sanity check, not as the source of certification.", "tab:reference-comparison", headers, rows, "c c c c c c c c")
    write_latex("table9_reference_comparison.tex", content)


def main() -> None:
    make_regime_selection_table("table1_regime_selection_a0.csv", "table1_regime_selection_a0.tex", r"Regime selection for $S_0(m)=\theta_3(e^{-m})$ with $D=1000$ and $\Delta=50$.", "tab:regime-selection-a0")
    make_regime_selection_table("table2_regime_selection_a_half.csv", "table2_regime_selection_a_half.tex", r"Regime selection for $S_{1/2}(m)=\theta_2(e^{-m})$ with $D=1000$ and $\Delta=50$.", "tab:regime-selection-a-half")
    make_componentwise_table("table1_regime_selection_a0.csv", "table3_components_a0.tex", r"Componentwise certified radii for $S_0(m)=\theta_3(e^{-m})$ with $D=1000$ and $\Delta=50$.", "tab:components-a0")
    make_componentwise_table("table2_regime_selection_a_half.csv", "table4_components_a_half.tex", r"Componentwise certified radii for $S_{1/2}(m)=\theta_2(e^{-m})$ with $D=1000$ and $\Delta=50$.", "tab:components-a-half")
    make_crossover_table()
    make_precision_scaling_table("table5_precision_scaling_a0.csv", "table5_precision_scaling_a0.tex", r"Precision scaling for $S_0(m)$; the column $N/\sqrt{D}$ displays the expected sublinear scaling except in minimum-depth saturation.", "tab:precision-a0")
    make_precision_scaling_table("table6_precision_scaling_a_half.csv", "table6_precision_scaling_a_half.tex", r"Precision scaling for $S_{1/2}(m)$; the column $N/\sqrt{D}$ displays the expected sublinear scaling except in minimum-depth saturation.", "tab:precision-a-half")
    make_rec_vs_exp_table()
    make_guard_stability_table()
    make_reference_comparison_table()


if __name__ == "__main__":
    main()
