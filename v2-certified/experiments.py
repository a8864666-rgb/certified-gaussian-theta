"""
Experiment runner for the v2-certified Gaussian theta framework.
Generates CSV tables used by the manuscript.
"""

from __future__ import annotations

import csv
import os
import platform
import sys
if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)
from datetime import datetime
from typing import Dict, Any, List

import mpmath as mp

from evaluator import eval_theta

OUTPUT_DIR = "outputs"


def ensure_output_dir() -> None:
    os.makedirs(OUTPUT_DIR, exist_ok=True)


def environment_metadata() -> Dict[str, str]:
    return {
        "timestamp_utc": datetime.utcnow().isoformat() + "Z",
        "python_version": sys.version.replace("\n", " "),
        "platform": platform.platform(),
        "processor": platform.processor(),
        "mpmath_version": getattr(mp, "__version__", "unknown"),
    }


def write_csv(filename: str, rows: List[Dict[str, Any]]) -> None:
    ensure_output_dir()
    path = os.path.join(OUTPUT_DIR, filename)
    if not rows:
        raise ValueError(f"No rows to write for {filename}")
    fieldnames = []
    for row in rows:
        for key in row.keys():
            if key not in fieldnames:
                fieldnames.append(key)
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    print(f"Wrote {path}")


def _safe_minus_log10_local(x) -> mp.mpf:
    x = mp.mpf(x)
    if x == 0:
        return mp.inf
    return -mp.log10(abs(x))


def compact_row_from_eval(m_label, a_label, D, Delta, S_hat, R_total, regime, N, diagnostics):
    row = {
        "m_label": m_label,
        "a_label": a_label,
        "D": D,
        "Delta": Delta,
        "N_dir": diagnostics["N_dir"],
        "N_tr": diagnostics["N_tr"],
        "selected_regime": regime,
        "final_N": N,
        "minus_log10_R_total": diagnostics["minus_log10_R_total"],
        "S_hat_prefix": mp.nstr(S_hat, 30),
        "minus_log10_R_tail": diagnostics["minus_log10_R_tail"],
        "minus_log10_R_rec": diagnostics["minus_log10_R_rec"],
        "minus_log10_R_sum": diagnostics["minus_log10_R_sum"],
        "minus_log10_R_trans": diagnostics["minus_log10_R_trans"],
        "p_bits": diagnostics["p_bits"],
        "mp_dps": diagnostics["mp_dps"],
        "lambda_eff": diagnostics["lambda_eff"],
        "c_rec": diagnostics["c_rec"],
        "c_init": diagnostics["c_init"],
        "c_alpha": diagnostics["c_alpha"],
    }
    return row


def run_regime_selection_tables(D: int = 1000, Delta: int = 50) -> None:
    m_values = [("1e-6", mp.mpf("1e-6")), ("1e-4", mp.mpf("1e-4")), ("1e-2", mp.mpf("1e-2")),
                ("1e-1", mp.mpf("1e-1")), ("1", mp.mpf("1")), ("10", mp.mpf("10")), ("100", mp.mpf("100"))]
    metadata = environment_metadata()
    for a_label, a_value, filename in [("0", mp.mpf("0"), "table1_regime_selection_a0.csv"),
                                       ("1/2", mp.mpf("0.5"), "table2_regime_selection_a_half.csv")]:
        rows = []
        for m_label, m_value in m_values:
            S_hat, R_total, regime, N, diagnostics = eval_theta(m=m_value, a=a_value, D=D, Delta=Delta)
            row = compact_row_from_eval(m_label, a_label, D, Delta, S_hat, R_total, regime, N, diagnostics)
            row.update(metadata)
            rows.append(row)
            print(f"[regime] a={a_label}, m={m_label}, regime={regime}, N={N}, -log10 R_total={diagnostics['minus_log10_R_total']}")
        write_csv(filename, rows)


def run_crossover_plateau_table(D: int = 1000, Delta: int = 50) -> None:
    m_values = [("2.0", mp.mpf("2.0")), ("2.5", mp.mpf("2.5")), ("2.9", mp.mpf("2.9")),
                ("pi", mp.pi), ("3.2", mp.mpf("3.2")), ("4.0", mp.mpf("4.0"))]
    rows = []
    metadata = environment_metadata()
    a_value = mp.mpf("0")
    for m_label, m_value in m_values:
        S_hat, R_total, regime, N, diagnostics = eval_theta(m=m_value, a=a_value, D=D, Delta=Delta)
        row = {
            "m_label": m_label,
            "m_value": mp.nstr(m_value, 50),
            "D": D,
            "Delta": Delta,
            "N_dir": diagnostics["N_dir"],
            "N_tr": diagnostics["N_tr"],
            "selected_regime": regime,
            "final_N": N,
            "minus_log10_R_total": diagnostics["minus_log10_R_total"],
        }
        row.update(metadata)
        rows.append(row)
        print(f"[crossover] m={m_label}, N_dir={diagnostics['N_dir']}, N_tr={diagnostics['N_tr']}, selected={regime}")
    write_csv("tableX_crossover_plateau.csv", rows)


def run_precision_scaling_tables(Delta: int = 50) -> None:
    m_values = [("1e-4", mp.mpf("1e-4")), ("1", mp.mpf("1")), ("100", mp.mpf("100"))]
    D_values = [1000, 5000, 10000]
    metadata = environment_metadata()
    for a_label, a_value, filename in [("0", mp.mpf("0"), "table5_precision_scaling_a0.csv"),
                                       ("1/2", mp.mpf("0.5"), "table6_precision_scaling_a_half.csv")]:
        rows = []
        for m_label, m_value in m_values:
            for D in D_values:
                S_hat, R_total, regime, N, diagnostics = eval_theta(m=m_value, a=a_value, D=D, Delta=Delta)
                row = compact_row_from_eval(m_label, a_label, D, Delta, S_hat, R_total, regime, N, diagnostics)
                row["N_over_sqrt_D"] = mp.nstr(mp.mpf(N) / mp.sqrt(D), 20)
                row.update(metadata)
                rows.append(row)
                print(f"[precision] a={a_label}, m={m_label}, D={D}, regime={regime}, N={N}, N/sqrtD={row['N_over_sqrt_D']}")
        write_csv(filename, rows)


def direct_exponential_value(m, a, N, regime):
    m = mp.mpf(m)
    a = mp.mpf(a)
    if regime == "direct":
        if a == 0:
            terms = [mp.exp(-m * j * j) for j in range(1, N + 1)]
            return 1 + 2 * mp.fsum(terms)
        if a == mp.mpf("0.5"):
            terms = [mp.exp(-m * (j + mp.mpf("0.5")) ** 2) for j in range(0, N + 1)]
            return 2 * mp.fsum(terms)
    if regime == "transformed":
        mu = mp.pi ** 2 / m
        alpha = mp.sqrt(mp.pi / m)
        if a == 0:
            terms = [mp.exp(-mu * j * j) for j in range(1, N + 1)]
            return alpha * (1 + 2 * mp.fsum(terms))
        if a == mp.mpf("0.5"):
            terms = [((-1) ** j) * mp.exp(-mu * j * j) for j in range(1, N + 1)]
            return alpha * (1 + 2 * mp.fsum(terms))
    raise ValueError("Unsupported a or regime")


def reference_mpmath_value(m, a, ref_dps):
    old_dps = mp.mp.dps
    mp.mp.dps = ref_dps
    m = mp.mpf(m)
    a = mp.mpf(a)
    q = mp.exp(-m)
    if a == 0:
        val = mp.jtheta(3, 0, q)
    elif a == mp.mpf("0.5"):
        val = mp.jtheta(2, 0, q)
    else:
        raise ValueError("Only a=0 and a=0.5 are supported")
    mp.mp.dps = old_dps
    return val


def run_rec_vs_exp_table(D: int = 1000, Delta: int = 50) -> None:
    cases = [("1e-4", mp.mpf("1e-4"), "0", mp.mpf("0"), "transformed"),
             ("1", mp.mpf("1"), "0", mp.mpf("0"), "transition"),
             ("100", mp.mpf("100"), "1/2", mp.mpf("0.5"), "direct")]
    rows = []
    metadata = environment_metadata()
    eps_exp = mp.mpf(10) ** (-(D + Delta + 120))
    for m_label, m_value, a_label, a_value, expected_regime in cases:
        S_rec, R_total, regime, N, diagnostics = eval_theta(m=m_value, a=a_value, D=D, Delta=Delta)
        old_dps = mp.mp.dps
        mp.mp.dps = int(D + Delta + 250)
        S_exp = direct_exponential_value(m=m_value, a=a_value, N=N, regime=regime)
        diff = abs(S_rec - S_exp)
        mp.mp.dps = old_dps
        status = "pass" if diff <= R_total + eps_exp else "fail"
        row = {
            "m_label": m_label, "a_label": a_label, "D": D, "Delta": Delta,
            "expected_regime": expected_regime, "selected_regime": regime, "N": N,
            "Delta_rec_exp": mp.nstr(diff, 30), "eps_exp": mp.nstr(eps_exp, 30),
            "minus_log10_Delta_rec_exp": mp.nstr(_safe_minus_log10_local(diff), 30),
            "R_total": diagnostics["R_total"], "minus_log10_R_total": diagnostics["minus_log10_R_total"],
            "status": status,
        }
        row.update(metadata)
        rows.append(row)
        print(f"[rec-vs-exp] m={m_label}, a={a_label}, status={status}")
    write_csv("table7_rec_vs_exp.csv", rows)


def run_guard_stability_table(D: int = 1000, Delta: int = 50, Delta_high: int = 150) -> None:
    cases = [("1e-4", mp.mpf("1e-4"), "0", mp.mpf("0")),
             ("1", mp.mpf("1"), "0", mp.mpf("0")),
             ("100", mp.mpf("100"), "1/2", mp.mpf("0.5"))]
    rows = []
    metadata = environment_metadata()
    for m_label, m_value, a_label, a_value in cases:
        S_base, R_base, regime_base, N_base, diag_base = eval_theta(m=m_value, a=a_value, D=D, Delta=Delta)
        S_high, R_high, regime_high, N_high, diag_high = eval_theta(m=m_value, a=a_value, D=D, Delta=Delta_high)
        diff = abs(S_base - S_high)
        bound = R_base + R_high
        status = "pass" if diff <= bound else "fail"
        row = {
            "m_label": m_label, "a_label": a_label, "D": D, "Delta": Delta, "Delta_high": Delta_high,
            "base_regime": regime_base, "high_regime": regime_high, "N_base": N_base, "N_high": N_high,
            "Delta_S": mp.nstr(diff, 30), "bound": mp.nstr(bound, 30),
            "minus_log10_Delta_S": mp.nstr(_safe_minus_log10_local(diff), 30),
            "minus_log10_bound": mp.nstr(_safe_minus_log10_local(bound), 30), "status": status,
        }
        row.update(metadata)
        rows.append(row)
        print(f"[guard] m={m_label}, a={a_label}, status={status}")
    write_csv("table8_guard_stability.csv", rows)


def run_reference_comparison_table(D: int = 1000, Delta: int = 50) -> None:
    cases = [("1e-4", mp.mpf("1e-4"), "0", mp.mpf("0")),
             ("1", mp.mpf("1"), "0", mp.mpf("0")),
             ("100", mp.mpf("100"), "1/2", mp.mpf("0.5"))]
    rows = []
    metadata = environment_metadata()
    ref_dps = D + Delta + 120
    eps_ref = mp.mpf(10) ** (-(D + Delta + 60))
    for m_label, m_value, a_label, a_value in cases:
        S_hat, R_total, regime, N, diagnostics = eval_theta(m=m_value, a=a_value, D=D, Delta=Delta)
        S_ref = reference_mpmath_value(m=m_value, a=a_value, ref_dps=ref_dps)
        diff = abs(S_hat - S_ref)
        bound = R_total + eps_ref
        status = "pass" if diff <= bound else "fail"
        row = {
            "m_label": m_label, "a_label": a_label, "D": D, "Delta": Delta, "ref_dps": ref_dps,
            "selected_regime": regime, "N": N, "Delta_ref": mp.nstr(diff, 30), "eps_ref": mp.nstr(eps_ref, 30),
            "bound": mp.nstr(bound, 30), "minus_log10_Delta_ref": mp.nstr(_safe_minus_log10_local(diff), 30),
            "minus_log10_R_total": diagnostics["minus_log10_R_total"],
            "minus_log10_bound": mp.nstr(_safe_minus_log10_local(bound), 30), "status": status,
        }
        row.update(metadata)
        rows.append(row)
        print(f"[reference] m={m_label}, a={a_label}, status={status}")
    write_csv("table9_reference_comparison.csv", rows)


def main() -> None:
    ensure_output_dir()
    run_regime_selection_tables(D=1000, Delta=50)
    run_crossover_plateau_table(D=1000, Delta=50)
    run_precision_scaling_tables(Delta=50)
    run_rec_vs_exp_table(D=1000, Delta=50)
    run_guard_stability_table(D=1000, Delta=50, Delta_high=150)
    run_reference_comparison_table(D=1000, Delta=50)


if __name__ == "__main__":
    main()
