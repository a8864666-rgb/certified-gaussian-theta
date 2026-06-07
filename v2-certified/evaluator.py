"""
Certified recurrence evaluator for the two canonical real Gaussian theta sums:

    S_0(m)     = sum_{n in Z} exp(-m n^2)
    S_{1/2}(m) = sum_{n in Z} exp(-m (n + 1/2)^2)

corresponding to theta_3(exp(-m)) and theta_2(exp(-m)).

The returned value is interpreted as:

    S(m,a) in S_hat +/- R_total
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Optional, Tuple, Dict, Any
import sys
import mpmath as mp

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)


@dataclass
class EvalDiagnostics:
    m: str
    a: str
    D: int
    Delta: int
    p_bits: int
    mp_dps: int
    N_dir: int
    N_tr: int
    N: int
    regime: str
    lambda_eff: str
    c_rec: int
    c_init: int
    c_alpha: int
    R_tail: str
    R_rec: str
    R_sum: str
    R_trans: str
    R_total: str
    minus_log10_R_tail: str
    minus_log10_R_rec: str
    minus_log10_R_sum: str
    minus_log10_R_trans: str
    minus_log10_R_total: str


def _to_mpf(x) -> mp.mpf:
    if isinstance(x, mp.mpf):
        return x
    return mp.mpf(str(x))


def _ceil_sqrt(x: mp.mpf) -> int:
    if x <= 0:
        return 0
    return int(mp.ceil(mp.sqrt(x)))


def _gamma(k: int, u: mp.mpf) -> mp.mpf:
    k = int(k)
    if k <= 0:
        return mp.mpf("0")
    ku = mp.mpf(k) * u
    if ku >= 1:
        raise ValueError(f"gamma_{k} undefined because k*u >= 1")
    return ku / (1 - ku)


def _safe_minus_log10(x: mp.mpf) -> mp.mpf:
    if x == 0:
        return mp.inf
    return -mp.log10(abs(x))


def _sci(x: mp.mpf, n: int = 30) -> str:
    if x == mp.inf:
        return "inf"
    if x == -mp.inf:
        return "-inf"
    return mp.nstr(x, n)


def _tail_direct(m: mp.mpf, a: mp.mpf, N: int) -> mp.mpf:
    if a == 0:
        num = 2 * mp.exp(-m * (N + 1) ** 2)
        den = 1 - mp.exp(-m * (2 * N + 3))
        return num / den
    if a == mp.mpf("0.5"):
        num = 2 * mp.exp(-m * (N + mp.mpf("1.5")) ** 2)
        den = 1 - mp.exp(-m * (2 * N + 4))
        return num / den
    raise ValueError("Only a=0 and a=0.5 are supported.")


def _tail_transformed(mu: mp.mpf, alpha: mp.mpf, N: int) -> mp.mpf:
    num = 2 * abs(alpha) * mp.exp(-mu * (N + 1) ** 2)
    den = 1 - mp.exp(-mu * (2 * N + 3))
    return num / den


def _planned_depths(m: mp.mpf, D: int, Delta: int) -> Tuple[int, int, mp.mpf]:
    H = mp.mpf(D + Delta)
    log10 = mp.log(10)
    N_dir = _ceil_sqrt((H * log10) / m)
    mu = mp.pi ** 2 / m
    H_tr = H + max(mp.mpf("0"), mp.mpf("0.5") * mp.log10(mp.pi / m))
    N_tr = _ceil_sqrt((H_tr * log10) / mu)
    return N_dir, N_tr, mu


def _increase_until_tail_ok(m: mp.mpf, a: mp.mpf, N: int, target: mp.mpf, regime: str,
                            mu: Optional[mp.mpf] = None, alpha: Optional[mp.mpf] = None) -> int:
    while True:
        if regime == "direct":
            R_tail = _tail_direct(m, a, N)
        elif regime == "transformed":
            if mu is None or alpha is None:
                raise ValueError("mu and alpha required for transformed tail")
            R_tail = _tail_transformed(mu, alpha, N)
        else:
            raise ValueError("regime must be direct or transformed")
        if R_tail <= target:
            return N
        N += 1


def _balanced_sum(terms):
    terms = list(terms)
    if not terms:
        return mp.mpf("0")
    while len(terms) > 1:
        new_terms = []
        for i in range(0, len(terms), 2):
            if i + 1 < len(terms):
                new_terms.append(terms[i] + terms[i + 1])
            else:
                new_terms.append(terms[i])
        terms = new_terms
    return terms[0]


def _direct_terms(m: mp.mpf, a: mp.mpf, N: int):
    q = mp.exp(-m)
    q2 = q * q
    if a == 0:
        if N <= 0:
            return []
        t = q
        r = q ** 3
        terms = [t]
        for _ in range(1, N):
            t = t * r
            r = r * q2
            terms.append(t)
        return terms
    if a == mp.mpf("0.5"):
        t = mp.exp(-m / 4)
        r = q2
        terms = [t]
        for _ in range(0, N):
            t = t * r
            r = r * q2
            terms.append(t)
        return terms
    raise ValueError("Only a=0 and a=0.5 are supported.")


def _transformed_terms(mu: mp.mpf, N: int):
    if N <= 0:
        return []
    q = mp.exp(-mu)
    q2 = q * q
    t = q
    r = q ** 3
    terms = [t]
    for _ in range(1, N):
        t = t * r
        r = r * q2
        terms.append(t)
    return terms


def _direct_value(m: mp.mpf, a: mp.mpf, N: int):
    terms = _direct_terms(m, a, N)
    if a == 0:
        return 1 + 2 * _balanced_sum(terms), terms
    if a == mp.mpf("0.5"):
        return 2 * _balanced_sum(terms), terms
    raise ValueError("Only a=0 and a=0.5 are supported.")


def _transformed_value(m: mp.mpf, a: mp.mpf, N: int):
    mu = mp.pi ** 2 / m
    alpha = mp.sqrt(mp.pi / m)
    terms = _transformed_terms(mu, N)
    if a == 0:
        bracket = 1 + 2 * _balanced_sum(terms)
        return alpha * bracket, terms, mu, alpha
    if a == mp.mpf("0.5"):
        signed_terms = [((-1) ** (j + 1)) * terms[j] for j in range(len(terms))]
        bracket = 1 + 2 * _balanced_sum(signed_terms)
        return alpha * bracket, terms, mu, alpha
    raise ValueError("Only a=0 and a=0.5 are supported.")


def _R_rec_direct(m: mp.mpf, a: mp.mpf, terms: list, u: mp.mpf, c_init: int = 4) -> mp.mpf:
    g_init = _gamma(c_init, u)
    if a == 0:
        abs_weighted_sum = 2 * mp.fsum([abs(t) for t in terms])
        R_init = g_init * abs_weighted_sum
        R_prop = mp.mpf("0")
        for idx, t in enumerate(terms, start=1):
            R_prop += 2 * abs(t) * _gamma(2 * idx, u)
        return R_init + R_prop
    if a == mp.mpf("0.5"):
        abs_weighted_sum = 2 * mp.fsum([abs(t) for t in terms])
        R_init = g_init * abs_weighted_sum
        R_prop = mp.mpf("0")
        for j, t in enumerate(terms):
            if j == 0:
                continue
            R_prop += 2 * abs(t) * _gamma(2 * j, u)
        return R_init + R_prop
    raise ValueError("Only a=0 and a=0.5 are supported.")


def _R_sum_direct(a: mp.mpf, terms: list, u: mp.mpf) -> mp.mpf:
    count = len(terms) + 1
    L_N = int(mp.ceil(mp.log(count, 2))) if count > 1 else 0
    g = _gamma(L_N, u)
    if a == 0:
        return g * (1 + 2 * mp.fsum([abs(t) for t in terms]))
    if a == mp.mpf("0.5"):
        return g * (2 * mp.fsum([abs(t) for t in terms]))
    raise ValueError("Only a=0 and a=0.5 are supported.")


def _R_rec_transformed(alpha: mp.mpf, terms: list, u: mp.mpf, c_init: int = 4) -> mp.mpf:
    g_init = _gamma(c_init, u)
    abs_weighted_sum = 2 * mp.fsum([abs(t) for t in terms])
    R_init_B = g_init * abs_weighted_sum
    R_prop_B = mp.mpf("0")
    for idx, t in enumerate(terms, start=1):
        R_prop_B += 2 * abs(t) * _gamma(2 * idx, u)
    return abs(alpha) * (R_init_B + R_prop_B)


def _R_sum_transformed(alpha: mp.mpf, terms: list, u: mp.mpf) -> mp.mpf:
    count = len(terms) + 1
    L_N = int(mp.ceil(mp.log(count, 2))) if count > 1 else 0
    g = _gamma(L_N, u)
    abs_bracket = 1 + 2 * mp.fsum([abs(t) for t in terms])
    return abs(alpha) * g * abs_bracket


def _R_transformed_factor(mu: mp.mpf, alpha: mp.mpf, terms: list, N: int, u: mp.mpf, c_alpha: int = 4) -> mp.mpf:
    g_alpha = _gamma(c_alpha, u)
    bracket_abs = 1 + 2 * mp.fsum([abs(t) for t in terms])
    bracket_tail = 2 * mp.exp(-mu * (N + 1) ** 2) / (1 - mp.exp(-mu * (2 * N + 3)))
    return g_alpha * abs(alpha) * (bracket_abs + bracket_tail)


def eval_theta(m, a, D: int, Delta: int, force_regime: Optional[str] = None, C_tr: int = 0,
               p_extra: int = 32, c_init: int = 4, c_alpha: int = 4):
    m = _to_mpf(m)
    a = _to_mpf(a)
    if m <= 0:
        raise ValueError("m must be positive")
    if a not in (mp.mpf("0"), mp.mpf("0.5")):
        raise ValueError("Only a=0 and a=0.5 are supported")
    mp_dps = int(D + Delta + mp.ceil(p_extra / mp.log(10, 2)) + 20)
    mp.mp.dps = mp_dps
    p_bits = int(mp.ceil((D + Delta) * mp.log(10, 2))) + p_extra
    u = mp.mpf(2) ** (-p_bits)
    N_dir, N_tr, mu = _planned_depths(m, D, Delta)
    if force_regime is not None:
        if force_regime not in ("direct", "transformed"):
            raise ValueError("force_regime must be None, 'direct', or 'transformed'")
        regime = force_regime
    else:
        regime = "direct" if N_dir <= N_tr + C_tr else "transformed"
    target_tail = mp.mpf(10) ** (-(D + Delta))
    if regime == "direct":
        N = _increase_until_tail_ok(m, a, N_dir, target_tail, "direct")
        S_hat, terms = _direct_value(m, a, N)
        R_tail = _tail_direct(m, a, N)
        R_rec = _R_rec_direct(m, a, terms, u, c_init=c_init)
        R_sum = _R_sum_direct(a, terms, u)
        R_trans = mp.mpf("0")
        lambda_eff = m
    else:
        alpha = mp.sqrt(mp.pi / m)
        N = _increase_until_tail_ok(m, a, N_tr, target_tail, "transformed", mu=mu, alpha=alpha)
        S_hat, terms, mu, alpha = _transformed_value(m, a, N)
        R_tail = _tail_transformed(mu, alpha, N)
        R_rec = _R_rec_transformed(alpha, terms, u, c_init=c_init)
        R_sum = _R_sum_transformed(alpha, terms, u)
        R_trans = _R_transformed_factor(mu, alpha, terms, N, u, c_alpha=c_alpha)
        lambda_eff = mu
    R_total = R_tail + R_rec + R_sum + R_trans
    diagnostics = EvalDiagnostics(
        m=_sci(m), a=_sci(a), D=int(D), Delta=int(Delta), p_bits=int(p_bits), mp_dps=int(mp_dps),
        N_dir=int(N_dir), N_tr=int(N_tr), N=int(N), regime=regime, lambda_eff=_sci(lambda_eff),
        c_rec=2, c_init=int(c_init), c_alpha=int(c_alpha),
        R_tail=_sci(R_tail), R_rec=_sci(R_rec), R_sum=_sci(R_sum), R_trans=_sci(R_trans), R_total=_sci(R_total),
        minus_log10_R_tail=_sci(_safe_minus_log10(R_tail)),
        minus_log10_R_rec=_sci(_safe_minus_log10(R_rec)),
        minus_log10_R_sum=_sci(_safe_minus_log10(R_sum)),
        minus_log10_R_trans=("inf" if R_trans == 0 else _sci(_safe_minus_log10(R_trans))),
        minus_log10_R_total=_sci(_safe_minus_log10(R_total)),
    )
    return S_hat, R_total, regime, N, asdict(diagnostics)


def reference_mpmath(m, a, dps: int = 100):
    m = _to_mpf(m)
    a = _to_mpf(a)
    old_dps = mp.mp.dps
    mp.mp.dps = dps
    q = mp.exp(-m)
    if a == 0:
        val = mp.jtheta(3, 0, q)
    elif a == mp.mpf("0.5"):
        val = mp.jtheta(2, 0, q)
    else:
        raise ValueError("Only a=0 and a=0.5 are supported")
    mp.mp.dps = old_dps
    return val


if __name__ == "__main__":
    S_hat, R_total, regime, N, diag = eval_theta(m=1, a=0, D=100, Delta=50)
    print("S_hat =", mp.nstr(S_hat, 50))
    print("R_total =", mp.nstr(R_total, 20))
    print("regime =", regime)
    print("N =", N)
    for k, v in diag.items():
        print(f"{k}: {v}")
