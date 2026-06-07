"""Certified recurrence evaluator for real Gaussian theta sums.

Implements the v2.1 near-boundary artifact for
    S(m, a) = sum_{n in Z} exp(-m (n+a)^2)
with a in {0, 1/2, 1/3, 1/4, 1/6}.

The implementation is intentionally transparent rather than optimized: it reports
componentwise radii for tail, recurrence, summation, transformation, coefficient,
and total errors.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from fractions import Fraction
from typing import Dict, Iterable, List, Tuple
import math
import sys
import mpmath as mp

try:
    sys.set_int_max_str_digits(0)
except AttributeError:
    pass

SUPPORTED_SHIFTS = {
    Fraction(0, 1),
    Fraction(1, 2),
    Fraction(1, 3),
    Fraction(1, 4),
    Fraction(1, 6),
}


def _mp(x) -> mp.mpf:
    if isinstance(x, Fraction):
        return mp.mpf(x.numerator) / x.denominator
    return mp.mpf(x)


def gamma_k(k: int, u: mp.mpf) -> mp.mpf:
    ku = mp.mpf(k) * u
    if ku >= 1:
        raise ValueError("gamma_k undefined because k*u >= 1; increase precision")
    return ku / (1 - ku)


def normalize_shift(a) -> Fraction:
    if isinstance(a, Fraction):
        f = a
    elif isinstance(a, int):
        f = Fraction(a, 1)
    elif isinstance(a, str):
        f = Fraction(a)
    else:
        f = Fraction(a).limit_denominator()
    f = f % 1
    if f not in SUPPORTED_SHIFTS:
        raise ValueError(f"unsupported shift {f}; supported: {sorted(SUPPORTED_SHIFTS)}")
    return f


def cosine_weight(k: int, a: Fraction) -> mp.mpf:
    # Exact dyadic weights for the supported shifts.
    if a == 0:
        return mp.mpf(1)
    if a == Fraction(1, 2):
        return mp.mpf(1 if k % 2 == 0 else -1)
    if a == Fraction(1, 3):
        return mp.mpf(1 if k % 3 == 0 else mp.mpf('-0.5'))
    if a == Fraction(1, 4):
        r = k % 4
        return mp.mpf(1 if r == 0 else -1 if r == 2 else 0)
    if a == Fraction(1, 6):
        r = k % 6
        table = {0: 1, 1: mp.mpf('0.5'), 2: mp.mpf('-0.5'), 3: -1, 4: mp.mpf('-0.5'), 5: mp.mpf('0.5')}
        return mp.mpf(table[r])
    raise ValueError("unsupported shift")


@dataclass
class ThetaResult:
    m: str
    q: str
    a: str
    requested_digits: int
    guard_digits: int
    working_bits: int
    planned_direct_depth: int
    planned_transformed_depth: int
    selected_regime: str
    retained_depth: int
    value: str
    R_tail: str
    R_rec: str
    R_sum: str
    R_trans: str
    R_coef: str
    R_total: str
    certified_digits: float

    def as_dict(self) -> Dict[str, object]:
        return asdict(self)


def plan_depths(m: mp.mpf, D: int, guard: int) -> Tuple[int, int]:
    # Mirrors the manuscript planning formulas.
    H = mp.mpf(D + guard) * mp.log(10)
    N_dir = int(mp.ceil(mp.sqrt(H / m)))
    mu = mp.pi ** 2 / m
    Htr_dec = mp.mpf(D + guard) + max(mp.mpf(0), mp.mpf('0.5') * mp.log10(mp.pi / m))
    N_tr = int(mp.ceil(mp.sqrt((Htr_dec * mp.log(10)) / mu)))
    return max(0, N_dir), max(0, N_tr)


def direct_tail_bound(m: mp.mpf, N: int, a: Fraction) -> mp.mpf:
    if a == 0:
        return 2 * mp.e ** (-m * (N + 1) ** 2) / (1 - mp.e ** (-m * (2 * N + 3)))
    if a == Fraction(1, 2):
        return 2 * mp.e ** (-m * (N + mp.mpf('1.5')) ** 2) / (1 - mp.e ** (-m * (2 * N + 4)))
    raise ValueError("direct representation only implemented for a=0 and a=1/2")


def transformed_tail_bound(m: mp.mpf, N: int, alpha: mp.mpf | None = None) -> mp.mpf:
    mu = mp.pi ** 2 / m
    if alpha is None:
        alpha = mp.sqrt(mp.pi / m)
    return 2 * abs(alpha) * mp.e ** (-mu * (N + 1) ** 2) / (1 - mp.e ** (-mu * (2 * N + 3)))


def _direct_terms(m: mp.mpf, a: Fraction, N: int) -> Tuple[mp.mpf, List[mp.mpf]]:
    if a == 0:
        beta0 = mp.mpf(1)
        terms = [2 * mp.e ** (-m * j * j) for j in range(1, N + 1)]
        return beta0, terms
    if a == Fraction(1, 2):
        beta0 = mp.mpf(0)
        terms = [2 * mp.e ** (-m * (j + mp.mpf('0.5')) ** 2) for j in range(0, N + 1)]
        return beta0, terms
    raise ValueError("direct representation only implemented for a=0 and a=1/2")


def _transformed_terms(m: mp.mpf, a: Fraction, N: int) -> Tuple[mp.mpf, List[mp.mpf]]:
    mu = mp.pi ** 2 / m
    beta0 = mp.mpf(1)
    terms = [2 * cosine_weight(j, a) * mp.e ** (-mu * j * j) for j in range(1, N + 1)]
    return beta0, terms


def _sum_abs_terms(beta0: mp.mpf, terms: Iterable[mp.mpf]) -> mp.mpf:
    return abs(beta0) + mp.fsum([abs(t) for t in terms])


def _finite_radii_direct(m: mp.mpf, a: Fraction, N: int, p_bits: int, cinit: int = 4) -> Tuple[mp.mpf, mp.mpf, mp.mpf]:
    u = mp.mpf(2) ** (-p_bits)
    beta0, terms = _direct_terms(m, a, N)
    abs_sum = _sum_abs_terms(beta0, terms)
    R_init = gamma_k(cinit, u) * mp.fsum([abs(t) for t in terms])
    R_prop = mp.mpf(0)
    if a == 0:
        for j in range(1, N + 1):
            R_prop += abs(2 * mp.e ** (-m * j * j)) * gamma_k(2 * j, u)
        count = N + 1
    else:
        for j in range(1, N + 1):
            R_prop += abs(2 * mp.e ** (-m * (j + mp.mpf('0.5')) ** 2)) * gamma_k(2 * j, u)
        count = N + 1
    L = int(math.ceil(math.log2(max(1, count))))
    R_sum = gamma_k(L, u) * abs_sum
    return R_init + R_prop, R_sum, mp.mpf(0)


def _finite_radii_transformed(m: mp.mpf, a: Fraction, N: int, p_bits: int, cinit: int = 4, calpha: int = 4) -> Tuple[mp.mpf, mp.mpf, mp.mpf, mp.mpf]:
    u = mp.mpf(2) ** (-p_bits)
    alpha = mp.sqrt(mp.pi / m)
    beta0, terms = _transformed_terms(m, a, N)
    abs_sum = _sum_abs_terms(beta0, terms)
    R_init_B = gamma_k(cinit, u) * mp.fsum([abs(t) for t in terms])
    mu = mp.pi ** 2 / m
    R_prop_B = mp.mpf(0)
    for j in range(1, N + 1):
        R_prop_B += abs(2 * cosine_weight(j, a) * mp.e ** (-mu * j * j)) * gamma_k(2 * j, u)
    L = int(math.ceil(math.log2(max(1, N + 1))))
    R_sum_B = gamma_k(L, u) * abs_sum
    R_B_tail = 2 * mp.e ** (-mu * (N + 1) ** 2) / (1 - mp.e ** (-mu * (2 * N + 3)))
    R_rec = abs(alpha) * (R_init_B + R_prop_B)
    R_sum = abs(alpha) * R_sum_B
    R_trans = gamma_k(calpha, u) * abs(alpha) * (abs_sum + R_B_tail)
    R_coef = mp.mpf(0)  # binary-exact weights for supported shifts
    return R_rec, R_sum, R_trans, R_coef


def evaluate(m, a=0, D: int = 1000, guard: int = 50, pextra: int = 64, force_regime: str | None = None) -> ThetaResult:
    a_frac = normalize_shift(a)
    m_mp = _mp(m)
    if m_mp <= 0:
        raise ValueError("m must be positive")

    p_bits = int(math.ceil((D + guard) * math.log2(10))) + pextra
    mp.mp.prec = p_bits
    m_mp = _mp(m)
    q = mp.e ** (-m_mp)
    N_dir, N_tr = plan_depths(m_mp, D, guard)

    if force_regime is not None:
        regime = force_regime
    elif a_frac in (Fraction(0, 1), Fraction(1, 2)) and N_dir <= N_tr:
        regime = "direct"
    else:
        regime = "transformed"

    if regime == "direct":
        if a_frac not in (Fraction(0, 1), Fraction(1, 2)):
            raise ValueError("direct representation only supported for a=0 and a=1/2")
        N = N_dir
        # Ensure manuscript-style target tail.
        while direct_tail_bound(m_mp, N, a_frac) > mp.mpf(10) ** (-(D + guard)):
            N += 1
        beta0, terms = _direct_terms(m_mp, a_frac, N)
        value = beta0 + mp.fsum(terms)
        R_tail = direct_tail_bound(m_mp, N, a_frac)
        R_rec, R_sum, R_trans = _finite_radii_direct(m_mp, a_frac, N, p_bits)
        R_coef = mp.mpf(0)
    elif regime == "transformed":
        N = N_tr
        alpha = mp.sqrt(mp.pi / m_mp)
        while transformed_tail_bound(m_mp, N, alpha) > mp.mpf(10) ** (-(D + guard)):
            N += 1
        beta0, terms = _transformed_terms(m_mp, a_frac, N)
        bracket = beta0 + mp.fsum(terms)
        value = alpha * bracket
        R_tail = transformed_tail_bound(m_mp, N, alpha)
        R_rec, R_sum, R_trans, R_coef = _finite_radii_transformed(m_mp, a_frac, N, p_bits)
    else:
        raise ValueError("force_regime must be 'direct', 'transformed', or None")

    R_total = R_tail + R_rec + R_sum + R_trans + R_coef
    certified_digits = float(-mp.log10(R_total)) if R_total > 0 else float("inf")
    return ThetaResult(
        m=mp.nstr(m_mp, 20),
        q=mp.nstr(q, 20),
        a=str(a_frac),
        requested_digits=D,
        guard_digits=guard,
        working_bits=p_bits,
        planned_direct_depth=N_dir,
        planned_transformed_depth=N_tr,
        selected_regime=regime,
        retained_depth=N,
        value=mp.nstr(value, min(D + 20, 80)),
        R_tail=mp.nstr(R_tail, 20),
        R_rec=mp.nstr(R_rec, 20),
        R_sum=mp.nstr(R_sum, 20),
        R_trans=mp.nstr(R_trans, 20),
        R_coef=mp.nstr(R_coef, 20),
        R_total=mp.nstr(R_total, 20),
        certified_digits=round(certified_digits, 3),
    )
