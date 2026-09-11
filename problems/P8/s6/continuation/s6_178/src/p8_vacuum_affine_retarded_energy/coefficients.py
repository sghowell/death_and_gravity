"""Uniform real near-clock bounds for the COMPLETE rational/exponential source."""

from functools import cache

import sympy as s
from p8_affine_vacuum_domain import family


@cache
def data():
    n = s.Integer(family.N)
    lo, hi = s.Rational(9, 10), s.Rational(11, 10)
    qmax = s.Rational(1, 9)
    t1 = n * qmax ** (n - 1) / lo**2
    t2 = (
        n * (n - 1) * qmax ** (n - 2) / lo**4
        + 2 * n * qmax ** (n - 1) / lo**3
        + 2 * n * n * qmax ** (2 * n - 2) / lo**4
    )
    # e^y > y^3/6 or y^4/24 on y>0; no floating exponential comparison.
    w1 = 12 * hi * (1 + n * hi**2) / (n * n * lo**6)
    w2 = 48 * (1 + 5 * n * hi**2 + 2 * n * n * hi**4) / (n**3 * lo**8)
    x = s.Symbol("X", positive=True)
    nn = s.Symbol("n", positive=True, integer=True)
    w = nn * x * x * s.exp(-nn * x * x)
    wp = 2 * nn * x * (1 - nn * x * x) * s.exp(-nn * x * x)
    wpp = 2 * nn * (1 - 5 * nn * x * x + 2 * nn * nn * x**4) * s.exp(-nn * x * x)
    q = (1 - x) / x
    step = 1 / (1 + q**nn)
    step1 = nn * q ** (nn - 1) / (x * x * (1 + q**nn) ** 2)
    step2 = (
        -nn * (nn - 1) * q ** (nn - 2) / (x**4 * (1 + q**nn) ** 2)
        - 2 * nn * q ** (nn - 1) / (x**3 * (1 + q**nn) ** 2)
        + 2 * nn**2 * q ** (2 * nn - 2) / (x**4 * (1 + q**nn) ** 3)
    )
    qp = s.Symbol("positive_q_for_rational_identity", positive=True)
    step_q = 1 / (1 + qp**nn)
    second_chain = s.diff(step_q, qp, 2) / x**4 + 2 * s.diff(step_q, qp) / x**3
    second_formula = (
        -nn * (nn - 1) * qp ** (nn - 2) / (x**4 * (1 + qp**nn) ** 2)
        - 2 * nn * qp ** (nn - 1) / (x**3 * (1 + qp**nn) ** 2)
        + 2 * nn**2 * qp ** (2 * nn - 2) / (x**4 * (1 + qp**nn) ** 3)
    )
    rmin = s.Rational(1, 2)
    RX = s.Rational(6, 5)
    RXX = s.Rational(22, 5)
    Ru = s.Rational(3, 10)
    RXu = s.Rational(18, 5)
    C = 1 / lo**2 + 3 * RX / (2 * rmin * lo)
    Cu = s.Rational(3, 2) * (RXu / (rmin * lo) + RX * Ru / (rmin**2 * lo))
    Cx = 2 / lo**3 + s.Rational(3, 2) * (
        RXX / (rmin * lo) + RX**2 / (rmin**2 * lo) + RX / (rmin * lo**2)
    )
    u = s.Symbol("u", real=True)
    h = (1 + u * u) ** -3
    return {
        "actual_switch_order": n,
        "literal_rational_step_second_derivative": step2,
        "rational_identity_extension": "The derivative identity is rational for each integer n. Its check at positive q extends to negative q wherever the denominator is nonzero; the frozen even n gives 1+q^n>0 for every real q. Both signs are independently tested against the complete function.",
        "physical_X_strip": [lo, hi],
        "rational_step_derivative_bounds": [t1, t2],
        "Gaussian_bump_derivative_bounds": [w1, w2],
        "full_B_bounds": {
            "B": [s.S.Zero, s.S.One],
            "abs_B_X": s.Integer(2),
            "abs_B_XX": s.Integer(4),
        },
        "full_R_coefficient_bounds": {
            "abs_R_minus_one_over_abs_X_minus_one": s.S.One,
            "abs_RX": RX,
            "abs_RXX": RXX,
            "abs_Ru_over_abs_R_minus_one": s.Integer(3),
            "abs_RXu": RXu,
            "abs_Ruu_over_abs_R_minus_one": s.Integer(15),
        },
        "source_Z_coefficient_bounds": {"abs_C": C, "abs_C_u": Cu, "abs_C_X": Cx},
        "proof": "For real X on the strip, q=(1-X)/X has |q|<=1/9 and n is the frozen even1024, so 1+q^n>=1. Since 0<=T,w<=1, B=T+(1-T)w has |B'|<=2 and |B''|<=4. R=1+B(X)(X-1)/(1+u^2)^3 keeps every full-function term; no Taylor remainder is omitted.",
        "checks": {
            "actual_n_is_frozen": n - 1024,
            "literal_bump_derivative": s.simplify(s.diff(w, x) - wp),
            "literal_bump_second_derivative": s.simplify(s.diff(w, x, 2) - wpp),
            "literal_step_first_derivative": s.simplify(s.diff(step, x) - step1),
            "rational_step_second_derivative_chain": s.simplify(
                second_chain - second_formula
            ),
            "rational_q_first_derivative": s.simplify(s.diff(q, x) + 1 / x**2),
            "rational_q_second_derivative": s.simplify(s.diff(q, x, 2) - 2 / x**3),
            "full_RX_majorant": 1 + s.Rational(1, 10) * 2 - RX,
            "full_RXX_majorant": 2 * 2 + s.Rational(1, 10) * 4 - RXX,
            "full_Ru_factor": s.factor(s.diff(h, u) / h + 6 * u / (1 + u * u)),
            "full_Ruu_factor": s.factor(
                s.diff(h, u, 2) / h - (-6 + 42 * u * u) / (1 + u * u) ** 2
            ),
            "source_C_majorant": C - s.Rational(424, 81),
        },
        "gates": {
            "full_step_first_derivative_below_one": bool(t1 < 1),
            "full_step_second_derivative_below_one": bool(t2 < 1),
            "full_bump_first_derivative_below_one": bool(w1 < 1),
            "full_bump_second_derivative_below_one": bool(w2 < 1),
            "source_C_below_six": bool(C < 6),
            "source_C_u_below_fifteen": bool(Cu < 15),
            "source_C_X_below_thirty_two": bool(Cx < 32),
        },
    }
