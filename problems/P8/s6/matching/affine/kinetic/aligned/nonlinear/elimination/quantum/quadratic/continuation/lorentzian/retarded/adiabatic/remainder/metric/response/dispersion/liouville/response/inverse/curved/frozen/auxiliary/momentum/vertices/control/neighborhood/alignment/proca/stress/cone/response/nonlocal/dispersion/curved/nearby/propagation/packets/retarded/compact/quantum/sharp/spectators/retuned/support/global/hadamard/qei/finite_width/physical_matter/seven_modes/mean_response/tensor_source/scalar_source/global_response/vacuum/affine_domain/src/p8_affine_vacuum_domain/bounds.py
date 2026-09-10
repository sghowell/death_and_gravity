"""Exact clock-tube majorants and rank-positive connected field domain."""

from functools import cache

import sympy as sp
from p8_exceptional_vacuum import uniform as previous_uniform

from . import family


@cache
def exact_checks():
    n = sp.Symbol("positive_integer_order", positive=True, integer=True)
    X = family.X
    z = sp.Symbol("formal_derivative_variable")
    q = 1 / X - 1
    poly = sum((-1) ** j * z**j / X ** (j + 1) for j in range(1, 7))
    P = sp.Integer(1)
    checks = {}
    for j in range(7):
        closed = sum(
            sp.prod(n - i for i in range(l))
            * q ** (j - l)
            * sp.expand(poly**l).coeff(z, j)
            / sp.factorial(l)
            for l in range(j + 1)
        )
        checks["normalized_rational_power_derivative_" + str(j)] = sp.factor(P - closed)
        P = sp.factor((q * sp.diff(P, X) + (n - j) * sp.diff(q, X) * P) / (j + 1))
    P = n * X * X
    for j in range(7):
        checks["literal_gaussian_bump_derivative_" + str(j)] = sp.expand(
            sp.diff(n * X * X * sp.exp(-n * X * X), X, j) * sp.exp(n * X * X) - P
        )
        P = sp.expand(sp.diff(P, X) - 2 * n * X * P)
    return checks


@cache
def switch_jet_bounds():
    n = sp.Symbol("positive_even_order", positive=True, integer=True)
    z = sp.Symbol("jet_variable")
    X = family.X
    order = 6
    q0 = sp.Rational(1, 9)
    qpoly = sum(sp.Rational(10, 9) ** (j + 1) * z**j for j in range(1, order + 1))
    D = []
    for j in range(order + 1):
        D.append(
            sum(
                sp.binomial(n, l) * q0 ** (n - l) * sp.expand(qpoly**l).coeff(z, j)
                for l in range(j + 1)
            )
        )
    # The finite inverse Taylor sum has norm <2 once sum(D)<1/2.
    # Keep this rational ceiling instead of needlessly expanding products
    # of seven enormous exact denominators.
    P = n * X * X
    bump = []
    polynomials = []
    for j in range(order + 1):
        polynomials.append(P)
        poly = sp.Poly(P, X, n)
        majorant = sum(
            abs(c) * sp.Rational(11, 10) ** powers[0] * n ** powers[1]
            for powers, c in poly.terms()
        )
        bump.append(majorant / (sp.factorial(j) * 2**n))
        P = sp.expand(sp.diff(P, X) - 2 * n * X * P)
    # ||1-w||_j<2 and ||1/(1+q^n)||_j<2, so each normalized
    # derivative of their product with q^n is bounded by 4 sum_(i<=j)D_i.
    S = [4 * sum(D[: j + 1]) for j in range(order + 1)]
    selected = lambda vals: [v.subs(n, family.N) for v in vals]
    # Leave the general formulas factored; exact selected rational entries
    # are separately evaluated without expanding a degree-n polynomial.
    return {
        "order_symbol": n,
        "normalized_q_power_bounds": D,
        "inverse_denominator_six_jet_norm_upper": sp.Integer(2),
        "one_minus_bump_six_jet_norm_upper": sp.Integer(2),
        "stripped_bump_derivative_polynomials": polynomials,
        "normalized_bump_bounds": bump,
        "normalized_switch_bounds": S,
        "selected_normalized_switch_bounds": selected(S),
        "checks": {
            "rational_exponential_majorant": sp.Rational(1)
            + sp.Rational(81, 100)
            + sp.Rational(81, 100) ** 2 / 2
            - sp.Rational(42761, 20000)
        },
        "bounds": {
            "selected_q_power_six_jet_norm_below_one_half": sum(selected(D))
            < sp.Rational(1, 2),
            "selected_bump_six_jet_norm_below_one": sum(selected(bump)) < 1,
            "e_power_081_above_two_by_positive_Taylor_terms": 1
            + sp.Rational(81, 100)
            + sp.Rational(81, 100) ** 2 / 2
            > 2,
            "q_power_even_step_majorant_below_one_fiftieth": sp.Rational(1, 81)
            * (1 + sp.Rational(2, family.N - 5)) ** 7
            < sp.Rational(1, 50),
            "bump_even_step_majorant_below_one_third": sp.Rational(1, 4)
            * (1 + sp.Rational(2, family.N)) ** 7
            < sp.Rational(1, 3),
        },
    }


@cache
def tube():
    sj = switch_jet_bounds()["selected_normalized_switch_bounds"]
    S = sum(sj[:5])
    Sd = sum((j + 1) * sj[j + 1] for j in range(5))
    H = sp.Integer(210)
    Xm = sp.Rational(11, 10)
    Xi = sum(sp.Rational(10, 9) ** (j + 1) for j in range(5))
    E = S + Xm * Sd
    F = Xm * S
    ro = 1 + Xm * H
    rn = ro + H * F
    inverse = sp.Rational(5, 4) * sum((sp.Rational(5, 4) * 233) ** j for j in range(5))
    T = (2 * E + E * E) * ro + H * F
    errs = {
        "F2": H * F / 2,
        "A3": H * E * Xi,
        "A4": H * E * Xi + sp.Rational(7, 4) * H * H * T * inverse**2,
        "A5": H * H * T * inverse**2 * Xi,
    }
    gu = previous_uniform.data()
    rows = gu["linear_and_constant_envelopes"]
    G = family.N * sum(rows["linear"].values()) + sum(rows["constant"].values())
    errs["F"] = S * G
    return {
        "weighted_four_jet_error_bounds": errs,
        "switch_normalized_four_jet": S,
        "switch_prime_normalized_four_jet": Sd,
        "new_R_norm_upper": rn,
        "inverse_R_jet_upper": inverse,
        "bounds": {
            "R_jet_below_233": rn < 233,
            "clock_tensor_above_four_fifths": sp.Rational(9, 10) - S / 10
            > sp.Rational(4, 5),
            **{
                key + "_error_below_one_e_minus_400": value < sp.Rational(1, 10**400)
                for key, value in errs.items()
            },
        },
    }


@cache
def domain():
    n = sp.Integer(family.N)
    target = sp.Rational(1, 2) + 1 / (8 * n)
    floors = {
        "negative_X": 1 - sp.Rational(5, 32) / n,
        "first_quarter": sp.Rational(3, 5),
        "quarter_to_half": sp.Rational(497, 800),
        "short_above_half": sp.Rational(499, 800),
        "remaining_X": target,
    }
    return {
        "selected_region_R_floors": floors,
        "selected_R_floor": target,
        "selected_full_quotient_factor_floor": 1 / (4 * n),
        "bounds": {
            "e_above_eight_thirds_from_positive_series": sum(
                sp.Rational(1, sp.factorial(j)) for j in range(5)
            )
            > sp.Rational(8, 3),
            "first_quarter_step_tail_small": 3 ** (-n) < sp.Rational(1, 40),
            "negative_X_step_tail_small": (1 / (4 * n)) ** n < 1 / (16 * n),
            "away_from_zero_bump_below_one_hundredth": sp.Rational(64, 2**64)
            < sp.Rational(1, 100),
            "short_above_half_log_odds_below_one": 1 / (2 * (1 - 1 / (4 * n))) < 1,
            "all_region_floors_above_target": all(v >= target for v in floors.values()),
            "strict_auxiliary_quotient_margin": 1 / (4 * n) > 0,
        },
    }
