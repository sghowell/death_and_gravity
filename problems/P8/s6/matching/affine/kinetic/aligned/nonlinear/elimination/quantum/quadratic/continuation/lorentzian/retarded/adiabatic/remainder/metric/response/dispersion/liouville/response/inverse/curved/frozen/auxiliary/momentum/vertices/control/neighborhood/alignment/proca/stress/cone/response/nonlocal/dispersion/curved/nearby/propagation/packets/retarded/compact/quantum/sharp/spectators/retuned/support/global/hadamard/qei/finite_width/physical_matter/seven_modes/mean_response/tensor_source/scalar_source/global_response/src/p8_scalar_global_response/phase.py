"""Exact global weighted scalar phase estimates with the actual moving CCR."""

from fractions import Fraction
from functools import cache

import sympy as sp
from p8_proca_mean_response import tails
from p8_scalar_mean_source import geometry, model

u, k = model.u, model.k
MOMENTUM_WEIGHT = sp.Rational(2, 5)


def half_line_integral(exponent):
    if isinstance(exponent, bool) or not isinstance(
        exponent, (int, Fraction, sp.Rational)
    ):
        raise TypeError("Require an exact integer or half-integer exponent")
    exponent = sp.Rational(exponent)
    if exponent <= sp.Rational(1, 2):
        raise ValueError("The half-line power is not integrable")
    if exponent.q == 1:
        return tails.half_line_integral_upper(int(exponent))
    if exponent.q == 2:
        n = int(exponent - sp.Rational(1, 2))
        return sp.prod(sp.Rational(2 * j, 2 * j + 1) for j in range(1, n))
    raise ValueError("Require integer or half-integer decay")


def envelope(value, extra_half_order=0, *, require_integrable=True):
    if isinstance(value, bool) or not isinstance(value, (int, sp.Expr)):
        raise TypeError("Require an exact scalar rational kernel")
    value = sp.sympify(value)
    if value.has(
        sp.Float, sp.oo, -sp.oo, sp.nan, sp.zoo, sp.I
    ) or not value.free_symbols <= {u, k}:
        raise ValueError("Require finite exact real declared-variable coefficients")
    if (
        isinstance(extra_half_order, bool)
        or not isinstance(extra_half_order, int)
        or extra_half_order not in (0, 1, 2)
    ):
        raise TypeError("Extra half order counts zero, one or two momentum indices")
    if not isinstance(require_integrable, bool):
        raise TypeError("The integrability requirement must be Boolean")
    value = sp.cancel(value)
    if value == 0:
        return {
            "constant": sp.Integer(0),
            "parity": 0,
            "power": 1,
            "extra_half_order": extra_half_order,
            "integral_upper": sp.Integer(0) if require_integrable else None,
            "coefficientwise_remainder_is_nonnegative": True,
        }
    numerator, denominator = sp.fraction(value)
    try:
        den = sp.Poly(denominator, u, domain=sp.QQ)
        original = sp.Poly(numerator, u, k, domain=sp.QQ)
    except (sp.PolynomialError, sp.polys.polyerrors.CoercionFailed) as exc:
        raise ValueError(
            "Require rational polynomial numerator and positive time denominator"
        ) from exc
    if any(powers[0] % 2 or c < 0 for powers, c in den.terms()) or den.nth(0) <= 0:
        raise ValueError("Not a positive even denominator")
    num = sp.Poly(
        sum(abs(c) * 2 ** powers[1] * u ** powers[0] for powers, c in original.terms()),
        u,
    )
    parities = {powers[0] % 2 for powers, c in num.terms() if c}
    if len(parities) != 1:
        raise ValueError("Mixed time parity")
    parity = next(iter(parities))
    power = (den.degree() - num.degree() + parity) // 2
    exponent = power + sp.Rational(extra_half_order, 2)
    if power < 0 or (
        require_integrable and exponent <= (1 if parity else sp.Rational(1, 2))
    ):
        raise ValueError("No integrable coefficient envelope: " + str(value))
    weighted = sp.Poly(sp.cancel(num.as_expr() / u**parity) * (1 + u * u) ** power, u)
    ratios = []
    for powers, c in weighted.terms():
        if c and den.nth(powers[0]) <= 0:
            raise ValueError("Denominator coefficient support mismatch")
        if c:
            ratios.append(c / den.nth(powers[0]))
    C = max(ratios)
    remainder = sp.Poly(C * den.as_expr() - weighted.as_expr(), u)
    if not all(co >= 0 for co in remainder.all_coeffs()):
        raise ValueError("The claimed coefficientwise envelope failed")
    integral = (
        (C / (2 * (exponent - 1)) if parity else C * half_line_integral(exponent))
        if require_integrable
        else None
    )
    return {
        "constant": C,
        "parity": parity,
        "power": int(power),
        "extra_half_order": extra_half_order,
        "integral_upper": integral,
        "coefficientwise_remainder_is_nonnegative": True,
    }


def uniform_upper(row):
    exponent = row["power"] + sp.Rational(row["extra_half_order"], 2)
    if row["parity"] == 0 and exponent >= 0:
        return row["constant"]
    if row["parity"] == 1 and exponent >= 1:
        return row["constant"] / 2
    if row["parity"] == 1 and exponent >= sp.Rational(1, 2):
        return row["constant"]
    raise ValueError("The pointwise envelope is unbounded")


@cache
def data():
    t = 1 + u * u
    T = sp.diag(1, 1, MOMENTUM_WEIGHT * t**-3, MOMENTUM_WEIGHT * t**-3)
    old = model.phase.data()["regular_density_generator"].subs(
        model.global_model.substitution(), simultaneous=True
    )
    old = old.subs(model.global_model.old.q, k * k / t**4)
    K = (T * old * T.inv() + sp.diff(T, u) * T.inv()).applyfunc(sp.factor)
    remainder = (K + sp.diag(0, 0, 6 * u / t, 6 * u / t)).applyfunc(sp.factor)
    envelopes = [envelope(value) for value in remainder]
    total = sum(row["integral_upper"] for row in envelopes)
    B = geometry.canonical()["old_to_natural_phase"]
    initial = (T * B.inv()).subs(u, 0)
    initial_norm = max(sum(abs(value) for value in initial.row(i)) for i in range(4))
    momentum_rows = [
        sum(
            row["constant"] * (sp.Rational(1, 2) if row["parity"] else 1)
            for row in envelopes[4 * i : 4 * i + 4]
        )
        for i in (2, 3)
    ]
    return {
        "weighted_generator": K,
        "negative_damping_removed_generator": remainder,
        "old_to_weighted_phase": T,
        "old_generator": old,
        "entry_envelopes": envelopes,
        "integrable_entry_sum_upper": total,
        "natural_anchor_to_weighted_initial_map": initial,
        "natural_anchor_initial_infinity_norm": initial_norm,
        "momentum_remainder_row_bounds_per_one_over_t": momentum_rows,
    }
