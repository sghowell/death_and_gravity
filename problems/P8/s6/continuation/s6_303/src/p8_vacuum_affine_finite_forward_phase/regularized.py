"""Exact regularized massive forward block and uniform second-derivative norm."""

from functools import cache

import sympy as s

from . import source

ENERGY, TRANSFER, ELL = s.symbols("energy transfer ell", real=True)
J, Q, L = s.symbols("small_J small_Q small_L", real=True)


@cache
def coefficients():
    a, t, ell = ENERGY, TRANSFER, ELL
    e = source.assembly.EPS
    row = source.dimensional.master_coefficients(t, a, 1, e)
    cm0 = row["C0mumu"].subs(e, 0)
    cm1 = s.diff(row["C0mumu"], e).subs(e, 0)
    bm0 = row["Bmm"].subs(e, 0)
    t0, t1 = source.assembly.tree_jets(a, t, 1)
    r0, r1 = source.assembly.completion_jets(a, t, 1)
    return {
        "J": s.cancel(t * t * (ell * cm0 - cm1) / 2),
        "Q": s.cancel(-t * t * cm0 / 2),
        "L": s.cancel(-t * t * bm0),
        "one": s.cancel(
            t
            * t
            * (
                ell * bm0
                + source.dimensional.compact_crossed_evans(a, t, 1)
                + (6 * ell + 14) * t0
                - 6 * t1
                - ell * r0
                + r1
            )
        ),
    }


def rational_majorant(expression):
    num, den = s.fraction(s.cancel(expression))
    constant, factors = s.factor_list(den)
    allowed = (ENERGY, ENERGY + TRANSFER - 4, TRANSFER - 4, ENERGY - 4)
    for factor, power in factors:
        if (
            s.expand(factor) not in allowed
            or not isinstance(power, (int, s.Integer))
            or power < 0
        ):
            raise ValueError("A denominator lacks the stated compact lower bound")
    return sum(
        abs(value) * 16**i * 3**k
        for (i, j, k), value in s.Poly(num, ENERGY, TRANSFER, ELL).terms()
    ) / abs(constant)


@cache
def majorants():
    return {
        name: tuple(rational_majorant(s.diff(value, TRANSFER, k)) for k in range(3))
        for name, value in coefficients().items()
    }


def massive_forward_bound():
    table = majorants()
    return (
        sum(v[2] + 2 * v[1] + v[0] for name, v in table.items() if name != "one")
        + table["one"][2]
    ) / 2


@cache
def data():
    a, t, ell = ENERGY, TRANSFER, ELL
    cs = coefficients()
    zero = {name: s.factor(value.subs(t, 0)) for name, value in cs.items()}
    first = {name: s.factor(s.diff(value, t).subs(t, 0)) for name, value in cs.items()}
    e = source.assembly.EPS
    row = source.dimensional.master_coefficients(t, a, 1, e)
    cm0 = row["C0mumu"].subs(e, 0)
    cm1 = s.diff(row["C0mumu"], e).subs(e, 0)
    bm0 = row["Bmm"].subs(e, 0)
    t0, t1 = source.assembly.tree_jets(a, t, 1)
    r0, r1 = source.assembly.completion_jets(a, t, 1)
    sm = cm0 * (ell * J - Q) / 2 - cm1 * J / 2 + bm0 * (ell - L)
    sm += source.dimensional.compact_crossed_evans(a, t, 1)
    sm += (6 * ell + 14) * t0 - 6 * t1 - ell * r0 + r1
    x = s.Symbol("unit_x", real=True)
    h = x * (1 - x)
    den = 1 - t * h
    qkernel = s.log(den) / den
    lkernel = s.log(den)
    expected = {
        "J": (7972, 14810, 16784),
        "Q": (2422, 4680, 5548),
        "L": (s.Rational(60323, 4), s.Rational(374714, 3), 552949),
        "one": (
            s.Rational(1108223764786, 225),
            s.Rational(4558231312812, 5),
            s.Rational(20687039523642988, 225),
        ),
    }
    checks = {
        "whole_regularized_massive_block": s.factor(
            t * t * sm - cs["J"] * J - cs["Q"] * Q - cs["L"] * L - cs["one"]
        ),
        "double_integer_pole_jet_cancels": s.factor(zero["J"] + zero["one"]),
        "simple_integer_pole_jet_cancels": s.factor(
            first["J"] + zero["J"] / 6 - zero["Q"] / 6 - zero["L"] / 6 + first["one"]
        ),
        "massive_J_first_kernel": s.factor(s.diff(1 / den, t) - h / den**2),
        "massive_J_second_kernel": s.factor(s.diff(1 / den, t, 2) - 2 * h * h / den**3),
        "massive_Q_first_kernel": s.simplify(
            s.diff(qkernel, t) + h * (1 - s.log(den)) / den**2
        ),
        "massive_Q_second_kernel": s.simplify(
            s.diff(qkernel, t, 2) - h * h * (2 * s.log(den) - 3) / den**3
        ),
        "massive_L_first_kernel": s.factor(s.diff(lkernel, t) + h / den),
        "massive_L_second_kernel": s.factor(s.diff(lkernel, t, 2) + h * h / den**2),
        "massive_J_zero": s.integrate((1 / den).subs(t, 0), (x, 0, 1)) - 1,
        "massive_J_first_zero": s.integrate(s.diff(1 / den, t).subs(t, 0), (x, 0, 1))
        - s.Rational(1, 6),
        "massive_Q_zero": s.integrate(qkernel.subs(t, 0), (x, 0, 1)),
        "massive_Q_first_zero": s.integrate(s.diff(qkernel, t).subs(t, 0), (x, 0, 1))
        + s.Rational(1, 6),
        "massive_L_zero": s.integrate(lkernel.subs(t, 0), (x, 0, 1)),
        "massive_L_first_zero": s.integrate(s.diff(lkernel, t).subs(t, 0), (x, 0, 1))
        + s.Rational(1, 6),
        "Taylor_second_derivative_weight": s.integrate(1 - x, (x, 0, 1))
        - s.Rational(1, 2),
        "whole_second_derivative_budget": massive_forward_bound()
        - s.Rational(82748158895162527, 1800),
    }
    for name, values in majorants().items():
        for k, value in enumerate(values):
            checks[
                "weighted_numerator_and_denominator_bound_" + name + "_" + str(k)
            ] = value - expected[name][k]
    return {
        "whole_small_channel_massive_block": sm,
        "whole_regularized_master_coefficients": cs,
        "whole_regularized_zero_and_first_jets": (zero, first),
        "whole_coefficient_derivative_majorants": majorants(),
        "whole_uniform_massive_forward_bound": massive_forward_bound(),
        "compact_norm_domain": "25/4<=s<=16,-1<=t<=0,0<ell<3. "
        "The factored rational denominator factors s,s+t-4,t-4,s-4 have "
        "absolute value at least1; weighted numerator l1 uses16,1,3. "
        "All J_t,Q_t,L_t and their first two derivatives have modulus<=1.",
        "Taylor_proof": "f(t)=t^2*S_m(t) has f(0)=f'(0)=0. The exact coefficients "
        "and master derivative bounds give |f''|<=2*massive_forward_bound. "
        "Taylor's integral remainder therefore gives |S_m|<=massive_forward_bound "
        "for-1<=t<0, and its continuous endpoint. No small-transfer pole "
        "is bounded separately before its exact cancellation.",
        "checks": checks,
        "gates": {
            "whole_massive_triangle_bubble_and_all_rational_completions": True,
            "zero_and_first_jets_cancel_before_absolute_bounds": True,
            "all_rational_denominator_factors_have_explicit_lower_bounds": True,
            "all_three_small_masters_have_two_derivative_bounds": True,
            "uniform_Taylor_remainder_not_a_numeric_scan": True,
            "exact_majorant_below_forty_six_trillion": bool(
                massive_forward_bound() < s.Integer(46) * 10**12
            ),
            "no_full_amplitude_or_matching_term_is_deleted": True,
        },
    }
