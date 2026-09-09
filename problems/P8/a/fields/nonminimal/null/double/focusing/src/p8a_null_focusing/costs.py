"""Exact finite spatial-width costs and actual null index margin."""

from functools import cache

import sympy as sp

R = sp.Rational(1, 1000)
L = sp.Integer(2)
SPATIAL_FIRST = sp.Integer(3)
PAST = (sp.Rational(9, 10), sp.Integer(2), (3, 16, 1000, 2000000))
FUTURE = (sp.Rational(1, 2), sp.Integer(2), (50, 10000, 2000000, 400000000))


def norms(length):
    return (sp.Rational(13, 35) * length, sp.Rational(6, 5) / length, 12 / length**3)


def interval_costs(length, geometry):
    lo, hi, derivatives = geometry
    d1, d2, d3, d4 = map(sp.Integer, derivatives)
    n0, n1, n2 = norms(length)
    P = 3 * (hi * hi * n2 + 4 * d1 * d1 * n1 + 4 * (d1 * d1 / lo + d2) ** 2 * n0)
    M = 2 * (hi * hi * n1 + 4 * d1 * d1 * n0) / lo**4
    quantum = sp.Rational(5, 9) * P + SPATIAL_FIRST * M / 3
    state = 2 * n1 + 8 * d1 * d1 * n0 / lo**2 + SPATIAL_FIRST * n0 / lo**4
    ref0 = 4 * d1 * d1 * d2 / lo * n0 / 360
    refbeta = (
        (48 * d1 * d1 * d2 / lo + 84 * d2 * d2 + 72 * d1 * d3 + 12 * hi * d4) * n0 / 360
    )
    return {
        "P_norm_majorant": P,
        "mixed_norm_majorant": M,
        "quantum_cost": quantum,
        "state_cost": state,
        "reference_zero_cost": ref0,
        "reference_beta_cost": refbeta,
        "sampler_mass": n0,
    }


@cache
def data():
    past = interval_costs(R, PAST)
    future = interval_costs(L, FUTURE)
    quantum = sp.factor(past["quantum_cost"] + future["quantum_cost"])
    ref0 = sp.factor(past["reference_zero_cost"] + future["reference_zero_cost"])
    refbeta = sp.factor(past["reference_beta_cost"] + future["reference_beta_cost"])
    c0 = quantum + ref0
    cbeta = refbeta
    state = sp.factor(past["state_cost"] + future["state_cost"])
    mass = sp.factor(past["sampler_mass"] + future["sampler_mass"])
    past_ricci = 2 * PAST[2][1] / PAST[0] * past["sampler_mass"]
    kinetic = 2 * norms(L)[1]
    delta = sp.Rational(1, 10**12)
    zeta = sp.Rational(1, 10**6)
    sigma = sp.Rational(1, 10)
    upper = kinetic + past_ricci + delta * max(c0, cbeta) + zeta * state + sigma * mass
    x = sp.Symbol("x", real=True)
    g = 1 - 3 * (x / L) ** 2 + 2 * (x / L) ** 3
    chi = 3 * ((x + R) / R) ** 2 - 2 * ((x + R) / R) ** 3
    checks = {}
    for name, poly, left, right, length in (
        ("past", chi, -R, 0, R),
        ("future", g, 0, L, L),
    ):
        for j, n in enumerate(norms(length)):
            checks[f"actual_{name}_sampler_derivative_norm_{j}"] = (
                sp.integrate(sp.diff(poly, x, j) ** 2, (x, left, right)) - n
            )
    checks.update(
        {
            "actual_sampler_value_joins_at_one": chi.subs(x, 0) - g.subs(x, 0),
            "actual_sampler_first_derivative_joins_at_zero": sp.diff(chi, x).subs(x, 0)
            - sp.diff(g, x).subs(x, 0),
            "actual_sampler_past_outer_value_and_first_traces": chi.subs(x, -R) ** 2
            + sp.diff(chi, x).subs(x, -R) ** 2,
            "actual_sampler_future_outer_value_and_first_traces": g.subs(x, L) ** 2
            + sp.diff(g, x).subs(x, L) ** 2,
            "future_index_kinetic_cost": kinetic - sp.Rational(6, 5),
        }
    )
    return {
        "past": past,
        "future": future,
        "quantum_cost": quantum,
        "C0": c0,
        "Cbeta": cbeta,
        "Cstate": state,
        "Csource": mass,
        "past_Ricci_upper_cost": past_ricci,
        "future_index_kinetic_cost": kinetic,
        "worst_gate_index_upper": upper,
        "strict_index_margin": 2 - upper,
        "checks": checks,
    }
