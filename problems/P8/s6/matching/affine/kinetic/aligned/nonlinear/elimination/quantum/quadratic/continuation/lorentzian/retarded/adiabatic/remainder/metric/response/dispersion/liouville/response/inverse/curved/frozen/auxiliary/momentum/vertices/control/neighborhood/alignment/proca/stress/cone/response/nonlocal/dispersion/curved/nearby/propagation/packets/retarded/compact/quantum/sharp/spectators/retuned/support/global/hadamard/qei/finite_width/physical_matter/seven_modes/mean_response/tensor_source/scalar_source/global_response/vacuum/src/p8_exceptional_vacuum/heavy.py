"""Exact massive tree exchange matching, not a bounce parent or UV completion."""

from functools import cache

import sympy as sp
from p8_uv import vacuum as original

from . import analytic, vacuum

CHANNEL_RADIUS = sp.Integer(10) ** 144


@cache
def data():
    s, t, w = original.s, original.transfer, original.w
    lam, gamma = sp.symbols("positive_fixed_lambda positive_fixed_gamma", positive=True)
    D = 2 * lam / gamma
    mu2 = 2 + D
    g2 = 8 * lam * D
    xs = [z - 2 for z in (s, t, w)]
    exact = g2 * sum(x * x / (D - x) for x in xs) / 4
    first_two = 2 * lam * sum(x * x for x in xs) + gamma * sum(x**3 for x in xs)
    remainder = g2 * sum(x**4 / (1 - x / D) for x in xs) / (4 * D**3)
    # Reconstruct from the actual independent local jets and labelled
    # amplitude, then co-scale the common n/kappa before any numeric choice.
    v = vacuum.data()
    d = analytic.local_jets()
    n, kap = (
        d["canonical_switch_order_symbol"],
        d["canonical_action_normalization_symbol"],
    )
    actual = v["analytic_candidate_scalar_contact"].subs(
        {kap: n / gamma, d["canonical_coupling_symbol"]: lam, original.mass: 1},
        simultaneous=True,
    )
    sub = {lam: analytic.VACUUM_LAMBDA_BAR, gamma: analytic.FIXED_GAMMA}
    actual_D = D.subs(sub)
    actual_g2 = g2.subs(sub)
    actual_mu2 = mu2.subs(sub)
    radius = CHANNEL_RADIUS
    remainder_upper = (3 * g2 * radius**4 / (4 * D**3 * (1 - radius / D))).subs(sub)
    delta = sp.Symbol("forward_crossing_variable", real=True)
    forward = exact.subs({t: 0, s: 2 + delta, w: 2 - delta}, simultaneous=True)
    exact_b2 = sp.factor(sp.diff(forward, delta, 2).subs(delta, 0) / 2)
    x = sp.Symbol("single_channel_crossing_shift", real=True)
    residue = sp.factor((g2 * x * x / 4).subs(x, D))
    return {
        "fixed_light_mass": sp.Integer(1),
        "canonical_heavy_mass_squared": mu2,
        "canonical_derivative_coupling_squared": g2,
        "channel_expansion_gap": D,
        "exact_crossing_symmetric_tree_exchange": exact,
        "first_two_crossing_channel_terms": first_two,
        "actual_coscaled_light_contact": sp.factor(actual),
        "exact_tree_omitted_terms": remainder,
        "positive_heavy_s_channel_residue": residue,
        "exact_heavy_forward_b2": exact_b2,
        "actual_parameters": {
            "lambda": sub[lam],
            "gamma": sub[gamma],
            "D": actual_D,
            "mu_squared": actual_mu2,
            "g_squared": actual_g2,
        },
        "complex_channel_polydisc_radius": radius,
        "exact_uniform_tree_remainder_upper": remainder_upper,
        "tree_vacuum_parent": "L=(dH)^2/2 + exp(gH)*(dPhi)^2/2 - mu^2 H^2/2 - Phi^2/2",
        "matching_boundary": "On-shell massive four-light-particle tree amplitude through cubic Mandelstam order; not off-shell action, loops, gravity, Regge, or common bounce matching",
        "checks": {
            "actual_contact_equals_both_heavy_expansion_terms_on_shell": sp.factor(
                (actual - first_two).subs(w, 4 - s - t)
            ),
            "exact_channel_geometric_remainder": sp.factor(
                g2 * x * x / (4 * (D - x))
                - g2 * x * x / (4 * D)
                - g2 * x**3 / (4 * D**2)
                - g2 * x**4 / (4 * D**3 * (1 - x / D))
            ),
            "crossing_s_t": sp.factor(exact - exact.xreplace({s: t, t: s})),
            "crossing_s_w": sp.factor(exact - exact.xreplace({s: w, w: s})),
            "heavy_residue_is_positive_coupling_square": sp.factor(
                residue - g2 * D**2 / 4
            ),
            "exact_heavy_forward_b2_equals_actual_contact": sp.factor(
                exact_b2 - 4 * lam
            ),
            "massive_constant_in_channel_cube_sum_not_dropped": sp.factor(
                (sum(x**3 for x in xs) - 3 * s * t * w + 8).subs(w, 4 - s - t)
            ),
            "matched_cubic_channel_coefficient": sp.factor(g2 / (4 * D**2) - gamma),
            "matched_quadratic_channel_coefficient": sp.factor(g2 / (4 * D) - 2 * lam),
        },
        "bounds": {
            "heavy_mass_above_two_light_particle_threshold": actual_mu2 > 4,
            "channel_radius_strictly_inside_heavy_gap": radius < actual_D,
            "channel_radius_over_gap_below_one_e_minus_53": radius / actual_D
            < sp.Rational(1, 10**53),
            "exact_tree_remainder_below_one_e_minus_400": remainder_upper
            < sp.Rational(1, 10**400),
            "heavy_squared_mass_above_channel_radius_by_one_e_50": actual_mu2
            > 10**50 * radius,
            "heavy_mass_below_derivative_interaction_scale": actual_mu2 * actual_g2 < 1,
            "heavy_s_channel_residue_positive": bool(residue.is_positive),
            "fixed_derivative_quartic_positive": sub[lam] > 0,
            "fixed_DHOST_interaction_nonzero": sub[gamma] > 0,
        },
    }
