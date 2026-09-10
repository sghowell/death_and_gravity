"""One fixed set of local full-model counterterms, compatible with the potential."""

from functools import cache

import sympy as sp
from p8_polynomial_vacuum import model
from p8_uv import vacuum


@cache
def data():
    lam, g, M, I0 = sp.symbols(
        "quartic_L cubic_squared heavy_mass_squared regulated_zero_bubble", real=True
    )
    h = sp.symbols("h_s h_t h_u", real=True)
    factor = I0 / (32 * sp.pi**2)
    delta_lam = 3 * lam * lam * factor
    delta_g = 2 * lam * g * factor
    delta_M = g * factor
    full_uv = sum((-lam + g * v) ** 2 for v in h) * factor
    ct = -delta_lam + delta_g * sum(h) - g * delta_M * sum(v * v for v in h)
    potential_delta = delta_lam - 3 * delta_g / M + 3 * g * delta_M / M**2
    A = (lam - g / M) / 2
    s, t, u = sp.symbols("generic_s generic_t generic_u", real=True)
    tree = -lam + g * sum(1 / (M - z) for z in (s, t, u))
    differential = (
        delta_lam * sp.diff(tree, lam)
        + delta_g * sp.diff(tree, g)
        + delta_M * sp.diff(tree, M)
    )
    actual = model.data()
    forward = actual["exact_tree_amplitude"].subs(
        {vacuum.transfer: 0, vacuum.w: 4 - vacuum.s}
    )
    fixed_lambda = sp.Symbol("positive_fixed_lambda", positive=True)
    finite = sp.Symbol("finite_potential_quartic_counterterm", real=True)
    return {
        "regulated_zero_momentum_light_bubble": I0,
        "delta_polynomial_quartic_UV": delta_lam,
        "delta_cubic_squared_UV": delta_g,
        "delta_heavy_mass_squared_UV": delta_M,
        "full_UV_counterterm_amplitude": ct,
        "constant_field_effective_quartic_counterterm": potential_delta,
        "finite_potential_contact": finite,
        "renormalization_convention": "Subtract I0 C(z)^2/(32pi^2) in each channel once, with I0 the regulated mass-one zero-external-momentum light bubble. The displayed total full-model interaction counterterms implement it. Add one momentum-independent finite quartic contact to keep the previous constant-potential quartic condition. Light mass and residue use S6.112 on shell; heavy M and G here are fixed renormalized parameters, not a stable heavy on-shell scheme.",
        "wavefunction_bookkeeping": "These are TOTAL vertex counterterms in the chosen canonical renormalized fields. Independent bare coupling relations absorb field counterterms. Unit renormalized light LSZ residue means no additional tree-amplitude rescaling is added; inserting one-loop counterterms inside one-loop integrals would be two-loop order.",
        "checks": {
            "all_three_channel_UV_terms_cancel_with_local_full_model_counterterms": sp.expand(
                full_uv + ct
            ),
            "counterterm_is_literal_tree_parameter_variation": sp.factor(
                differential - ct.subs(dict(zip(h, [1 / (M - z) for z in (s, t, u)])))
            ),
            "potential_quartic_UV_counterterm_is_same_scheme": sp.factor(
                potential_delta - 3 * (lam - g / M) ** 2 * factor
            ),
            "full_Hessian_potential_UV_coefficient_cancels": sp.factor(
                potential_delta / 24 - A * A * I0 / (64 * sp.pi**2)
            ),
            "finite_potential_contact_has_zero_forward_second_derivative": sp.diff(
                -finite, s, 2
            ),
            "actual_full_tree_forward_coefficient_is_four_fixed_lambda": sp.factor(
                sp.diff(forward, vacuum.s, 2).subs(vacuum.s, 2) / 2 - 4 * fixed_lambda
            ),
        },
    }
