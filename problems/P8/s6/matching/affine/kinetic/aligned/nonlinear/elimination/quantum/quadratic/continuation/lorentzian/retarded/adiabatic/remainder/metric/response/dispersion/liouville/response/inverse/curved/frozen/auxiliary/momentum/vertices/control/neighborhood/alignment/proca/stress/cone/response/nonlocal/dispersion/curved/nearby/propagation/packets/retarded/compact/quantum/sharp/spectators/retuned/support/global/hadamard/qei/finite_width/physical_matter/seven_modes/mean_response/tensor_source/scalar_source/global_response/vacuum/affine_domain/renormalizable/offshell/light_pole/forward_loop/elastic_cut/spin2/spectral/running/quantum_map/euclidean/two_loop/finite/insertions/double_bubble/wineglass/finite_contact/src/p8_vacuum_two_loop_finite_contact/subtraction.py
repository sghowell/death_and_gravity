"""The linear inherited local subtraction for a fixed finite quartic insertion."""

from functools import cache

import sympy as sp
from p8_vacuum_forward_loop import subtraction as previous


@cache
def data():
    old = previous.data()
    L, g = sp.symbols("quartic_L cubic_squared", real=True)
    sigma = sp.Symbol("once_fixed_finite_contact", real=True)
    hs = sp.symbols("h_s h_t h_u", real=True)
    I0 = old["regulated_zero_momentum_light_bubble"]
    factor = I0 / (32 * sp.pi**2)
    dL, dg, dM = 6 * sigma * L * factor, 2 * sigma * g * factor, sp.Integer(0)
    local = -dL + dg * sum(hs)
    raw_UV = factor * sum((-L + g * h) ** 2 for h in hs)
    s = sp.Symbol("channel_invariant", real=True)
    finite2 = sp.Symbol("fixed_two_loop_potential_contact", real=True)
    return {
        "same_entire_regulated_reference": I0,
        "fixed_finite_contact": sigma,
        "local_reference_continuation_counterterms": {
            "delta_polynomial_quartic": dL,
            "delta_cubic_squared": dg,
            "delta_heavy_mass_squared": dM,
        },
        "total_canonical_counterterm_amplitude_for_this_group": local,
        "scope": "The finite contact is already fixed by the constant-field quartic condition. Its one-loop insertion uses the linear variation of the same entire I0 reference counterterms. These are total canonical-field interaction counterterms, not an additional external tree LSZ factor or a new momentum-dependent contact chosen to tune b2.",
        "checks": {
            "local_counterterms_cancel_variation_of_all_channel_UV": sp.expand(
                local + sigma * sp.diff(raw_UV, L)
            ),
            "quartic_reference_counterterm_is_same_linear_variation": sp.expand(
                dL - sigma * sp.diff(old["delta_polynomial_quartic_UV"], L)
            ),
            "cubic_reference_counterterm_is_same_linear_variation": sp.expand(
                dg - sigma * sp.diff(old["delta_cubic_squared_UV"], L)
            ),
            "heavy_mass_reference_has_no_L_variation": dM
            - sigma * sp.diff(old["delta_heavy_mass_squared_UV"], L),
            "full_local_counterterm_is_same_linear_variation": sp.expand(
                local - sigma * sp.diff(old["full_UV_counterterm_amplitude"], L)
            ),
            "new_order_constant_potential_contact_cannot_tune_b2": sp.diff(
                finite2, s, 2
            ),
        },
    }
