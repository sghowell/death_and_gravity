"""Literal two-line insertion count and local-reference forest ownership."""

from functools import cache

import sympy as sp
from p8_vacuum_forward_loop import normalization


@cache
def data():
    h, L, Q = sp.symbols("loop_order L Q")
    D1, D2, P1, P2 = sp.symbols("D1 D2 insertion1 insertion2")
    channel = L * L * (D1 + h * P1) * (D2 + h * P2) / (2 * Q)
    derivative = sp.diff(channel, h).subs(h, 0)
    placements = [
        {
            "channel": ch,
            "inserted_light_line": line,
            "outer_symmetry_weight": sp.Rational(1, 2),
        }
        for ch in ("s", "t", "u")
        for line in (1, 2)
    ]
    z, c0, c1 = sp.symbols("z c0 c1")
    affine = c0 + c1 * z
    paired = affine - affine.subs(z, 1) - (z - 1) * sp.diff(affine, z).subs(z, 1)
    v = sp.symbols("forward_displacement")
    contact = sp.symbols("fixed_outer_Phi4_constant")
    f = sp.Function("whole_fermion_inverse")
    fR = f(z) - f(1) - (z - 1) * sp.diff(f(z), z).subs(z, 1)
    checks = {
        "independent_product_derivative_two_positions": sp.factor(
            derivative - L * L * (P1 * D2 + D1 * P2) / (2 * Q)
        ),
        "two_positions_times_outer_symmetry": 2 * sp.Rational(1, 2) - 1,
        "same_frozen_outer_bubble_symmetry": normalization.data()[
            "three_channel_symmetry_factor"
        ]
        - sp.Rational(1, 2),
        "same_frozen_radial_Q_normalization": normalization.data()[
            "four_dimensional_radial_channel_prefactor"
        ]
        * 2
        * (16 * sp.pi**2)
        - 1,
        "three_channels_two_positions": len(placements) - 6,
        "all_s_channel_weights": sum(
            p["outer_symmetry_weight"] for p in placements if p["channel"] == "s"
        )
        - 1,
        "all_t_channel_weights": sum(
            p["outer_symmetry_weight"] for p in placements if p["channel"] == "t"
        )
        - 1,
        "all_u_channel_weights": sum(
            p["outer_symmetry_weight"] for p in placements if p["channel"] == "u"
        )
        - 1,
        "whole_inner_affine_forest_cancellation": sp.expand(paired),
        "inner_finite_mass_reference_retained": fR.subs(z, 1),
        "inner_finite_residue_reference_retained": sp.diff(fR, z).subs(z, 1),
        "overall_local_contact_cannot_tune_b2": sp.diff(contact, v, 2),
        "parameter_shift_inside_two_loop_family_is_later_order": sp.diff(
            h * h * (L + h * c0) ** 2, h, 2
        ).subs(h, 0)
        / 2
        - L * L,
    }
    return {
        "placements": placements,
        "literal_outer_channel_before_one_insertion": channel,
        "literal_one_insertion_sum": derivative,
        "proper_subgraph": "Pair the complete one-loop fermion self-energy with its mass and field counterterms in the already fixed physical on-shell scheme, obtaining f_R.",
        "overall_subtraction": "After that pairing, subtract the outer Phi^4 constant A_channel(2) before removing the common regulator or integrating the spectral mass. This is a local quartic term; it has no forward b2.",
        "reference_conversion": "A different prescribed finite outer Phi^4 constant changes only the local contact for this literal quartic-vertex family. No independent b2 counterterm is added. The full two-loop MS conversion in other families remains to be computed.",
        "formal_order": "Use leading reference L,Y and mass-one free scalar lines. Exact powers of the finite one-loop kinetic normalization are not resummed into this formal order-two family.",
        "scope": "Six named local-quartic bubble contributions plus their paired inner on-shell references and overall quartic subtraction. Fermion boxes, heavy-exchange vertices, other primitive graphs, wave/coupling conversions of other lower-order amplitudes and later loops are not included.",
        "checks": checks,
    }
