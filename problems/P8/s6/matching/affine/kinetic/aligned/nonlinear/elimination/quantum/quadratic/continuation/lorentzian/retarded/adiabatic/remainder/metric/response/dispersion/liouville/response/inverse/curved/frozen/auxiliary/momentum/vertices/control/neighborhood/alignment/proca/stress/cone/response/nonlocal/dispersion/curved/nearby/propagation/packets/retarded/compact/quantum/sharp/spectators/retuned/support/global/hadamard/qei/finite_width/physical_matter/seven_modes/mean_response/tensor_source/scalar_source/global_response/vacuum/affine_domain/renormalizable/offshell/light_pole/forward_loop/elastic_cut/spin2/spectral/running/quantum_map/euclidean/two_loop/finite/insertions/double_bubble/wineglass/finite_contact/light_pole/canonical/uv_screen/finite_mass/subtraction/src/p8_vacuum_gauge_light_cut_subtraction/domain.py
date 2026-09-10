"""Complex-s extension of the inherited entire finite-mass box remainder."""

from functools import cache

import sympy as sp
from p8_vacuum_finite_mass_gauge_cut import kinematics, series


@cache
def data():
    m, a, Y = sp.symbols(
        "positive_fermion_mass positive_gauge_squared positive_Yukawa_squared",
        positive=True,
    )
    eps, z = sp.symbols("positive_box_epsilon positive_box_scale", positive=True)
    q = 16 * sp.pi**2
    color = 8
    eta = sp.Rational(5, 2) * eps + eps * eps / 8
    rho_error = color * z * z * eta / sp.pi
    checks = {
        "same_six_box_routes": len(kinematics.data()["six_cyclic_box_routes"]) - 6,
        "same_all_degree_spin_order_radial_normalization": series.data()[
            "unrounded_box_remainder_upper"
        ]
        / (a * Y / (q * m**3))
        - 96 * 6 * 12**3,
        "complex_invariant_absolute_upper": 2 + 5 - 7,
        "complex_crossed_invariant_absolute_upper": abs(2 - 4) + 5 - 7,
        "energy_and_spatial_square_margin": sp.Rational(3, 2) ** 2
        - sp.Rational(7, 4)
        - sp.Rational(1, 2),
        "scalar_component_norm_bound": 2 * sp.Rational(3, 2) - 3,
        "gauge_component_norm_bound": 3 * sp.Rational(3, 2) - sp.Rational(9, 2),
        "strict_routing_margin": 18 - 3 * sp.Rational(9, 2) - sp.Rational(9, 2),
        "new_unrounded_box_constant": 96 * 6 * 18**3 - 3359232,
        "rounded_box_constant_strict_margin": 4000000 - 3359232 - 640768,
        "complex_disc_leading_amplitude_margin": 10
        - sp.Rational(4, 3) * 7
        - sp.Rational(2, 3),
        "four_pair_cut_product_bound": 4 * (2 * 10 * eps + eps * eps)
        - (80 * eps + 4 * eps * eps),
        "optical_cut_remainder_dictionary": sp.factor(
            color * (80 * eps + 4 * eps * eps) * z * z / (32 * sp.pi) - rho_error
        ),
    }
    return {
        "complex_s_disc_center": 2,
        "complex_s_disc_radius": 5,
        "absolute_s_and_s_minus_four_upper": 7,
        "absolute_E_and_p_upper": sp.Rational(3, 2),
        "strict_cumulative_component_norm_upper": 18,
        "minimum_fermion_mass": 36,
        "uniform_finite_mass_box_error_upper": 4000000 * a * Y / (q * m**3),
        "box_epsilon": 4000000 / m,
        "normalized_single_channel_rho_error_eta": eta,
        "holomorphic_single_channel_rho_error_upper": rho_error,
        "symmetry_descent": "The complete transition is even in p by scalar Bose symmetry and even under (E,p)->(-E,-p) by Lorentz tensor parity. It is therefore separately even in E and p and descends holomorphically to s through s=0 and s=4.",
        "soft_endpoint": "Analyticity and both off-shell gauge Ward identities supply one soft gauge momentum per leg: B=O(E^2), hence rho=O(s^2) at s=0.",
        "scope": "The same first two-gauge channel and same loop order, now uniformly controlled on the complex disc needed to differentiate the finite-cut subtraction. No new full-amplitude or all-loop claim.",
        "checks": checks,
    }
