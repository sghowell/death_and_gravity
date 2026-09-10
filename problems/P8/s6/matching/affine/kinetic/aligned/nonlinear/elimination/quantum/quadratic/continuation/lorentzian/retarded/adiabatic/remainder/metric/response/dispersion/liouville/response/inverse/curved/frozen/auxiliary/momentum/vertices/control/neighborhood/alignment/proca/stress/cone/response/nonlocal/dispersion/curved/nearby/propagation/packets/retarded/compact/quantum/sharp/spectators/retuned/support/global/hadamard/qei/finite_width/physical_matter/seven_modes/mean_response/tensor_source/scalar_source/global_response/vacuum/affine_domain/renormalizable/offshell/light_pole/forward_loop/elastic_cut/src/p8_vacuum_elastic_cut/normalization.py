"""Independent identical-particle phase space and bubble-cut normalization."""

from functools import cache

import sympy as sp


@cache
def data():
    beta, p, E, s = sp.symbols(
        "positive_beta momentum_magnitude particle_energy physical_s", positive=True
    )
    x = sp.Symbol("unit_Feynman_parameter", real=True)
    C = sp.Symbol("constant_tree_vertex", real=True)
    radial_delta = p * p / E**2 / (2 * p / E)
    full_phase_solid_angle = beta / (32 * sp.pi**2)
    full_phase_cosine = 2 * sp.pi * full_phase_solid_angle
    cut_cosine = full_phase_cosine / 4
    roots = [(1 - beta) / 2, (1 + beta) / 2]
    kernel = 1 - s * x * (1 - x)
    t = sp.Symbol("unit_arctangent_parameter", nonnegative=True)
    return {
        "two_body_phase_space_per_solid_angle": full_phase_solid_angle,
        "two_body_phase_space_per_cosine": full_phase_cosine,
        "one_loop_imaginary_part_per_cosine_squared_vertex": cut_cosine,
        "even_angular_average_imaginary_prefactor": beta / (32 * sp.pi),
        "light_bubble_cut_parameter_endpoints": roots,
        "normalization_convention": "S=1+iT; 2 Im A=(1/2! identical intermediate particles) integral dPhi2 |A_tree|^2. The extra factor 1/2 converting 2 Im to Im is retained.",
        "checks": {
            "radial_energy_delta_Jacobian": sp.factor(radial_delta - p / (2 * E)),
            "radial_phase_space_beta_conversion": sp.factor(
                radial_delta.subs(p, beta * E) - beta / 2
            ),
            "literal_four_dimensional_phase_space_normalization": sp.factor(
                (beta / 2) / (4 * (2 * sp.pi) ** 2) - full_phase_solid_angle
            ),
            "azimuthal_phase_integration": sp.factor(
                full_phase_cosine - beta / (16 * sp.pi)
            ),
            "identical_intermediate_and_two_Im_factors": sp.factor(
                cut_cosine - beta / (64 * sp.pi)
            ),
            "independent_constant_vertex_bubble_cut": sp.factor(
                C * C * sp.pi * beta / (32 * sp.pi**2) - 2 * C * C * cut_cosine
            ),
            "lower_Feynman_cut_endpoint": sp.factor(
                kernel.subs({s: 4 / (1 - beta * beta), x: roots[0]})
            ),
            "upper_Feynman_cut_endpoint": sp.factor(
                kernel.subs({s: 4 / (1 - beta * beta), x: roots[1]})
            ),
            "bubble_cut_parameter_width": sp.expand(roots[1] - roots[0] - beta),
            "arctangent_primitive_derivative": sp.factor(
                sp.diff(sp.atan(t), t) - 1 / (1 + t * t)
            ),
            "arctangent_endpoint_pi_over_four": sp.atan(1) - sp.pi / 4,
            "strict_pi_upper_integrand_gap": sp.factor(
                1 - 1 / (1 + t * t) - t * t / (1 + t * t)
            ),
        },
    }
