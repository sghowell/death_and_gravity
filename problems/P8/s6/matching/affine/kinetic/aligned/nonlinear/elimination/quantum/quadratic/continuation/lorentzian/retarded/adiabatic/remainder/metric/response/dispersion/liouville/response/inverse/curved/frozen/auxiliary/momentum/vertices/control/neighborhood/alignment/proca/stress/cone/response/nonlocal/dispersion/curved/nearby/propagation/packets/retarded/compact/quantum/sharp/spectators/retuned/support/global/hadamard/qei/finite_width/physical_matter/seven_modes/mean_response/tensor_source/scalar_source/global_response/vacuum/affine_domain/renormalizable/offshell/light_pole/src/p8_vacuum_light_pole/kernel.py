"""Actual one-loop light two-point kernel and explicit on-shell subtractions."""

from functools import cache

import sympy as sp
from p8_polynomial_vacuum import model, potential


@cache
def data():
    d = model.data()
    M = d["heavy_mass_squared"]
    G2 = d["cubic_coupling_squared"]
    x = sp.Symbol("unit_Feynman_parameter", real=True)
    s = sp.Symbol("Minkowski_light_momentum_squared", real=True)
    a = x * (1 - x)
    Delta = x * M + 1 - x - a * s
    d1 = x * M + (1 - x) ** 2
    b = a / d1
    z = sp.Symbol("complex_Taylor_increment")
    h = sp.Symbol("positive_real_Taylor_increment", nonnegative=True)
    remainder = -sp.log(1 - z) - z
    p_kernel = h - sp.log(1 + h)
    # A one-loop local affine polynomial includes all momentum-independent
    # tadpoles and mass terms plus any finite kinetic counterterm.
    local_a, local_b = sp.symbols("local_constant local_kinetic", real=True)
    loc = local_a + local_b * s
    qv = sp.Matrix(sp.symbols("loop_q0:4", real=True))
    pv = sp.Matrix(sp.symbols("external_Euclidean_p0:4", real=True))
    A = (qv.T * qv)[0] + 1
    B = ((qv + pv).T * (qv + pv))[0] + M
    shifted = ((qv + x * pv).T * (qv + x * pv))[0] + x * M + 1 - x + a * (pv.T * pv)[0]
    y, L = sp.symbols("positive_radial_y positive_radial_Delta", positive=True)
    radial = sp.integrate(y / (y + L) ** 3, (y, 0, sp.oo))
    momentum_kernel = -sp.log(Delta)
    on_shell = (
        momentum_kernel
        - momentum_kernel.subs(s, 1)
        - (s - 1) * sp.diff(momentum_kernel, s).subs(s, 1)
    )
    # Log branches are NOT algebraically combined globally. Their equal
    # derivatives and anchored value identify the analytic branch locally.
    target = -sp.log(1 - (s - 1) * b) - (s - 1) * b
    Mp = d["actual_parameters"]["heavy_mass_squared"]
    Gp = d["actual_parameters"]["cubic_coupling_squared"]
    residue_upper = Gp / (288 * Mp)
    curvature_upper = Gp / (864 * Mp**2)
    epsilon = Gp / (864 * Mp**2 * (1 - 1 / Mp))
    p = sp.Symbol("positive_once_subtracted_curvature_shift", positive=True)
    Phi = sp.Symbol("constant_canonical_light_field", real=True)
    mu = sp.Symbol("positive_heavy_mass", positive=True)
    return {
        "heavy_mass_squared": M,
        "cubic_squared": G2,
        "actual_heavy_mass_squared": Mp,
        "actual_cubic_squared": Gp,
        "Feynman_parameter": x,
        "Minkowski_invariant": s,
        "Delta": Delta,
        "on_shell_Delta": d1,
        "normalized_parameter_weight": b,
        "self_energy_momentum_kernel": momentum_kernel,
        "self_energy_prefactor": G2 / (16 * sp.pi**2),
        "first_two_particle_bubble_threshold": (sp.sqrt(M) + 1) ** 2,
        "inverse_propagator_convention": "D(s)=s-1+Pi(s); Pi(s)=G^2/(16pi^2)[local UV constant - integral_0^1 log Delta(s,x) dx]. Imaginary Pi is positive above its cut in the Feynman prescription.",
        "on_shell_kernel_before_log_combination": on_shell,
        "local_analytic_on_shell_kernel": target,
        "complex_log_remainder": remainder,
        "positive_curvature_conversion_integrand": p_kernel,
        "actual_rational_finite_kinetic_counterterm_upper": residue_upper,
        "actual_rational_curvature_shift_upper": curvature_upper,
        "actual_unit_disc_quadratic_self_energy_upper": epsilon,
        "on_shell_potential_quadratic": (1 - p) * Phi**2 / 2,
        "unchanged_positive_one_loop_Phi6_upper": potential.data()[
            "actual_mass_weighted_rational_Phi6_upper_coefficient"
        ],
        "scheme": "For this new polynomial vacuum choose Pi_R(s)=Pi(s)-Pi(1)-(s-1)Pi'(1) once; the one-loop inverse has pole mass one and residue one. Potential curvature becomes 1-Pi_R(0), not one. S6.110 fixed potential quartic subtraction is retained to the same order.",
        "scope": "Actual one-loop two-point function, local analytic disc and explicit finite subtraction conversion only. No higher-loop pole/cut control, full four-point matching, old cosmological state/counterterm transfer or V/G/B closure.",
        "checks": {
            "literal_four_dimensional_Feynman_parameter_shift": sp.expand(
                (1 - x) * A + x * B - shifted
            ),
            "actual_on_shell_parameter_denominator": sp.factor(Delta.subs(s, 1) - d1),
            "parameter_weight_upper_gap": sp.factor(
                (1 - x) / M - b - (1 - x) ** 3 / (M * d1)
            ),
            "first_momentum_derivative_of_bubble_kernel": sp.factor(
                sp.diff(momentum_kernel, s) - a / Delta
            ),
            "second_momentum_derivative_of_bubble_kernel": sp.factor(
                sp.diff(momentum_kernel, s, 2) - a * a / Delta**2
            ),
            "local_analytic_log_remainder_derivative": sp.factor(
                sp.diff(on_shell - target, s)
            ),
            "local_analytic_log_remainder_anchor": sp.expand(
                (on_shell - target).subs(s, 1)
            ),
            "arbitrary_local_affine_counterterm_cancels": sp.expand(
                loc - loc.subs(s, 1) - (s - 1) * sp.diff(loc, s)
            ),
            "normalized_log_remainder_zero": remainder.subs(z, 0),
            "normalized_log_remainder_linear_zero": sp.diff(remainder, z).subs(z, 0),
            "normalized_log_remainder_derivative": sp.factor(
                sp.diff(remainder, z) - z / (1 - z)
            ),
            "positive_curvature_kernel_derivative": sp.factor(
                sp.diff(p_kernel, h) - h / (1 + h)
            ),
            "positive_curvature_kernel_anchor": p_kernel.subs(h, 0),
            "finite_radial_bubble_derivative_integral": radial - 1 / (2 * L),
            "parameter_weight_linear_integral": sp.integrate(1 - x, (x, 0, 1))
            - sp.Rational(1, 2),
            "parameter_weight_square_integral": sp.integrate((1 - x) ** 2, (x, 0, 1))
            - sp.Rational(1, 3),
            "on_shell_potential_curvature_is_not_pole_mass": sp.diff(
                (1 - p) * Phi**2 / 2, Phi, 2
            )
            - (1 - p),
            "actual_two_particle_threshold_square": sp.expand(
                x * mu**2 + 1 - x - a * (mu + 1) ** 2 - ((mu + 1) * x - 1) ** 2
            ),
        },
        "bounds": {
            "actual_heavy_mass_positive_and_above_two_light_threshold": Mp > 4,
            "actual_unit_disc_below_parameter_analyticity_radius": Mp > 1,
            "finite_light_kinetic_counterterm_below_three_e_minus_208": 0
            < residue_upper
            < sp.Rational(3, 10**208),
            "positive_curvature_shift_below_one_e_minus_405": 0
            < curvature_upper
            < sp.Rational(1, 10**405),
            "unit_disc_quadratic_remainder_below_one_e_minus_405": 0
            < epsilon
            < sp.Rational(1, 10**405),
            "unit_disc_factorization_excludes_extra_one_loop_poles": epsilon < 1,
            "converted_constant_field_curvature_positive": curvature_upper < 1,
            "finite_local_light_kinetic_coefficient_positive": residue_upper < 1,
            "same_positive_mass_weighted_Phi6_bound_retained": potential.data()[
                "actual_mass_weighted_rational_Phi6_upper_coefficient"
            ]
            < sp.Rational(1, 10**815),
        },
    }
