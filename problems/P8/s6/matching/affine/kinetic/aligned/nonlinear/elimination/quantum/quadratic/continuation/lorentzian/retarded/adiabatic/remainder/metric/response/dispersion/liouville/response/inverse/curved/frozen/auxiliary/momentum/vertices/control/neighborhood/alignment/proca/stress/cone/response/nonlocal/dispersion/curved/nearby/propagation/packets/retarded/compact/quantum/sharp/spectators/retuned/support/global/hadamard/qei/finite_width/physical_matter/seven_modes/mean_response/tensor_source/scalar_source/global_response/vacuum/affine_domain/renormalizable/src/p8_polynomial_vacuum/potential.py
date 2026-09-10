"""Positive one-loop constant-background remainder with fixed local subtractions."""

from functools import cache

import sympy as sp

from . import model


@cache
def data():
    d = model.data()
    mu2 = d["heavy_mass_squared"]
    G2 = d["cubic_coupling_squared"]
    lam4 = d["polynomial_quartic"]
    A = sp.factor((lam4 - G2 / mu2) / 2)
    delta = d["positive_quartic_margin"]
    y = sp.Symbol("Euclidean_momentum_squared", nonnegative=True)
    phi = sp.Symbol("constant_canonical_light_field", real=True)
    G = sp.Symbol("positive_cubic_coupling", positive=True)
    F = A - G2 / (y + mu2)
    z = sp.Symbol("positive_log_argument_increment", nonnegative=True)
    remainder = sp.log(1 + z) - z + z * z / 2
    H = sp.Symbol("canonical_heavy_field", real=True)
    V = phi * phi / 2 + mu2 * H * H / 2 + G * H * phi * phi / 2 + lam4 * phi**4 / 24
    Hstar = -G * phi * phi / (2 * mu2)
    full = sp.hessian(V, (phi, H)) + sp.diag(y, y)
    schur = (
        (full[0, 0] - full[0, 1] * full[1, 0] / full[1, 1])
        .subs(H, Hstar)
        .subs(G, sp.sqrt(G2))
    )
    lower = delta**3 * phi**6 / (1536 * sp.pi**2 * (1 + A * phi * phi))
    upper = A**3 * phi**6 / (192 * sp.pi**2)
    radial = sp.integrate(y / (y + 1) ** 3, (y, 0, sp.oo))
    actual = {
        sp.Symbol("positive_fixed_lambda", positive=True): d["actual_parameters"][
            "lambda"
        ],
        sp.Symbol("positive_fixed_gamma", positive=True): d["actual_parameters"][
            "gamma"
        ],
    }
    actualA = sp.factor(A.subs(actual))
    actualdelta = d["actual_parameters"]["positive_completed_square_quartic_margin"]
    # pi^2>9 gives a deliberately rational global ceiling.
    rational_upper = actualA**3 / 1728
    # Retain the momentum dependence F=delta/2+b*y/(y+mu^2).
    # (a+b)^3<=4(a^3+b^3), followed by y^3/(y+1)^3<=1,
    # controls the high-momentum region with the actual heavy mass.
    mass_weighted = ((delta / 2) ** 3 + (G2 / mu2) ** 3 / mu2) / (48 * sp.pi**2)
    actual_weighted = sp.factor(
        (((delta / 2) ** 3 + (G2 / mu2) ** 3 / mu2) / 432).subs(actual)
    )
    a0, b0 = sp.symbols("nonnegative_a nonnegative_b", nonnegative=True)
    mass2 = sp.Symbol("positive_radial_mass_squared", positive=True)
    # A counterterm-free third derivative in r=Phi^2 is positive; the
    # zeroth, first and second Taylor coefficients are chosen only once.
    r = sp.Symbol("field_square", nonnegative=True)
    q = sp.Symbol("positive_kernel_ratio", positive=True)
    t = sp.Symbol("unit_interval_variable", nonnegative=True)
    pi_lower = 4 * sum(sp.Rational((-1) ** j, 2 * j + 1) for j in range(8))
    return {
        "constant_stationary_heavy": Hstar,
        "light_Schur_kernel": y + 1 + F * phi * phi,
        "vertex_kernel_F": F,
        "positive_F_lower": delta / 2,
        "F_upper": A,
        "once_subtracted_log_remainder": remainder,
        "one_loop_integral_convention": "V1_ren=one half integral d^4q/(2pi)^4 [log(1+z)-z+z^2/2], z=F(q^2)Phi^2/(q^2+1), subtract constant, Phi^2 and Phi^4 at the same constant-field vacuum",
        "global_positive_remainder_lower": lower,
        "global_remainder_upper": upper,
        "actual_kernel_upper": actualA,
        "actual_quartic_margin": actualdelta,
        "actual_rational_Phi6_upper_coefficient": rational_upper,
        "mass_weighted_Phi6_upper_coefficient": mass_weighted,
        "actual_mass_weighted_rational_Phi6_upper_coefficient": actual_weighted,
        "elementary_pi_lower": pi_lower,
        "renormalization_scope": "Constant-background effective-potential curvature and quartic renormalization, once for this new model. These are not on-shell pole-residue/mass conditions, a b2 matching choice, a state-dependent subtraction, or the frozen P8 counterterms.",
        "loop_scope": "Exact one-loop constant-background determinant and its convergent Taylor-subtracted remainder. No higher-loop bound, inhomogeneous gradient estimate, spectral/dispersion reconstruction or full V/G/B claim.",
        "checks": {
            "literal_two_field_Hessian_Schur_kernel": sp.factor(
                schur - (y + 1 + F * phi * phi)
            ),
            "kernel_lower_margin": sp.factor(
                F - delta / 2 - G2 * y / (mu2 * (y + mu2))
            ),
            "kernel_upper_gap": sp.factor(A - F - G2 / (y + mu2)),
            "kernel_monotone": sp.factor(sp.diff(F, y) - G2 / (y + mu2) ** 2),
            "subtracted_log_derivative": sp.factor(
                sp.diff(remainder, z) - z * z / (1 + z)
            ),
            "subtracted_log_zero": remainder.subs(z, 0),
            "log_remainder_integral": sp.integrate(z * z / (1 + z), (z, 0, z))
            - remainder,
            "four_dimensional_radial_integral": radial - sp.Rational(1, 2),
            "counterterm_free_third_field_square_derivative": sp.factor(
                sp.diff(sp.log(1 + q * r) / 2, r, 3) - q**3 / (1 + q * r) ** 3
            ),
            "stationary_tree_light_potential": sp.factor(
                V.subs(H, Hstar).subs(G, sp.sqrt(G2))
                - phi * phi / 2
                - delta * phi**4 / 24
            ),
            "rational_pi_upper_normalization": sp.Rational(1, 192 * 9)
            - sp.Rational(1, 1728),
            "finite_geometric_pi_integrand_lower": sp.factor(
                1 / (1 + t * t)
                - sum((-t * t) ** j for j in range(8))
                - t**16 / (1 + t * t)
            ),
            "positive_two_term_cube_majorant": sp.expand(
                4 * (a0**3 + b0**3) - (a0 + b0) ** 3 - 3 * (a0 - b0) ** 2 * (a0 + b0)
            ),
            "mass_weighted_radial_majorant": sp.factor(
                y - y**4 / (y + 1) ** 3 - y * (3 * y**2 + 3 * y + 1) / (y + 1) ** 3
            ),
            "actual_heavy_radial_scale_integral": sp.integrate(
                y / (y + mass2) ** 3, (y, 0, sp.oo)
            )
            - 1 / (2 * mass2),
        },
        "bounds": {
            "actual_kernel_lower_positive": actualdelta > 0,
            "actual_kernel_upper_below_one_e_minus_205": actualA
            < sp.Rational(1, 10**205),
            "actual_positive_remainder_below_one_e_minus_618_times_Phi6": 0
            < rational_upper
            < sp.Rational(1, 10**618),
            "finite_log_kernel_lower_and_upper": actualA > actualdelta / 2 > 0,
            "positive_integral_proves_pi_greater_than_three": pi_lower > 3,
            "actual_mass_weighted_remainder_below_one_e_minus_815_times_Phi6": 0
            < actual_weighted
            < sp.Rational(1, 10**815),
            "mass_weighted_bound_improves_coarse_bound": actual_weighted
            < rational_upper,
        },
    }
