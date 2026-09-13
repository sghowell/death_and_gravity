"""Complete new-parent matching margins and explicit finite clock extensions."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_scalar_four_point_loop import symmetric
from p8_vacuum_affine_heavy_scalar_loop_coefficients import (
    coefficients as old_coefficients,
)
from p8_vacuum_affine_heavy_scalar_loop_remainder import bounds as old_angle
from p8_vacuum_affine_heavy_scalar_parent import heavy
from p8_vacuum_affine_heavy_scalar_tree_matching import model

from . import germs, loops

S_MAX = s.Integer(10) ** 196
EXTENSION_RADIUS = s.Rational(1, 64)


@cache
def data():
    ka, n, g, la, ga = germs.KAPPA, germs.MASS2, germs.G, germs.LAMBDA, germs.GAMMA
    D = n - 2
    _, (k0, k2, k3) = loops.local_polynomial()
    q = heavy.data()["actual_quartic_q"]
    beta = heavy.data()["actual_coercive_quartic_beta"]
    s0 = s.Rational(4, 3)
    mixed = 32 * g * g * 1386 / ka
    vertex = 64 * g * g / (ka * (n - s0))
    local = abs(k0) + s.Rational(16, 3) * abs(k2) + s.Rational(64, 27) * abs(k3)
    extra_value = (mixed + vertex + local) / 144
    poly_upper = symmetric.ABS_UPPER
    poly_lower = 5985 * g**4 / (640 * D * D)
    full_contact_upper = poly_upper + extra_value
    full_contact_lower = poly_lower - extra_value
    eps_vertex = 8 * (4 - n) / (16 * s.pi**2 * ka)
    eps_upper = 8 * (n - 4) / (144 * ka)
    local_angle = (16 * abs(k2) + 4 * abs(k3) * S_MAX) / (144 * la)
    extra_angle = 2 * eps_upper + local_angle
    old_angle_upper = old_angle.data()["actual_rational_relative_upper"]
    extra_coeff = {
        "b20": eps_upper + abs(k2) / (288 * la),
        "b21": eps_upper + abs(k3) / (432 * ga),
        "b40": eps_upper,
    }
    old = old_coefficients.data()["actual_relative_rational_bounds"]
    old20 = old["delta_b20_over_4lambda"]
    old21 = old["delta_b21_over_minus3gamma"]
    old40 = old["abs_delta_b40_over_gamma2_over_lambda"]
    old20_lower = (1000 * g**4 / (160 * n**4)) / (4 * la)
    old21_lower = (1000 * g**4 / (160 * n**5)) / (3 * ga)
    tube = 2 * s.Integer(10) ** 403 * s.Rational(9, 55) ** germs.N
    u, X = germs.u, germs.X
    V = 1 - germs.parent.original.data()["T"]
    dc, dm, dz, j = s.symbols("deltaC_finite delta_mass_squared deltaZ j1", real=True)
    finiteF = V * (dc * ka * u**4 / 24 - dm * u * u / 2 + dz * X / 2)
    finiteJ = -j * V / s.sqrt(ka)
    P, P1 = s.symbols("Pi_MS_at_1 Pi_MS_prime_at_1")
    spoint = s.Symbol("s")
    checks = {
        "same_old_complete_symmetric_upper": s.cancel(
            poly_upper - 14793 * model.G2**2 / (192 * model.GAP**2)
        ),
        "same_old_angular_bound": old_angle_upper - 10**7 * g * g / (9 * n),
        "same_old_b20_bound": old20 - 1500 * g * g / (288 * n),
        "same_old_b21_bound": old21 - 1500 * g * g / (432 * n),
        "same_old_b40_bound": old40 - 8 * 10**12 * g * g / (9 * n),
        "new_mass_and_residue_counterterms_complete": s.expand(
            (-P1) * spoint - (P - P1) - (-P - (spoint - 1) * P1)
        ),
        "new_first_b20_shift": s.cancel(
            (eps_vertex * 4 * la - 2 * k2 / (16 * s.pi**2)) / (4 * la)
            - eps_vertex
            + k2 / (32 * s.pi**2 * la)
        ),
        "new_first_b21_shift": s.cancel(
            (eps_vertex * (-3 * ga) + k3 / (16 * s.pi**2)) / (-3 * ga)
            - eps_vertex
            + k3 / (48 * s.pi**2 * ga)
        ),
        "new_first_b40_shift_retained": s.cancel(
            eps_vertex * (ga * ga / la) / (ga * ga / la) - eps_vertex
        ),
        "finite_value_lower_retains_extra_sign_uncertainty": full_contact_lower
        - (poly_lower - extra_value),
        "outer_counterterm_switch_ratio": s.Rational(9, 64) / s.Rational(55, 64)
        - s.Rational(9, 55),
        "all_four_jet_Cauchy_factor": s.factorial(4) * 64**4 - 402653184,
        "normalized_contact_coefficient_scaling": dc * ka * u**4 / 24
        - (dc * (s.sqrt(ka) * u) ** 4 / 24) / ka,
        "normalized_heavy_onepoint_source_scaling": -j / s.sqrt(ka)
        - (-j * s.sqrt(ka)) / ka,
        "complete_contact_only_coercive_margin": beta
        - full_contact_upper / 24
        - q / 2
        - (beta - q / 2 - full_contact_upper / 24),
    }
    # Positive real parameter denominator for the complete mixed bubble at p²1.
    x, nn = s.symbols("x n", real=True)
    F = (1 - x) ** 2 + nn * x
    checks["mixed_B_denominator_above_one"] = s.expand(F - 1 - x * (nn - 2 + x))
    checks["mixed_B_denominator_below_n"] = s.expand(nn - F - (1 - x) * (nn - 1 + x))
    return {
        "new_prescription_name": "H8A420-VAC-OS4",
        "full_new_value_condition": "deltaC_new=-[A1_polynomial_base+A1_new_extra](4/3,4/3,4/3), with the same first-order light mass/residue and zero-heavy-onepoint conditions. This is a new full limiting-action prescription, not borrowing the old OS4 contact unchanged.",
        "complete_mixed_integral_bound": "For0<=x<=1 and n>2,1<=(1-x)^2+n*x<=n. The independently shifted complete mixed bracket is integral2(x+1)log F dx, so it lies between0 and3 log n<1386. Integration by parts recovers the full A0/B0 expression, without a heavy-endpoint expansion or a cancellation-sensitive numerical subtraction.",
        "new_symmetric_value_absolute_upper": extra_value,
        "new_contact_positive_lower": full_contact_lower,
        "new_contact_absolute_upper": full_contact_upper,
        "new_contact_only_global_potential": "The new positive finite contact is below10^-408. At Y0 the chosen switch is1, and the S238 completed-square quartic beta-deltaC_new/24 remains above q/2. The contact-only classical comparison therefore retains V>=H²/6+q Phi4/2. This is not the quantum effective potential or the bare vacuum after all finite onepoint and quadratic counterterms.",
        "new_heavy_vertex_factor": eps_vertex,
        "new_extra_angular_relative_upper": extra_angle,
        "combined_full_first_loop_angular_relative_upper": old_angle_upper
        + extra_angle,
        "physical_window": "Every4<=s=S<=10^196 and physical angle -1<=x<=1, with the same original complete massive tree and all old loop parts. The new extra contribution is below10^-600 relative to the original tree; the combined complete first loop stays below10^-199. Tree plus this loop stays within1/59, explicitly as a truncation.",
        "full_angle_bound_proof": "The centered heavy exchange is nonnegative by its exact second-order divided remainder around s0=4/3 and the sum of centered channels0. It is at most the complete separate tree, bounded by61/60 of the original tree. The looser factor2 is used. Also |sigma2-16/3|<=4S² and |sigma3-64/27|<=S³, while A_original>=lambda S²/4. These bounds retain the full rational exchange, not its mass expansion.",
        "extra_relative_coefficient_bounds": extra_coeff,
        "combined_first_loop_coefficient_bounds": {
            "b20": old20 + extra_coeff["b20"],
            "b21": old21 + extra_coeff["b21"],
            "b40": old40 + extra_coeff["b40"],
        },
        "complete_coefficient_signs": "The extra relative changes are below10^-600. The old strictly positive b20 and negative b21 intervals dominate them, so0<delta_b20/(4lambda)<10^-203 and0<delta_b21/(-3gamma)<10^-203 for the complete new first loop. The b40 shift includes the nonzero heavy-vertex correction; its sign remains unassigned and its magnitude is below10^-192 of gamma²/lambda.",
        "finite_counterterm_extension": "Use V(X)=1-T(X) for every stated FINITE counterterm: V[deltaC_new Phi4/24-delta_m² Phi²/2+deltaZ Y/2-j1 H]. Here deltaZ=-Pi_MS'(1),delta_m²=Pi_MS(1)-Pi_MS'(1),j1=-g/(32pi²) is the finite mass1 coincident contribution in the source convention. The relevant vacuum jets are unchanged, and clock jets through order1023 vanish. This defines a finite off-vacuum choice, not the full UV functional or a quantum bounce.",
        "normalized_finite_scalar_extension": finiteF,
        "normalized_finite_heavy_source_extension": finiteJ,
        "finite_extension_complex_radius": EXTENSION_RADIUS,
        "finite_extension_complete_four_jet_upper": tube,
        "finite_extension_tube_bound": "On real |u|<=1,7/8<=X<=9/8, a joint complex radius1/64 gives |(1-X)/X|<=9/55 and |V|<2(9/55)^1024. All finite normalized polynomial prefactors are below10^394, using |deltaC_new|<10^-408, |delta_m²|<2*10^-7, |deltaZ|<10^-207 and |j1|<g/288. Every mixed four-jet factor is below10^9; hence the complete finite F and j extensions are below10^-400. This is a NEW finite-counterterm budget, not S238's classical10^-2700 bound or a bound on the entire loop effective action.",
        "not_inferred": "No omitted-loop control, regulator-removed quantum gravitational decoupling, full curved counterterm functional/state/response, nonlinear same-state bounce, physical UV or finite-gravity Regge theorem follows.",
        "checks": {key: s.cancel(value) for key, value in checks.items()},
        "gates": {
            "actual_mixed_integral_parameter_domain": n > 4,
            "new_local_sigma2_coefficient_positive": k2 > 0,
            "new_local_sigma3_coefficient_positive": k3 > 0,
            "extra_value_smaller_than_old_positive_contact_lower": extra_value
            < poly_lower,
            "new_extra_value_fits_quartic_margin": extra_value / (24 * q)
            < s.Rational(1, 10**200),
            "new_full_positive_contact_below_one_e_minus_408": 0
            < full_contact_lower
            < full_contact_upper
            < s.Rational(1, 10**408),
            "new_contact_only_coercive_quartic_above_half_q": beta
            - full_contact_upper / 24
            > q / 2,
            "extra_local_physical_shape_below_one_e_minus_795": local_angle
            < s.Rational(1, 10**795),
            "new_extra_full_physical_shape_below_one_e_minus_600": extra_angle
            < s.Rational(1, 10**600),
            "combined_full_first_loop_below_one_e_minus_199": old_angle_upper
            + extra_angle
            < s.Rational(1, 10**199),
            "complete_tree_plus_new_first_loop_below_one_over_59": s.Rational(1, 60)
            + old_angle_upper
            + extra_angle
            < s.Rational(1, 59),
            "new_b20_extra_below_one_e_minus_600": extra_coeff["b20"]
            < s.Rational(1, 10**600),
            "new_b21_extra_below_one_e_minus_600": extra_coeff["b21"]
            < s.Rational(1, 10**600),
            "new_b40_extra_below_one_e_minus_600": extra_coeff["b40"]
            < s.Rational(1, 10**600),
            "complete_new_b20_remains_positive": old20_lower > extra_coeff["b20"],
            "complete_new_b21_remains_negative": old21_lower > extra_coeff["b21"],
            "complete_new_b20_relative_error": old20 + extra_coeff["b20"]
            < s.Rational(1, 10**203),
            "complete_new_b21_relative_error": old21 + extra_coeff["b21"]
            < s.Rational(1, 10**203),
            "complete_new_b40_relative_error": old40 + extra_coeff["b40"]
            < s.Rational(1, 10**192),
            "finite_extension_denominator_modulus_margin": s.Rational(9, 55) ** germs.N
            < s.Rational(1, 2),
            "finite_normalized_prefactor_majorant": (
                s.Rational(1, 10**408) * ka * s.Rational(65, 64) ** 4 / 24
                + s.Rational(2, 10**7) * s.Rational(65, 64) ** 2 / 2
                + s.Rational(1, 10**207) * s.Rational(73, 64) / 2
            )
            < 10**394,
            "finite_heavy_source_prefactor_majorant": g / (288 * s.sqrt(ka)) < 10**394,
            "all_mixed_four_jet_factors_bounded": s.factorial(4) * 64**4 < 10**9,
            "complete_finite_extension_four_jet_budget": tube < s.Rational(1, 10**400),
        },
    }
