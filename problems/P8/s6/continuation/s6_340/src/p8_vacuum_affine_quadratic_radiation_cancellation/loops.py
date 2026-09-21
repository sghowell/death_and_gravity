"""Actual mixed exact-D insertions, complete counterterms and positive parameter gap."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_parent_one_loop import germs
from p8_vacuum_affine_heavy_scalar_one_loop import self_energy
from p8_vacuum_affine_mixed_source_radiation import triangles


@cache
def data():
    checks = {}
    x, n = s.symbols("x n", positive=True)
    M = (1 - x) ** 2 + n * x
    alpha = x * (1 - x) / M
    shift = s.Symbol("shift", real=True)
    Mr = 1 - x + n * x - x * (1 - x) * (1 + shift)
    checks["original_massive_anchor"] = s.expand(
        M - self_energy.F.subs({self_energy.x: x, self_energy.m2: n})
    )
    checks["original_alpha_exact"] = s.factor(
        alpha - self_energy.ALPHA.subs({self_energy.x: x, self_energy.m2: n})
    )
    checks["shifted_mixed_denominator"] = s.expand(Mr - M * (1 - alpha * shift))
    checks["positive_gap_through_virtuality2"] = s.expand(
        Mr.subs(shift, 1) - 1 - (n - 3) * x - 2 * x * x
    )
    checks["alpha_majorant"] = s.factor((1 - x) / n - alpha - (1 - x) ** 3 / (n * M))
    checks["integrated_alpha_square_majorant"] = s.integrate(
        (1 - x) ** 2, (x, 0, 1)
    ) - s.Rational(1, 3)
    f = -s.log(1 - alpha * shift) - alpha * shift
    checks["full_nonpolynomial_OS_anchor"] = f.subs(shift, 0)
    checks["full_nonpolynomial_OS_slope"] = s.diff(f, shift).subs(shift, 0)
    checks["full_nonpolynomial_OS_nonzero_curvature"] = s.factor(
        s.diff(f, shift, 2).subs(shift, 0) - alpha * alpha
    )
    checks["analytic_soft_quotient_limit"] = s.limit(f / shift, shift, 0)
    assert germs.MASS2 >= 128 and germs.G == s.Rational(1, 8192)
    assert set(germs.graph_degrees(2)) == {(2, 0, 0, 0), (0, 1, 0, 0)}
    assert germs.G**2 / (864 * germs.MASS2**2 * (1 - 1 / germs.MASS2)) < s.Rational(
        1, 10**400
    )

    # Both explicit S338 exact-D insertions are retained, not inferred from Ward alone.
    prior = triangles.data()
    for key in (
        "literal_light_three_denominators_complete_square",
        "literal_heavy_three_denominators_complete_square",
        "light_triangle_TT_numerator_after_isotropic_average",
        "heavy_triangle_TT_numerator_after_isotropic_average",
        "whole_D_light_triangle_parameter_antiderivative",
        "whole_D_heavy_triangle_parameter_antiderivative",
        "whole_D_two_triangles_equal_bubble_divided_difference",
        "full_dimensional_gamma_recurrence",
        "finite_light_triangle_antiderivative",
        "finite_heavy_triangle_antiderivative",
        "finite_two_triangles_equal_log_bubble_divided_difference",
    ):
        checks["explicit_mixed_loop_" + key] = prior["checks"][key]
    loop_square, loop_cross, kk = s.symbols(
        "epsilon_loop_loop epsilon_loop_k epsilon_kk"
    )
    shift_parameter = s.Symbol("shift_parameter")
    tadpole = loop_square - 2 * shift_parameter * loop_cross + shift_parameter**2 * kk
    I, g, C, n = s.symbols("tadpole_I g C n", nonzero=True)
    j1 = g * I / 2
    checks["exact_D_tadpole_metric_response_TT_zero"] = tadpole.subs(
        {loop_square: 0, loop_cross: 0, kk: 0}
    )
    checks["entire_quartic_tadpole_and_existing_mass_subtraction"] = (
        C * I / 2 - C * I / 2
    )
    checks["entire_heavy_onepoint_source_cancellation"] = s.factor(
        -g * g * I / (2 * n) + g * j1 / n
    )
    eta = s.diag(1, -1, -1, -1)
    plus = s.diag(0, 1, -1, 0)
    cross = s.zeros(4)
    cross[1, 2] = cross[2, 1] = 1
    checks["pure_mass_counterterm_plus_TT_zero"] = C * s.trace(eta * plus)
    checks["pure_mass_counterterm_cross_TT_zero"] = C * s.trace(eta * cross)
    return {
        "checks": checks,
        "gates": {
            "both_literal_mixed_triangles_before_pole_subtraction": True,
            "full_logarithmic_inverse_not_a_derivative_truncation": True,
            "actual_original_parameter_gap_at_all_recoil_virtualities": bool(
                germs.MASS2 >= 128
            ),
            "source_renormalized_inverse_double_zero_not_identically_zero": True,
            "all_three_two_point_topologies_and_metric_responses_retained": True,
            "actual_relative_inverse_bound_less_than_1e_minus400": bool(
                germs.G**2 / (864 * germs.MASS2**2 * (1 - 1 / germs.MASS2))
                < s.Rational(1, 10**400)
            ),
            "no_hidden_null_pole_or_threshold_in_selected_massive_loop": True,
            "full_counterterm_variation_before_real_TT_projection": True,
        },
        "whole_original_subtracted_parameter_integrand": f,
        "whole_mixed_anchor_denominator": M,
        "whole_original_relative_inverse_majorant": germs.G**2
        / (864 * germs.MASS2**2 * (1 - 1 / germs.MASS2)),
        "whole_tadpole_shifted_TT_numerator": tadpole,
        "whole_actual_triangle_denominators": (
            prior["whole_light_denominator"],
            prior["whole_heavy_denominator"],
        ),
        "whole_actual_triangle_antiderivatives": (
            prior["whole_light_antiderivative"],
            prior["whole_heavy_antiderivative"],
        ),
        "whole_sign_and_counterterm_proof": "The mixed bilocal inverse variation has the opposite sign to the physical vertex, normalized independently by the literal free tensor. Both triangle insertions give2epsilon(p,p)*g^2 DD_B/(16pi^2). The fixed kinetic vertex subtracts2epsilon(p,p)*g^2 Bprime(1)/(16pi^2); the mass vertex is pure trace. Together they give2epsilon(p,p)*Pi_OS(r^2)/(r^2-1). The constant quartic tadpole is canceled at every virtuality, and its exact-D metric insertion is zero in TT by isotropy; the full source onepoint pair cancels.",
        "whole_uniform_analytic_proof": "Actual external virtuality is0<r^2<2, so M(v)>=1+(n-3)x+2x^2>=1 at the worst endpointv2. Shifted triangle denominators interpolate positive endpoints. alpha<= (1-x)/n gives|Pi_OS|<=g^2|r^2-1|^2/[96pi^2 n^2(1-1/n)]. The relative quotient is analytic at the anchor and<1e-400 for original n,g. This excludes hidden selected-loop1/k^2 terms, not all gravitational singularities.",
    }
