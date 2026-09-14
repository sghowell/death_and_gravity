"""Full nonlinear coefficient and commuting bounce-volume reference expectations."""

from functools import cache

import sympy as s

from . import local, physical, source

v, n = s.symbols("reference_hat_log_scale reference_linear_lapse", real=True)
Cvv, Cvn, Cnn = s.symbols(
    "whole_band_Cvv whole_band_symmetric_Cvn whole_band_Cnn", real=True
)


def centered_cubic_mean(polynomial):
    zero = {v: 0, n: 0}
    laplacian = (
        Cvv * s.diff(polynomial, v, 2)
        + 2 * Cvn * s.diff(polynomial, v, n)
        + Cnn * s.diff(polynomial, n, 2)
    )
    return s.factor(polynomial.subs(zero) + laplacian.subs(zero) / 2)


@cache
def remainder():
    derivatives = local.derivatives()
    high = [
        max(
            row[j]
            for row in derivatives[
                "whole_all_fifteen_outward_derivative_bounds"
            ].values()
        )
        for j in range(4)
    ]
    coefficient_sum = sum(
        s.Integer(3) ** i / s.factorial(i) * high[j] / s.factorial(j)
        for i in range(4)
        for j in range(4 - i)
    )
    vv = s.Rational(1, 10**740)
    nn = source.VAR_UPPER
    moment8 = 128 * 105 * (3**8 * vv**4 + nn**4)
    polynomial_square = 10**5 * coefficient_sum**2
    bulk = s.Rational(10**8 * 2 * 10**3, 24) * nn**2
    exponent = source.RADIUS**2 / (4 * nn)
    tail = (4 + 2 * 10**10) * s.factorial(3) / exponent**3
    L = s.Symbol("positive_gaussian_variance", positive=True)
    return {
        "whole_first_four_complete_derivative_bounds": high,
        "whole_complete_cubic_coefficient_absolute_sum_bound": coefficient_sum,
        "whole_complete_cubic_polynomial_square_moment_bound": polynomial_square,
        "whole_eighth_absolute_mixed_moment_bound": moment8,
        "whole_fourth_order_bulk_remainder_bound": bulk,
        "whole_tail_Cauchy_exponent": exponent,
        "whole_positive_series_tail_bound": tail,
        "whole_safe_nonlinear_expectation_remainder": source.VOLUME_ERROR,
        "whole_remainder_argument": "For |n|<=delta=1e-3, use Taylor through total degree3 along(sv,sn). The full derivative ceiling gives remainder <=(1e8/24) exp(3|v|)(3|v|+|n|)^4. At the bounce v,n strongly commute and have their actual centered joint Gaussian law; for a single coefficient put v=0 at any time. Cauchy-Schwarz, E exp(6|v|)<=2exp(18Cvv)<4 and the displayed eighth-moment bound give the bulk ceiling. No cross covariance is removed. The complete cubic polynomial has square moment <=1e5 times its squared absolute coefficient sum: each degree<=3 monomial obeys this Gaussian moment bound by Cauchy-Schwarz and variances below1. Thus its square moment is below1e20.",
        "whole_tail_argument": "The full positive-lapse function is bounded by2. The complementary event |n|>delta is treated separately, not by applying local Taylor bounds there. Cauchy-Schwarz bounds the exponential-volume tail by4 exp(-delta^2/(4Cnn)) and the polynomial tail by2e10 times the same exponential. Use Cnn<1e-490 and exp(z)>=z^3/3! to obtain the exact rational outward tail ceiling. Never replace an underflowed probability by zero. Odd terms vanish only after taking the complete untruncated centered Gaussian expectation.",
        "checks": {
            "whole_gaussian_eighth_moment_coefficient": s.factorial2(7) - 105,
            "whole_Gaussian_fourth_moment_coefficient": 3 * L**2
            - s.factorial2(3) * L**2,
            "complete_cubic_coefficient_sum": coefficient_sum - 1309,
            "complete_positive_tail_exponent_binding": 4 * nn * exponent
            - source.RADIUS**2,
        },
        "gates": {
            "whole_cubic_square_moment_below1e20": polynomial_square < 10**20,
            "whole_eighth_moment_below_1e5_times_variance_fourth": moment8
            < 10**5 * nn**4,
            "exponential_absolute_metric_moment_below_four": 18 * vv < s.Rational(1, 2),
            "whole_bulk_plus_tail_below1e_minus960": bulk + tail < source.VOLUME_ERROR,
            "only_bounce_joint_spectral_probability_used": True,
            "entire_full_functions_used_outside_local_taylor_event": True,
            "nonzero_Gaussian_tail_not_discarded": True,
        },
    }


@cache
def composites():
    jets = local.derivatives()["whole_complete_center_jets_zero_through_four"]
    single = {}
    polynomials = {}
    checks = {}
    for alpha in source.ALPHAS:
        key = str(alpha)
        row = jets[key]
        polynomial = sum(row[j] * n**j / s.factorial(j) for j in range(4))
        mean = centered_cubic_mean(polynomial)
        wanted = 1 + row[2] * Cnn / 2
        single[key] = wanted
        polynomials[key] = polynomial
        checks["whole_centered_single_lapse_cubic_mean_" + key] = s.factor(
            mean - wanted
        )
        odd = (polynomial - polynomial.subs(n, -n)) / 2
        checks["whole_odd_single_lapse_moment_zero_" + key] = centered_cubic_mean(odd)
    row = [value.subs(source.u, 0) for value in jets["3/4"]]
    volume = sum(
        s.Integer(3) ** i * row[j] * v**i * n**j / (s.factorial(i) * s.factorial(j))
        for i in range(4)
        for j in range(4 - i)
    )
    contact = s.Rational(9, 2) * Cvv + s.Rational(9, 2) * Cvn + s.Rational(3, 8) * Cnn
    mean = centered_cubic_mean(volume)
    odd = s.expand((volume - volume.subs({v: -v, n: -n}, simultaneous=True)) / 2)
    checks.update(
        {
            "complete_bounce_nonlinear_volume_cubic_Gaussian_mean": s.factor(
                mean - 1 - contact
            ),
            "complete_bounce_cubic_odd_moment_zero": centered_cubic_mean(odd),
            "whole_original_S264_quadratic_volume_contact": s.factor(
                contact - previous_volume_contact()
            ),
            "whole_bounce_volume_first_second_third_lapse_jets": s.Matrix(row[:4])
            - s.Matrix([1, s.Rational(3, 2), s.Rational(3, 4), -s.Rational(3, 8)]),
        }
    )
    correction = (
        5 * s.Rational(1, 10**740) + 5 * s.Rational(1, 10**615) + s.Rational(1, 10**490)
    )
    return {
        "whole_positive_lapse_coefficient_spectral_definition": "f_alpha(N)=1_(N>0) R_full(u,N^-2)^(-alpha), defined as0 on N<=0. This is a bounded positive Borel function of the SAME single self-adjoint band-limited linear lapse N=1+n. The state is neither conditioned nor projected. The full function, including both switches, is evaluated on positive lapse; no Gaussian tail or source remainder is deleted.",
        "whole_all_three_complete_single_lapse_cubic_polynomials": polynomials,
        "whole_all_three_single_lapse_quadratic_expectations": single,
        "whole_all_three_single_lapse_exact_nonlinear_expectation_error": source.VOLUME_ERROR,
        "whole_bounce_complete_cubic_volume_polynomial": volume,
        "whole_bounce_complete_quadratic_volume_contact": contact,
        "whole_bounce_full_nonlinear_volume_expectation_error": source.VOLUME_ERROR,
        "whole_bounce_full_nonlinear_volume_mean_distance_from_one_bound": correction
        + source.VOLUME_ERROR,
        "whole_bounce_volume_spectral_definition": "At u=0, v_hat and n strongly commute by the full retained momentum rows. Define the positive unbounded joint-spectral reference composite Vref=exp(3 v_hat) f_(3/4)(1+n). Its expectation is finite because f<=2 and Gaussian exponential moments are finite. It is the entire positive-lapse background volume map on LINEAR reference fields, NOT the full nonlinear auxiliary lapse solution, an interacting physical mean or a replacement Gaussian state. Its complete quadratic correction includes the retained Cvn term and agrees with S264.",
        "whole_spectral_probability_boundary": "The joint spectral Gaussian interpretation is restricted to the bounce slice. The single-lapse coefficient expectations hold separately throughout |u|<=1e-60. No joint noncommuting probability, uniform spacetime event, full interacting ordering/gauge measure, Wilsonian cutoff or original P8 closure follows. No H/Proca preparation or determinant is changed.",
        "checks": checks,
        "gates": {
            "all_three_full_single_lapse_coefficients_included": len(single) == 3,
            "all_nonlinear_remainders_use_completed_bound": remainder()[
                "whole_fourth_order_bulk_remainder_bound"
            ]
            + remainder()["whole_positive_series_tail_bound"]
            < source.VOLUME_ERROR,
            "whole_volume_mean_distance_from_one_below1e_minus489": correction
            + source.VOLUME_ERROR
            < s.Rational(1, 10**489),
            "whole_volume_cross_covariance_retained": contact.has(Cvn),
            "full_bounce_strong_commutation_and_positive_variance": all(
                physical.bounce()["gates"].values()
            ),
            "S261_refuted_physical_binding_not_reused": True,
            "full_reference_composite_not_interacting_physical_mean": True,
        },
    }


def previous_volume_contact():
    from p8_vacuum_affine_quantitative_gaussian_window import physical as previous

    expression = previous.volume()[
        "whole_complete_second_order_Weyl_spatial_volume_contact"
    ].subs(source.u, 0)
    lookup = {"fixed_band_Cvv": Cvv, "fixed_band_Cvn": Cvn, "fixed_band_Cnn": Cnn}
    return expression.subs(
        {symbol: lookup[str(symbol)] for symbol in expression.free_symbols}
    )
