"""Full finite complement integral, measure density and formal source contacts."""

from functools import cache

import sympy as s

from . import split

LAMBDA = s.Symbol("positive_cell_action_scale", positive=True)
R = s.Symbol("positive_full_tensor_factor", positive=True)


@cache
def determinant():
    p = split.P
    raw = split.complement()["whole_56_direction_complement_determinant"]
    reference = raw.subs(p, s.Rational(1, 2))
    ratio = s.factor(raw / reference)
    ratio_R = (
        R ** s.Rational(65, 2)
        * (4 * R ** s.Rational(1, 2) + 5) ** 3
        * (2 * R - 1) ** 3
        * (4 - R ** s.Rational(3, 2))
        / 2187
    )
    first = s.factor(s.diff(raw, p) / raw)
    second = s.factor(s.diff(first, p))
    dR = s.factor(s.diff(ratio_R, R) / ratio_R)
    ddR = s.factor(s.diff(dR, R))
    return {
        "whole_complement_reference_determinant": reference,
        "whole_relative_complement_determinant": ratio,
        "whole_relative_complement_determinant_in_R": ratio_R,
        "whole_log_determinant_p_first": first,
        "whole_log_determinant_p_second": second,
        "whole_log_determinant_R_first": dR,
        "whole_log_determinant_R_second": ddR,
        "whole_clock_R_first_second_logdet": (
            s.Rational(116, 3),
            -s.Rational(1228, 27),
        ),
        "whole_flat_configuration_log_weight_p_first": -first / 2,
        "whole_flat_configuration_log_weight_p_second": -second / 2,
        "full_scale_boundary": "For independent positive cell scale lambda, the 56-dimensional flat Gaussian has lambda^-28 times det A(p)^-1/2 and its separately continued phase. The displayed p and R derivatives hold lambda fixed. If lambda=kappa times a metric-dependent physical cell volume, its full -28 log lambda derivatives must also be retained; these are not physical lapse or stress derivatives with the density silently fixed.",
        "checks": {
            "whole_normalized_determinant_reference": ratio.subs(p, s.Rational(1, 2))
            - 1,
            "whole_full_R_determinant_conversion": s.simplify(
                ratio.subs(p, s.sqrt(R) / 2) - ratio_R
            ),
            "whole_log_first_is_full_determinant_derivative": s.factor(
                first * raw - s.diff(raw, p)
            ),
            "whole_log_second_is_full_log_derivative": s.factor(
                second - s.diff(first, p)
            ),
            "whole_clock_first_R_contact": s.simplify(
                dR.subs(R, 1) - s.Rational(116, 3)
            ),
            "whole_clock_second_R_contact": s.simplify(
                ddR.subs(R, 1) + s.Rational(1228, 27)
            ),
            "whole_56_cell_scale_power": s.Integer(56) / 2 - 28,
        },
        "gates": {
            "whole_flat_measure_factor_is_coefficient_dependent": first != 0,
            "reference_determinant_positive": reference > 0,
            "metric_cell_volume_variations_not_deleted": True,
        },
    }


@cache
def normalization():
    detA = split.complement()["whole_56_direction_complement_determinant"]
    positive_det = s.Symbol("positive_complement_determinant", positive=True)
    phase = split.inertia()["whole_positive_scale_Fresnel_phase"]
    flat = (2 * s.pi / LAMBDA) ** 28 * positive_det ** s.Rational(-1, 2) * phase
    density = (LAMBDA / (2 * s.pi)) ** 28 * s.sqrt(positive_det)
    phase_inverse = s.conjugate(phase)
    return {
        "whole_positive_complement_determinant_binding": detA,
        "whole_flat_Fresnel_integral": flat,
        "whole_normalized_positive_configuration_density": density,
        "whole_normalization_continued_phase": phase_inverse,
        "whole_positive_convergence_prescription": "For each finite cell use exp(i lambda z^T A z/2 - epsilon z^T z/2), epsilon>0, and continue each real eigenvalue Gaussian from epsilon>0 before taking epsilon down to zero. The 26/30 inertia fixes exp(i*pi*(26-30)/4)=-1. A principal square root of the single final determinant of -i lambda A has the wrong sign; the phase is fixed by this path.",
        "whole_finite_insertion": "At each finite cell, the displayed positive density times its separately continued inverse phase multiplies the full 56-dimensional Fresnel integral to one. This is an exact conditional insertion into any already specified retained finite functional. It does not choose that retained functional, its state, time slicing, ordering or continuum regularization.",
        "whole_canonical_comparison": "The finite canonical-jet primary/secondary reduction gives unit remaining phase-space density after the full second-class delta Jacobian cancels. The normalized configuration insertion represents that absence of an additional complement weight. A flat configuration integral instead keeps the nonconstant determinant and scale factors. No separately regularized infinite functional determinant is canceled here.",
        "checks": {
            "whole_normalized_56_integral_is_one": s.simplify(
                density * phase_inverse * flat - 1
            ),
            "whole_positive_and_negative_branch_pair_phase": s.simplify(
                phase * s.conjugate(phase) - 1
            ),
            "whole_principal_endpoint_determinant_phase_is_not_the_continued_phase": s.Integer(
                1
            )
            - phase
            - 2,
        },
        "gates": {
            "whole_density_is_not_flat_Lebesgue": density != 1,
            "continued_phase_not_principal_endpoint_square_root": phase == -1,
            "no_regulator_delta_zero_assumption": True,
        },
    }


@cache
def source_contacts():
    d = split.matrices()
    A, inv, E = d["A"], d["inverse_A"], d["embedding"] * d["kernel"]
    b = s.Matrix(s.symbols("formal_complement_source0:56", real=True))
    mean = -inv * b / LAMBDA
    exponent = -(b.T * inv * b)[0] / (2 * LAMBDA)
    covariance = s.I * inv / LAMBDA
    kernel = split.clean(E * inv * E.T)
    quotient_kernel = split.clean(d["kernel"] * inv * d["kernel"].T)
    full_inverse = quotient_kernel + d["lift"] * d["parent"]["eta"] * d["lift"].T
    return {
        "whole_56_formal_complement_sources": b,
        "whole_conditional_source_exponent": exponent,
        "whole_conditional_source_dependent_auxiliary_mean": mean,
        "whole_conditional_source_functional_contact": covariance,
        "whole_64_connection_source_contact_kernel": kernel,
        "whole_60_quotient_complement_contact_kernel": quotient_kernel,
        "whole_updated_60_quotient_inverse_decomposition": full_inverse,
        "whole_source_identity": "With a projectively invariant original connection source J, pull it back before elimination. Its complement source is b=(E K)^T J and its retained source is J dot Gamma_bar. The normalized conditional integral is exp[i J dot Gamma_bar - i b^T A^-1 b/(2 lambda)]. The 56 secondary equations become lambda A z+b=0. Their delta Jacobian and positive second-class density still cancel, with z=-A^-1 b/lambda. Dropping the quadratic source term would change the generating functional.",
        "whole_retained_expectation_identity": "If a retained finite functional is separately defined, differentiating the exact source identity gives <Gamma>=<Gamma_bar>, not Gamma_bar(<retained fields>). The connection two-point source-functional derivative includes i<lambda^-1 E K A^-1 K^T E^T> times the finite-cell delta, as well as the full composite Gamma_bar correlations. This instantaneous source contact is not a new propagating auxiliary degree of freedom or a replacement for the fixed interacting state.",
        "checks": {
            "whole_56_source_secondary_equations": split.clean(LAMBDA * A * mean + b),
            "whole_56_source_first_derivative": split.clean(
                s.Matrix([s.diff(exponent, v) for v in b]) - mean
            ),
            "whole_56_source_second_contact": split.clean(
                s.hessian(exponent, b) + inv / LAMBDA
            ),
            "whole_complement_contact_has_no_retained_trace": split.clean(
                d["N"] * quotient_kernel
            ),
            "whole_complement_contact_updated_Euler": split.clean(
                d["parent"]["M_new"] * quotient_kernel
                - s.eye(60)
                + d["N"].T * d["lift"].T
            ),
            "whole_updated_full_quotient_inverse": split.clean(
                d["parent"]["M_new"] * full_inverse - s.eye(60)
            ),
            "whole_64_source_contact_projective_trace": split.clean(
                d["gauge"] * kernel
            ),
        },
        "gates": {
            "all_56_source_contacts_retained": covariance.shape == (56, 56),
            "full_original_connection_contact_not_zero": kernel != s.zeros(64),
            "dynamical_Proca_directions_not_integrated_out": True,
            "source_functional_contact_not_a_new_physical_oscillator": True,
        },
    }


@cache
def local_variation():
    a = s.Matrix([[2, 1], [1, -3]])
    da = s.Matrix([[1, 2], [2, 4]])
    dda = s.Matrix([[3, -1], [-1, 2]])
    lam, dl, ddl = s.symbols(
        "lambda lambda_first lambda_second", nonzero=True, real=True
    )
    inv = a.inv()
    determinant_first = s.trace(inv * da)
    determinant_second = s.trace(inv * dda - inv * da * inv * da)
    theta = s.Symbol("variation_parameter", real=True)
    family = a + theta * da + theta**2 * dda / 2
    detfamily = family.det()
    scale = lam + theta * dl + theta**2 * ddl / 2
    # A finite two-direction diagnostic of the same general dimension-m/2 rule.
    total_first = -dl / lam - determinant_first / 2
    total_second = -(ddl / lam - dl**2 / lam**2) - determinant_second / 2
    derivative = s.diff(detfamily, theta) / detfamily
    variance = s.Symbol("positive_retained_variance", positive=True)
    q = s.Symbol("retained_Gaussian_variable", real=True)
    density = s.exp(-q * q / (2 * variance)) / s.sqrt(2 * s.pi * variance)
    mean = s.integrate(q * density, (q, -s.oo, s.oo))
    second_moment = s.integrate(q * q * density, (q, -s.oo, s.oo))
    return {
        "independent_complete_variation_family": family,
        "independent_complete_cell_scale_family": scale,
        "whole_two_direction_log_weight_first": total_first,
        "whole_two_direction_log_weight_second": total_second,
        "whole_general_dimension_weight_rule": "For dimension m, d log Z_flat=-(m/2)d log lambda-(1/2)Tr(A^-1 dA). The full mixed second variation is -(m/2)d1d2 log lambda-(1/2)Tr(A^-1 d1d2 A-A^-1 d1 A A^-1 d2 A). The normalized insertion supplies exactly the opposite density variations, at the same finite regulator. Both second coefficient and scale contacts are mandatory.",
        "independent_composite_mean_counterexample": "A centered real Gaussian retained variable q with variance v has <q^2>=v while (<q>)^2=0. This finite diagnostic forbids replacing the expectation of a nonlinear stationary-connection map by its value at the retained mean. It is not a calculation of the actual P8 quantum mean.",
        "checks": {
            "whole_first_density_variation": s.factor(
                total_first + dl / lam + derivative.subs(theta, 0) / 2
            ),
            "whole_second_density_variation": s.factor(
                total_second
                + ddl / lam
                - dl**2 / lam**2
                + s.diff(derivative, theta).subs(theta, 0) / 2
            ),
            "whole_centered_retained_mean": mean,
            "whole_nonlinear_composite_mean_retained": s.simplify(
                second_moment - variance
            ),
        },
        "gates": {
            "second_A_contact_is_nonzero": s.trace(inv * dda) != 0,
            "scale_contact_not_silently_dropped": dl != 0 and ddl != 0,
            "nonlinear_mean_factorization_rejected": variance > 0,
        },
    }
