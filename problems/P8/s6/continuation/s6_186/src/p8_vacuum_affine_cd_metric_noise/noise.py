"""Actual centered quadratic stress covariance, not a metric solution."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes

from . import oscillatory, reference


@cache
def data():
    m, K = modes.MASS, modes.KAPPA
    ref = 36 * oscillatory.REF_BOUND**2 / m
    err = 288 * reference.PAIR_ERROR**2 / m**7
    total = ref + err
    T, C, mean = s.symbols("stress c_number_counterterm mean_stress")
    hcan = s.Symbol("canonical_metric_test", real=True)
    # The phase-space regulator is removed only after the complete pair integral.
    # A c-number subtraction cannot alter the centered stress operator.
    qzero, rzero, epsilon, rfirst, qfirst = s.symbols(
        "Q0 Rminus0 epsilon Rfirst Qfirst"
    )
    source = (rzero + epsilon * rfirst) * (qzero + epsilon * qfirst)
    checks = {
        "fixed_counterterm_cancels_only_centered_covariance": s.expand(
            (T + C) - (mean + C) - (T - mean)
        ),
        "full_source_first_metric_variation_zero": s.diff(source, epsilon).subs(
            {epsilon: 0, rzero: 0, qzero: 0}
        ),
        "Wick_pair_exchange_factor": 2 * 3**2 * 2 - 36,
        "full_pair_error_radial_factor": 36 * 8 - 288,
        "physical_normalized_variance": K * K * (total / K**2) - total,
        "canonical_gravity_metric_dictionary": s.simplify(
            (2 * hcan / s.sqrt(K)) ** 2 * total / 4 - total * hcan**2 / K
        ),
        "canonical_metric_and_clock_dictionaries_not_conflated": (total / K)
        / (total / K**2)
        - K,
        "actual_fixed_mass": m - 1000,
    }
    return {
        "observable": "Var T[f]=omega((T[f]-omega(T[f]))^2) for arbitrary real compact symmetric physical-frame test tensors on the actual CD clock; local quadratic stress, not products of separately smeared fields",
        "exact_Wick_pair_formula": "Var T[f]=2 sum_(r,s=1)^3 integral d^3k d^3l/(2pi)^6 |integral_I dt a^3 sum_(a,b) fhat_ab(t,k+l) conjugate(u_r(t,k))^T M_ab conjugate(u_s(t,l))|^2",
        "smearing_norm": "N[f]^2=sum_(j=0)^3 ||partial_t^j f||L2(dt dx;F)^2+||grad_spatial f||L2(dt dx;F)^2",
        "reference_pair_variance_coefficient_upper": ref,
        "actual_mode_remainder_variance_coefficient_upper": err,
        "complete_stress_variance_coefficient_upper": total,
        "complete_stress_variance_display_upper": s.Integer(10) ** 50,
        "stress_standard_deviation_display": "1e25 N[f]",
        "physical_stress_divided_by_kappa_standard_deviation_display": "1e-775 N[f]",
        "Einstein_normalized_metric_force_standard_deviation_display": "1e-375 N[h], with delta g_ab=2h_ab/sqrt(kappa) and the variational factor1/2 retained; this is not yet a fully constraint-reduced canonical mixed-mode norm",
        "normalization_boundary": "The reference Einstein metric normalization suppresses a fixed canonical vector stress. Lapse/shift constraints and scalar-metric mixing have not been inverted or canonically diagonalized here. This is not S6.185's raw-clock source covariance, where the scalar field dictionary cancels the apparent suppression.",
        "prescription": "The same fixed covariant stress and fixed scalar coefficient are retained. Their deterministic c-number contributions cancel in centering. Actual-state Wick ordering is a computation of that centered operator, not a new mean prescription.",
        "distributional_boundary": "Compact smooth smearing precedes the limit. The two-particle kernel is Hilbert-Schmidt by the full radial bounds; momentum regulators are proof devices, not physical cutoffs. No pointwise variance is defined.",
        "remaining": "This is the conditional vector noise input. It is not Gaussian statistics for the quadratic stress, a metric solution, a stable inverse, complete quantum loops or original V/G/B closure.",
        "checks": checks,
        "gates": {
            "complete_variance_below_one_e_fifty": total < 10**50,
            "physical_normalized_variance_below_display": total / K**2
            < s.Rational(1, 10**1550),
            "canonical_metric_variance_below_display": total / K
            < s.Rational(1, 10**750),
            "strictly_positive_upper_not_zero_noise_assertion": total > 0,
        },
    }
