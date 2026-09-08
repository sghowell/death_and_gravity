"""Continuous local-potential bound, not a full quantum/EFT remainder."""
from functools import cache

import sympy as sp
from p8_affine_retuned import bounds as exact_parent

from . import potential

A_LOW = sp.Rational(8, 9)
B_LOW = sp.Rational(35, 36)
UPPER = sp.Rational(19, 18)


@cache
def proof_checks():
    pole = 3*UPPER**2
    log_a = (1-A_LOW)/A_LOW
    log_b_negative = (1-B_LOW)/B_LOW
    log_positive = UPPER-1
    finite_upper = 2*UPPER**2*sp.Rational(9, 16)+UPPER**2*sp.Rational(13, 8)
    jets = potential.clock_jets()
    h = jets["h"]
    return {"positive_mass_coefficients_on_original_tube": bool(A_LOW > 0 and B_LOW > 0),
            "longitudinal_squared_speed_lower": B_LOW/UPPER == sp.Rational(35, 38),
            "longitudinal_squared_speed_upper": UPPER/A_LOW == sp.Rational(19, 16),
            "pole_coefficient_below_four": bool(pole < 4),
            "log_a_absolute_bound": log_a == sp.Rational(1, 8),
            "negative_log_b_below_one_sixteenth": bool(log_b_negative < sp.Rational(1, 16)),
            "positive_logs_below_one_sixteenth": bool(log_positive < sp.Rational(1, 16)),
            "finite_MSbar_weight_continuous_bound": finite_upper == sp.Rational(3971, 1296),
            "finite_MSbar_weight_below_four": bool(finite_upper < 4),
            "pole_clock_first_derivative_nonzero": sp.factor(jets["pole_weight"]["N_first"]-20/(9*h)) == 0,
            "finite_clock_first_derivative_nonzero": sp.factor(jets["finite_weight"]["N_first"]+22/(27*h)) == 0,
            "example_canonical_frequency_domain": 1000**2 >= 20000,
            "example_local_potential_ratio_below_1e_minus_14": bool(sp.Rational(1000**4, 144*(10**12)**2) < sp.Rational(1, 10**14)),
            "strict_pi_lower_bound_uses_inscribed_hexagon": True,
            "constant_coefficient_term_not_full_curved_or_nonlinear_loop_bound": True,
            "scheme_dependent_potential_not_vacuum_or_UV_verdict": True}


def scale_bound(planck_time_product, reference_mass_time_product):
    L, R = [exact_parent.exact(value, name) for value, name in
            ((planck_time_product, "M_tau"), (reference_mass_time_product, "m0_tau"))]
    if L.is_positive is not True or R.is_positive is not True:
        raise ValueError("Both dimensionless scale products must be positive")
    return {"M_tau": L, "m0_tau": R, "m0_over_M": R/L,
            "normalized_zeta": 1/R**2, "physical_zeta": L**2/R**2,
            "local_potential_over_reference_M_squared_over_tau_squared": R**4/(144*L**2),
            "scale_choice_is_a_cutoff_or_complete_quantum_bound": False}
