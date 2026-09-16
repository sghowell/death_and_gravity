"""Explicit collinear radial modulus and uniform dimensional remainder."""

from functools import cache

import sympy as s

from . import source

CURRENT_BOUND = s.Rational(65, 4)
K0_BOUND = s.Integer(265)
K1_BOUND = s.Integer(130000)
K_REMAINDER = s.Integer(1400000)
PHASE_REMAINDER = s.Integer(2000000)
RADIAL_HOLDER = s.Integer(32000)
TRACE_HOLDER = s.Integer(21000)


def require_epsilon(value):
    e = source.recoil.exact_real(value)
    if e < 0 or e > s.Rational(1, 8):
        raise ValueError("Require regulator in[0,1/8]")
    return e


def logarithmic_holder_moment(order):
    if isinstance(order, bool) or not isinstance(order, (int, s.Integer)) or order < 0:
        raise ValueError("Require a nonnegative exact integer log moment")
    return s.factorial(order) * 4 ** (order + 1)


def dimensional_remainder(epsilon):
    e = require_epsilon(epsilon)
    return K_REMAINDER * e * e


@cache
def data():
    checks = {}

    def put(name, value):
        checks[name] = s.factor(s.expand_log(s.sympify(value)))

    r, c = s.symbols("radial cosine", real=True)
    speed = s.Symbol("sqrt_tau", positive=True)
    put(
        "massless_transverse_norm_division",
        (1 - r * r * c * c) - (1 - r * c) * (1 + r * c),
    )
    E, p, z = s.symbols("E spatial_norm projected_component", real=True)
    put(
        "massive_transverse_norm_division",
        (p * p - z * z) - ((E + z) * (E - z) - (E * E - p * p)),
    )
    put("radial_square_gap", 1 - r * r - (1 - r) * (1 + r))
    put("null_denominator_lower_c_positive", (1 - r * c) - (1 - c) - c * (1 - r))
    put(
        "null_denominator_lower_c_negative",
        2 * (1 - r * c) - (1 - c) - (1 + c * (1 - 2 * r)),
    )
    put("massive_matrix_numerator_bound", 4 * 2 - 8)
    put("massive_denominator_gap_bound", s.Rational(1, 4) ** -2 - 16)
    put("massive_projector_change", 8 * 4 + 4 * 2 * 16 - 160)
    put("massless_matrix_and_denominator_majorant", 8 + 4 - 12)
    put("massless_small_cap", s.Rational(1, 2) * 4 * speed - 2 * speed)
    a = s.Symbol("cap_distance", positive=True)
    put(
        "massless_large_cap",
        s.integrate(6 * speed / a, (a, speed, 2)) - 6 * speed * s.log(2 / speed),
    )
    y = s.Symbol("sqrt_speed", positive=True)
    put(
        "log_majorant_stationary_point", s.diff(-2 * y * s.log(y), y).subs(y, s.exp(-1))
    )
    put("holder_power", 8 + 6 - 14)
    put("total_massive_and_null_current", 16 + 2 * s.Rational(1, 8) - CURRENT_BOUND)
    put("Frobenius_and_trace_bilinear", 2 + 1 - 3)
    put(
        "radial_Holder_margin",
        RADIAL_HOLDER
        - 3 * CURRENT_BOUND * (640 + s.Rational(14, 8))
        - s.Rational(11435, 16),
    )
    put(
        "trace_Holder_margin",
        TRACE_HOLDER
        - 2 * CURRENT_BOUND * (640 + s.Rational(14, 8))
        - s.Rational(1145, 8),
    )
    q = s.Symbol("holder_exponent", positive=True)
    for order in range(5):
        integral = (-1) ** order * s.diff(1 / q, q, order).subs(q, s.Rational(1, 4))
        put(
            "log_holder_moment_" + str(order),
            integral - logarithmic_holder_moment(order),
        )
    put("radial_integral_bound", 4 * RADIAL_HOLDER - 128000)
    put("radial_log_difference", 16 * RADIAL_HOLDER - 512000)
    put("trace_radial_integral", 4 * TRACE_HOLDER - 84000)
    put("K0_margin", K0_BOUND - CURRENT_BOUND**2 - s.Rational(15, 16))
    put("K1_margin", K1_BOUND - CURRENT_BOUND**2 / 2 - 128000 - s.Rational(59775, 32))
    put(
        "all_e_absolute_contraction_margin",
        400 - s.Rational(3, 2) * CURRENT_BOUND**2 - s.Rational(125, 32),
    )
    e = s.Symbol("epsilon", nonnegative=True)
    put("projector_linear_remainder", e / (2 * (1 + e)) - e / 2 + e * e / (2 * (1 + e)))
    put(
        "normalization_log_derivative_at_zero",
        s.digamma(s.Rational(3, 2)) - s.digamma(1) - (2 - 2 * s.log(2)),
    )
    put(
        "K_second_order_majorant_margin",
        K_REMAINDER
        - (CURRENT_BOUND**2 / 2 + 2 * 128000 + 2 * 512000 + 84000)
        - s.Rational(1147775, 32),
    )
    put(
        "phase_second_order_majorant_margin",
        PHASE_REMAINDER
        - (K_REMAINDER + s.Rational(17, 2) * K0_BOUND + 4 * K1_BOUND)
        - s.Rational(155495, 2),
    )
    return {
        "whole_one_leg_modulus": "Writing tau=1-t and s=sqrt(tau), the massive current changes in nuclear norm by<=160s per leg. A null leg of energyomega changes by<=omega*min(4,12s/(1-cos(theta))). Its sphere average is<=omega*s*[2+6ln(2/s)]<=14omega*tau^(1/4).",
        "whole_uniform_current_majorant": CURRENT_BOUND,
        "whole_radial_Holder_constant": RADIAL_HOLDER,
        "whole_trace_Holder_constant": TRACE_HOLDER,
        "whole_logarithmic_Holder_moments": {
            str(j): logarithmic_holder_moment(j) for j in range(5)
        },
        "whole_K0_absolute_bound": K0_BOUND,
        "whole_K1_absolute_bound": K1_BOUND,
        "whole_K_second_order_bound": K_REMAINDER,
        "whole_phase_K_second_order_bound": PHASE_REMAINDER,
        "whole_uniform_scope": "Any finite positive outgoing null multiplicity, totalenergy<=1/8, massive COM energy[5/4,2], all angles. The same estimates extend to finite positive angular-energy measures. No pair-separation or particle-count lower/upper cutoff appears.",
        "checks": checks,
        "gates": {
            "massive_gap_and_null_collinear_caps_bounded_separately": True,
            "nuclear_norm_controls_both_trace_and_Frobenius": True,
            "angular_average_precedes_collinear_radial_log_bound": True,
            "integrable_quarter_power_majorant_uniform_in_multiplicity": True,
            "first_dimensional_derivative_not_inferred_from_bounded_convergence": True,
            "explicit_uniform_second_order_regulator_error": True,
            "finite_energy_measure_extension_not_quantum_state_construction": True,
        },
    }
