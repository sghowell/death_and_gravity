"""One real conformal scalar: actual occupation state and exact named SEE history."""

from functools import cache

import sympy as sp

from . import cosmology, stress

PARENT_SHA = "d77f1a4ef293c53c70bb7fd540b473c557c36018dfbea558ea1af9afc147375c"
ERROR_CAPS = tuple(map(sp.Rational, ("1/100", "1/2", "3", "16")))
RESPONSE = tuple(map(sp.Rational, ("16/25", "1728/25", "155136/25", "14770176/25")))


def coupling(delta):
    delta = cosmology.exact_nonnegative(delta)
    if not 0 < delta <= cosmology.DELTA_MAX:
        raise ValueError("Require 0 < delta <= 1e-8 in this named scalar SEE example")
    return delta / 360


def velocity(y, lam):
    return 2 * y * y * (1 - lam * y * y) / (1 - 2 * lam * y * y)


def jets(y, lam):
    z = sp.Dummy("y", positive=True)
    v = velocity(z, lam)
    second = sp.factor(v * sp.diff(v, z))
    third = sp.factor(v * sp.diff(second, z))
    return [y, v.subs(z, y), second.subs(z, y), third.subs(z, y)]


def clock(y, lam):
    return 1 / (2 * y) + sp.sqrt(lam) * sp.atanh(sp.sqrt(lam) * y) / 2


def branch_point(y, delta):
    y = cosmology.exact_nonnegative(y)
    lam = coupling(delta)
    if y <= 0 or 1 - 2 * lam * y * y <= 0:
        raise ValueError("Require an interior point of the actual positive low branch")
    a4 = 4 * (1 - 4 * lam) / (y * y * (1 - lam * y * y))
    return {
        "lambda": lam,
        "y": y,
        "branch_denominator": 1 - 2 * lam * y * y,
        "a_fourth": a4,
        "normal_Hubble_jets": [-v for v in jets(y, lam)],
        "zeta_squared": sp.Rational(20, 9) * sp.Rational(delta) * (1 - 4 * lam) / a4,
        "kappa_EED_times_tau_squared": 3 * (velocity(y, lam) - y * y),
    }


def history(delta):
    lam = coupling(delta)
    errors = [v * lam for v in RESPONSE]
    margins = [cap - error for cap, error in zip(ERROR_CAPS, errors, strict=True)]
    zeta2 = sp.Rational(20, 9) * sp.Rational(delta) * (1 - 4 * lam)
    if min(margins) <= 0 or cosmology.ZETA_MAX**2 - zeta2 <= 0:
        raise ValueError("The actual scalar state/history fails its named strict gates")
    return {
        "delta": sp.Rational(delta),
        "lambda": lam,
        "past_x": [-cosmology.RATIO, sp.Integer(0)],
        "backward_y_bounds": [sp.Rational(50, 27), sp.Integer(2)],
        "a_past_lower": sp.Integer(1),
        "C3_errors": errors,
        "strict_C3_margins": margins,
        "anchor_H0_times_tau": sp.Integer(2),
        "initial_and_past_zeta_squared_upper": zeta2,
        "strict_squared_field_gate_margin": cosmology.ZETA_MAX**2 - zeta2,
        "sigma": sp.Integer(0),
        "beta_S": sp.Integer(0),
        "actual_scalar_state_and_SEE_past": True,
        "future_state_cap_proved_on_every_shorter_segment": False,
    }


@cache
def checks():
    hb, bT, a, kappa, tau, y, lam, b = sp.symbols(
        "hbar beta_T a kappa tau y lambda b", positive=True
    )
    h, hd, q = sp.symbols("H Hdot Q", real=True)
    # Bose moments come from the positive exponential series and Tonelli.
    energy = hb / (2 * sp.pi**2) * sp.gamma(4) * sp.zeta(4) / bT**4
    wick = hb / (2 * sp.pi**2) * sp.gamma(2) * sp.zeta(2) / bT**2
    r = hb / (960 * sp.pi**2)
    rho = q / a**4 + r * h**4
    pressure = q / (3 * a**4) - r * h**4 - 4 * r * h * h * hd / 3
    b_actual = kappa * r / 3
    delta = kappa * hb / (8 * sp.pi**2 * tau * tau)
    constraint_q = a**4 * (3 * h * h / kappa - r * h**4)
    reduced_pressure = (2 * hd + 3 * h * h + kappa * pressure).subs(q, constraint_q)
    F = velocity(y, lam)
    a4 = 4 * (1 - 4 * lam) / (y * y * (1 - lam * y * y))
    d = 1 - 2 * lam * y * y
    z = lam * y * y
    actual_jets = jets(y, lam)
    deviations = [
        2 * lam * y**4 / d,
        8 * lam * y**5 * (6 * z * z - 8 * z + 3) / d**3,
        16 * lam * y**6 * (6 * z * z - 5 * z + 2) * (14 * z * z - 22 * z + 9) / d**5,
    ]
    Q0 = 12 * (1 - 4 * lam) / (kappa * tau * tau)
    # The occupation-state relation is squared, avoiding inexact square roots.
    initial_zeta2 = (kappa / 3) ** 2 * Q0 * 5 * hb / (24 * sp.pi**2)
    result = {
        "one_real_scalar_thermal_Bose_energy": sp.simplify(
            energy - hb * sp.pi**2 / (30 * bT**4)
        ),
        "actual_relative_Wick_square_Bose_moment": sp.simplify(
            wick - hb / (12 * bT * bT)
        ),
        "same_state_energy_Wick_square_relation": sp.simplify(
            energy - 24 * sp.pi**2 * wick**2 / (5 * hb)
        ),
        "positive_oscillator_CCR_difference": (1 + sp.Symbol("n")) - sp.Symbol("n") - 1,
        "own_scalar_reference_density_coefficient": sp.simplify(
            stress.reference_jets(h, hd, 0, 0, 0)["rho"] * hb / (2880 * sp.pi**2)
            - r * h**4
        ),
        "actual_scalar_SEE_density": sp.factor(
            3 * h * h
            - kappa * rho
            - 3 * (h * h - b_actual * h**4 - kappa * q / (3 * a**4))
        ),
        "independent_scalar_SEE_pressure": sp.factor(
            reduced_pressure
            - 2 * (1 - 2 * b_actual * h * h) * hd
            - 4 * h * h * (1 - b_actual * h * h)
        ),
        "actual_scalar_delta_dictionary_not_photon": sp.simplify(
            b_actual / tau**2 - delta / 360
        ),
        "full_scalar_stress_conservation": sp.factor(
            -4 * h * q / a**4 + 4 * r * h**3 * hd + 3 * h * (rho + pressure)
        ),
        "exact_low_branch_pressure_solution": sp.factor(
            (2 * (1 - 2 * b * h * h) * hd + 4 * h * h * (1 - b * h * h)).subs(
                hd, -2 * h * h * (1 - b * h * h) / (1 - 2 * b * h * h)
            )
        ),
        "exact_scalar_proper_clock": sp.simplify(sp.diff(clock(y, lam), y) * F + 1),
        "exact_scalar_scale_factor_clock": sp.simplify(sp.diff(a4, y) * F + 4 * y * a4),
        "exact_scalar_density_first_integral": sp.factor(
            y * y * (1 - lam * y * y) * a4 - 4 * (1 - 4 * lam)
        ),
        "anchored_scale_factor": sp.factor(a4.subs(y, 2) - 1),
        "reciprocal_clock_history_identity": sp.factor(
            -F / y**2 + 2 + 2 * lam * y * y / d
        ),
        "initial_actual_field_amplitude_from_same_thermal_state": sp.simplify(
            initial_zeta2 - sp.Rational(20, 9) * delta * (1 - 4 * lam)
        ),
        "endpoint_field_budget_not_uniformly_small": sp.factor(
            (sp.Rational(20, 9) * (360 * lam) * (1 - 4 * lam) / a4).subs(
                y, 1 / sp.sqrt(2 * lam)
            )
            - 50
        ),
        "actual_low_branch_positive_EED": sp.factor(3 * (F - y * y) - 3 * y * y / d),
    }
    result.update(
        {
            f"continuous_scalar_C3_response_identity_{j + 1}": sp.factor(
                actual - base - diff
            )
            for j, (actual, base, diff) in enumerate(
                zip(
                    actual_jets[1:],
                    [2 * y * y, 8 * y**3, 48 * y**4],
                    deviations,
                    strict=True,
                )
            )
        }
    )
    return result


def calibration():
    maximum = cosmology.DELTA_MAX
    return {
        "actual_thermal_past_at_delta_max": history(maximum),
        "lambda_to_delta": sp.Rational(1, 360),
        "generic_ODE_comparison_coupling_margin": sp.Rational(1, 16)
        - coupling(maximum),
        "state_Wick_to_energy_coefficient_without_hbar_over_pi_squared": sp.Rational(
            5, 24
        ),
        "endpoint_zeta_squared": sp.Integer(50),
        "endpoint_x_strict_bounds": [sp.Rational(1, 16), sp.Rational(1, 4)],
        "thermal_branch_is_not_a_new_QEI_only_endpoint_argument": True,
        "future_Wick_bound_is_a_conditional_theorem_hypothesis": True,
    }
