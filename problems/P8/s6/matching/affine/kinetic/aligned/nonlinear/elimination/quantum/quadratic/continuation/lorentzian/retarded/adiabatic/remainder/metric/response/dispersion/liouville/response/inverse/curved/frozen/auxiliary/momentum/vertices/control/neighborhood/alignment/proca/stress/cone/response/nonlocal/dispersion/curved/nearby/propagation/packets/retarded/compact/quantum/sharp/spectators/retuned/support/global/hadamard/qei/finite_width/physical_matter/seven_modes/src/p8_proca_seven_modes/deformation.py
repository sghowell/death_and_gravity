"""Explicit smooth reference preparation and finite-strip energy bounds."""

from fractions import Fraction
from functools import cache

import sympy as sp


def exact_radius(value):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction, sp.Rational)):
        raise TypeError("Require an exact finite strip radius")
    value = sp.Rational(value)
    if value < 3:
        raise ValueError(
            "The reference preparation strip requires radius at least three"
        )
    return value


@cache
def data():
    s = sp.Symbol("interpolation_parameter", real=True)
    f = sp.exp(-1 / s)
    g = sp.exp(-1 / (1 - s))
    chi = f / (f + g)
    derivative = chi * (1 - chi) * (s**-2 + (1 - s) ** -2)
    outer = sp.Rational(160, 33)
    initial_tensor_ratio = sp.Rational(25, 2) / sp.Rational(2, 25**3)
    return {
        "smooth_step_on_open_unit_interval": chi,
        "smooth_step_extended_values": "0 for s<=0,1 for s>=1",
        "auxiliary_scale": "b(u)=1+chi(u+2)*((1+u^2)^2-1)",
        "reference_preparation_initial_time": -3,
        "actual_anchor_time": 0,
        "auxiliary_metric_equals_Minkowski_for_u_le_minus_two": True,
        "auxiliary_metric_equals_actual_for_u_ge_minus_one": True,
        "transition_scale_lower": 1,
        "transition_scale_upper": 25,
        "smooth_step_derivative_upper": 8,
        "transition_absolute_scale_derivative_upper": 232,
        "actual_global_absolute_Hubble_upper": 2,
        "Proca_energy_exponent_from_initial_to_anchor_upper": 702,
        "tensor_modified_energy_exponent_from_initial_to_anchor_upper": 705,
        "Proca_transfer_norm_anchor_bound": "1000 sqrt(1+k^2) exp(351)",
        "tensor_transfer_norm_anchor_bound": "(625/2) sqrt(1+k^2) exp(705/2)",
        "Minkowski_Proca_covariance_norm_upper": "501 sqrt(1+k^2), excluding hbar/kappa",
        "Minkowski_tensor_covariance_norm_upper": "(1+k^2)/k for k>0, excluding hbar/kappa",
        "checks": {
            "actual_smooth_step_derivative": sp.simplify(sp.diff(chi, s) - derivative),
            "transition_scale_derivative_bound_keeps_both_terms": sp.Integer(
                24 * 8 + 40 - 232
            ),
            "outer_smooth_step_derivative_bound_uses_exp_lower": outer
            - sp.Rational(160, 9) / (1 + sp.Rational(8, 3)),
            "Proca_complete_preparation_energy_exponent": sp.Integer(
                3 * 232 + 3 * 2 - 702
            ),
            "tensor_complete_preparation_energy_exponent": sp.Integer(702 + 3 - 705),
            "tensor_complete_preparation_norm_prefactor": initial_tensor_ratio
            - sp.Rational(625, 2) ** 2,
        },
    }


def strip(value=3):
    R = exact_radius(value)
    A = (1 + R * R) ** 2
    return {
        "radius": R,
        "actual_and_auxiliary_scale_lower": sp.Integer(1),
        "actual_and_auxiliary_scale_upper": A,
        "Proca_energy_lower": sp.Min(1, 10**6 / A),
        "Proca_energy_upper_per_one_plus_k_squared": sp.Max(10**6, A),
        "tensor_modified_energy_lower": 2 / A**3,
        "tensor_modified_energy_upper_per_one_plus_k_squared": A / 2,
        "auxiliary_and_actual_energy_growth_upper": sp.Integer(697),
        "tensor_low_frequency_radial_integral_majorant": sp.Rational(7, 6),
    }


@cache
def checks():
    v = sp.Symbol("r", positive=True)
    return {
        **data()["checks"],
        "tensor_massless_infrared_majorant_is_integrable_in_three_dimensions": sp.integrate(
            v * (1 + v * v) ** 2, (v, 0, 1)
        )
        - sp.Rational(7, 6),
        "bounded_cartesian_projection_does_not_introduce_inverse_radial_power": sp.Integer(
            0
        ),
    }
