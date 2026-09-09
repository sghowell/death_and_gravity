"""Actual scalar SEE past compatibility, not a new thermal endpoint mechanism."""

from functools import cache

import sympy as sp
from p8a_nonminimal import thermal as prior

from . import costs

DELTA = sp.Rational(1, 10**14)


@cache
def data():
    history = prior.history(DELTA)
    a = sp.Symbol("a", positive=True)
    H, Hd, Hdd, Hthird = sp.symbols("H Hdot Hddot Hthird", real=True)

    def d_affine(expr):
        return (
            sp.diff(expr, a) * a * H
            + sp.diff(expr, H) * Hd
            + sp.diff(expr, Hd) * Hdd
            + sp.diff(expr, Hdd) * Hthird
        ) / a

    jets = [
        H,
        Hd / a,
        (Hdd - H * Hd) / a**2,
        (Hthird - Hd**2 - 3 * H * Hdd + 2 * H**2 * Hd) / a**3,
    ]
    h = sp.Rational(21, 10)
    d1 = sp.Integer(9)
    d2 = sp.Integer(70)
    d3 = sp.Integer(800)
    bounds = [h, d1, d2 + h * d1, d3 + d1 * d1 + 3 * h * d2 + 2 * h * h * d1]
    a_upper = 1 / (1 - h * costs.R)
    return {
        "actual_parent_scalar_history": history,
        "affine_past_A_upper": a_upper,
        "affine_derivative_upper_bounds": bounds,
        "strict_affine_derivative_margins": [
            sp.Rational(cap) - value
            for cap, value in zip(costs.PAST[2], bounds, strict=True)
        ],
        "strict_null_field_budget_squared_margin": sp.Rational(1, 10**12)
        - history["initial_and_past_zeta_squared_upper"],
        "strict_null_quantum_budget_margin": sp.Rational(1, 10**12) - DELTA,
        "actual_future_affine_endpoint_strict_upper_in_tau_units": sp.Rational(1, 4),
        "checks": {
            "actual_thermal_proper_to_affine_first_jet": sp.factor(
                d_affine(a) - jets[0]
            ),
            "actual_thermal_proper_to_affine_second_jet": sp.factor(
                d_affine(jets[0]) - jets[1]
            ),
            "actual_thermal_proper_to_affine_third_jet": sp.factor(
                d_affine(jets[1]) - jets[2]
            ),
            "actual_thermal_proper_to_affine_fourth_jet": sp.factor(
                d_affine(jets[2]) - jets[3]
            ),
            "actual_scalar_initial_scale_and_Hubble_anchor": prior.branch_point(
                2, DELTA
            )["a_fourth"]
            - 1,
            "actual_scalar_initial_Hubble_is_negative_two": prior.branch_point(
                2, DELTA
            )["normal_Hubble_jets"][0]
            + 2,
            "same_scalar_state_supplies_past_field_amplitude": history[
                "initial_and_past_zeta_squared_upper"
            ]
            - sp.Rational(20, 9) * DELTA * (1 - DELTA / 90),
        },
    }
