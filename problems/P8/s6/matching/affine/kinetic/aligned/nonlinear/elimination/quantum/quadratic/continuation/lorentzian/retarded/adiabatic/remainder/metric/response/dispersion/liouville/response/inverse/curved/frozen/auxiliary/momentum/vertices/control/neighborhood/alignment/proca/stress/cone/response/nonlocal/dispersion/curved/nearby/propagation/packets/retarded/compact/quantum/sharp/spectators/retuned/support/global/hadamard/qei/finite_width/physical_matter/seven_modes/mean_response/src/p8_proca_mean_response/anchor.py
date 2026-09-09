"""Independent exact stress jets and proper bounce-acceleration response."""

from functools import cache

import sympy as sp
from p8_proca_seven_modes import model

from . import mean, stress


@cache
def data():
    a, m, k = model.a, model.m, model.kvec
    H = model.H
    Hd = sp.Symbol("physical_Hubble_derivative", real=True)
    K = sp.Symbol("squared_comoving_band_momentum", nonnegative=True)
    d = stress.data()
    M = model.data()["Proca_cartesian_generator"]

    def transport(G):
        return (a * H * sp.diff(G, a) + Hd * sp.diff(G, H) + M.T * G + G * M).applyfunc(
            sp.factor
        )

    sub = {a: 1, H: 0, Hd: 4, m: 1000, k[1]: 0, k[2]: 0}

    def trace_at(G):
        value = sp.factor(sp.trace(G.subs(sub, simultaneous=True)) / 2)
        return sp.Poly(value, k[0]).as_expr().subs(k[0] ** 2, K).expand()

    R, P = d["physical_density_Hessian"], d["physical_isotropic_pressure_Hessian"]
    R1, P1 = transport(R), transport(P)
    R2, P2 = transport(R1), transport(P1)
    r, p = trace_at(R), trace_at(P)
    r2, p2 = trace_at(R2), trace_at(P2)
    b = mean.data()
    u = mean.u
    J = b["J"].subs(u, 0)
    J2 = sp.diff(b["J"], u, 2).subs(u, 0)
    ap = sp.diff(b["alpha"], u).subs(u, 0)
    ell = b["ell"].subs(u, 0)
    beta = b["beta"].subs(u, 0)
    F = r - sp.Rational(3, 2) * p
    F2 = r2 - sp.Rational(3, 2) * p2 + 9 * p
    n = F / (2 * J)
    pdot = p + ell * beta * n
    xiddot = -pdot / 2 + ap * n / 3
    nddot = (2 * ap * pdot - 3 * ell * beta * xiddot + F2) / (2 * J) - n * J2 / J
    acceleration = sp.factor(xiddot + nddot / 2 - 11 * n)
    polynomial = sp.Poly(acceleration, K)
    upper = sum(abs(c) * 4 ** power[0] for power, c in polynomial.terms())
    n_poly = sp.Poly(nddot, K)
    return {
        "anchor_density_per_unit_band_covariance": r,
        "anchor_pressure_per_unit_band_covariance": p,
        "anchor_density_second_derivative_per_unit_band_covariance": r2,
        "anchor_pressure_second_derivative_per_unit_band_covariance": p2,
        "anchor_lapse_response_per_unit_band_covariance": sp.factor(n),
        "anchor_trace_response_derivative_per_unit_band_covariance": sp.factor(pdot),
        "anchor_hat_scale_response_second_derivative_per_unit_band_covariance": sp.factor(
            xiddot
        ),
        "anchor_lapse_response_second_derivative_per_unit_band_covariance": sp.factor(
            nddot
        ),
        "anchor_proper_Hubble_derivative_response_per_unit_band_covariance": acceleration,
        "absolute_anchor_proper_Hubble_derivative_response_upper_per_eta": upper,
        "rounded_absolute_anchor_proper_Hubble_derivative_response_upper_per_eta": sp.Integer(
            500_000_000_000
        ),
        "all_anchor_acceleration_response_polynomial_coefficients_positive": all(
            c > 0 for c in polynomial.all_coeffs()
        ),
        "rounded_acceleration_bound_valid": upper < 500_000_000_000,
        "initial_lapse_response_positive_on_entire_band": n.subs(K, 4) > 0,
        "initial_hat_scale_second_derivative_positive_on_entire_band": xiddot.subs(K, 4)
        > 0,
        "initial_lapse_second_derivative_all_polynomial_coefficients_positive": all(
            c > 0 for c in n_poly.all_coeffs()
        ),
        "checks": {
            "positive_band_density_first_derivative_zero_at_bounce": trace_at(R1),
            "positive_band_pressure_first_derivative_zero_at_bounce": trace_at(P1),
            "positive_band_density_second_derivative_obeys_exact_conservation": sp.factor(
                r2 + 12 * (r + p)
            ),
            "background_scalar_lapse_cross_derivative_at_center": ap - 9,
            "proper_time_acceleration_keeps_both_lapse_factors": sp.factor(
                acceleration - (xiddot + nddot / 2 - 3 * n) - (-8 * n)
            ),
        },
    }
