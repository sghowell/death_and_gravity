"""Massless-exchange MS mass anchors from the full dimensionally regulated forest."""

from functools import cache

import sympy as s
from p8_vacuum_fermion_ms_slopes import masters
from p8_vacuum_fermion_ms_slopes import tensors as t


@cache
def data():
    e = t.e
    A = s.exp(s.EulerGamma * e) * s.gamma(1 + e)
    vs = 4 / ((1 - e) * (2 * e - 1)) + 1 / (1 - e) ** 2
    vg = -4 / ((1 - e) * (2 * e - 1)) + (2 - 2 * e) / (1 - e) ** 2
    Rs = s.factor((4 - 4 * e) * (3 - 4 * e) * vs)
    Rg = s.factor((4 - 4 * e) * (3 - 4 * e) * vg)
    H = 4 - 2 / (e - 1)
    cs = 2 * (6 - 3 * e) * H
    cg = -4 * (6 - 3 * e) * H
    ns = A * A * Rs + A * cs
    ng = A * A * Rg + A * cg
    js = tuple(
        s.simplify(s.diff(ns, e, n).subs(e, 0) / s.factorial(n)) for n in range(3)
    )
    jg = tuple(
        s.simplify(s.diff(ng, e, n).subs(e, 0) / s.factorial(n)) for n in range(3)
    )
    raw_s, _ = masters.integrate(
        2 * t.trace(t.Sk, t.Sk, t.Sk, t.Sl) + t.trace(t.Sk, t.Sk, t.Sl, t.Sl)
    )
    raw_g, _ = masters.integrate(
        -(
            2 * t.trace(t.Sk, t.mu, t.Sl, t.mu, t.Sk, t.Sk)
            + t.trace(t.mu, t.Sk, t.Sk, t.mu, t.Sl, t.Sl)
        )
    )
    return {
        "scalar_raw_mass_rational_Gamma_factor": Rs,
        "gauge_raw_mass_rational_Gamma_factor": Rg,
        "scalar_paired_mass_reference_in_NY_squared_m_squared_over_Q_squared_units": ns
        / e**2,
        "gauge_paired_mass_reference_in_NYaCf_m_squared_over_Q_squared_units": ng
        / e**2,
        "scalar_double_pole_simple_pole_finite": js,
        "gauge_double_pole_simple_pole_finite": jg,
        "finite_massless_reference": "f_MS(0)=N m^2/Q^2 [-56Y^2+40Ya C_F]. Scalar exchange is a reference with b=0; gauge exchange is actually massless.",
        "checks": {
            "scalar_raw_mass_tensor_matches_vacuum_derivative": s.factor(raw_s - Rs),
            "gauge_raw_mass_tensor_matches_vacuum_derivative": s.factor(raw_g - Rg),
            "scalar_double_pole": js[0] - 36,
            "scalar_simple_pole": js[1] + 48,
            "scalar_finite": js[2] + 56,
            "gauge_double_pole": jg[0] + 72,
            "gauge_simple_pole": jg[1] - 24,
            "gauge_finite": jg[2] - 40,
            "scalar_pi_squared_cancels": Rs.subs(e, 0) / 6 + cs.subs(e, 0) / 12,
            "gauge_pi_squared_cancels": Rg.subs(e, 0) / 6 + cg.subs(e, 0) / 12,
            "fixed_mu_vacuum_mass_homogeneity": (4 - 4 * e) * (3 - 4 * e)
            - (12 - 28 * e + 16 * e * e),
            "actual_massless_reference_coupling_dictionary": 40 * s.Rational(4, 3)
            - 56 * s.Rational(19, 45)
            - s.Rational(1336, 45),
        },
    }
