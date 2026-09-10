"""Dimension-symbolic vacuum tensors and their complete proper fermion MS forest."""

from functools import cache

import sympy as s
from p8_vacuum_fermion_ms_slopes import masters
from p8_vacuum_fermion_ms_slopes import tensors as t
from p8_vacuum_fermion_proper_references import conversion


@cache
def data():
    e = t.e
    A = s.exp(s.EulerGamma * e) * s.gamma(1 + e)
    vs = 4 / ((1 - e) * (2 * e - 1)) + 1 / (1 - e) ** 2
    vg = -4 / ((1 - e) * (2 * e - 1)) + 2 / (1 - e)
    cs, cg = 6 / (1 - e), -12 / (1 - e)
    ns, ng = A * A * vs + A * cs, A * A * vg + A * cg
    js = tuple(
        s.simplify(s.diff(ns, e, n).subs(e, 0) / s.factorial(n)) for n in range(3)
    )
    jg = tuple(
        s.simplify(s.diff(ng, e, n).subs(e, 0) / s.factorial(n)) for n in range(3)
    )
    rs, _ = masters.integrate(s.Rational(1, 2) * t.trace(t.Sk, t.Sl))
    rg, _ = masters.integrate(-s.Rational(1, 2) * t.trace(t.mu, t.Sk, t.mu, t.Sl))
    p = conversion.data()
    Y, a, Cf, Q = s.symbols("Y a Cf Q", positive=True)
    effective = (
        p["fermion_total_mass_UV_counterterm_over_m"]
        - p["fermion_total_kinetic_UV_counterterm"]
    )
    fullct = s.factor(2 * (4 - 2 * e) * effective * Q / ((e - 1) * (e - 2)))
    return {
        "symbols": {"epsilon": e},
        "scalar_raw_vacuum_Gamma_rational_factor": vs,
        "gauge_raw_vacuum_Gamma_rational_factor": vg,
        "scalar_proper_fermion_CT_rational_factor": cs,
        "gauge_proper_fermion_CT_rational_factor": cg,
        "scalar_paired_massless_reference_in_NY_m_four_over_Q_squared_units": ns / e**2,
        "gauge_paired_massless_reference_in_NaCf_m_four_over_Q_squared_units": ng
        / e**2,
        "scalar_double_pole_simple_pole_finite": js,
        "gauge_double_pole_simple_pole_finite": jg,
        "one_loop_vacuum": "V1,F=2N e^(gamma epsilon)mu^(2epsilon)Gamma(epsilon-2)m^(4-2epsilon)/Q. Its proper fermion insertion is (eta-z)m partial_m V1 at fixed mu. The constant kinetic Jacobian is a scaleless trace.",
        "scope": "The scalar massless chord is an auxiliary anchor. The gauge chord is actually massless. N=6 for scalar Yukawa exchange; N=42 for the gauge vacuum.",
        "checks": {
            "scalar_direct_vacuum_tensor": s.factor(rs - vs),
            "gauge_direct_vacuum_tensor": s.factor(rg - vg),
            "complete_fermion_mass_kinetic_CT_dictionary": s.factor(
                fullct - Y * cs - a * Cf * cg
            ),
            "scalar_double_pole": js[0] - 3,
            "scalar_simple_pole": js[1] + 4,
            "scalar_finite": js[2] + 19,
            "gauge_double_pole": jg[0] + 6,
            "gauge_simple_pole": jg[1] - 2,
            "gauge_finite": jg[2] - 18,
            "scalar_pi_squared_cancellation": vs.subs(e, 0) / 6 + cs.subs(e, 0) / 12,
            "gauge_pi_squared_cancellation": vg.subs(e, 0) / 6 + cg.subs(e, 0) / 12,
            "scalar_combined_rational_factor": s.factor(
                vs - (3 - 2 * e) / ((1 - e) ** 2 * (2 * e - 1))
            ),
        },
    }
