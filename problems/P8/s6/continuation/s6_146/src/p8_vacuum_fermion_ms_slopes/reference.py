"""Complete primitive MS slope forests, exact Laurent coefficients and signs."""

from functools import cache

import sympy as s
from p8_vacuum_fermion_proper_references import conversion as old

from . import masters
from . import tensors as t


@cache
def data():
    e = t.e
    Rs, scalar_rows = masters.integrate(t.scalar_self + t.scalar_vertex)
    Rg, gauge_rows = masters.integrate(-(t.gauge_self + t.gauge_vertex))
    A = s.exp(s.EulerGamma * e) * s.gamma(1 + e)
    cs = 6 * (1 - e) * (1 - 2 * e / 3)
    cg = -12 * (1 - e) * (1 - 2 * e / 3)
    ns = A * A * Rs + A * cs
    ng = A * A * Rg + A * cg
    js = tuple(
        s.simplify(s.diff(ns, e, n).subs(e, 0) / s.factorial(n)) for n in range(3)
    )
    jg = tuple(
        s.simplify(s.diff(ng, e, n).subs(e, 0) / s.factorial(n)) for n in range(3)
    )
    Y, a, Cf, Q = s.symbols("Y a Cf Q", positive=True)
    prev = old.data()
    z = prev["fermion_total_kinetic_UV_counterterm"] / e
    eta = prev["fermion_total_mass_UV_counterterm_over_m"] / e
    nu = prev["total_Yukawa_UV_counterterm_over_y"] / e
    coefficient = s.simplify(2 * (nu - z) - 2 * e * (eta - z))
    ms, mg = (
        t.trace(t.Sk, t.Sk, t.Sk, t.Sl) * 2 + t.trace(t.Sk, t.Sk, t.Sl, t.Sl),
        -(
            2 * t.trace(t.Sk, t.mu, t.Sl, t.mu, t.Sk, t.Sk)
            + t.trace(t.mu, t.Sk, t.Sk, t.mu, t.Sl, t.Sl)
        ),
    )
    masss, _ = masters.integrate(ms)
    massg, _ = masters.integrate(mg)
    vs = 4 / ((1 - e) * (2 * e - 1)) + 1 / (1 - e) ** 2
    vg = -4 / ((1 - e) * (2 * e - 1)) + (2 - 2 * e) / (1 - e) ** 2
    return {
        "scalar_raw_rational_Gamma_factor": Rs,
        "gauge_raw_rational_Gamma_factor_with_vertex_sign": Rg,
        "scalar_master_terms": [
            {"powers": k, "coefficient": v} for k, v in sorted(scalar_rows.items())
        ],
        "gauge_master_terms": [
            {"powers": k, "coefficient": v} for k, v in sorted(gauge_rows.items())
        ],
        "common_MSbar_Gamma_normalization": A,
        "scalar_paired_slope_in_NY_squared_over_Q_squared_units": ns / e**2,
        "gauge_paired_slope_in_NYaCf_over_Q_squared_units": ng / e**2,
        "scalar_double_pole_simple_pole_finite": js,
        "gauge_double_pole_simple_pole_finite": jg,
        "zero_boson_mass_Euclidean_slope": "N/Q^2 [49Y^2/12-37Ya C_F/6] at mu=m. The Minkowski s slope has the opposite sign.",
        "proper_counterterm_scope": "The two marked fermion kinetic/mass insertions and two Yukawa counterterm insertions are paired exactly once. The whole-fermion-cycle local tadpole affects only mass, not slope. The remaining overall MS kinetic poles are subtracted, retaining their finite products.",
        "massless_reference_not_actual_scalar": "Gauge exchange is actually massless. Scalar exchange has mass one; its full nonzero-mass difference is bounded separately.",
        "checks": {
            "same_parent_proper_slope_counterterm": s.simplify(
                coefficient - 3 * (1 - e) * (Y - 2 * a * Cf) / (Q * e)
            ),
            "scalar_double_pole": js[0] - 3,
            "scalar_simple_pole": js[1] + s.Rational(5, 2),
            "scalar_finite": js[2] - s.Rational(49, 12),
            "gauge_double_pole": jg[0] + 6,
            "gauge_simple_pole": jg[1] - 5,
            "gauge_finite": jg[2] + s.Rational(37, 6),
            "vacuum_mass_derivative_scalar_normalization": s.factor(
                masss - (4 - 4 * e) * (3 - 4 * e) * vs
            ),
            "vacuum_mass_derivative_gauge_normalization": s.factor(
                massg - (4 - 4 * e) * (3 - 4 * e) * vg
            ),
            "trace_contracted_pair": s.simplify(t.traceword(("mu", "mu")) - 4 * t.d),
            "trace_contracted_vector": s.simplify(
                t.traceword(("mu", "k", "mu", "l")) - 4 * (2 - t.d) * t.z
            ),
            "trace_four_vectors": s.simplify(
                t.traceword(("k", "l", "k", "l")) - 4 * (2 * t.z**2 - t.x * t.y)
            ),
            "scalar_pi_squared_cancels": Rs.subs(e, 0) / 6 + cs.subs(e, 0) / 12,
            "gauge_pi_squared_cancels": Rg.subs(e, 0) / 6 + cg.subs(e, 0) / 12,
            "actual_positive_Minkowski_coupling_coefficient": s.Rational(37, 6)
            * s.Rational(4, 3)
            - s.Rational(49, 12) * s.Rational(19, 45)
            - s.Rational(3509, 540),
        },
    }
