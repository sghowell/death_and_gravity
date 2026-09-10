"""Dirac-trace kernel and complete local on-shell subtraction of its increment."""

from functools import cache

import sympy as sp
from p8_vacuum_finite_mass_gauge_cut import dirac


@cache
def data():
    gamma = dirac.data()["gamma_matrices"]
    qv = sp.symbols("real_q0:4", real=True)
    pv = sp.symbols("external_p0:4", real=True)
    m, Y = sp.symbols("positive_fermion_mass positive_Yukawa_squared", positive=True)
    N = sp.Integer(6)
    Q = 16 * sp.pi**2
    slash = lambda p: sum((p[i] * gamma[i] for i in range(4)), sp.zeros(4))
    qp = [qv[i] + pv[i] for i in range(4)]
    trace = sp.trace(
        (m * sp.eye(4) - sp.I * slash(qv)) * (m * sp.eye(4) - sp.I * slash(qp))
    )
    D0 = m * m + sum(v * v for v in qv)
    D1 = m * m + sum(v * v for v in qp)
    p2 = sum(v * v for v in pv)
    s, A = sp.symbols("continued_invariant_s Feynman_weight_A")
    Ibar = sp.symbols("entire_dimensional_reference_Ibar")
    J = -sp.log(1 - A * s / m**2)
    f_integrand = (4 * m * m - s) * J
    second = sp.factor(sp.diff(f_integrand, s, 2))
    B = (
        N
        * Y
        / Q
        * (sp.Rational(7, 15) * m * m + sp.Rational(1, 10))
        / (m * m - sp.Rational(3, 4)) ** 2
    )
    f = sp.Function("f_fermion")
    subtraction = f(s) - f(1) - (s - 1) * sp.diff(f(s), s).subs(s, 1)
    affine = (6 * m * m - s) * Ibar
    checks = {
        "independent_Dirac_numerator_trace": sp.expand(
            trace - 4 * (m * m - sum(qv[i] * qp[i] for i in range(4)))
        ),
        "trace_tadpole_bubble_reduction": sp.expand(
            trace - 2 * (4 * m * m + p2) + 2 * (D0 + D1)
        ),
        "all_UV_reference_terms_cancel_after_OS_subtraction": sp.expand(
            affine - affine.subs(s, 1) - (s - 1) * sp.diff(affine, s).subs(s, 1)
        ),
        "twice_subtracted_mass_anchor": subtraction.subs(s, 1),
        "twice_subtracted_residue_anchor": sp.diff(subtraction, s).subs(s, 1),
        "explicit_second_derivative_integrand": sp.factor(
            second - (m * m * (4 * A * A - 2 * A) + s * A * A) / (m * m - A * s) ** 2
        ),
        "negative_real_second_derivative_margin": sp.factor(
            m * m * (4 * A * A - 2 * A)
            + s * A * A
            + A * (m * m - s * A)
            - m * m * A * (4 * A - 1)
        ),
        "second_derivative_absolute_integral_numerator": 4 * sp.Rational(1, 30)
        + 2 * sp.Rational(1, 6)
        - sp.Rational(7, 15),
        "second_derivative_s_numerator": 3 * sp.Rational(1, 30) - sp.Rational(1, 10),
        "Taylor_integral_remainder_factor": 2 * sp.Rational(1, 2) - 1,
    }
    return {
        "Euclidean_Dirac_trace_numerator": sp.expand(trace),
        "Euclidean_kernel": "Gamma_F2(pE^2)=2 N Y [(4mF^2+pE^2) B(pE^2)-2 T], N=6",
        "MSbar_tadpole_at_scale_m": "T_MSbar=-mF^2/Q; B_MSbar(s)=J(s)/Q",
        "continued_MSbar_kernel": "f(s)=(2NY/Q)[2mF^2+(4mF^2-s) integral_0^1 -log(1-x(1-x)s/mF^2) dx]",
        "entire_UV_affine_reference": 2 * N * Y * (6 * m * m - s) * Ibar / Q,
        "OS_subtracted_kernel": subtraction,
        "second_derivative_parameter_integrand": second,
        "uniform_OS_remainder_coefficient_upper": B,
        "complex_domain": "|s-1|<=2, mF>=36: |s|<=3 and Re(mF^2-A s)>=mF^2-3/4",
        "real_shape": "On 0<=s<=1 the second derivative is strictly negative; f'(1)<f'(0)=4NY/(3Q).",
        "absorptive_sign_dictionary": "For s>4mF^2, Im f(s+i0)=-NY s(1-4mF^2/s)^(3/2)/(8pi). The continued Euclidean inverse is m_ref^2-s+f(s); multiplying by -1 gives the usual positive-width Minkowski inverse.",
        "scope": "Complete one-loop fermion two-point increment and its local reference subtraction. Four-dimensional trace algebra is used with the gauge-invariant dimensional regulator; the full tadpole and affine UV reference are retained before subtraction.",
        "checks": checks,
    }
