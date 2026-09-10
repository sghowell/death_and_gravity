"""Full regulated first scalar and fermion slopes, including epsilon terms."""

from functools import cache

import sympy as s
from p8_vacuum_scalar_insertion_ms import inner


@cache
def data():
    x, M, m, g, Y, N, Q, e, ell, u = s.symbols(
        "x M m g Y N Q epsilon ell invariant", positive=True
    )
    a = x * (1 - x)
    Ds = x * M + (1 - x) ** 2
    b = a / Ds
    Delta = m**2 - a * u
    ellx = s.log(m**2 / (m**2 - a))
    J1 = s.Integral(a * ellx, (x, 0, 1))
    J2 = s.Integral(a * ellx**2, (x, 0, 1))
    r0 = g * s.Integral(b, (x, 0, 1)) / Q
    r1 = g * s.Integral(b * (ell - s.log(Ds)), (x, 0, 1)) / Q
    fp0 = 4 * N * Y * (s.Rational(1, 3) - 3 * J1) / Q
    fp1 = 4 * N * Y * (2 * J1 - s.Rational(3, 2) * J2 - s.pi**2 / 24) / Q
    j1, j2 = s.symbols("first_log_moment second_log_moment")
    regulator = (
        (1 + s.zeta(2) * e**2 / 2)
        * (3 - 2 * e)
        * (s.Rational(1, 6) + e * j1 + e**2 * j2 / 2)
    )
    fp = (
        s.series(-4 * N * Y * (regulator - s.Rational(1, 2)) / (Q * e), e, 0, 2)
        .removeO()
        .expand()
    )
    F = a * (1 - 2 * x) * Delta ** (-e)
    trace_derivative = -(Delta ** (-e)) + e * (4 * m**2 - u) * a * Delta ** (-e - 1)
    reduced = -2 * (3 - 2 * e) * a * Delta ** (-e)
    old = inner.data()
    sx = old["symbols"]
    old_first = old["alpha_first_epsilon"].subs(
        {sx["x"]: x, sx["M"]: M, sx["g"]: g, sx["Q"]: Q, sx["ell"]: ell}
    )
    return {
        "scalar_slope_zero": r0,
        "scalar_slope_first_epsilon": r1,
        "fermion_slope_zero": fp0,
        "fermion_slope_first_epsilon": fp1,
        "complete_Phi_k_zero": r0 - fp0,
        "complete_Phi_k_first_epsilon": r1 - fp1,
        "normalized_fermion_MS_slope_dimensional": "-4NY/(Q epsilon) [exp(gamma epsilon) Gamma(1+epsilon) (3-2epsilon) integral a exp(epsilon log(mu^2/Delta)) - 1/2]",
        "checks": {
            "same_full_scalar_first_epsilon_coefficient": old_first - r1,
            "trace_derivative_reduction_by_parameter_IBP": s.factor(
                trace_derivative - reduced + s.diff(F, x)
            ),
            "IBP_lower_endpoint": F.subs(x, 0),
            "IBP_upper_endpoint": F.subs(x, 1),
            "fermion_pole_residue": -4 * N * Y * 3 * s.integrate(a, (x, 0, 1)) / Q
            + 2 * N * Y / Q,
            "fermion_zero_coefficient": s.expand(
                fp.coeff(e, 0) - 4 * N * Y * (s.Rational(1, 3) - 3 * j1) / Q
            ),
            "fermion_first_epsilon_coefficient": s.expand(
                fp.coeff(e, 1)
                - 4 * N * Y * (2 * j1 - s.Rational(3, 2) * j2 - s.pi**2 / 24) / Q
            ),
            "required_second_Gamma_coefficient": s.diff(
                s.exp(s.EulerGamma * e) * s.gamma(1 + e), e, 2
            ).subs(e, 0)
            / 2
            - s.pi**2 / 12,
            "scalar_positive_log_weight_upper": s.integrate((1 - x) / M, (x, 0, 1))
            - 1 / (2 * M),
            "fermion_log_moment_weight": s.integrate(a, (x, 0, 1)) - s.Rational(1, 6),
        },
        "scope": "k_D=r_D-fp_MS,D is holomorphic after the complete first-order MS affine fermion subtraction. Keep its first epsilon coefficient. The scale is fixed while differentiating; mu=m is imposed afterward. The second-loop normalization coefficient remains separate.",
    }
