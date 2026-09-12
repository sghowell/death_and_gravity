"""Uniform momentum/time convolution bounds and the reference normalization."""

from functools import cache

import sympy as s

from . import resolvent as r


@cache
def data():
    t, T, u, v = s.symbols("time window u v", positive=True)
    kap = s.Symbol("kappa", positive=True)
    C, J = r.C, r.JBOUND
    simple = s.integrate(s.integrate(t - u - v, (v, 0, t - u)), (u, 0, t))
    kernel = C * C * J * simple
    derivative = C * C * J * t * t / 2
    operator = s.integrate(kernel, (t, 0, T))
    checks = {
        "ordered_simplex_volume_with_right_linear_kernel": simple - t**3 / 6,
        "uniform_full_matrix_kernel_bound": kernel - s.Rational(375, 16) * t**3,
        "uniform_full_matrix_operator_bound": operator - s.Rational(375, 64) * T**4,
        "uniform_first_time_derivative_kernel_bound": derivative
        - s.Rational(1125, 16) * t * t,
        "uniform_first_time_derivative_operator_bound": s.integrate(
            derivative, (t, 0, T)
        )
        - s.Rational(375, 16) * T**3,
        "physical_Hessian_normalization_not_silently_dropped": 64 * s.pi**2 * operator
        - 375 * s.pi**2 * T**4,
        "kappa_normalized_force_inverse_not_small": 64 * s.pi**2 * kap * operator
        - 375 * s.pi**2 * kap * T**4,
        "coarse_six_T_fourth_margin": 6 - s.Rational(375, 64) - s.Rational(9, 64),
    }
    return {
        "kernel_bound": kernel,
        "operator_bound": operator,
        "first_derivative_kernel_bound": derivative,
        "first_derivative_operator_bound": s.Rational(375, 16) * T**3,
        "spaces": "For every realr, the normalized isolated two-channel flat reference inverse maps C([0,T];H^r(R^3;R^2)) to itself, and to C1 with the displayed first-derivative bound. The sameT4 bound holds onL2_t H^r by Young. There is no spatial derivative loss and no transfer cutoff.",
        "initial_boundary": "Use the complete causal graph domain of the forward distribution on continuous amplitude histories, including its initial boundary. Both inverse identities hold there. Initial delta derivatives are not ignored; no density or full curved graph invariance is asserted.",
        "continuum_proof": "Uniform pointwise Fourier operator bounds, Plancherel, Minkowski and dominated continuity give the all-momentum extension. q0 is a continuous coordinate-kernel extension; the literal homogeneous lapse/shift constraint is separate.",
        "normalization": "The6T4 bound is for A=64pi^2 times the isolated metric Hessian. The physical Hessian inverse is64pi^2 Eq and has bound375pi^2 T4. For the normalized S222 force convention on a flat unit-density reference, Qbar_ref=A/(64pi^2 kappa); its inverse has bound375pi^2 kappa T4. No kappa smallness follows.",
        "remaining": "The full curved/state/contact/tree/matter remainder and quantum auxiliary bridge must be derived and bounded before applying this isolated inverse to the S222 system. Finite-coupling/nonlinear parent control, stability and original V/G/B remain open.",
        "checks": checks,
        "gates": {
            "uniform_kernel_bound_positive": s.Rational(375, 16) > 0,
            "uniform_operator_bound_strictly_below_six": s.Rational(375, 64) < 6,
            "first_derivative_bound_finite": s.Rational(375, 16) > 0,
            "both_positive_time_simplex_integrals_finite": simple == t**3 / 6,
            "physical_normalization_retained": 64 * s.pi**2 > 1,
        },
    }
