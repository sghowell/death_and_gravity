"""Uniform on-shell fermion insertion bound without a momentum expansion."""

from fractions import Fraction
from functools import cache

import sympy as sp
from p8_vacuum_fermion_local_matching import twopoint

N = sp.Integer(6)


def exact(value):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction, sp.Rational)):
        raise TypeError("Require a finite exact rational")
    return sp.Rational(value)


def enclosure(mass, yukawa_upper, Q_lower, sigma=1):
    m, Y, Q, sig = map(exact, (mass, yukawa_upper, Q_lower, sigma))
    if m <= 0 or Y < 0 or Q <= 0 or sig < 1 or sig >= 4 * m * m:
        raise ValueError("Require m>0, Y>=0, Q>0 and 1<=sigma<4m^2")
    gap = m * m - sig / 4
    B = N * Y / Q * (1 / (6 * gap) + m * m / (30 * gap * gap))
    return {
        "fermion_mass": m,
        "Y_upper": Y,
        "Q_lower": Q,
        "halfplane_real_part_upper": sig,
        "parameter_gap_lower": gap,
        "unscaled_OS_insertion_uniform_upper": B,
        "scope": "One-loop f_R(s)/(s-1)^2 on Re(s)<=sigma; no bound on a complete outer two-loop integral.",
    }


@cache
def data():
    m, Y, Q = sp.symbols("m Y Q", positive=True)
    s, A, sigma = sp.symbols("s A sigma")
    x, u = sp.symbols("x u", real=True)
    D = m * m - A * s
    kernel = -(4 * m * m - s) * sp.log(1 - A * s / (m * m))
    second = -A / D - m * m * A * (1 - 4 * A) / D**2
    gap = m * m - sigma / 4
    bound = N * Y / Q * (1 / (6 * gap) + m * m / (30 * gap**2))
    old = twopoint.data()["second_derivative_parameter_integrand"]
    replacement = {
        v: {
            "positive_fermion_mass": m,
            "Feynman_weight_A": A,
            "continued_invariant_s": s,
        }[v.name]
        for v in old.free_symbols
    }
    c0, c1 = sp.symbols("whole_local_constant whole_local_slope")
    affine = c0 + c1 * s
    sub = affine - affine.subs(s, 1) - (s - 1) * sp.diff(affine, s).subs(s, 1)
    f = sp.Function("complete_fermion_kernel")
    fR = f(s) - f(1) - (s - 1) * sp.diff(f(s), s).subs(s, 1)
    checks = {
        "independent_second_derivative": sp.factor(sp.diff(kernel, s, 2) - second),
        "same_frozen_whole_fermion_kernel": sp.factor(
            old.xreplace(replacement) - second
        ),
        "parameter_A_integral": sp.integrate(x * (1 - x), (x, 0, 1))
        - sp.Rational(1, 6),
        "positive_second_weight_integral": sp.integrate(
            x * (1 - x) * (1 - 4 * x * (1 - x)), (x, 0, 1)
        )
        - sp.Rational(1, 30),
        "second_weight_square": 1 - 4 * x * (1 - x) - (2 * x - 1) ** 2,
        "whole_affine_local_reference_cancels": sp.expand(sub),
        "finite_mass_anchor_fixed": fR.subs(s, 1),
        "finite_residue_anchor_fixed": sp.diff(fR, s).subs(s, 1),
        "Taylor_integral_half": sp.integrate(1 - u, (u, 0, 1)) - sp.Rational(1, 2),
        "bound_single_rational_form": sp.factor(
            bound - N * Y * (6 * m * m - 5 * sigma / 4) / (30 * Q * gap**2)
        ),
        "segment_real_part_margin": sp.expand(
            sigma - (1 + u * (s - 1)) - (1 - u) * (sigma - 1) - u * (sigma - s)
        ),
        "parameter_gap_margin": sp.expand(
            m * m - A * s - gap - A * (sigma - s) - sigma * (sp.Rational(1, 4) - A)
        ),
    }
    return {
        "symbols": {"m": m, "Y": Y, "Q": Q, "s": s, "A": A, "sigma": sigma},
        "same_full_finite_kernel_before_OS_subtraction": kernel,
        "second_derivative_resolvent_decomposition": second,
        "on_shell_subtraction": fR,
        "whole_affine_subtraction": sub,
        "unscaled_uniform_insertion_bound": bound,
        "integral_representation": "f_R(s)/(s-1)^2=(2NY/Q) integral dx integral_0^1 (1-u) g_x''(1+u(s-1)) du, extended removably at s=1.",
        "domain": "Re(s)<=sigma, 1<=sigma<4m^2. A=x(1-x) in [0,1/4]. No restriction on Im(s) or |s|; no expansion in momentum/m.",
        "scope": "An exact one-loop insertion bound after pairing the complete affine counterterm and finite on-shell anchors. It is not a full two-loop bound.",
        "checks": checks,
    }
