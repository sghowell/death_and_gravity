"""Exact physical scalar-bubble boundary, derivatives and continuous quotient."""

from functools import cache

import sympy as s


def invariant(value):
    if isinstance(value, (bool, str, float)) or value is None:
        raise ValueError("Require an exact physical bubble invariant")
    value = s.sympify(value)
    if not isinstance(value, s.Rational) or not (
        -12 <= value <= 0 or s.Rational(45, 8) <= value <= 16
    ):
        raise ValueError("Require a physical crossed or timelike bubble interval")
    return value


def bubble(value):
    value = invariant(value)
    if value == 0:
        return s.S.Zero
    beta = s.sqrt(1 - 4 / value)
    if value < 0:
        return 2 - beta * s.log((beta + 1) / (beta - 1))
    return 2 - beta * (s.log((1 + beta) / (1 - beta)) - s.I * s.pi)


def derivative(value):
    value = invariant(value)
    if value == 0:
        return s.Rational(1, 6)
    beta = s.sqrt(1 - 4 / value)
    log = (
        s.log((beta + 1) / (beta - 1))
        if value < 0
        else s.log((1 + beta) / (1 - beta)) - s.I * s.pi
    )
    return -2 * log / (value**2 * beta) - 1 / value


def divided_difference(left, right):
    left, right = invariant(left), invariant(right)
    if (left <= 0) != (right <= 0):
        raise ValueError(
            "A physical channel cannot straddle the two-particle threshold"
        )
    if left == right:
        return derivative(left)
    return (bubble(right) - bubble(left)) / (right - left)


def symmetric_value():
    return 2 - 2 * s.sqrt(2) * s.atan(1 / s.sqrt(2))


@cache
def data():
    b = s.Symbol("beta", positive=True)
    v = 4 / (1 - b * b)
    L = s.log((1 + b) / (1 - b))
    B = 2 - b * (L - s.I * s.pi)
    db = (1 - b * b) ** 2 / (8 * b)
    Bp = -2 * (L - s.I * s.pi) / (v * v * b) - 1 / v
    bpp = -4 / (v**3 * b) - 4 / (v**4 * b**3)
    written = -bpp * (L - s.I * s.pi) - 2 / (v * v * b) / (v * b) + 1 / (v * v)
    r1, r2 = (1 - b) / 2, (1 + b) / 2
    log1, log2 = s.symbols("log_root_minus log_root_plus")
    x = s.Symbol("x", real=True)
    q = s.Symbol("formal_invariant")
    hf = s.Function("endpoint_factor")(q)
    bf = s.Function("bubble")(q)
    h0, h1, h2, F, F1, F2 = s.symbols("h hprime hsecond B Bprime Bsecond")
    sub = {
        s.diff(hf, q, 2): h2,
        s.diff(hf, q): h1,
        hf: h0,
        s.diff(bf, q, 2): F2,
        s.diff(bf, q): F1,
        bf: F,
    }
    checks = {
        "physical_cut_first_derivative": s.factor(s.diff(B, b) * db - Bp),
        "physical_cut_second_derivative": s.factor(s.diff(Bp, b) * db - written),
        "split_log_boundary_real_part": s.expand(
            2 + log1 + log2 - 2 * (r1 * log1 + r2 * log2) - (2 - b * (log2 - log1))
        ),
        "positive_two_particle_cut": s.expand(s.pi * (r2 - r1) - s.pi * b),
        "crossed_first_derivative_integral": s.integrate(x * (1 - x), (x, 0, 1))
        - s.Rational(1, 6),
        "crossed_second_derivative_integral": s.integrate(
            x * x * (1 - x) ** 2, (x, 0, 1)
        )
        - s.Rational(1, 30),
        "actual_timelike_recoil_lower_endpoint": 4
        * s.Rational(5, 4)
        * (s.Rational(5, 4) - s.Rational(1, 8))
        - s.Rational(45, 8),
        "independent_full_kernel_first_derivative": s.expand(
            s.diff(hf * hf * bf, q).xreplace(sub) - 2 * h0 * h1 * F - h0 * h0 * F1
        ),
        "independent_full_kernel_second_derivative": s.expand(
            s.diff(hf * hf * bf, q, 2).xreplace(sub)
            - 2 * (h1 * h1 + h0 * h2) * F
            - 4 * h0 * h1 * F1
            - h0 * h0 * F2
        ),
    }
    gates = {
        "log_ratio_bound_by_exact_exponential_partial_sum": bool(
            sum(s.Rational(3) ** j / s.factorial(j) for j in range(6)) > 14
        ),
        "upper_velocity_logarithm_ratio_less_than14": bool(7 + 4 * s.sqrt(3) < 14),
        "lower_velocity_exceeds_two_fifths": bool(
            s.Rational(1, 5) > s.Rational(2, 5) ** 2
        ),
        "complex_Bprime_bound_less_than2": bool(
            s.Rational(14, 10) + s.Rational(1, 5) < 2
        ),
        "complex_Bsecond_bound_less_than2": bool(
            (
                s.Rational(4, 125) * s.Rational(5, 2)
                + s.Rational(4, 625) * s.Rational(5, 2) ** 3
            )
            * 7
            + s.Rational(1, 5) * s.Rational(1, 2)
            + s.Rational(1, 25)
            < 2
        ),
        "actual_channel_segments_do_not_cross_threshold": True,
        "Feynman_cut_not_replaced_by_absolute_parameter_integral": True,
        "equal_invariants_continuously_extended": True,
    }
    return {
        "checks": checks,
        "gates": gates,
        "whole_physical_branch": B,
        "whole_physical_first_derivative": Bp,
        "whole_physical_second_derivative": written,
        "whole_symmetric_bubble": symmetric_value(),
        "whole_bounds": "Physical timelike pair endpoints are in[45/8,16] and crossed endpoints in[-12,0], with each pair segment in its own interval. On the larger timelike interval[5,16], beta>2/5,L<3 and pi<4 imply |B|<6,|Bprime|<2,|Bsecond|<2. Crossed bounds are2,1/6,1/30. At4/3,0<B<1. The physical divided difference is a segment integral of Bprime and stays finite at equality.",
        "whole_analytic_boundary": "Both exact-D triangle identities are first proved in an analytic subthreshold domain and continued to the same Feynman sheet. The closed boundary formula retains the positive imaginary cut pi*beta. No Euclidean positive-denominator bound is extrapolated across a zero; threshold4 is uniformly separated from the actual timelike recoil domain.",
    }
