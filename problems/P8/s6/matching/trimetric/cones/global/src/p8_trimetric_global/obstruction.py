"""Exact algebraic inputs to the written concavity and tail theorem.

The all-time statement is calculus on one connected regular interval,
not finite sampling. The finite CD test is a sufficient strict three-slice
condition; equality or a failed test is inconclusive, not existence.
"""

from fractions import Fraction

import sympy as sp


def exact_rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction, sp.Rational)):
        raise TypeError("Exact finite rational inputs are required")
    return sp.Rational(value)


def negative_link_three_slice(r, length, scale_error=0):
    """r=R_i(center)/Rcap in (0,1], L>0, 0<=relative scale error<1.

    CD a=(1+(T/tau)²)² at the three physical times -Ltau,0,+Ltau.
    Optional pointwise relative errors affect all three physical scale
    values. No frame-derivative or equal-parent-proper-time premise enters.
    """
    r, length, error = map(exact_rational, (r, length, scale_error))
    if not (0 < r <= 1 and length > 0 and 0 <= error < 1):
        raise ValueError("Require 0<r<=1, L>0 and 0<=scale_error<1")
    growth = (1+length**2)**2
    margin = sp.factor(r*growth*(1-error)-(1+error))
    threshold = sp.factor((r*growth-1)/(r*growth+1))
    return {"r": r, "L": length, "scale_error": error, "CD_scale_growth": growth,
            "strict_chord_margin": margin, "necessary_error_threshold_when_positive": threshold,
            "status": "EXCLUDED_BY_NEGATIVE_LINK_CONCAVITY" if margin > 0 else "INCONCLUSIVE"}


def checks():
    g, magnitude, ratio, h = sp.symbols("G s R H", positive=True)
    # For p=-s, the constraint implies b=s R³+3G H²/2>0.
    b_on_shell = magnitude*ratio**3+3*g*h**2/2
    cap_cube = b_on_shell/magnitude
    t = sp.Symbol("t", real=True)
    a, lapse = sp.Function("a", positive=True)(t), sp.Function("n", positive=True)(t)
    hi = sp.diff(a, t)/(lapse*a)
    derivative = sp.diff(hi, t)/lapse
    expanded = sp.diff(a, t, 2)/(lapse**2*a)-sp.diff(a, t)**2/(lapse**2*a**2)-sp.diff(a, t)*sp.diff(lapse, t)/(lapse**3*a)
    left, right = sp.symbols("left_duration right_duration", positive=True)
    ym, y0, yp = sp.symbols("log_a_minus log_a_zero log_a_plus", real=True)
    slopes = (y0-ym)/left-(yp-y0)/right
    chord_excess = right*ym+left*yp-(left+right)*y0
    r, length, error = sp.symbols("r L epsilon_a", real=True)
    growth = (1+length**2)**2
    threshold = (r*growth-1)/(r*growth+1)
    pe, pf, de, df, nh = sp.symbols("P_g P_f deficit_g deficit_f n_h", nonnegative=True)
    # c_i=1-deficit_i; positive weights and nonnegative deficits force zero.
    weighted = pe*((1-de)-1)+pf*((1-df)-1)-nh/2
    return {
        "negative_link_endpoint_forced_positive": sp.expand(3*g*h**2-2*(-magnitude*ratio**3+b_on_shell)),
        "negative_link_cap_cube_gap": sp.expand(cap_cube-ratio**3-3*g*h**2/(2*magnitude)),
        "proper_clock_log_scale_concavity": sp.simplify(derivative-expanded),
        "arbitrary_positive_duration_secant_chord_identity": sp.expand(left*right*slopes+chord_excess),
        "nonnegative_links_NEC_deficit_sum": sp.expand(weighted+pe*de+pf*df+nh/2),
        "strict_three_slice_geometry_error_threshold": sp.factor((r*growth*(1-error)-(1+error)).subs(error, threshold)),
    }


def controls():
    return {"exact_CD_r_half_L1": negative_link_three_slice(sp.Rational(1, 2), 1),
            "strict_equality_inconclusive": negative_link_three_slice(sp.Rational(1, 4), 1),
            "smaller_window_inconclusive": negative_link_three_slice(sp.Rational(1, 10), 1),
            "relative_geometry_error_control": negative_link_three_slice(sp.Rational(1, 2), 1, sp.Rational(1, 10)),
            "geometry_error_threshold_inconclusive": negative_link_three_slice(sp.Rational(1, 2), 1, sp.Rational(1, 3))}
