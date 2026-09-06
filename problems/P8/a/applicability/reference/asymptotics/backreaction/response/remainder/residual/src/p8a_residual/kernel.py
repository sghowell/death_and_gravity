"""Correct-clock logarithmic history and derived finite-family norm bounds."""

import sympy as sp
from p8a_remainder import geometry as prior_geometry
from p8a_remainder import kernel as prior_kernel
from p8a_remainder.mode_bounds import exact_nonnegative


def retarded(a, eta, eta_start, *, span, length):
    """Use the full scale factor and actual eta; no inverse-clock expansion."""
    return prior_kernel.retarded(a, eta, eta_start, span=span, length=length)


def norm_bounds(potential_bounds, active_span, local_log_abs, hubble_cap, hubble_derivative_cap):
    """Conditional generic K0..2 bounds; flat past U jets are required.

    The named calibration below derives every input from the same A.7 metric.
    This generic assembly does not independently validate a supplied geometry.
    """
    if len(potential_bounds) != 4:
        raise ValueError("Need four actual conformal potential norms through U'''")
    b0, b1, b2, b3 = map(exact_nonnegative, potential_bounds)
    t, local, h, hp = map(exact_nonnegative,
                          (active_span, local_log_abs, hubble_cap, hubble_derivative_cap))
    return [t*b1+local*b0,
            t*b2+local*b1+h*b0,
            t*b3+local*b2+2*h*b1+hp*b0]


def calibration():
    data = prior_geometry.calibration()
    b = list(map(sp.Rational, data["potential_bounds"]))
    cap = sp.Rational(data["delta_bar"])
    # Strip delta/eta_star^(j+2), using Hc<=2/eta_star and
    # |Hc'|<=U_bound+Hc^2. All caps hold on the active history, not just target.
    return norm_bounds(b[:4], 3, 3, 2, 4+cap*b[0])


def named_length(normalization, eta_star):
    normalization = exact_nonnegative(normalization, positive=True)
    eta_star = exact_nonnegative(eta_star, positive=True)
    return 2*sp.sqrt(2)*normalization*eta_star**2


def elementary_margins():
    return {"local_log_argument_fourth_power_lower_minus_one": sp.Rational(3, 2)**4/2-1,
            "exp_two_series_lower_minus_log_argument_upper": sp.Integer(1)+2+2-sp.Rational(9, 2),
            "local_factor_cap_margin": sp.Integer(3)-(2+sp.Rational(5, 6))}


def identities():
    eta, span, duration = sp.symbols("eta T D", positive=True)
    a = sp.Function("a", positive=True)(eta)
    local = sp.log(a*span)+sp.Rational(5, 6)
    h = sp.diff(a, eta)/a
    log_integral = duration*(1+sp.log(span/duration))
    return {"local_log_first_derivative": sp.simplify(sp.diff(local, eta)-h),
            "local_log_second_derivative": sp.simplify(sp.diff(local, eta, 2)-sp.diff(h, eta)),
            "absolute_log_integral_monotonic_derivative": sp.simplify(
                sp.diff(log_integral, duration)-sp.log(span/duration)),
            "absolute_log_integral_at_full_span": log_integral.subs(duration, span)-span,
            "absolute_log_integral_at_zero_span_limit": sp.limit(log_integral, duration, 0, dir="+")}
