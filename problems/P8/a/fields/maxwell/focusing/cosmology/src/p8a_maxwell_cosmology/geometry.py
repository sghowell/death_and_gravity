"""Actual short-history bounds for a normalized C3 tube of smooth metrics."""

import sympy as sp
from p8a_maxwell.domain import nonnegative, rational

RATIO = sp.Rational(1, 100)
CONTRACTION = sp.Rational(19, 10)
REFERENCE_CAPS = tuple(map(sp.Integer, (2, 8, 64, 768)))
EXPANDED_CAPS = (sp.Rational(21, 10), sp.Integer(9), sp.Integer(70), sp.Integer(800))
FUTURE_CAPS = tuple(map(sp.Integer, (4, 128, 16384, 1048576)))
ERROR_CAPS = (sp.Rational(1, 100), sp.Rational(1, 2), sp.Integer(3), sp.Integer(16))


def reference_jets(power, clock):
    """Formal dimensionless jets: H_p*tau=-p/(p/2-x), x=s/tau."""
    p, x = map(sp.sympify, (power, clock))
    return [-sp.factorial(j)*p/(p/2-x)**(j+1) for j in range(4)]


def history_point(power, clock):
    """Exact evaluation only in p in[1/2,2/3], x in[-1/100,0]."""
    p, x = map(rational, (power, clock))
    if not sp.Rational(1, 2) <= p <= sp.Rational(2, 3) or not -RATIO <= x <= 0:
        raise ValueError("the point is outside the certified short-history rectangle")
    return reference_jets(p, x)


def neighborhood(error_caps=ERROR_CAPS, *, anchored=False):
    """Bounds for actual smooth H with the displayed scaled C3 errors.

The optional anchored slice additionally imposes H(0)*tau=-2. The full
unanchored tube does NOT preserve an exact observed H0*tau=2 identity.
The logarithmic scale-factor bound compares a(s)/a(0) with the likewise
normalized reference, not arbitrary absolute normalizations of a.
    """
    if not isinstance(anchored, bool):
        raise TypeError("anchored must be boolean")
    if not isinstance(error_caps, (tuple, list)) or len(error_caps) != 4:
        raise TypeError("four positive exact dimensionless jet-error bounds are required")
    errors = tuple(nonnegative(value, positive=True) for value in error_caps)
    gaps = [cap-base-error for cap, base, error in zip(EXPANDED_CAPS, REFERENCE_CAPS, errors, strict=True)]
    contraction_lower = sp.Rational(25, 13)-errors[0]
    if min(*gaps, contraction_lower-CONTRACTION) <= 0:
        raise ValueError("the proposed tube has no strict cap/history margin")
    observer = [sp.Integer(2), sp.Integer(2)] if anchored else [2-errors[0], 2+errors[0]]
    return {"error_caps": list(errors), "reference_jet_caps": list(REFERENCE_CAPS),
            "expanded_past_caps": list(EXPANDED_CAPS), "strict_cap_margins": gaps,
            "baseline_contraction_lower": sp.Rational(25, 13),
            "actual_contraction_lower": contraction_lower,
            "strict_contraction_margin": contraction_lower-CONTRACTION,
            "observed_H0_times_tau_range": observer, "anchored_at_observer": anchored,
            "log_scale_factor_error_upper": RATIO*errors[0],
            "scale_factor_comparison_is_normalized_at_observer": True,
            "future_relative_caps_verified_by_tube": False,
            "actual_quantum_SEE_solution_asserted": False}


def identities():
    x = sp.Symbol("x", real=True)
    p = sp.Symbol("p", positive=True)
    h = -p/(p/2-x)
    data = reference_jets(p, x)
    result = {f"proper_clock_jet_{j}": sp.simplify(sp.diff(h, x, j)-data[j]) for j in range(4)}
    result.update({f"observer_jet_{j}": sp.simplify(-data[j].subs(x, 0)
                                                   -sp.factorial(j)*2**(j+1)/p**j) for j in range(4)})
    result.update({"scale_factor_time_orientation": sp.simplify(
                       sp.diff(sp.log((p/2-x)/(p/2)), x)*p-h),
                   "reference_Hstar_times_tau": sp.simplify(-h.subs(x, 0)-2),
                   "minimum_history_at_radiation_edge": sp.simplify(-h.subs({p: sp.Rational(1, 2), x: -RATIO})-sp.Rational(25, 13)),
                   "history_increases_with_power": sp.simplify(sp.diff(p/(p/2+RATIO), p)-RATIO/(p/2+RATIO)**2)})
    return result


def calibration():
    return {"power_range": [sp.Rational(1, 2), sp.Rational(2, 3)],
            "dimensionless_history_interval": [-RATIO, sp.Integer(0)],
            "ratio": RATIO, "contraction": CONTRACTION,
            "reference_Hstar_times_tau": sp.Integer(2),
            "reference_age_over_tau_range": [sp.Rational(1, 4), sp.Rational(1, 3)],
            "unanchored_C3_tube": neighborhood(), "anchored_C3_slice": neighborhood(anchored=True),
            "this_tube_Ricci_upper_times_tau_squared": 3*(-sp.Rational(60000, 10609)+ERROR_CAPS[1]+(2+ERROR_CAPS[0])**2),
            "this_tube_demonstrates_initial_SEC_violation": False,
            "power_law_geometries_are_quantum_SEE_solutions": False}
