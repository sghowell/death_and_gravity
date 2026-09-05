"""Exact interval bookkeeping for a supplied conditional dispersion bound.

This does not derive a gravitational remainder or check theorem hypotheses.
"""

import sympy as sp


def _rational(value):
    result = sp.sympify(value)
    if result.is_Rational is not True:
        raise ValueError("An exact rational input is required")
    return result


def necessary_bound(b2, *, error, gravity_lower_radius):
    """Given R_grav>=-Delta, test only the necessary b2+Delta>=0.

    error bounds the total error in b2, including any omitted low-energy
    pieces. None is an unknown bound, never an assumed zero correction.
    Compatibility with this weak necessary inequality is not UV admissibility.
    """
    b2 = _rational(b2)
    bounds = {"error": error, "gravity_lower_radius": gravity_lower_radius}
    for name, value in bounds.items():
        if value is not None:
            bounds[name] = _rational(value)
            if bounds[name] < 0:
                raise ValueError("Error and gravitational radii must be nonnegative")
    if any(value is None for value in bounds.values()):
        return {"status": "UNTESTED_MISSING_BOUND", "necessary_margin": None}
    error, delta = bounds["error"], bounds["gravity_lower_radius"]
    lower, upper = b2-error+delta, b2+error+delta
    if upper < 0:
        status = "EXCLUDED_BY_SUPPLIED_CONDITIONAL_BOUND"
    elif lower > 0:
        status = "NOT_EXCLUDED_BY_THIS_NECESSARY_BOUND"
    else:
        status = "BOUNDARY_OR_UNRESOLVED"
    return {"status": status, "necessary_margin": [str(lower), str(upper)]}
