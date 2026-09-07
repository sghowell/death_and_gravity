"""Complete smooth geometric comparators and a rejected tight-future envelope."""

import sympy as sp
from p8a_maxwell.domain import nonnegative

from . import geometry

WIDTH = sp.Rational(1, 10**12)


def switch(variable):
    z = sp.sympify(variable)
    return 35*z**4-84*z**5+70*z**6-20*z**7


def raw_global_jet_caps():
    """Uniform all-p j0..4 caps for the C3, W4-infinity compact-future H."""
    # On 0<=x<=1/8, p>=1/2: |d^j[-p/(p/2-x)]|<=j! 2^(3j+2).
    base = [sp.factorial(j)*2**(3*j+2) for j in range(5)]
    unit_switch = [sp.Integer(1), sp.Rational(35, 16), sp.Rational(105, 4),
                   sp.Integer(210), sp.Integer(52920)]
    cutoff = [value*8**j for j, value in enumerate(unit_switch)]
    return [sum(sp.binomial(j, k)*base[j-k]*cutoff[k] for k in range(j+1)) for j in range(5)]


def complete_geometry(width=WIDTH):
    """Bounds for positive unit-L1 smoothing of the specified raw profile.

This is a geometric family for every p in [1/2,2/3], NOT an actual Maxwell
SEE/state construction. Its future completeness prevents a geometry-only
reading of the calibrated incompleteness theorem.
    """
    width = nonnegative(width, positive=True)
    raw = raw_global_jet_caps()
    errors = [width*raw[j+1] for j in range(4)]
    tube_gaps = [cap-error for cap, error in zip(geometry.ERROR_CAPS, errors, strict=True)]
    future_gaps = [cap-value for cap, value in zip(geometry.FUTURE_CAPS, raw[:4], strict=True)]
    if min(*tube_gaps) <= 0 or min(*future_gaps) < 0:
        raise ValueError("the proposed smoothing does not fit the robust history/future calibration")
    return {"power_range": [sp.Rational(1, 2), sp.Rational(2, 3)],
            "cutoff_ends_at_dimensionless_time": sp.Rational(1, 8),
            "mollifier_width_over_tau": width,
            "raw_global_jet_caps_through_four": raw,
            "past_C3_error_upper": errors,
            "strict_past_tube_margins": tube_gaps,
            "future_cap_margins_from_global_bounds": future_gaps,
            "future_static_by_time_over_tau": sp.Rational(1, 8)+width,
            "future_log_scale_factor_loss_upper": 4*(sp.Rational(1, 8)+width),
            "covers_every_reference_power_in_interval": True,
            "future_timelike_and_null_complete": True,
            "actual_quantum_SEE_solution_asserted": False,
            "geometric_assumptions_alone_force_incompleteness": False}


def rejected_tight_future():
    """The old candidate c=d forces a contradiction before any QSEI is used."""
    x = sp.Rational(1, 20)
    hubble0, first0, second0 = -sp.Rational(199, 100), -sp.Rational(11, 2), -sp.Integer(33)
    third = sp.Integer(800)/(1-x)**4
    upper = hubble0+first0*x+second0*x*x/2+third*x**3/6
    required_lower = -sp.Rational(21, 10)/(1-x)
    return {"test_time_over_tau": x,
            "initial_upper_jets": [hubble0, first0, second0],
            "third_jet_upper_on_test_interval": third,
            "Taylor_Hubble_upper": upper,
            "tight_future_required_Hubble_lower": required_lower,
            "strict_incompatibility_gap": required_lower-upper,
            "QSEI_used_in_this_exclusion": False,
            "rejected_as_nontrivial_cosmological_calibration": True}


def identities():
    z, x = sp.symbols("z x", real=True)
    p = sp.Symbol("p", positive=True)
    polynomial = switch(z)
    raw_past = -p/(p/2-x)
    raw_future = raw_past*(1-switch(8*x))
    result = {"switch_first_factor": sp.expand(sp.diff(polynomial, z)-140*z**3*(1-z)**3),
              "switch_second_factor": sp.expand(sp.diff(polynomial, z, 2)-420*z**2*(1-z)**2*(1-2*z)),
              "switch_third_factor": sp.expand(sp.diff(polynomial, z, 3)-840*z*(1-z)*(1-5*z+5*z*z)),
              "switch_fourth_coefficient_majorant": sum(abs(value) for value in sp.Poly(sp.diff(polynomial, z, 4), z).all_coeffs())-52920,
              "rejected_tight_future_gap": rejected_tight_future()["strict_incompatibility_gap"]-sp.Rational(4707907, 62554080)}
    for j in range(4):
        result[f"raw_C3_join_{j}"] = sp.simplify(sp.diff(raw_future-raw_past, x, j).subs(x, 0))
        result[f"raw_C3_flat_future_join_{j}"] = sp.simplify(sp.diff(raw_future, x, j).subs(x, sp.Rational(1, 8)))
    return result


def calibration():
    return {"all_p_complete_geometric_family": complete_geometry(),
            "rejected_past_caps_as_future_caps": rejected_tight_future()}
