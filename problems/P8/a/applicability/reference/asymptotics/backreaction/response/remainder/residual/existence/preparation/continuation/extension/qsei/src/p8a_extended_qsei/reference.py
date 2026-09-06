"""Signed reference bound from the actual longer unforced SEE solution."""

from functools import cache

import sympy as sp
from p8a_extension import bounds as actual
from p8a_see_qsei.scattering import nonnegative


def energy_dimensionless(scale, hubble, potential, delta):
    """E_ref in hbar/(pi² T0⁴) units, conditional on the fixed unforced SEE.

    This is not an EED formula for arbitrary target states, and it is not
    valid during the external preparation source without its source term.
    """
    a, h, u, delta = map(sp.sympify, (scale, hubble, potential, delta))
    return (a**2*(u+h**2)-1)/(960*delta*a**4)


@cache
def calibration():
    data = actual.calibration()
    geo, gains = data["geometry"], data["metric_gains"]
    delta, weight = data["delta"], data["weight_cap"]
    distance = data["gate"]["fixed_point_weighted"]
    cap = sp.Rational(9, 10**8)
    # N-Nbar=a²(Delta u+(h+hbar)Delta h)
    #         +Delta(a²)(ubar+hbar²), with the *weighted* metric gains.
    numerator_pair = geo["a_max"]**2*(gains["u"]+2*geo["h_max"]*gains["h"])
    numerator_pair += (geo["u_cap"]+geo["h_max"]**2)*gains["a_squared"]
    loss = weight*cap*numerator_pair
    negative = loss/(960*delta*geo["a_min"]**4)
    rounded = sp.Rational(1, 40)
    margins = {"sharp_weighted_distance": cap-distance,
               "signed_reference_coarsening": rounded-negative}
    if any(value <= 0 for value in margins.values()):
        raise ValueError("the extended actual-SEE signed reference bound failed")
    return {"actual_weighted_distance": distance, "distance_cap": cap,
            "weight_cap": weight, "numerator_weighted_pair": numerator_pair,
            "numerator_loss_upper": loss, "baseline_numerator_is_positive": True,
            "positive_reference_credit_asserted": False,
            "negative_EED_magnitude_upper": negative,
            "rounded_negative_EED_magnitude": rounded,
            "reference_EED_lower_dimensionless": -rounded,
            "strict_margins": margins}


def lower_bound(*, time_scale=1, hbar=1):
    """Pointwise lower E_ref on the actual A.14 source-free domain only."""
    time_scale = nonnegative(time_scale, positive=True)
    hbar = nonnegative(hbar)
    return -hbar*calibration()["rounded_negative_EED_magnitude"]/(sp.pi**2*time_scale**4)


def identities():
    x = sp.Symbol("x", real=True)
    a = sp.Function("a", positive=True)(x)
    delta, y = sp.symbols("delta y", positive=True)
    h, u = sp.diff(a, x)/a, -sp.diff(a, x, 2)/a
    hh = h/a
    geometric = -3*(sp.diff(hh, x)/a+hh**2)/(2880*delta)
    ordinary_radiation = 3/(2880*delta*a**4)
    a2, b2, u1, u0, h1, h0 = sp.symbols("a2 b2 u1 u0 h1 h0", real=True)
    fraction = 1-delta/y**4
    abar, hbar = y*fraction**sp.Rational(1, 4), 1/(y*fraction**sp.Rational(3, 4))
    ubar = 2*delta/(y**6*fraction**sp.Rational(3, 2))
    return {
        "actual_SEE_reference_not_old_metric_transfer": sp.simplify(geometric-ordinary_radiation-energy_dimensionless(a, h, u, delta)),
        "exact_weighted_numerator_difference": sp.expand(a2*(u1+h1**2)-b2*(u0+h0**2)
                                                          -a2*((u1-u0)+(h1+h0)*(h1-h0))-(a2-b2)*(u0+h0**2)),
        "unchanged_baseline_numerator_positive": sp.simplify(abar**2*(ubar+hbar**2)-1-3*delta/(y**4-delta)),
    }
