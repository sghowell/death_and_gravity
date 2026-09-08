"""Locally finite Cauchy preparation in the separately named all-order state."""
import sympy as sp
from p8_affine_retuned.bounds import exact
from p8_vector_state import comparison, wkb

from . import cutoffs, series


def initial_data(kind, mass_time_product, initial_frequency):
    m, nu = exact(mass_time_product, "m0_tau"), exact(initial_frequency, "initial_frequency")
    if not bool(m >= wkb.MASS_TIME_MIN) or not bool(nu >= m):
        raise ValueError("Require exact initial frequency >= m0*tau >=1000")
    u0 = -sp.Rational(1, 2)
    z0 = 1-m**2/nu**2
    point = {wkb.u: u0, wkb.z: z0}
    lam = wkb.background()["lambda"].subs(point)
    frequency, slope = nu, lam*nu
    active, terms = [], {}
    for order in (1, 2):
        P = series.coefficient(kind, order)
        frequency += P.subs(point)*nu**(1-2*order)
        slope += (series.D0(P)+(1-2*order)*wkb.background()["lambda"]*P).subs(point)*nu**(1-2*order)
    old_frequency, old_slope = sp.factor(frequency), sp.factor(slope)
    order = 3
    while True:
        threshold = cutoffs.threshold(kind, order, m)
        if bool(nu <= threshold):
            break
        P = series.coefficient(kind, order)
        weight = cutoffs.turn_on(nu/threshold)
        frequency_term = weight*P.subs(point)*nu**(1-2*order)
        slope_term = weight*(series.D0(P)+(1-2*order)*wkb.background()["lambda"]*P).subs(point)*nu**(1-2*order)
        frequency += frequency_term
        slope += slope_term
        active.append(2*order)
        terms[2*order] = {"cutoff": threshold, "weight": weight,
                         "frequency_correction": frequency_term, "slope_correction": slope_term}
        order += 1
    return {"m0_tau": m, "initial_frequency": nu,
            "comoving_momentum_squared": comparison.AMAX**2*(nu**2-m**2),
            "initial_u": u0, "initial_z": z0,
            "frozen_frequency": old_frequency, "frozen_frequency_slope": old_slope,
            "all_order_frequency": sp.factor(frequency),
            "all_order_frequency_slope": sp.factor(slope),
            "all_order_half_log_rate": sp.factor(slope/(2*frequency)),
            "active_higher_derivative_orders": active,
            "included_correction_terms": terms,
            "first_inactive_threshold": threshold}
