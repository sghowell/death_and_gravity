"""Continuous homogeneous local physical-stress response bounds."""
from functools import cache

import sympy as sp
from p8_affine_retuned.bounds import exact
from p8_vector_dimensional import local
from p8_vector_mass_adiabatic.bounds import rational_bound

from . import consistency, counterterms, jets, radial


def component(value):
    if type(value) is not str or value not in ("energy", "pressure"):
        raise ValueError("Require physical energy or pressure output")
    return value


def physical_operator(output, order):
    output, order = component(output), jets.order(order)
    return _physical_operator(output, order)


@cache
def _physical_operator(output, order):
    old = counterterms.first.finite_coefficients()[order][output].subs(local.ell, 0)
    background = counterterms.mass**(4-2*order)*jets.actual(old)
    if output == "energy":
        value = -radial.matched("N", order)-3*jets.v[0]*background
    else:
        value = radial.matched("Z", order)/3-(jets.n[0]+3*jets.v[0])*background
    return jets.linear_clean(value)


@cache
def envelopes():
    result = {}
    for output in ("energy", "pressure"):
        result[output] = {}
        for order in (0, 1, 2):
            value = sp.expand(physical_operator(output, order)/counterterms.mass**(4-2*order))
            result[output][order] = {
                source: {j: rational_bound(sp.factor(value.coeff(row[j]))) for j in range(5)}
                for source, row in (("N", jets.n), ("Z", jets.v))}
    return result


def operator_bound(planck_time_product, mass_time_product):
    L, m = exact(planck_time_product, "M_tau"), exact(mass_time_product, "m0_tau")
    if not bool(L > 0) or not bool(m >= 1000):
        raise ValueError("Require M*tau>0 and m0*tau>=1000")
    pieces = {output: {j: sp.factor(m**(4-2*j)*sum(sum(row.values()) for row in sources.values())/(576*L**2))
                       for j, sources in table.items()}
              for output, table in envelopes().items()}
    totals = {output: sum(table.values()) for output, table in pieces.items()}
    return {"by_component_and_adiabatic_half_order": pieces,
            "physical_C4_to_C0_component_bounds": totals,
            "joint_physical_C4_to_C0_bound": max(totals.values()),
            "source_norm": "maximum of both physical lapse and log-scale source time derivatives 0..4 on I",
            "scope": "Finite homogeneous local component only; no state remainder or no-loss feedback claim"}


@cache
def checks():
    result = {}
    for order in (0, 1, 2):
        for output in ("energy", "pressure"):
            result[output+"_physical_normalization_cancels_constant_spatial_rescaling_"+str(2*order)] = sp.factor(
                sp.expand(physical_operator(output, order)).coeff(jets.v[0]))
    # Matrix self-adjointness belongs to the coordinate-density Euler currents,
    # not separately to the volume-normalized rho and p components.
    for output in ("N", "Z"):
        result[output+"_top_source_order_is_four"] = sp.expand(
            radial.matched(output, 2)-sum(value*field for source, row in (("N", jets.n), ("Z", jets.v))
                                         for value, field in zip(consistency.coefficients(output, 2, source), row[:5])))
    return result
