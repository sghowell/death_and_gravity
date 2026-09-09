"""Recomputed constant-mass prepared mode transport and physical integral bounds."""
from functools import cache

import sympy as sp
from p8_proca_local_response import bounds as new_local
from p8_vector_mass_remainder.evolution import algebra_checks
from p8_vector_regularity import preparation
from p8_vector_state import comparison, wkb

from . import envelopes, source, tail

AMAX = comparison.AMAX
__all__ = ["AMAX", "algebra_checks", "bound", "constants", "weight_checks"]


@cache
def weight_checks():
    z = source.z
    target = {("T", "energy"): (sp.Integer(1), sp.Integer(1)),
              ("L", "energy"): (sp.Integer(1), sp.Integer(1)),
              ("T", "pressure"): (sp.Rational(1, 3), (2*z-1)/3),
              ("L", "pressure"): ((1+2*z)/3, -sp.Rational(1, 3))}
    return {sector+"_"+component+"_closed_weight_"+key: sp.factor(source.readout(sector, component)[key]-expected)
            for (sector, component), row in target.items() for key, expected in zip(("A", "B"), row)}


@cache
def constants():
    by_sector = {}
    for sector in ("T", "L"):
        data = envelopes.reference(sector)
        C = data["base_residual_over_inverse_frequency_eighth_upper"]
        dC = data["delta_residual_over_inverse_frequency_eighth_upper"]
        B = data["delta_W_over_frequency_upper"]
        rb = data["source_log_frequency_variation_upper"]
        m = wkb.MASS_TIME_MIN
        E = B+data["delta_c_upper"]/m
        outputs = {}
        for component in ("energy", "pressure"):
            weights = source.readout(sector, component)
            da = source.linear_bound(weights["delta_A"])["absolute_upper"]
            db = source.linear_bound(weights["delta_B"])["absolute_upper"]
            # On the exact baseline, |A|<=1 and |B_weight|<=2. Thus
            # |J_ref|,|Z_ref|<=3 omega. Differentiate p_ref and the
            # normalization/phase before replacing each factor by a bound.
            K = AMAX*((2*da+2*E+db/2+2*rb+6*B)/m+6*B*AMAX)
            outputs[component] = {
                "delta_A_upper": da, "delta_B_upper": db,
                "delta_reference_bilinear_over_nu_squared_upper": sp.ceiling(K),
                "varied_reference_tail_over_inverse_frequency_fifth_upper":
                    tail.reference_tail(sector, component)["varied_reference_tail_over_inverse_frequency_fifth_upper"]}
        by_sector[sector] = {
            "delta_mixing_ninth_inverse_frequency_coefficient": sp.ceiling(16*(dC+2*C*B)),
            "delta_mixing_eighth_inverse_frequency_coefficient": sp.ceiling(16*C*B*AMAX),
            "physical_readouts": outputs}
    common = {key: max(row[key] for row in by_sector.values()) for key in
              ("delta_mixing_ninth_inverse_frequency_coefficient", "delta_mixing_eighth_inverse_frequency_coefficient")}
    common["physical_readouts"] = {
        component: {key: max(row["physical_readouts"][component][key] for row in by_sector.values())
                    for key in by_sector["T"]["physical_readouts"][component]}
        for component in ("energy", "pressure")}
    return {"by_sector": by_sector, "common": common,
            "frozen_initial_and_evolved_mixing": preparation.constants()["common"]}


def bound(planck_time_product):
    return _bound(new_local.exact_scale(planck_time_product))


@cache
def _bound(L):
    m = new_local.MASS
    data, old = constants()["common"], constants()["frozen_initial_and_evolved_mixing"]
    beta6 = old["initial_mixing_envelopes"][6]+old["evolution_mixing_envelope"]/m**4
    I4, I7, I8 = AMAX**3/(24*m), AMAX**3/(135*m**4), AMAX**3/(192*m**5)
    D9 = data["delta_mixing_ninth_inverse_frequency_coefficient"]
    D8 = data["delta_mixing_eighth_inverse_frequency_coefficient"]
    parts, nonlocal_bounds, complete = {}, {}, {}
    local = {key: row["C4_to_C0_upper"] for key, row in new_local.at_scale(L)["response_bounds"]["physical_finite_local_response"].items()}
    for component, readout in data["physical_readouts"].items():
        K = readout["delta_reference_bilinear_over_nu_squared_upper"]
        T = readout["varied_reference_tail_over_inverse_frequency_fifth_upper"]
        parts[component] = {
            "exact_minus_reference_variation": 3*(30*AMAX*(D9*I8+D8*I7)+6*beta6*K*I4)/L**2,
            "varied_reference_minus_adiabatic": 3*T/(54*m**2*L**2)}
        nonlocal_bounds[component] = sum(parts[component].values())
        complete[component] = nonlocal_bounds[component]+local[component]
    return {"nonlocal_components": parts,
            "nonlocal_subtracted_C10_to_C0_component_bounds": nonlocal_bounds,
            "finite_local_component_bounds": local,
            "complete_metric_response_C10_to_C0_component_bounds": complete,
            "joint_nonlocal_C10_to_C0_upper_bound": max(nonlocal_bounds.values()),
            "joint_complete_metric_C10_to_C0_upper_bound": max(complete.values()),
            "source_norm": "Maximum of both physical sources and time derivatives zero through ten on I, zero on an initial neighborhood",
            "scope": "Homogeneous prepared physical-metric Gaussian response; not a coupled feedback contraction or arbitrary spatial response"}

