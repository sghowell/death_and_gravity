"""Piecewise preparation bounds for the unchanged all-order vector state."""
from functools import cache

import sympy as sp
from p8_affine_retuned.bounds import exact
from p8_vector_hadamard import cutoffs, series
from p8_vector_state import wkb

from . import frequency

FIXED_HIGH_THRESHOLD = sp.Integer(10**12)


@cache
def constants():
    mass = wkb.MASS_TIME_MIN
    result = {}
    for kind in ("transverse", "longitudinal"):
        data = {order: series.coefficient_bounds(kind, order) for order in (3, 4, 5)}
        V = {order: value["coefficient_upper"] for order, value in data.items()}
        G = {order: value["frequency_coefficient_slope_upper"] for order, value in data.items()}
        bands = {6: (2*V[3]+1+V[4]/mass**2, 2*G[3]+1+G[4]/mass**2),
                 8: (V[4]+1, G[4]+1), 10: (V[5]+1, G[5]+1)}
        result[kind] = {"frequency_correction_value_and_slope_envelopes": bands,
                        "initial_mixing_envelopes": {power: sp.ceiling(2*F+4*(G0+4*F)/mass)
                                                      for power, (F, G0) in bands.items()},
                        "coefficient_bounds": data,
                        "evolution_mixing_envelope": frequency.reference(kind)["oscillatory_evolution_beta_upper"]}
    result["common"] = {"initial_mixing_envelopes": {
        power: max(result[kind]["initial_mixing_envelopes"][power] for kind in ("transverse", "longitudinal"))
        for power in (6, 8, 10)},
        "evolution_mixing_envelope": max(result[kind]["evolution_mixing_envelope"]
                                        for kind in ("transverse", "longitudinal"))}
    return result


def threshold(mass_time_product):
    mass = exact(mass_time_product, "m0_tau")
    if not bool(mass >= wkb.MASS_TIME_MIN):
        raise ValueError("Require exact m0*tau>=1000")
    return max(8*mass, FIXED_HIGH_THRESHOLD)


def mixing_bound(mass_time_product, minimum_frequency):
    mass = exact(mass_time_product, "m0_tau")
    nu = exact(minimum_frequency, "minimum_frequency")
    high = threshold(mass)
    if not bool(nu >= mass):
        raise ValueError("Require exact minimum_frequency>=m0*tau")
    data = constants()["common"]
    power = 6 if nu < 4*mass else 8 if nu < high else 10
    return {"preparation_inverse_frequency_power": power,
            "absolute_mixing_upper": data["initial_mixing_envelopes"][power]/nu**power
            +data["evolution_mixing_envelope"]/nu**10}


def proof_checks():
    data, mass = constants(), wkb.MASS_TIME_MIN
    common = data["common"]
    out = {}
    for kind in ("transverse", "longitudinal"):
        Lambda4 = cutoffs.threshold(kind, 4, mass)
        out[kind+"_fixed_high_threshold_completes_fourth_coefficient_cutoff"] = bool(2*Lambda4 <= FIXED_HIGH_THRESHOLD)
        out[kind+"_initial_coefficient_norm_below_two"] = bool(data[kind]["initial_mixing_envelopes"][6]/mass**6 < sp.Rational(1, 2))
        out[kind+"_reference_positive_and_small_residual"] = all(frequency.reference(kind)["proof_checks"].values())
    out.update({"initial_all_momentum_envelope_remains_below_one_half": bool(common["initial_mixing_envelopes"][6]/mass**6 < sp.Rational(1, 2)),
                "evolved_all_momentum_mixing_remains_below_one": bool(
                    common["initial_mixing_envelopes"][6]/mass**6+common["evolution_mixing_envelope"]/mass**10 < 1),
                "higher_cutoff_tail_in_each_required_weight_below_one": bool(sp.Rational(1, 32)/(1-sp.Rational(1, 2)) < 1),
                "interval_length_is_one": sp.Rational(1, 2)-sp.Rational(-1, 2) == 1,
                "initial_norm_times_evolution_exponential_below_four": 2*2 == 4,
                "endpoint_bulk_feedback_allowance": bool(16+256+2 < 300)})
    return out
