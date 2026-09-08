"""Quantified transfer from the frozen Gaussian to the new Cauchy state."""
from functools import cache

import sympy as sp
from p8_vector_clock_matching import bounds as matched
from p8_vector_evolution import estimates as differentiated
from p8_vector_state import comparison, wkb

from . import cutoffs


@cache
def proof_checks():
    A, m = comparison.AMAX, wkb.MASS_TIME_MIN
    return {"old_exact_mode_mixed_product_below_seven": bool(4*12 < 7**2),
            "old_exact_absolute_reference_density_constant": (12+sp.Rational(3, 2)*4)/2 == 9,
            "state_change_value_integral_constant": bool(3*A**4 < 20),
            "state_change_derivative_readout_constant": bool(sp.Rational(816, 1)/m+84 < 85),
            "state_change_derivative_integral_constant": bool(255*A**5/24 < 100),
            "complex_old_exact_bilinear_readout_constant": (sp.Rational(5, 4)*64+sp.Rational(3, 2)*16)/2 == 52,
            "complex_state_change_all_polarization_constant": bool(3*52*sp.Rational(13, 4) < 1000),
            "complex_state_change_radial_high_tail_integrable": bool(sp.Rational(9, 4)-5 < -1),
            "complex_state_change_radial_low_tail_integrable": bool(sp.Rational(7, 4) > -1)}


def scale_checks():
    data = physical_bounds(10**12, 1000)
    keys = [name for name in data if "new_state" in name]
    return {name+"_below_stated_bound": bool(data[name] < sp.Rational(1, 10**14)) for name in keys}


def physical_bounds(planck_time_product, reference_mass_time_product):
    old = matched.physical_bounds(planck_time_product, reference_mass_time_product)
    rates = differentiated.physical_bounds(planck_time_product, reference_mass_time_product)
    L, R, C = old["M_tau"], old["m0_tau"], cutoffs.INITIAL_MIXING_CONSTANT
    value, derivative = 20*C/(R**2*L**2), 100*C/(R*L**2)
    source = derivative+12*value
    result = {"M_tau": L, "m0_tau": R,
              "state_change_energy_or_pressure_over_reference_density": value,
              "state_change_density_derivative_over_reference_density_rate": derivative,
              "state_change_clock_source_over_reference_scalar_equation": source}
    for name in ("energy", "pressure"):
        result[name+"_new_state_total_over_reference_density"] = old[name+"_candidate_matched_total_over_reference_density"]+value
        result[name+"_new_state_derivative_over_reference_density_rate"] = rates[name+"_matched_derivative_over_reference_density_rate"]+derivative
    result["clock_source_new_state_over_reference_scalar_equation"] = rates["clock_source_matched_over_reference_scalar_equation"]+source
    result["full_quantum_solution_or_VGB_claim"] = False
    return result
