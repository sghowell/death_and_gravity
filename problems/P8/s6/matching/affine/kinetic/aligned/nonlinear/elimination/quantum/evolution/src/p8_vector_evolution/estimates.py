"""Physical C1 first-variation bounds in the unchanged matching prescription."""
from functools import cache

import sympy as sp
from p8_vector_clock_matching import bounds as matched
from p8_vector_subtraction import tail

from . import majorants, mixing, ward


def physical_bounds(planck_time_product, reference_mass_time_product):
    prior = matched.physical_bounds(planck_time_product, reference_mass_time_product)
    L, R = prior["M_tau"], prior["m0_tau"]
    state = tail.physical_bounds(L, R)
    evolved_derivative = 40*mixing.MIXING_CONSTANT/(R*L**2)
    tail_derivative = majorants.TAIL_DERIVATIVE_CONSTANT/(72*R**2*L**2)
    derivative = evolved_derivative+tail_derivative
    result = {"M_tau": L, "m0_tau": R,
              "exact_reference_derivative_over_reference_density_rate": evolved_derivative,
              "reference_tail_derivative_over_reference_density_rate": tail_derivative,
              "full_subtracted_state_derivative_over_reference_density_rate": derivative}
    envelopes = ward.local_envelopes()
    for name in ("energy", "pressure"):
        key = name+"_derivative"
        local = sum(envelopes[n][key]*R**(4-2*n) for n in (0, 1, 2))/(576*L**2)
        result[name+"_local_derivative_over_reference_density_rate"] = local
        result[name+"_matched_derivative_over_reference_density_rate"] = derivative+local
    source_local = sum(envelopes[n]["negative_clock_source"]*R**(4-2*n) for n in (0, 1, 2))/(576*L**2)
    source_state = derivative+6*(state["full_subtracted_energy_integral_over_reference_density"]
                                +state["full_subtracted_pressure_integral_over_reference_density"])
    result.update({"clock_source_local_over_reference_scalar_equation": source_local,
                   "clock_source_state_over_reference_scalar_equation": source_state,
                   "clock_source_matched_over_reference_scalar_equation": source_local+source_state,
                   "all_order_Hadamard_or_quantum_corrected_solution_claim": False})
    return result


@cache
def proof_checks():
    d = physical_bounds(10**12, 1000)
    return {"energy_derivative_example_below_stated_bound": bool(d["energy_matched_derivative_over_reference_density_rate"] < sp.Rational(1, 10**14)),
            "pressure_derivative_example_below_stated_bound": bool(d["pressure_matched_derivative_over_reference_density_rate"] < sp.Rational(1, 10**14)),
            "clock_source_example_below_stated_bound": bool(d["clock_source_matched_over_reference_scalar_equation"] < sp.Rational(1, 10**14))}
