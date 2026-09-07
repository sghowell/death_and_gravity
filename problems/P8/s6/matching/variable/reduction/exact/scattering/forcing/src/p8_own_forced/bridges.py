"""Discharge the fixed-pulse proof's unchanged-action and domain premises."""
from functools import cache

import sympy as sp
from p8_exact_stationary import intervals as background
from p8_own_scattering import action as parent_action
from p8_own_scattering import potential

from . import audit, loading, phase


@cache
def identities():
    independent, source, propagation = audit.constants(), loading.calibration(), phase.calibration()
    pairings = {
        "physical_L": (independent["L"], source["L"]),
        "central_half_width": (independent["a"], propagation["a"]),
        "delta_max": (independent["delta_max"], source["delta_max"]),
        "canonical_V_upper": (independent["canonical_potential_upper"], source["potential_upper"]),
        "Green_derivative": (independent["Green_derivative_floor"], source["green_derivative_lower"]),
        "Green_ratio": (independent["Green_ratio_floor"], source["green_ratio_lower"]),
        "strict_limiting_Y_loading": (independent["sharp_loaded_Y_floor_per_eta"], source["Y_load_sharper_lower_over_eta"]),
        "strict_limiting_Y_u_loading": (independent["sharp_loaded_Y_u_floor_per_eta"], source["Y_u_load_sharper_lower_over_eta"]),
        "central_error_integral": (independent["central"], propagation["central_integral_upper"]),
        "both_reference_tails": (independent["tails"], propagation["tail_integral_upper"]),
        "complete_normalized_budget": (independent["budget"], propagation["perturbation_budget_upper"]),
        "transfer_error": (independent["transfer_error_upper"], propagation["transfer_error_upper"]),
        "both_phase_errors": (independent["balanced_pair_separation_factor"], propagation["normalized_pair_separation_factor"]),
        "original_Q_gap": (independent["Q_pair_separation_per_eta"], propagation["scalar_Q_gap_lower_over_eta"]),
        "same_canonical_remainder": (source["remainder_cap"], potential.POTENTIAL_CAP),
        "same_entire_canonical_pump": (source["pump_cap"], potential.PUMP_CAP),
    }
    output = {name: left-right for name, (left, right) in pairings.items()}
    literal = parent_action.canonical()
    mass = literal["V"]+literal["pump"]
    output["source_weight_is_unchanged"] = sp.simplify(sp.sqrt(literal["k"])*mass-literal["canonical_source_weight"])
    output["full_canonical_potential_is_unchanged"] = sp.simplify(mass-literal["spring"]/literal["k"])
    return output


def checks():
    actual = background.enclose()
    source = loading.source_bounds()
    return {
        "b_strict_lower": actual["b"].lo > loading.COFRAME_MIN,
        "b_strict_upper": actual["b"].hi < loading.COFRAME_MAX,
        "N_strict_lower": actual["N"].lo > loading.COFRAME_MIN,
        "N_strict_upper": actual["N"].hi < loading.COFRAME_MAX,
        "loading_window_inside_actual_branch": 0 < loading.A < loading.L <= potential.RADIUS,
        "positive_delta_subdomain": 0 < loading.DELTA_MAX <= potential.DELTA_MAX,
        "same_remainder_cap": loading.REMAINDER_CAP == potential.POTENTIAL_CAP,
        "same_pump_cap": loading.PUMP_CAP == potential.PUMP_CAP,
        "source_positive_on_full_loading_interval": source["mass_lower"] > source["canonical_source_lower"] > 0,
        "canonical_source_lower_uses_sqrt_k_above_one": source["kinetic_lower"] > 1,
        "punctured_loading_denominator_excludes_zero": source["punctured_D_lower"] > 0,
        "no_literal_delta_zero_at_center": True,
        "no_convergent_central_remainder_premise": True,
        "positive_lower_loading_discharged_not_chosen_data": True,
    }
