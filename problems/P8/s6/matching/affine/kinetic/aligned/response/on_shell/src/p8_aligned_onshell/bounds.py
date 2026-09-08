"""Continuous actual phase/readout bounds; no nonlinear evolution assertion."""
from functools import cache

import sympy as sp
from p8_affine_retuned.bounds import exact
from p8_aligned_response import bounds as response
from p8_aligned_response import leading

from . import phase


def matrix_infinity_bound(matrix):
    return max(sum(leading.compact_rational_bound(value) for value in matrix.row(i))
               for i in range(matrix.rows))


@cache
def continuous():
    data = phase.system()
    A = [matrix_infinity_bound(data["A"].diff(phase.u, j)) for j in range(3)]
    lapse = [matrix_infinity_bound(data["lapse_row"].diff(phase.u, j)) for j in range(4)]
    trace = [matrix_infinity_bound(data["trace_row"].diff(phase.u, j)) for j in range(4)]
    state = [sp.Integer(1), A[0], A[1]+A[0]**2, A[2]+3*A[0]*A[1]+A[0]**3]
    n_jets = [sum(sp.binomial(j, i)*lapse[i]*state[j-i] for i in range(j+1)) for j in range(4)]
    k_jets = [sum(sp.binomial(j, i)*trace[i]*state[j-i] for i in range(j+1)) for j in range(4)]
    BN, BK = max(n_jets), max(k_jets)
    return {"A_time_jet_bounds": A, "lapse_observer_time_jet_bounds": lapse,
            "trace_observer_time_jet_bounds": trace, "state_time_jet_multipliers": state,
            "lapse_time_jet_multipliers": n_jets, "trace_time_jet_multipliers": k_jets,
            "B_N": BN, "B_K": BK, "normalization_C": 1+BN+BK,
            "growth_exponent": A[0], "growth_integer_ceiling": sp.ceiling(A[0]),
            "growth_bound": sp.Integer(3)**sp.ceiling(A[0])}


def prepared_family(light_jet_bound):
    epsilon = exact(light_jet_bound, "epsilon")
    if epsilon.is_positive is not True:
        raise ValueError("Require positive epsilon for the nontrivial prepared family")
    data, jets = continuous(), phase.initial_jets()
    scale = epsilon/(data["growth_bound"]*data["normalization_C"])
    return {"epsilon": epsilon, "incoming_comoving_squared": phase.K_IN2,
            "quadratic_output_comoving_squared": 4*phase.K_IN2,
            "initial_scale": scale, "initial_state": scale*jets["initial_vector"],
            "lapse_jet_upper": epsilon*data["B_N"]/data["normalization_C"],
            "trace_jet_upper": epsilon*data["B_K"]/data["normalization_C"],
            "source_jet_E": leading.source_envelope(epsilon, epsilon)["source_jet_E"],
            "nonzero_source_third_cosine_squared_coefficient": scale**2*jets["source_third_initial"],
            "linear_light_and_leading_heavy_only": True,
            "nonlinear_solution_or_higher_order_remainder_claim": False}


def transfer(light_jet_bound, coupling):
    family = prepared_family(light_jet_bound)
    result = response.bound_values(coupling, phase.LEFT, phase.RIGHT, 4*phase.K_IN2, family["source_jet_E"])
    return {"family": family, "response": result,
            "source_class_realized_by_nonzero_linear_light_solution": True,
            "higher_order_nonlinear_or_quantum_error_claim": False}


@cache
def proof_checks():
    data, jets = continuous(), phase.initial_jets()
    vector = jets["initial_vector"]
    return {"strict_positive_continuous_matrix_bounds": all(value > 0 for value in data["A_time_jet_bounds"]),
            "exact_growth_integer_ceiling": data["growth_integer_ceiling"] == 238,
            "rational_growth_majorant": bool(data["growth_exponent"] <= data["growth_integer_ceiling"]),
            "strict_positive_readout_bounds": bool(data["B_N"] > 0 and data["B_K"] > 0),
            "unit_initial_direction_norm": max(abs(value) for value in vector) == 1,
            "rank_three_preparation": jets["prepared_matrix"].rank() == 3,
            "rank_four_lapse_observability": jets["observability_determinant"] != 0,
            "nonzero_initial_trace": jets["trace_initial"] != 0,
            "nonzero_third_lapse_derivative": jets["lapse_third"] != 0,
            "actual_nonzero_third_source_derivative": jets["source_third_initial"] != 0,
            "prepared_lapse_jet_bound_strict": bool(data["B_N"] < data["normalization_C"]),
            "prepared_trace_jet_bound_strict": bool(data["B_K"] < data["normalization_C"]),
            "output_quadratic_momentum_inside_response_window": 4*phase.K_IN2 == sp.Rational(1, 4),
            "time_window_covers_the_crossing": phase.LEFT < 0 < phase.RIGHT and phase.RIGHT-phase.LEFT == 1}
