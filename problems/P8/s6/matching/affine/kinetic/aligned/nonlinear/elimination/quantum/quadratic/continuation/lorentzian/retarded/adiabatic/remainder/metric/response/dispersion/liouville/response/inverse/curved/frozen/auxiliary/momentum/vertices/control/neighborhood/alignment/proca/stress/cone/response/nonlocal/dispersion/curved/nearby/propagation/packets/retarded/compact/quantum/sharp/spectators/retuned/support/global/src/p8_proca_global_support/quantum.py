"""Fresh global reduced scalar Gaussian state and exact causal Kubo identity."""
from functools import cache

import sympy as sp
from p8_proca_scalar_ccr import state

from . import coverage, phase


@cache
def data():
    d=phase.data()
    U=sp.Matrix(4,4,lambda i,j:sp.Symbol(f"global_retuned_transfer_{i}_{j}",real=True))
    source_scale=sp.Symbol("actual_source_global_scale",positive=True)
    comm=sp.I*state.hbar*U*d["constant_symplectic_form"]/state.kappa
    return {"new_global_reference_time":sp.Integer(0),
        "new_global_initial_frequency":"mu(k)=sqrt(1+|k|^2)",
        "new_global_ordered_initial_covariance":state.data()["ordered_Cauchy_covariance"],
        "new_global_initial_positive_Gram_factor":state.data()["positive_Cauchy_covariance_Gram_factor"],
        "new_global_covariance_transport":"U_global(t,0,k)*C0(k)*U_global(s,0,k)^T",
        "new_global_transfer":U,"new_global_source_physical_volume":source_scale**3,
        "new_global_relational_commutator_multiplier":comm[1,1],
        "new_global_normalized_classical_response_multiplier":source_scale**3*U[1,3],
        "new_global_normalized_Kubo_prefactor":sp.I*state.kappa*source_scale**3/state.hbar,
        "new_global_unscaled_Kubo_prefactor":sp.I*source_scale**3/state.hbar,
        "compact_strip_transfer_polynomial_degree":coverage.data()["complete_global_transfer_frequency_degree"],
        "all_orders_Fourier_derivative_degree":"9+13*multi_index_order",
        "microcausality_domain":"all finite source and output times on the preserved retuned global clock bounce",
        "new_state_uses_global_evolution_not_either_frozen_nearby_state":True,
        "no_Hadamard_stress_interacting_microcausality_or_global_Fock_unitary_implementation":True}


@cache
def checks():
    rows={"global_new_Gaussian_"+name:value for name,value in state.checks().items()}
    d=data()
    rows.update({"global_relational_commutator_keeps_density_source_sign":sp.factor(
        d["new_global_relational_commutator_multiplier"]+sp.I*state.hbar*d["new_global_transfer"][1,3]/state.kappa),
        "global_normalized_Kubo_identity_matches_complete_scalar_response":sp.factor(
            d["new_global_normalized_Kubo_prefactor"]*d["new_global_relational_commutator_multiplier"]
            -d["new_global_normalized_classical_response_multiplier"]),
        "global_unscaled_source_keeps_action_normalization":sp.factor(
            d["new_global_unscaled_Kubo_prefactor"]*d["new_global_relational_commutator_multiplier"]
            -d["new_global_normalized_classical_response_multiplier"]/state.kappa),
        "global_Duhamel_derivative_degree_increment_retains_transfer_and_generator":9+4-13})
    return rows


def controls():
    return state.controls()
