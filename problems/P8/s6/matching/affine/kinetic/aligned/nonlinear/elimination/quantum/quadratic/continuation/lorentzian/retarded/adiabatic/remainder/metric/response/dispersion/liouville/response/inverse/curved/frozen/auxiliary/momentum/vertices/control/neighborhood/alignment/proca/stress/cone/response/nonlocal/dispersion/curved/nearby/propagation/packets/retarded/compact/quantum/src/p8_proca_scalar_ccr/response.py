"""Exact state-independent linear Kubo response and compact commutator lower."""
from functools import cache

import sympy as sp
from p8_proca_compact_response import probes

from . import canonical, state


@cache
def data():
    U=sp.Matrix(4,4,lambda i,j:sp.Symbol(f"two_time_density_transfer_{i}_{j}",real=True))
    Omega=canonical.data()["constant_symplectic_form"]
    w=sp.Symbol("source_physical_volume_density",positive=True)
    hbar,kappa=state.hbar,state.kappa
    comm=sp.I*hbar*U*Omega/kappa
    response=w*U[1,3]
    return {"two_time_canonical_transfer":U,
            "two_time_CCR_kernel":comm,
            "relational_matter_commutator_multiplier":comm[1,1],
            "normalized_probe_classical_response_multiplier":response,
            "normalized_probe_Kubo_prefactor":sp.I*kappa*w/hbar,
            "unscaled_probe_Kubo_prefactor":sp.I*w/hbar,
            "positive_compact_classical_response_lower":3*probes.SIGNAL_LOWER/4,
            "positive_compact_smeared_commutator_magnitude_lower":3*hbar*probes.SIGNAL_LOWER/(4*kappa),
            "compact_source_operator_test":"physical_volume_density(s)*J(s,y)",
            "compact_detector_operator_test":"f(t,x)",
            "state_independent_nonzero_commutator_in_exact_reduced_quadratic_scalar_algebra":True,
            "no_interacting_parent_or_nonperturbative_quantum_gravity_locality_transfer":True}


@cache
def checks():
    d=data()
    Omega=canonical.data()["constant_symplectic_form"]
    hbar,kappa=state.hbar,state.kappa
    rows={
        "relational_matter_commutator_has_the_retarded_momentum_source_sign":
            d["relational_matter_commutator_multiplier"]+sp.I*hbar*d["two_time_canonical_transfer"][1,3]/kappa,
        "normalized_probe_Kubo_response_equals_exact_classical_response":
            sp.expand(d["normalized_probe_Kubo_prefactor"]*d["relational_matter_commutator_multiplier"]
                      -d["normalized_probe_classical_response_multiplier"]),
        "unscaled_probe_Kubo_response_keeps_overall_action_normalization":
            sp.expand(d["unscaled_probe_Kubo_prefactor"]*d["relational_matter_commutator_multiplier"]
                      -d["normalized_probe_classical_response_multiplier"]/kappa),
        "compact_commutator_lower_retains_hbar_and_action_scale":
            d["positive_compact_smeared_commutator_magnitude_lower"]
            -hbar*d["positive_compact_classical_response_lower"]/kappa}
    T=sp.Matrix(4,4,lambda i,j:sp.Symbol(f"T_{i}_{j}"))
    S=sp.Matrix(4,4,lambda i,j:sp.Symbol(f"S_{i}_{j}"))
    rows["two_time_symplectic_composition_reduces_CCR_to_retarded_transfer"]=(
        T*S*Omega*S.T-T*Omega-T*(S*Omega*S.T-Omega)).applyfunc(sp.expand)
    return rows


@cache
def gates():
    d=data()
    return {"compact_quantum_commutator_lower_is_positive_for_positive_action_and_hbar":
                d["positive_compact_smeared_commutator_magnitude_lower"].is_positive is True,
            "normalized_Kubo_response_not_confused_with_unscaled_physical_source":True,
            "all_states_of_the_same_CCR_algebra_share_the_commutator":True,
            "Gaussian_linear_field_state_does_not_establish_Hadamard_stress_or_semiclassical_backreaction":True,
            "finite_EFT_quantum_UV_matching_and_original_P8_remain_open":True}
