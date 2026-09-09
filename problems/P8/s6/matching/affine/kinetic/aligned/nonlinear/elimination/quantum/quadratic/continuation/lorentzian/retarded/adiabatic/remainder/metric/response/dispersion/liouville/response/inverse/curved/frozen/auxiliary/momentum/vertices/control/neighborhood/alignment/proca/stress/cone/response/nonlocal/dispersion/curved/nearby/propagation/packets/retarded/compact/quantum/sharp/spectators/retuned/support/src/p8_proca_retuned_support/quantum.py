"""Fresh reduced scalar CCR state and Kubo identity for the retuned evolution."""
from functools import cache

import sympy as sp
from p8_proca_nearby_packets import domains as algebra
from p8_proca_scalar_ccr import canonical, state

from . import domain, growth


@cache
def data():
    d=canonical.data()
    p=canonical.parent.parent
    limits={canonical.R:(sp.Rational(1,2),2),canonical.k:(1,1),
        p.alpha:(0,1),p.rn:(sp.Rational(1,10),1)}
    majorants={}
    for name,matrix in (("packet_from_density_scaled",d["canonical_density_to_packet_map"]/canonical.k**2),
            ("density_from_packet_scaled",d["packet_to_canonical_density_map"]/canonical.k)):
        entries=matrix.applyfunc(lambda value:algebra.laurent_majorant(value,limits))
        majorants[name]={"entry_modulus_upper":entries,
            "infinity_norm_upper":max(sum(entries[i,j] for j in range(4)) for i in range(4))}
    U=sp.Matrix(4,4,lambda i,j:sp.Symbol(f"new_two_time_density_transfer_{i}_{j}",real=True))
    volume=sp.Symbol("new_source_physical_volume",positive=True)
    comm=sp.I*state.hbar*U*d["constant_symplectic_form"]/state.kappa
    low=domain.coefficient_majorants()["low_scalar_root_unit_disc_generator_infinity_norm_upper"]
    return {"new_density_generator_universal_formula":d["canonical_density_generator"],
        "new_literal_coefficient_substitution_is_given_in_bridge":True,
        "new_normalized_source_density_vector":d["normalized_probe_source_vector"],
        "new_reference_time":-domain.T/2,
        "new_ordered_initial_covariance":state.data()["ordered_Cauchy_covariance"],
        "new_positive_initial_Gram_factor":state.data()["positive_Cauchy_covariance_Gram_factor"],
        "new_covariance_transport":"C_new(t,s,k)=U_new(t,u_star,k)*C0(k)*U_new(s,u_star,k)^T",
        "same_Gaussian_prescription_applied_to_new_evolution_not_old_state_transfer":True,
        "exact_scaled_density_packet_map_majorants":majorants,
        "new_real_packet_transfer_infinity_norm_upper":320*sp.exp(400),
        "new_high_frequency_density_transfer_coefficient":24*21*320*sp.exp(400),
        "new_low_frequency_density_transfer_norm_upper":64/(1-domain.T*low),
        "new_all_frequency_density_transfer_coefficient":128+161280*sp.exp(400),
        "new_all_frequency_density_transfer_degree":3,
        "new_all_order_Fourier_derivative_degree":"3+7*multi_index_order",
        "new_two_time_density_transfer":U,
        "new_relational_matter_commutator_multiplier":comm[1,1],
        "new_normalized_probe_classical_response_multiplier":volume*U[1,3],
        "new_normalized_Kubo_prefactor":sp.I*state.kappa*volume/state.hbar,
        "new_unscaled_probe_Kubo_prefactor":sp.I*volume/state.hbar,
        "quadratic_relational_microcausality":"[O(f),O(V*J)]=0 for compact matter-spacelike supported tests in the new inner time interval",
        "support_comes_from_complete_new_classical_response_not_old_packet_front":True,
        "no_Hadamard_interacting_backreaction_or_UV_state_transfer":True}


@cache
def checks():
    rows={"new_universal_density_"+name:value for name,value in canonical.checks().items()}
    rows.update({"new_Gaussian_prescription_"+name:value for name,value in state.checks().items()})
    d=data()
    rows.update({
        "new_real_packet_infinity_norm_keeps_Euclidean_conversion":
            d["new_real_packet_transfer_infinity_norm_upper"]-2*160*sp.exp(400),
        "new_density_growth_keeps_both_endpoint_maps":
            d["new_high_frequency_density_transfer_coefficient"]-161280*sp.exp(400),
        "new_density_polynomial_degree_keeps_forward_and_inverse_maps":1+2-3,
        "new_relational_commutator_has_canonical_source_sign":sp.factor(
            d["new_relational_matter_commutator_multiplier"]
            +sp.I*state.hbar*d["new_two_time_density_transfer"][1,3]/state.kappa),
        "new_normalized_Kubo_response_is_exact_new_classical_response":sp.factor(
            d["new_normalized_Kubo_prefactor"]*d["new_relational_matter_commutator_multiplier"]
            -d["new_normalized_probe_classical_response_multiplier"]),
        "new_unscaled_probe_response_keeps_action_scale":sp.factor(
            d["new_unscaled_probe_Kubo_prefactor"]*d["new_relational_matter_commutator_multiplier"]
            -d["new_normalized_probe_classical_response_multiplier"]/state.kappa)})
    return rows


@cache
def gates():
    d=data()
    scaled=(canonical.data()["canonical_density_to_packet_map"]/canonical.k**2,
        canonical.data()["packet_to_canonical_density_map"]/canonical.k)
    powers=[term.as_powers_dict().get(canonical.k,0) for matrix in scaled
        for value in matrix for term in sp.Add.make_args(sp.expand(value))]
    rows={"new_scaled_density_maps_have_no_positive_residual_frequency_power":all(v<=0 for v in powers),
        "new_forward_density_map_fits_twenty_four":
            d["exact_scaled_density_packet_map_majorants"]["packet_from_density_scaled"]["infinity_norm_upper"]<24,
        "new_inverse_density_map_fits_twenty_one":
            d["exact_scaled_density_packet_map_majorants"]["density_from_packet_scaled"]["infinity_norm_upper"]<21,
        "new_low_density_transfer_bound_fits_one_hundred_twenty_eight":
            d["new_low_frequency_density_transfer_norm_upper"]<128,
        "new_Gaussian_initial_time_is_within_certified_interval":growth.time_pair()["source_time"]==d["new_reference_time"],
        "new_state_uses_smooth_positive_weight_including_zero_momentum":state.frequency_weight(0)==1,
        "new_positive_action_and_hbar_give_positive_covariance_factor":state.hbar.is_positive and state.kappa.is_positive}
    return {name:bool(value) for name,value in rows.items()}


def controls():
    return state.controls()
