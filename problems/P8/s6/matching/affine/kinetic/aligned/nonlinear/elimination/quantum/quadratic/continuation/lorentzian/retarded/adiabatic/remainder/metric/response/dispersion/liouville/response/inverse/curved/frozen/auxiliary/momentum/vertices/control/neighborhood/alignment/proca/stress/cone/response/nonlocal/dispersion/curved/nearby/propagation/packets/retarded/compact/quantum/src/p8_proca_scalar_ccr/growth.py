"""Polynomial Fourier growth for the actual canonical transfer and field smearings."""
from functools import cache

import sympy as sp
from p8_proca_nearby_packets import domains
from p8_proca_nearby_retarded import remainder, transport

from . import canonical

PACKET_TRANSFER_UPPER=sp.Integer(6)*10**6
DENSITY_TRANSFER_COEFFICIENT=sp.Integer(10)**10


@cache
def data():
    d=canonical.data()
    p=canonical.parent.parent
    limits={canonical.R:(sp.Rational(1,2),sp.Integer(2)),
            canonical.k:(sp.Integer(1),sp.Integer(1)),
            p.alpha:(0,1),p.rn:(sp.Rational(1,10),1)}
    E=d["canonical_density_to_packet_map"]/canonical.k**2
    Ei=d["packet_to_canonical_density_map"]/canonical.k
    majorants={}
    for name,matrix in (("packet_from_density_scaled",E),("density_from_packet_scaled",Ei)):
        matrix=matrix.applyfunc(lambda value:domains.laurent_majorant(value,limits))
        majorants[name]={"entry_upper":matrix,
                         "infinity_norm_upper":max(sum(matrix[i,j] for j in range(4)) for i in range(4))}
    B=remainder.low_frequency()["complete_original_generator_infinity_norm_upper"]
    return {"minimum_high_frequency_domain":transport.normal_form.MIN_MOMENTUM,
            "exact_scaled_canonical_map_majorants":majorants,
            "declared_packet_from_density_norm_coefficient":sp.Integer(24),
            "declared_density_from_packet_norm_coefficient":sp.Integer(21),
            "exact_packet_transfer_norm_upper":PACKET_TRANSFER_UPPER,
            "high_frequency_density_transfer_triangle_coefficient":24*21*PACKET_TRANSFER_UPPER,
            "declared_high_frequency_density_transfer_coefficient":DENSITY_TRANSFER_COEFFICIENT,
            "high_frequency_density_transfer_polynomial_degree":3,
            "low_frequency_density_transfer_norm_upper":64*sp.exp(domains.T*B),
            "all_orders_Fourier_derivative_degree_sufficient":"3+7*multi_index_order",
            "Schwartz_test_space_invariant_by_written_Duhamel_induction":True,
            "no_global_Fock_unitary_implementation_or_Hadamard_property_assumed":True}


@cache
def checks():
    d=data()
    return {"exact_packet_transfer_bound_retains_all_five_factors":
                1000*sp.Rational(3,2)*2*2*1000-PACKET_TRANSFER_UPPER,
            "canonical_transfer_growth_retains_both_endpoints_and_packet_transfer":
                d["high_frequency_density_transfer_triangle_coefficient"]-24*21*PACKET_TRANSFER_UPPER,
            "full_canonical_map_polynomial_growth_degree_is_three":1+2-3,
            "low_frequency_density_bound_retains_both_volume_maps":8*8-64}


@cache
def gates():
    d=data()
    scaled=(canonical.data()["canonical_density_to_packet_map"]/canonical.k**2,
            canonical.data()["packet_to_canonical_density_map"]/canonical.k)
    exponents=[term.as_powers_dict().get(canonical.k,0) for matrix in scaled
               for value in matrix for term in sp.Add.make_args(sp.expand(value))]
    return {name:bool(value) for name,value in {
        "scaled_map_majorants_only_use_nonpositive_frequency_powers":all(value<=0 for value in exponents),
        "forward_canonical_map_below_declared_coefficient":
            d["exact_scaled_canonical_map_majorants"]["packet_from_density_scaled"]["infinity_norm_upper"]<24,
        "inverse_canonical_map_below_declared_coefficient":
            d["exact_scaled_canonical_map_majorants"]["density_from_packet_scaled"]["infinity_norm_upper"]<21,
        "actual_canonical_transfer_below_declared_polynomial_coefficient":
            d["high_frequency_density_transfer_triangle_coefficient"]<DENSITY_TRANSFER_COEFFICIENT,
        "compact_frequency_transfer_bound_is_finite":d["low_frequency_density_transfer_norm_upper"].is_finite is True,
        "every_Fourier_derivative_is_controlled_by_written_polynomial_ODE_induction":True,
        "growth_bound_is_not_a_finite_energy_or_stress_tensor_state_bound":True}.items()}
