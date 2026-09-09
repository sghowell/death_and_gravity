"""Actual canonical density phase and scalar CCR symplectic evolution."""
from functools import cache

import sympy as sp
from p8_proca_nearby_retarded import source as parent

R,k=parent.R,parent.k
piv,pic=sp.symbols("metric_scalar_canonical_density matter_canonical_density",real=True)


@cache
def data():
    d=parent.data()
    Omega=d["constant_symplectic_form"]
    density=sp.diag(1,1,R**3,R**3)
    density_dot=density.diff(R)*parent.parent.H*R
    generator=(density_dot+density*d["weighted_first_order_generator"])*density.inv()
    ham=R**3*d["Hamiltonian"].subs({parent.p:piv/R**3,parent.parent.P:pic/R**3},simultaneous=True)
    phase=sp.Matrix([parent.v,parent.parent.chi,piv,pic])
    packet=(d["old_to_packet_canonical_map"]*density.inv()).subs(parent.parent.q,k*k/R**2)
    return {"canonical_density_phase":phase,"normalized_to_density_phase":density,
            "canonical_density_Hamiltonian":ham,"canonical_density_generator":generator,
            "constant_symplectic_form":Omega,
            "canonical_density_to_packet_map":packet,
            "packet_to_canonical_density_map":packet.inv().applyfunc(sp.factor),
            "normalized_probe_source_vector":sp.Matrix([0,0,0,d["physical_volume_density"]]),
            "relational_observable_row":sp.Matrix([[0,1,0,0]]),
            "overall_positive_action_normalization":"kappa; normalized J multiplies kappa in the probe action",
            "unscaled_physical_probe":"j=kappa*J; its response is normalized classical response divided by kappa"}


@cache
def checks():
    d=data()
    M,Omega=d["canonical_density_generator"],d["constant_symplectic_form"]
    E=d["canonical_density_to_packet_map"]
    source=d["normalized_probe_source_vector"]
    rows={
        "actual_density_generator_is_the_literal_canonical_Hamiltonian_generator":
            M-Omega*sp.hessian(d["canonical_density_Hamiltonian"],tuple(d["canonical_density_phase"])),
        "actual_density_generator_preserves_the_equal_time_CCR_form":M*Omega+Omega*M.T,
        "packet_map_has_exact_k_scaled_symplectic_form":E*Omega*E.T-k*Omega,
        "canonical_density_probe_force_matches_actual_packet_force":
            E*source-parent.data()["packet_chart_source_vector"],
        "density_phase_observable_is_exact_original_relational_scalar":
            d["relational_observable_row"]*d["canonical_density_phase"]-sp.Matrix([parent.parent.chi]),
        "both_density_and_packet_maps_are_mutual_inverses":
            E*d["packet_to_canonical_density_map"]-sp.eye(4)}
    return {name:value.applyfunc(sp.factor) for name,value in rows.items()}
