"""Local relational forcing and a polynomial low-frequency Hamiltonian chart."""
from functools import cache

import sympy as sp
from p8_proca_nearby_cones import principal as parent
from p8_proca_nearby_packets import finite

v,p=sp.symbols("metric_scalar trace_scalar_momentum",real=True)
R,k=finite.R,finite.k
N,e,probe=sp.symbols("physical_lapse conformal_scale local_probe",positive=True)


@cache
def data():
    a,m,g,r,rn,h,ell=parent.a,parent.m,parent.g,parent.r,parent.rn,parent.h,parent.ell
    q,H,alpha,beta=parent.q,parent.H,parent.alpha,parent.beta
    chi,P=parent.chi,parent.P
    force=alpha*p+beta*(P-3*ell*v)+4*rn*q*v
    ham=m*(P-3*ell*v)**2+sp.Rational(2,3)*a*ell*p*chi
    ham+=(g*q-a*ell**2)*chi**2+2*r*q*v**2-force**2/(2*h)
    phase=sp.Matrix([v,chi,p,P])
    zero=sp.zeros(2)
    omega=zero.row_join(sp.eye(2)).col_join((-sp.eye(2)).row_join(zero))
    hessian=sp.hessian(ham,tuple(phase))
    generator=omega*hessian-sp.diag(0,0,3*H,3*H)
    B=parent.formula()["B"]
    change=R**sp.Rational(3,2)*sp.Matrix([
        [0,0,-k/(2*q),0],[0,k,0,0],[2*q,0,B/2,0],[0,0,0,1]])
    Bsymbol=sp.Symbol("momentum_shift_coefficient",real=True)
    Bdot=parent.formula()["B_dot"]
    generic_change=change.subs(B,Bsymbol)
    change_dot=generic_change.diff(R)*H*R+generic_change.diff(q)*(-2*H*q)
    change_dot+=generic_change.diff(Bsymbol)*Bdot
    change_dot=change_dot.subs(Bsymbol,B)
    volume=R**3*N*e**3
    forcing=sp.Matrix([0,0,0,N*e**3])
    return {"phase":phase,"Hamiltonian":ham,"Hessian":hessian,
            "weighted_first_order_generator":generator,
            "polynomial_q_coefficients":[generator.applyfunc(lambda x,j=j:sp.expand(x).coeff(q,j)) for j in range(3)],
            "constant_symplectic_form":omega,
            "old_to_packet_canonical_map":change,"old_to_packet_map_dot":change_dot,
            "local_relational_source_action_density":volume*probe*chi,
            "physical_volume_density":volume,
            "original_chart_source_vector":forcing,
            "packet_chart_source_vector":R**sp.Rational(3,2)*forcing,
            "packet_chart_observable_row":sp.Matrix([[0,1/(k*R**sp.Rational(3,2)),0,0]]),
            "no_small_frequency_Legendre_inverse_or_division_by_Hubble":True,
            "q_zero_extension_used_only_for_bounded_frequency_Fourier_estimate":True}


@cache
def checks():
    d=data()
    q,H=parent.q,parent.H
    old=parent.data()["Hamiltonian"]+H*parent.b*parent.Pb
    recovered=old.subs({parent.b:-p/(2*q),parent.Pb:2*q*v},simultaneous=True)
    D=d["old_to_packet_canonical_map"]
    transformed=(d["old_to_packet_map_dot"]+D*d["weighted_first_order_generator"])*D.inv()
    fd=finite.data()
    expected=k*fd["J"]+fd["L0"]+fd["L1"]/k+fd["L2"]/k**2+fd["L3"]/k**3
    polynomial=sum((entry*q**j for j,entry in enumerate(d["polynomial_q_coefficients"])),sp.zeros(4))
    out={
        "literal_original_scalar_Hamiltonian_recovered_before_momentum_chart":recovered-d["Hamiltonian"],
        "low_frequency_generator_is_polynomial_degree_at_most_two_in_q":
            d["weighted_first_order_generator"]-polynomial,
        "full_time_dependent_canonical_map_reproduces_exact_packet_equations":
            transformed.subs(q,k*k/R**2)-expected,
        "covariant_local_source_reproduces_weighted_canonical_matter_force":
            sp.diff(d["local_relational_source_action_density"],parent.chi)/R**3-probe*N*e**3,
        "momentum_shift_and_volume_change_retain_exact_local_source_vector":
            D*d["original_chart_source_vector"]-d["packet_chart_source_vector"],
        "packet_reconstruction_is_exact_original_relational_coordinate":
            d["packet_chart_observable_row"]*D*d["phase"]-sp.Matrix([parent.chi]),
        "physical_spacetime_volume_matches_hat_coordinates":d["physical_volume_density"]-N*(e*R)**3,
        "weighted_original_equations_keep_canonical_volume_damping":
            d["weighted_first_order_generator"].T*d["constant_symplectic_form"]
            +d["constant_symplectic_form"]*d["weighted_first_order_generator"]
            +3*H*d["constant_symplectic_form"]}
    return {name:value.applyfunc(sp.factor) if isinstance(value,sp.MatrixBase) else sp.factor(value)
            for name,value in out.items()}
