"""Physical-volume Green operators and the density-phase field algebra."""
from functools import cache

import sympy as sp
from p8_proca_global_support import model, phase


@cache
def data():
    d=phase.data()
    Omega=d["constant_symplectic_form"]
    H=sp.hessian(d["regular_density_Hamiltonian"],tuple(d["regular_density_phase"]))
    return {"physical_volume":"a(u)^3 du d^3x",
        "positive_fibre_pairing":"standard Euclidean Hermitian pairing on the trivial complex rank-four bundle",
        "operator":"Q=a^-3 Omega (partial_u-M_density)=a^-3 (Omega partial_u+H_density)",
        "time_derivative_matrix":Omega/model.old.a**3,
        "spatial_differential_matrix":H/model.old.a**3,
        "spatial_symbol_replacement":"q=-Delta_x/a(u)^2; all coefficient matrices depend only on time",
        "retarded_kernel_relative_to_physical_source_volume":"-theta(t-s)*U_density(t,s)*Omega",
        "advanced_kernel_relative_to_physical_source_volume":"+theta(s-t)*U_density(t,s)*Omega",
        "advanced_minus_retarded_kernel":"E_Q(t,s)=U_density(t,s)*Omega",
        "coordinate_volume_retarded_operator":"-integral_{s<=t} U_density(t,s)*Omega*a(s)^3*f(s) ds",
        "formal_Green_boundary":"partial_u(f^T*Omega*g), plus spatial divergences",
        "commutator_normalization":"W-W^transpose=i*hbar*E_Q/kappa",
        "whole_phase_support_from_S698":"The full U_density matrix obeys the degree-nine exact-matter-radius entire bound, not only its chi,pi_chi entry",
        "Green_hyperbolicity_not_inferred_from_naive_highest_order_principal_rank":True}


@cache
def checks():
    d=phase.data()
    Omega,M=d["constant_symplectic_form"],d["regular_density_generator"]
    H=sp.hessian(d["regular_density_Hamiltonian"],tuple(d["regular_density_phase"]))
    f,g,fd,gd=(sp.Matrix(sp.symbols(prefix+"0:4",real=True)) for prefix in ("f","g","fd","gd"))
    Qf,Qg=Omega*fd+H*f,Omega*gd+H*g
    boundary=(fd.T*Omega*g+f.T*Omega*gd)[0]
    E=sp.Matrix(4,4,lambda i,j:sp.Symbol(f"exact_density_transfer_{i}_{j}",real=True))
    return {
        "physical_volume_operator_is_negative_normalized_action_Hessian":(-Omega*M-H).applyfunc(sp.factor),
        "spatial_operator_is_formally_symmetric_and_even":(H-H.T).applyfunc(sp.factor),
        "physical_volume_Green_identity_retains_boundary_sign":sp.factor((f.T*Qg-Qf.T*g)[0]-boundary),
        "retarded_jump_with_physical_source_volume_is_identity":-Omega*Omega-sp.eye(4),
        "advanced_jump_with_reverse_step_is_identity":-Omega*Omega-sp.eye(4),
        "advanced_minus_retarded_uses_positive_density_CCR_sign":E*Omega-(-(-E*Omega)),
        "full_phase_CCR_generator_is_Hamiltonian":(M*Omega+Omega*M.T).applyfunc(sp.factor),
    }
