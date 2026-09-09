"""Actual relational-matter residues of the normalized global scalar modes."""
from functools import cache

import sympy as sp
from p8_proca_global_hadamard import canonical
from p8_proca_global_support import model

clean=canonical.clean


@cache
def data():
    mass,a,c,ell,k=sp.symbols("positive_clock_mass positive_scale positive_clock_speed positive_matter_density positive_momentum",positive=True)
    V=sp.Matrix([[1,0],[-ell,1]])
    K=V.inv().T*sp.diag(mass,1)*V.inv()
    omega=sp.diag(c/a,1/a)
    Q=V*sp.diag(sp.sqrt(a/(2*mass*c)),sp.sqrt(a/2))
    P=V.inv().T*sp.diag(sp.sqrt(mass*c/(2*a)),1/sp.sqrt(2*a))
    S=Q.col_join(-sp.I*P)
    Omega=sp.zeros(4)
    Omega[:2,2:]=sp.eye(2)
    Omega[2:,:2]=-sp.eye(2)
    E=sp.Matrix(2,2,lambda i,j:sp.Symbol(f"real_order_zero_transport_{i}_{j}",real=True))
    L=E.row_join(sp.zeros(2)).col_join(sp.zeros(2).row_join(-E.T))
    dots={v:sp.Symbol(str(v)+"_time_derivative",real=True) for v in (mass,a,c,ell)}
    Sd=sum((S.diff(v)*jet for v,jet in dots.items()),sp.zeros(4,2))
    B0=clean(S.conjugate().T*sp.I*Omega*(L*S-Sd))
    chi=clean(S[1,:]/(sp.sqrt(k)*a**sp.Rational(3,2)))
    weights=[sp.factor(2*k*a*a*v*sp.conjugate(v)) for v in chi]
    bg=model.background()
    actual_mass=sp.factor(2*bg["J_new"]/bg["lam"]**2)
    actual_weight=bg["ell"]**2/(actual_mass*sp.sqrt(bg["clock_physical_squared_speed"]))
    lower=1/(100*(sp.Rational(17,16))**12)*sp.Rational(1,25)/(2*(36+sp.Rational(1,50)))
    return {"clock_mass":mass,"scale":a,"clock_speed":c,"matter_density":ell,"momentum":k,
        "leading_configuration_basis":V,"leading_kinetic":K,"leading_positive_frequencies":omega,
        "normalized_negative_frequency_modes":S,"constant_symplectic_form":Omega,
        "full_generic_order_zero_Hamiltonian_generator":L,
        "actual_time_chain_of_mode_basis":Sd,"leading_negative_mode_transport":B0,
        "density_relational_matter_mode_row":chi,"relative_clock_and_matter_spectral_weights":weights,
        "actual_clock_kinetic":actual_mass,"actual_clock_relational_spectral_weight":actual_weight,
        "uniform_clock_relational_weight_lower_on_gamma_slab":lower,
        "center_clock_speed":sp.sqrt(sp.Rational(1199,1215)),
        "center_clock_relative_weight":1/(1215*sp.sqrt(sp.Rational(1199,1215))),
        "leading_density_chi_two_point":"hbar/(2*kappa*a^2*k) times [r_c*clock_phase + matter_phase]",
        "Gamma_clock_principal_amplitude_is_nonzero_not_a_free_matter_projection":True}


@cache
def checks():
    d=data()
    S,Omega=d["normalized_negative_frequency_modes"],d["constant_symplectic_form"]
    K,V,omega=d["leading_kinetic"],d["leading_configuration_basis"],d["leading_positive_frequencies"]
    J=sp.zeros(4)
    J[:2,2:]=K.inv()
    J[2:,:2]=-V.inv().T*sp.diag(d["clock_mass"],1)*omega**2*V.inv()
    bg=model.background()
    actual_K=canonical.data()["principal_kinetic"].subs(model.substitution(),simultaneous=True)
    expected_K=K.subs({d["clock_mass"]:d["actual_clock_kinetic"],d["matter_density"]:bg["ell"]},simultaneous=True)
    r=d["relative_clock_and_matter_spectral_weights"][0]
    return {
        "action_normalized_negative_modes_have_positive_symplectic_norm":clean(S.conjugate().T*sp.I*Omega*S-sp.eye(2)),
        "opposite_frequency_modes_are_symplectically_orthogonal":clean(S.conjugate().T*sp.I*Omega*sp.conjugate(S)),
        "normalized_modes_solve_both_actual_type_principal_frequencies":clean(J*S+sp.I*S*omega),
        "full_order_zero_and_all_basis_jets_have_zero_diagonal_mode_transport":clean(sp.diag(
            *(d["leading_negative_mode_transport"][j,j] for j in range(2)))),
        "actual_global_kinetic_is_the_normalized_mode_kinetic":clean(actual_K-expected_K),
        "clock_relational_spectral_weight_retains_matter_mixing":sp.factor(r-d["matter_density"]**2/(d["clock_mass"]*d["clock_speed"])),
        "matter_relational_spectral_weight_is_one":d["relative_clock_and_matter_spectral_weights"][1]-1,
        "actual_bounce_clock_residue_is_nonzero_exact_expression":sp.simplify(
            d["actual_clock_relational_spectral_weight"].subs(model.u,0)-d["center_clock_relative_weight"]),
        "actual_bounce_relative_clock_weight_times_1215c_is_one":sp.simplify(
            1215*d["center_clock_speed"]*d["center_clock_relative_weight"]-1),
    }


@cache
def gates():
    d=data()
    return {"Gamma_clock_relational_spectral_lower_is_strictly_positive":bool(d["uniform_clock_relational_weight_lower_on_gamma_slab"]>0),
        "actual_center_clock_weight_is_strictly_positive":bool(d["center_clock_relative_weight"]>0),
        "actual_center_clock_sheet_differs_from_matter_metric_null_sheet":bool(0<d["center_clock_speed"]<1)}
