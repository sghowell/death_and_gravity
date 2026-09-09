"""Literal new source/contact term and exact nonlinear ADM Hamiltonian."""
from functools import cache

import sympy as sp
from p8_affine_aligned import alignment
from p8_auxiliary_neighborhood import model as old
from p8_coupled_vertices import coefficients, taylor


@cache
def data():
    g=old.generic()
    before=g["Hamiltonian_over_hat_volume"]
    after=before.subs({g["d"]:0,g["c"]:0},simultaneous=True)
    p=g["p0"]+old.dp
    difference=old.N*((p-g["b"])*g["d"]*old.j/(2*g["a"])
                       -g["d"]**2*old.j**2/(4*g["a"])+g["c"]*old.j)
    return {"old_Hamiltonian":before,"new_Hamiltonian":after,
            "new_minus_old_Hamiltonian":difference,
            "new_trace_reconstruction":(p-g["b"])/(2*g["a"]),
            "new_temporal_reconstruction":-g["gt"]*old.j/g["U"],
            "unchanged_scalar_and_mass_coefficients":old.coefficients(),
            "changed_terms":"Only the normal source and its completing contact are changed. The full curled one-form kinetic operator, mass functions, physical matter metric, scalar boundary primitive and margin remain fixed."}


@cache
def source_jets():
    eps,n,k=sp.symbols("field_degree lapse_amplitude trace_amplitude",real=True)
    bg=old.coefficients()["background"]
    h,H=bg["h"],bg["H"]
    QNN=taylor.derivatives(coefficients.boundary_jets()["Q"])[2]
    N=1+eps*n
    S=(N**-2-1)/h*(3*H+eps*k-3*H/N)+sp.Rational(3,2)*QNN*eps**2*n*n/(2*N)
    expansion=sp.series(S,eps,0,3).removeO().expand()
    return {"actual_source_expansion":expansion,"Q_second_lapse_derivative":QNN,
            "zero":sp.factor(expansion.coeff(eps,0)),
            "first":sp.factor(expansion.coeff(eps,1)),
            "second":sp.factor(expansion.coeff(eps,2)),
            "expected_second":sp.factor(-2*n*(k+3*H*n)/h-sp.Rational(9,4)*sp.diff(h,old.u)*n*n/h**3)}


@cache
def checks():
    d=data()
    T,K,a,U,gt,p,b,f=sp.symbols("normal_vector trace trace_a positive_volume positive_gamma trace_p trace_b scalar_f",real=True)
    delta,c=sp.symbols("source_delta source_c",real=True)
    before=a*K*K+b*K+f+U*(T-delta*K-c)**2/(2*gt)
    source=delta*K+c
    added=U*(T*source-source**2/2)/gt
    after=a*K*K+b*K+f+U*T*T/(2*gt)
    Knew=(p-b)/(2*a)
    Hnew=sp.factor(p*Knew-after.subs(K,Knew))
    expected=(p-b)**2/(4*a)-f-U*T*T/(2*gt)
    eps=sp.Symbol("field_degree",real=True)
    c0,c1,t1,t2,s2,s3=sp.symbols("coefficient0 coefficient1 vector1 vector2 source2 source3",real=True)
    change=sp.expand((c0+eps*c1)*((eps*t1+eps**2*t2)*(eps**2*s2+eps**3*s3)
                                -(eps**2*s2+eps**3*s3)**2/2))
    out={"literal_normal_mass_source_and_contact_replacement":sp.factor(before+added-after),
         "new_trace_Legendre_map":sp.factor(Hnew-expected),
         "new_nonlinear_Hamiltonian_difference":sp.factor(d["new_Hamiltonian"]-d["old_Hamiltonian"]-d["new_minus_old_Hamiltonian"]),
         "full_unchanged_source_background":alignment.rolling()["zero_order"],
         "full_unchanged_source_first_variation":alignment.rolling()["first_variation"],
         "actual_Q_fixed_basepoint_second_derivative":sp.factor(source_jets()["Q_second_lapse_derivative"]+3*sp.diff(old.coefficients()["background"]["h"],old.u)/old.coefficients()["background"]["h"]**3),
         "source_second_order_retains_original_Q_jet":sp.factor(source_jets()["second"]-source_jets()["expected_second"]),
         "source_zero_order":source_jets()["zero"],"source_first_order":source_jets()["first"],
         "new_density_cubic_change":sp.factor(change.coeff(eps,3)-c0*t1*s2),
         "new_density_quartic_change":sp.factor(change.coeff(eps,4)-c0*(t1*s3+t2*s2-s2*s2/2)-c1*t1*s2)}
    for degree in range(3):
        out[f"exact_action_agreement_through_field_degree_{degree}"]=change.coeff(eps,degree)
    g=old.generic()
    force=sp.diff(Hnew,T)-sp.Symbol("divergence",real=True)
    solution=-gt*sp.Symbol("divergence",real=True)/U
    out["new_temporal_constraint_at_fixed_canonical_trace"]=sp.factor(force.subs(T,solution))
    out["old_and_new_homogeneous_light_Hamiltonians_identical"]=sp.factor(d["new_minus_old_Hamiltonian"].subs(old.j,0))
    out["new_source_free_temporal_pivot"]=sp.factor(sp.diff(Hnew,T,2)+U/gt)
    out["unchanged_background_scalar_Hamiltonian"]=sp.factor(
        (d["new_Hamiltonian"]-d["old_Hamiltonian"]).subs({g["d"]:0,g["c"]:0}))
    out["new_complete_Hamiltonian_has_exact_vector_parity"]=sp.factor(
        d["new_Hamiltonian"].subs(old.j,-old.j)-d["new_Hamiltonian"])
    return out
