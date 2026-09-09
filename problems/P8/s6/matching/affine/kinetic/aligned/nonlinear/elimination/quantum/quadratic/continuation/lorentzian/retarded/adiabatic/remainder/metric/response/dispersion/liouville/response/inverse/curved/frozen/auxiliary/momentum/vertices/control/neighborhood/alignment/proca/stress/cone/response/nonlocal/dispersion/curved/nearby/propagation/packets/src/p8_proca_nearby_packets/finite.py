"""Exact finite-q Hamiltonian normalization, without a Legendre inverse."""
from functools import cache

import sympy as sp
from p8_proca_nearby_cones import principal as parent

q,H=parent.q,parent.H
R,k=sp.symbols("hat_scale comoving_frequency",positive=True)


@cache
def data():
    ham=parent.data()["Hamiltonian"]
    mat=sp.hessian(ham,(parent.b,parent.chi,parent.Pb,parent.P))
    A,E,C=mat[2:,2:],mat[2:,:2],mat[:2,:2]
    B=parent.formula()["B"]
    Bd=parent.formula()["B_dot"]
    shift=sp.diag(q*B,0)
    Enew=(E+A*shift).applyfunc(sp.factor)
    Cnew=(C+E.T*shift+shift*E+shift*A*shift+sp.diag(q*(Bd+H*B),0)).applyfunc(sp.factor)
    A0=A.applyfunc(parent.infinity_limit)
    A1=((A-A0)*q).applyfunc(parent.infinity_limit)
    A2=((A-A0-A1/q)*q*q).applyfunc(parent.infinity_limit)
    E0=Enew.applyfunc(parent.infinity_limit)
    E1=((Enew-E0)*q).applyfunc(parent.infinity_limit)
    G=(Cnew/q).applyfunc(parent.infinity_limit)
    C0=(Cnew-q*G).applyfunc(parent.infinity_limit)
    zero=sp.zeros(2)
    block=lambda a,b,c,d:a.row_join(b).col_join(c.row_join(d))
    J=block(zero,A0,-G/R**2,zero)
    L0=block(E0+3*H*sp.eye(2)/2,zero,zero,-E0.T-3*H*sp.eye(2)/2)
    L1=block(zero,R**2*A1,-C0,zero)
    L2=block(R**2*E1,zero,zero,-R**2*E1.T)
    L3=block(zero,R**4*A2,zero,zero)
    return {"original_A":A,"shifted_E":Enew,"shifted_C":Cnew,"momentum_shift":shift,
            "A0":A0,"A1":A1,"A2":A2,"E0":E0,"E1":E1,"G":G,"C0":C0,
            "J":J,"L0":L0,"L1":L1,"L2":L2,"L3":L3,
            "energy_symmetrizer":block(G/R**2,zero,zero,A0),
            "symplectic_form":block(zero,sp.eye(2),-sp.eye(2),zero)}


@cache
def checks():
    d=data()
    out={
        "exact_finite_q_momentum_Hessian":d["original_A"]-d["A0"]-d["A1"]/q-d["A2"]/q**2,
        "exact_finite_q_momentum_coordinate_block":d["shifted_E"]-d["E0"]-d["E1"]/q,
        "exact_finite_q_coordinate_block":d["shifted_C"]-q*d["G"]-d["C0"],
        "leading_kinetic_inverse_is_actual_parent":d["A0"]*parent.formula()["kinetic"]-sp.eye(2),
        "shift_keeps_actual_Euler_gradient":d["G"]-parent.formula()["Euler_gradient"],
        "uniform_energy_principal_cancellation":d["energy_symmetrizer"]*d["J"]+d["J"].T*d["energy_symmetrizer"]}
    zero=sp.zeros(2)
    E=d["shifted_E"]+3*H*sp.eye(2)/2
    exact=E.row_join(k*d["original_A"]).col_join((-d["shifted_C"]/k).row_join(-E.T))
    out["all_inverse_frequency_orders_after_volume_rescaling"]=(
        exact.subs(q,k*k/R**2)-k*d["J"]-d["L0"]-d["L1"]/k-d["L2"]/k**2-d["L3"]/k**3)
    for name in ("J","L0","L1","L2","L3"):
        out["constant_symplectic_structure_"+name]=d[name].T*d["symplectic_form"]+d["symplectic_form"]*d[name]
    out["volume_rescaling_removes_both_weighted_damping_halves"]=zero+sp.eye(2)*(-3*H+3*H/2+3*H/2)
    return {name:value.applyfunc(sp.factor) for name,value in out.items()}
