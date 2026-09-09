"""Action-normalized eigenbasis and exact off-diagonal adiabatic transport."""
from functools import cache

import sympy as sp

ell,kc,km,wc,wm=sp.symbols("ell K_clock K_matter omega_clock omega_matter",positive=True)
variables=(ell,kc,km,wc,wm)
velocities=sp.symbols("ell_dot K_clock_dot K_matter_dot omega_clock_dot omega_matter_dot",real=True)


@cache
def data():
    ec,em=sp.Matrix([1,-ell]),sp.Matrix([0,1])
    T=sp.Matrix([[1,0],[-ell,1]])
    K=T.T.inv()*sp.diag(kc,km)*T.inv()
    G=T.T.inv()*sp.diag(kc*wc*wc,km*wm*wm)*T.inv()
    zero=sp.zeros(2)
    Om=zero.row_join(sp.eye(2)).col_join((-sp.eye(2)).row_join(zero))
    J=zero.row_join(K.inv()).col_join((-G).row_join(zero))
    columns=[]
    for sign in (1,-1):
        for v,kinetic,frequency in ((ec,kc,wc),(em,km,wm)):
            columns.append(v.col_join(sign*sp.I*frequency*K*v)/sp.sqrt(2*frequency*kinetic))
    S=sp.Matrix.hstack(*columns)
    Sinv=S.inv().applyfunc(sp.factor)
    Sd=sum((S.diff(x)*xd for x,xd in zip(variables,velocities,strict=True)),sp.zeros(4))
    e00,e01,e10,e11=sp.symbols("E00 E01 E10 E11",real=True)
    E=sp.Matrix([[e00,e01],[e10,e11]])
    L0=E.row_join(zero).col_join(zero.row_join(-E.T))
    B0=(Sinv*L0*S-Sinv*Sd).applyfunc(sp.factor)
    return {"kinetic":K,"coordinate_gradient":G,"symplectic_form":Om,"J":J,
            "S":S,"S_inverse":Sinv,"S_dot":Sd,"L0":L0,"B0":B0,
            "Lambda":sp.diag(wc,wm,-wc,-wm),"Krein_signature":sp.diag(-1,-1,1,1)}


@cache
def checks():
    d=data()
    return {name:value.applyfunc(sp.factor) for name,value in {
        "complete_mode_eigenbasis":d["J"]*d["S"]-sp.I*d["S"]*d["Lambda"],
        "exact_basis_inverse":d["S_inverse"]*d["S"]-sp.eye(4),
        "positive_and_negative_action_normalization":
            sp.I*d["S"].conjugate().T*d["symplectic_form"]*d["S"]-d["Krein_signature"],
        "leading_diagonal_transport_identically_zero":sp.diag(*(d["B0"][i,i] for i in range(4))),
        "leading_transport_Krein_symmetry":
            d["B0"].conjugate().T*d["Krein_signature"]+d["Krein_signature"]*d["B0"]}.items()}
