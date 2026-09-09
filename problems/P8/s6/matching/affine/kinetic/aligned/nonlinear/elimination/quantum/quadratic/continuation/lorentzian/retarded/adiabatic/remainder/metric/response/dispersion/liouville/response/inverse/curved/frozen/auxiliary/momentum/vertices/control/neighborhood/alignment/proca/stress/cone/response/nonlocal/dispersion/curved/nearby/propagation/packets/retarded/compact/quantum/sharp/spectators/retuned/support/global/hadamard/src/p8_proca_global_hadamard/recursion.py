"""Finite native tests of the all-order induction, including derivative controls."""
from functools import cache

import sympy as sp

from . import riccati

clean=riccati.clean


@cache
def data():
    t=sp.Symbol("formal_recursion_time",real=True)
    omega=sp.diag(2,3)
    H0=sp.diag(4,9,1,1)
    E=sp.Matrix([[t,sp.Rational(1,3)],[sp.Rational(-1,4),-t/2]])
    zero=sp.zeros(2)
    H1=zero.row_join(E.T).col_join(E.row_join(zero))
    H2=sp.Matrix([[sp.Rational(1,2),sp.Rational(1,5),0,0],
        [sp.Rational(1,5),sp.Rational(-1,3),0,0],
        [0,0,sp.Rational(-1,2),sp.Rational(1,7)],
        [0,0,sp.Rational(1,7),sp.Rational(2,3)]])
    hs=[H0,H1,H2]
    rs=[-sp.I*omega]
    for n in range(1,5):
        rs.append(riccati.next_coefficient(n,rs,hs,lambda value:value.diff(t),sp.eye(2),omega))
    R=sum((riccati.rho**j*v for j,v in enumerate(rs)),sp.zeros(2))
    H=sum((riccati.rho**j*v for j,v in enumerate(hs)),sp.zeros(4))
    defect=riccati.rho*R.diff(t)+R*H[2:,2:]*R+R*H[2:,:2]+H[:2,2:]*R+H[:2,:2]
    wrong2=riccati.next_coefficient(2,rs[:2],hs,lambda value:sp.zeros(2),sp.eye(2),omega)
    return {"time":t,"complete_scaled_Hamiltonian":hs,"four_computed_corrections":rs[1:],
        "complete_defect_coefficients_through_four":[clean(defect.applyfunc(
            lambda value,n=n:sp.expand(value).coeff(riccati.rho,n))) for n in range(5)],
        "missing_time_derivative_second_correction":wrong2,
        "missing_derivative_control_difference":clean(wrong2-rs[2])}


@cache
def checks():
    d=data()
    out={f"universal_Riccati_induction_native_order_{n}":v
        for n,v in enumerate(d["complete_defect_coefficients_through_four"])}
    for n,value in enumerate(d["four_computed_corrections"],1):
        out[f"native_Riccati_coefficient_{n}_is_symmetric"]=clean(value-value.T)
    return out


@cache
def gates():
    return {"deleting_time_derivative_changes_the_second_coefficient":any(
        v!=0 for v in data()["missing_derivative_control_difference"])}


@cache
def mode_data():
    c,a=sp.symbols("positive_cone_speed positive_conformal_scale",positive=True)
    lam=sp.diag(c/a,1/a,-c/a,-1/a)
    B=sp.Matrix(4,4,lambda i,j:sp.Symbol(f"transport_{i}_{j}"))
    T=sp.zeros(4)
    for i in range(4):
        for j in range(4):
            if i!=j:
                T[i,j]=-sp.I*B[i,j]/(lam[i,i]-lam[j,j])
    diagonal=sp.diag(*(B[i,i] for i in range(4)))
    return {"signed_positive_frequency_matrix":lam,"complete_generic_order_zero_transport":B,
        "first_four_mode_change":T,"retained_diagonal_transport":diagonal,
        "opposite_sign_frequency_gaps":[2*c/a,(c+1)/a,2/a],
        "same_sign_gap_magnitude":(1-c)/a,
        "initial_slab_compact_speed_margin_lower":sp.Rational(1,2000),
        "initial_slab_scale_upper":sp.Rational(17,16)**2,
        "initial_slab_same_sign_frequency_gap_lower":sp.Rational(1,4000)/(sp.Rational(17,16)**2)}


@cache
def mode_checks():
    d=mode_data()
    lam,B,T,D=d["signed_positive_frequency_matrix"],d["complete_generic_order_zero_transport"],d["first_four_mode_change"],d["retained_diagonal_transport"]
    return {"full_four_mode_change_removes_only_off_diagonal_transport":
        clean(-sp.I*(lam*T-T*lam)+B-D),
        "first_mode_change_has_no_assigned_diagonal_terms":sp.diag(*(T[i,i] for i in range(4)))}
