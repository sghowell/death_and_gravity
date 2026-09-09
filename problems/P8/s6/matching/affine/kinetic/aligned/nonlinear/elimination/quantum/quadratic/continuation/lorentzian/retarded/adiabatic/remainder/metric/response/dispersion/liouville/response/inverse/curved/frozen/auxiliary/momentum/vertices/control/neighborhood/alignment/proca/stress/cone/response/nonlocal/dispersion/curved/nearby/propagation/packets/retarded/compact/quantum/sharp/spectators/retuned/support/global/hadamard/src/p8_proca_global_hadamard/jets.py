"""Actual initial-slab bridge and nontrivial retuned center Riccati jets."""
from functools import cache

import sympy as sp

from . import canonical, riccati

clean=riccati.clean


@cache
def data():
    d=canonical.data()
    model=canonical.model
    u=model.u
    bg=model.background()
    V=sp.Matrix([[1,0],[bg["w"]/bg["lam"],1]])
    mass=2*bg["J_new"]/bg["lam"]**2
    omega=sp.diag(sp.sqrt(bg["J"]/bg["J_new"])/bg["a"],1/bg["a"])
    B=clean(V.inv().T*sp.diag(mass,1)*omega*V.inv())
    R0=-sp.I*B
    sub=model.substitution()
    H=[clean((-d["constant_symplectic_form"]*d["complete_Laurent_coefficients"][j]).subs(
        sub,simultaneous=True)) for j in (1,0,-1)]
    V0,omega0=V.subs(u,0),omega.subs(u,0)
    R00=R0.subs(u,0)
    R0d=R0.diff(u).subs(u,0)
    R0dd=clean(R0.diff(u,2).subs(u,0))
    H0=H[0].subs(u,0)
    H1=H[1].subs(u,0)
    H1d=clean(H[1].diff(u).subs(u,0))
    H2=H[2].subs(u,0)
    force1=R0d+R00*H1[2:,2:]*R00+R00*H1[2:,:2]+H1[:2,2:]*R00+H1[:2,:2]
    R1=riccati.solve_sylvester(-force1,V0,omega0)
    force1d=R0dd+R00*H1d[2:,:2]+H1d[:2,2:]*R00
    R1d=riccati.solve_sylvester(-force1d,V0,omega0)
    force2=R1d+R00*H2[2:,2:]*R00+R00*H2[2:,:2]+H2[:2,2:]*R00+H2[:2,:2]
    R2=riccati.solve_sylvester(-force2,V0,omega0)
    return {"initial_time":sp.Integer(0),"actual_global_negative_graph_leading_symbol":R0,
        "actual_global_configuration_eigenbasis":V,"actual_global_positive_frequencies":omega,
        "actual_global_clock_kinetic":mass,"actual_complete_scaled_Hamiltonian_coefficients":H,
        "center_leading_graph":R00,"center_first_time_jet_leading_graph":R0d,
        "center_second_time_jet_leading_graph":R0dd,
        "center_Hamiltonian_leading":H0,"center_Hamiltonian_order_zero":H1,
        "center_Hamiltonian_order_zero_first_time_jet":H1d,"center_Hamiltonian_inverse_order":H2,
        "center_Riccati_first_coefficient":R1,"center_Riccati_first_coefficient_time_jet":R1d,
        "center_Riccati_second_coefficient":R2,
        "center_first_recursion_forcing":force1,"center_second_recursion_forcing":force2,
        "center_naive_leading_graph_second_defect":clean(force2-R1d),
        "center_positive_frequencies":omega0,"center_configuration_eigenbasis":V0}


@cache
def checks():
    d=data()
    R0=d["actual_global_negative_graph_leading_symbol"]
    H=d["actual_complete_scaled_Hamiltonian_coefficients"]
    V=d["actual_global_configuration_eigenbasis"]
    R00=d["center_leading_graph"]
    H0=d["center_Hamiltonian_leading"]
    R1=d["center_Riccati_first_coefficient"]
    R2=d["center_Riccati_second_coefficient"]
    c0=sp.sqrt(sp.Rational(1199,1215))
    expectedB=sp.Matrix([[sp.Rational(243,20)*c0+sp.Rational(1,100),sp.Rational(1,10)],
        [sp.Rational(1,10),1]])
    return {
        "fresh_complete_global_leading_graph_solves_Riccati":clean(R0*H[0][2:,2:]*R0+H[0][:2,:2]),
        "fresh_global_negative_graph_evolution_kinetic_bridge":clean(
            H[0][2:,2:]*V.inv().T*sp.diag(d["actual_global_clock_kinetic"],1)-V),
        "actual_center_graph_keeps_coupled_kinetic_and_speed":clean(R00+sp.I*expectedB),
        "actual_center_leading_graph_first_jet_is_zero":d["center_first_time_jet_leading_graph"],
        "actual_center_complete_order_zero_Hamiltonian_is_zero":d["center_Hamiltonian_order_zero"],
        "global_complete_order_zero_Hamiltonian_has_only_mixed_blocks":H[1][:2,:2].row_join(H[1][2:,2:]),
        "actual_center_first_Riccati_coefficient_vanishes":R1,
        "actual_center_first_Riccati_recursion":clean(R00*H0[2:,2:]*R1+R1*H0[2:,2:]*R00+d["center_first_recursion_forcing"]),
        "actual_center_second_Riccati_recursion_retains_derivative":clean(R00*H0[2:,2:]*R2+R2*H0[2:,2:]*R00+d["center_second_recursion_forcing"]),
        "actual_second_Riccati_coefficient_is_symmetric":clean(R2-R2.T),
        "actual_second_Riccati_coefficient_is_pure_imaginary":clean(R2+sp.conjugate(R2)),
        "global_exact_Laurent_orders_below_inverse_one_vanish":canonical.data()["complete_Laurent_coefficients"][-2].row_join(
            canonical.data()["complete_Laurent_coefficients"][-3]),
    }


@cache
def gates():
    d=data()
    return {"actual_first_coefficient_time_derivative_is_not_zero":any(v!=0 for v in d["center_Riccati_first_coefficient_time_jet"]),
        "actual_second_coefficient_is_not_zero":any(v!=0 for v in d["center_Riccati_second_coefficient"]),
        "naive_leading_graph_alone_has_nonzero_inverse_order_defect":any(v!=0 for v in d["center_naive_leading_graph_second_defect"])}
