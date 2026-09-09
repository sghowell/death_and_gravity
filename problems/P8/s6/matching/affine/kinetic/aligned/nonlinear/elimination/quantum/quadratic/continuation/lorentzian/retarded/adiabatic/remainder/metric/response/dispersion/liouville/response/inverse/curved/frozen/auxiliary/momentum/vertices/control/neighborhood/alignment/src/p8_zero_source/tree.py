"""Finite hard-tree change budget at a separately named scale, not a cutoff."""
from functools import cache

import sympy as sp
from p8_coupled_energy import tree as old

from . import bounds

SCALE=sp.Integer(10**400)


def at_scale(value):
    if isinstance(value,bool) or not isinstance(value,(int,sp.Rational)):
        raise TypeError("Require a positive exact rational M*tau scale")
    scale=sp.Rational(value)
    if scale<=0:
        raise ValueError("Require a positive exact rational M*tau scale")
    d=data()
    cubic=SCALE/scale
    quartic=cubic*cubic
    return {"common_M_tau":scale,
            "new_cubic_transition_bound":cubic*d["new_cubic_transition_bound"],
            "new_connected_quartic_tree_bound":quartic*d["new_connected_quartic_tree_bound"],
            "new_minus_old_cubic_transition_norm_upper":cubic*d["new_minus_old_cubic_transition_norm_upper"],
            "new_minus_old_connected_quartic_tree_norm_upper":quartic*d["new_minus_old_connected_quartic_tree_norm_upper"]}


@cache
def data():
    previous=old.data()
    R=bounds.RATIO
    C3=previous["cubic_transition_numerator"]
    C4=previous["quartic_connected_tree_numerator"]
    return {"same_seven_free_modes_and_hard_window":previous["domain"],
            "new_named_M_tau":SCALE,"new_named_M_tau_power10":400,
            "new_cubic_transition_bound":R**4*C3/SCALE,
            "new_connected_quartic_tree_bound":R**8*C4/SCALE**2,
            "new_minus_old_cubic_transition_norm_upper":(R**4+1)*C3/SCALE,
            "new_minus_old_connected_quartic_tree_norm_upper":(R**8+1)*C4/SCALE**2,
            "comparison_target":old.TARGET,
            "comparison_convention":"Both models are compared at the same new M*tau=10^400 using their identical free canonical phase coordinates, seven free mode columns, original vector covariance and fixed total-momentum hard masks. Nonlinear canonical momenta change; their tangent Darboux identification fixes this finite-order interaction comparison.",
            "not_the_old_M_tau_10_to_24_or_10_to_353_examples":True,
            "not_an_all_orders_or_Wilsonian_cutoff":True}


@cache
def gates():
    d=data()
    return {"new_finite_cubic_block_below_one_e_minus_three":bool(d["new_cubic_transition_bound"]<old.TARGET),
            "new_finite_connected_quartic_block_below_one_e_minus_three":bool(d["new_connected_quartic_tree_bound"]<old.TARGET),
            "finite_cubic_model_change_below_one_e_minus_fourteen":bool(d["new_minus_old_cubic_transition_norm_upper"]<sp.Rational(1,10**14)),
            "finite_connected_quartic_model_change_below_one_e_minus_one_hundred_fifty":bool(d["new_minus_old_connected_quartic_tree_norm_upper"]<sp.Rational(1,10**150))}
