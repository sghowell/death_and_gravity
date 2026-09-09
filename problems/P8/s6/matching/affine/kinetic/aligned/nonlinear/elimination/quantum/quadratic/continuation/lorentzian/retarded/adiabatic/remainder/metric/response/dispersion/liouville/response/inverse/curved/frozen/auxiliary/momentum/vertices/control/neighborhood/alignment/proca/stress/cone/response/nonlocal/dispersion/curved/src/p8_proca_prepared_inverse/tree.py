"""Actual unchanged quadratic tree, fixed matter charge and adapted derivative orders."""
from functools import cache

import sympy as sp
from p8_prepared_volterra import tree as parent

from . import coordinates

u,H=parent.u,parent.H
eta,scale=coordinates.eta,coordinates.scale
delta=parent.delta


def weighted(value):
    return parent.time(value)+3*H*value


def euler(value,field):
    return sp.expand(sum((-1)**j*iterate(weighted,sp.diff(value,sp.diff(field,u,j)),j) for j in range(3)))


def iterate(function,value,count):
    for _ in range(count):
        value=function(value)
    return value


@cache
def data():
    original=parent.data()
    ed=sp.diff(eta,u)
    hatv=scale+H*eta-delta*ed
    density=sp.expand(original["effective_two_metric_density"].subs({parent.n:ed,parent.v:hatv},simultaneous=True).doit())
    currents=sp.ImmutableMatrix([euler(density,field) for field in (eta,scale)])
    physical=original["physical_metric_Euler_currents"].subs(
        {parent.n:ed,parent.zeta:scale+H*eta},simultaneous=True).doit()
    transformed=sp.ImmutableMatrix([-weighted(physical[0])+H*physical[1],physical[1]])
    coefficients={str(i)+str(j):{order:sp.factor(sp.diff(currents[i],sp.diff(field,u,order))) for order in range(5)}
                  for i in range(2) for j,field in enumerate((eta,scale))}
    return {"actual_matter_reduced_adapted_density":density,"actual_adapted_tree_currents":currents,
            "independent_transformed_physical_tree_currents":transformed,
            "complete_local_derivative_coefficients":coefficients,
            "highest_time_channel_coefficient":-6*delta**2,
            "clock_delta":1/(2*(1+u*u)**3),"clock_delta_lower":sp.Rational(32,125),
            "clock_delta_upper":sp.Rational(1,2),
            "highest_time_channel_absolute_lower":sp.Rational(6144,15625),
            "instantaneous_time_channel_inverse_absolute_upper":sp.Rational(15625,6144),
            "maximum_local_derivative_orders":sp.ImmutableMatrix([[4,3],[3,2]]),
            "matter_reconstruction_rate":-3*parent.ell*hatv-parent.w*ed,
            "prepared_reconstruction":"n=eta'; physical v=w+H eta; clock vhat=w+H eta-delta eta'",
            "background_scale_L_squared_factored_out_of_tree":True}


@cache
def checks():
    item=data()
    currents=item["actual_adapted_tree_currents"]
    out={"actual_tree_charge_reduction_replayed":parent.checks()["effective_two_current_variation"],
         "actual_literal_tree_action_replayed":parent.checks()["literal_density_replays_S6_57"],
         "independent_adapted_Euler_equals_transformed_actual_tree":sp.ImmutableMatrix(
             (currents-item["independent_transformed_physical_tree_currents"]).applyfunc(sp.simplify)),
         "actual_fourth_time_channel_sign":sp.factor(sp.diff(currents[0],sp.diff(eta,u,4))+6*delta**2),
         "actual_highest_time_scale_cross":sp.factor(sp.diff(currents[0],sp.diff(scale,u,3))-6*delta),
         "actual_highest_scale_time_cross":sp.factor(sp.diff(currents[1],sp.diff(eta,u,3))+6*delta),
         "actual_second_scale_channel":sp.factor(sp.diff(currents[1],sp.diff(scale,u,2))-6),
         "time_scale_channel_has_no_fourth_derivative":sp.diff(currents[0],sp.diff(scale,u,4)),
         "scale_time_channel_has_no_fourth_derivative":sp.diff(currents[1],sp.diff(eta,u,4)),
         "scale_scale_channel_has_no_third_or_fourth_derivative":sp.ImmutableMatrix(
             [sp.diff(currents[1],sp.diff(scale,u,j)) for j in (3,4)]),
         "highest_time_channel_bound_from_actual_clock_map":6*item["clock_delta_lower"]**2-item["highest_time_channel_absolute_lower"],
         "instantaneous_time_inverse_bound_reciprocal":item["instantaneous_time_channel_inverse_absolute_upper"]
             *item["highest_time_channel_absolute_lower"]-1,
         "actual_clock_delta_maximum_at_center":item["clock_delta"].subs(u,0)-item["clock_delta_upper"],
         "actual_clock_delta_minimum_at_both_endpoints":sp.ImmutableMatrix(
             [item["clock_delta"].subs(u,t)-item["clock_delta_lower"] for t in (-sp.Rational(1,2),sp.Rational(1,2))]),
         "actual_clock_delta_monotone_on_positive_half":sp.factor(
             sp.diff(item["clock_delta"],u)+3*u/(1+u*u)**4)}
    for key in ("01","10","11"):
        top=3 if key in ("01","10") else 2
        for order in range(top+1,5):
            out["complete_derivative_inventory_zero_"+key+"_"+str(order)]=item["complete_local_derivative_coefficients"][key][order]
    return out
