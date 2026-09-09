"""Exact full chart-to-phase intertwining and polynomial endpoint weights."""
from functools import cache

import sympy as sp

from . import charts, model, phase


def rational_degree(value,variable):
    num,den=sp.fraction(sp.cancel(value))
    if num==0:
        return -sp.oo
    return sp.degree(num,variable)-sp.degree(den,variable)


@cache
def data():
    d=phase.data()
    out={}
    for name,E in d["chart_velocity_maps_before_q_to_comoving_substitution"].items():
        c=charts.data()[name]
        predicted=sp.zeros(4)
        predicted[:2,2:]=model.k*sp.eye(2)
        predicted[2:,:2]=-c["complete_velocity_force"]/model.k
        predicted[2:,2:]=-c["complete_velocity_friction"]
        Ek=E.subs(model.q,model.k**2/model.a**2).applyfunc(sp.factor)
        inverse=Ek.inv().applyfunc(sp.factor)
        o=model.old
        Dk=model.k**2*o.lam**2-model.a**2*(o.J+model.margin.delta_J+o.w**2/2)
        allowed=(model.a,model.k,o.J+model.margin.delta_J,o.theta) if name=="unitary" else (
            model.a,model.k,o.J+model.margin.delta_J,o.lam,Dk)
        matrices=(Ek,inverse,c["complete_velocity_friction"],c["complete_bounded_force_remainder"])
        factors={factor for matrix in matrices for value in matrix
            for factor,power in sp.factor_list(sp.denom(sp.cancel(value.subs(model.q,model.k**2/model.a**2))))[1]}
        unexpected=[factor for factor in factors if not any(sp.cancel(factor/pivot).is_number for pivot in allowed)]
        out[name]={"exact_scaled_velocity_generator":predicted,
            "original_phase_to_scaled_velocity":Ek,"scaled_velocity_to_original_phase":inverse,
            "forward_map_maximum_frequency_degree":max(rational_degree(v,model.k) for v in Ek),
            "inverse_map_maximum_frequency_degree":max(rational_degree(v,model.k) for v in inverse),
            "complete_friction_maximum_q_degree":max(rational_degree(v,model.q) for v in c["complete_velocity_friction"]),
            "complete_force_remainder_maximum_q_degree":max(rational_degree(v,model.q) for v in c["complete_bounded_force_remainder"]),
            "all_endpoint_and_remainder_denominator_factors":sorted(map(str,factors)),
            "declared_regular_chart_pivot_factors":allowed,
            "unexpected_denominator_factors":sorted(map(str,unexpected)),
            "complete_intertwining_residual":(model.dt(E)+E*d["regular_density_generator"]-predicted*E).applyfunc(sp.factor)}
    return out


@cache
def checks():
    out={}
    for name,d in data().items():
        out[name+"_full_Euler_and_original_phase_generators_exactly_intertwine"]=d["complete_intertwining_residual"]
        out[name+"_exact_scaled_velocity_map_has_literal_inverse"]=(
            d["original_phase_to_scaled_velocity"]*d["scaled_velocity_to_original_phase"]-sp.eye(4)).applyfunc(sp.factor)
    return out


@cache
def gates():
    rows={}
    for name,d in data().items():
        rows[name+"_complete_friction_is_uniformly_rationally_bounded_at_large_complex_q"]=d["complete_friction_maximum_q_degree"]<=0
        rows[name+"_complete_force_remainder_is_uniformly_rationally_bounded_at_large_complex_q"]=d["complete_force_remainder_maximum_q_degree"]<=0
        rows[name+"_forward_scaled_velocity_map_grows_at_most_quadratically"]=d["forward_map_maximum_frequency_degree"]<=2
        rows[name+"_inverse_scaled_velocity_map_grows_at_most_linearly"]=d["inverse_map_maximum_frequency_degree"]<=1
        rows[name+"_all_endpoint_and_remainder_poles_are_exactly_declared_regular_pivots"]=not d["unexpected_denominator_factors"]
    return {name:bool(value) for name,value in rows.items()}
