"""Literal matter-charge elimination, with no lapse or bounce denominator."""
from functools import cache

import sympy as sp
from p8_margin_response import system

u=sp.Symbol("time",real=True)
n,v,matter=(sp.Function(name)(u) for name in ("n","v","matter"))
zeta,delta=(sp.Function(name)(u) for name in ("zeta","delta"))
H,theta,ell,w,Je=(sp.Function(name)(u) for name in ("H","theta","ell","w","J_e"))


def time(value):
    return sp.diff(value,u).subs(sp.diff(ell,u),-3*H*ell)


def euler(L,field):
    return sp.expand(sp.diff(L,field)-time(sp.diff(L,sp.diff(field,u)))-3*H*sp.diff(L,sp.diff(field,u)))


@cache
def data():
    vd,sd=sp.diff(v,u),sp.diff(matter,u)
    L=-3*vd**2+(Je+w*w/2-3*theta**2)*n*n+6*theta*n*vd+sd*sd/2+w*n*sd-3*ell*vd*matter
    charge=sd+w*n+3*ell*v
    effective=-3*vd**2+6*theta*n*vd+(Je-3*theta**2)*n*n-3*ell*w*n*v-sp.Rational(9,2)*ell**2*v**2
    substitution={sd:-3*ell*v-w*n}
    actual=sp.ImmutableMatrix([sp.simplify(euler(L,field).subs(substitution)) for field in (n,v)])
    physical=actual.subs(v,zeta-delta*n).doit()
    physical=sp.ImmutableMatrix([sp.expand(physical[0]-delta*physical[1]),sp.expand(physical[1])])
    return {"literal_homogeneous_density":L,"zero_prepared_charge":charge,
            "matter_reconstruction_rate":-3*ell*v-w*n,
            "effective_two_metric_density":effective,"metric_Euler_currents":actual,
            "physical_metric_Euler_currents":physical,
            "physical_effective_density":effective.subs(v,zeta-delta*n).doit(),
            "normalized_matter_Euler":euler(L,matter),
            "highest_metric_derivative_order":2}


@cache
def checks():
    item=data()
    old=system.old
    mapping={old.v:v,old.n:n,old.matter:matter,old.vd:sp.diff(v,u),old.sd:sp.diff(matter,u),
             old.theta:theta,old.ell:ell,old.w:w,old.J:Je,system.margin.delta_J:0,
             system.Fn:0,system.Fv:0}
    # The literal action contains the symbolic delta_J, distinct from its clock profile.
    imported=system.action()["L"].subs({old.J:Je-system.margin.delta_J},simultaneous=True)
    imported=imported.subs(mapping,simultaneous=True)
    out={"literal_density_replays_S6_57":sp.simplify(imported-item["literal_homogeneous_density"]),
         "matter_charge_evolution":sp.simplify(
             item["normalized_matter_Euler"]+time(item["zero_prepared_charge"])+3*H*item["zero_prepared_charge"]),
         "effective_two_current_variation":sp.ImmutableMatrix([
             sp.simplify(item["metric_Euler_currents"][j]-euler(item["effective_two_metric_density"],field))
             for j,field in enumerate((n,v))])}
    expected=sp.Matrix([
        2*(Je-3*theta**2)*n+6*theta*sp.diff(v,u)-3*ell*w*v,
        6*sp.diff(v,u,2)+18*H*sp.diff(v,u)-6*theta*sp.diff(n,u)
        -(6*sp.diff(theta,u)+18*H*theta+3*ell*w)*n-9*ell**2*v])
    out["explicit_local_two_metric_equations"]=sp.ImmutableMatrix([
        sp.simplify(value) for value in item["metric_Euler_currents"]-expected])
    out["zero_charge_matter_sign"]=sp.simplify(
        item["zero_prepared_charge"].subs(sp.diff(matter,u),item["matter_reconstruction_rate"]))
    physical=item["physical_metric_Euler_currents"]
    top=sp.ImmutableMatrix([[sp.diff(row,sp.diff(field,u,2)) for field in (n,zeta)] for row in physical])
    out["physical_local_principal_matrix"]=sp.ImmutableMatrix((top-6*sp.Matrix([[delta**2,-delta],[-delta,1]])).applyfunc(sp.simplify))
    out["physical_point_chart_preserves_both_Euler_currents"]=sp.ImmutableMatrix([
        sp.simplify(physical[j]-euler(item["physical_effective_density"],field))
        for j,field in enumerate((n,zeta))])
    out["physical_point_chart_has_no_bounce_denominator"]=sp.Matrix([[1,0],[delta,1]]).det()-1
    return out
