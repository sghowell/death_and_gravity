"""Actual clock-chart pullback plus the fixed local-profile part of S6.82."""
from functools import cache

import sympy as sp
from p8_constant_proca import quantum

from . import local

u,n=local.u,local.n
vhat=sp.symbols("hat_logscale_source0:10",real=True)
h=(1+u*u)**3


@cache
def mapping():
    w1=1/(2*h)
    w2=-sp.Rational(3,2)/h+1/h**2
    row=[vhat[0]+w1*n[0]]
    for _ in range(9):
        row.append(local.time(row[-1],extra=(vhat,)))
    return {"omega_N":w1,"omega_NN":w2,"linear_physical_logscale_jet":tuple(row)}


def density(order):
    order=local.jets.order(order)
    return _density(order)


@cache
def _density(order):
    d=mapping()
    coefficients=quantum.local_coefficients()[order]
    rho,p=coefficients["energy"],coefficients["pressure"]
    linear_map=dict(zip(local.v,d["linear_physical_logscale_jet"]))
    metric=sp.expand(local.density(order).subs(linear_map,simultaneous=True))
    metric+=3*p*d["omega_NN"]*n[0]**2/2
    # This local part is fixed by linear decomposition of the already
    # chosen S6.82 profile, not by reselecting its state or counterterms.
    N,V=sp.symbols("literal_physical_lapse literal_hat_logscale",real=True)
    omega=-sp.log((h-1+N**-2)/h)/4
    P=-p+(rho+p)*(1-N**-2)/2
    expression=N*sp.exp(3*omega+3*V)*P
    point={N:1,V:0}
    tad=sp.diff(expression,N,2).subs(point)*n[0]**2/2
    tad+=sp.diff(expression,N,V).subs(point)*n[0]*vhat[0]
    tad+=sp.diff(expression,V,2).subs(point)*vhat[0]**2/2
    metric=local.clean(metric,extra=(vhat,))
    tad=local.clean(tad,extra=(vhat,))
    return {"metric_linear_pullback_and_second_map_contact":metric,
            "actual_fixed_local_profile_quadratic_density":tad,
            "combined_local_quadratic_density":local.clean(metric+tad,extra=(vhat,))}


def operator(output,order):
    if type(output) is not str or output not in ("N","V"):
        raise ValueError("Require clock-chart lapse N or hat log-scale V current")
    order=local.jets.order(order)
    return local.euler(density(order)["combined_local_quadratic_density"],
                       n if output=="N" else vhat,extra=(vhat,))


def direct_metric_density(order):
    return _direct_metric_density(local.jets.order(order))


@cache
def _direct_metric_density(order):
    e,D=local.counterterms.e,local.jets.D
    geometry=local.counterterms.geometry()
    invariant={0:sp.Rational(5,2),1:5*geometry["R"]/3,2:-4*geometry["a4sc"]}[order]
    volume=1+e*(n[0]+D*local.v[0])+e**2*(D*n[0]*local.v[0]+D**2*local.v[0]**2/2)
    second=[mapping()["omega_NN"]*n[0]**2/2]
    # The finite covariant invariants need log-scale jets only through
    # order two. No higher second-map derivative enters this action.
    for _ in range(2):
        second.append(local.time(second[-1],extra=(vhat,)))
    # Compose the exact quadratic Taylor polynomial before introducing
    # rational clock coefficients. The first covariant variation acts
    # on the second metric-map jet; the second variation acts on two
    # first map jets. Higher parameter powers cannot contribute.
    first=local.actual(local.counterterms.keep(volume*invariant).coeff(e,1).subs(D,3))
    linear_map=dict(zip(local.v,mapping()["linear_physical_logscale_jet"]))
    value=local.density(order).subs(linear_map,simultaneous=True)
    value+=sum(sp.diff(first,field)*second[j] for j,field in enumerate(local.v[:3]))
    return local.clean(value,extra=(vhat,))


@cache
def top():
    matrix=sp.Matrix([[sp.factor(operator(output,2).coeff(row[4])) for row in (n,vhat)] for output in ("N","V")])
    e=sp.Matrix([mapping()["omega_N"],1])
    return {"actual_fourth_derivative_coefficient":sp.ImmutableMatrix(matrix),
            "expected_rank_one_coefficient":sp.ImmutableMatrix(-4*e*e.T),
            "kernel_vector":sp.ImmutableMatrix([1,-mapping()["omega_N"]])}


@cache
def checks():
    N=sp.Symbol("literal_positive_lapse",positive=True)
    omega=-sp.log((h-1+N**-2)/h)/4
    out={"actual_clock_chart_first_logscale_derivative":sp.factor(sp.diff(omega,N).subs(N,1)-mapping()["omega_N"]),
         "actual_clock_chart_second_logscale_derivative":sp.factor(sp.diff(omega,N,2).subs(N,1)-mapping()["omega_NN"]),
         "constant_local_potential_cancels_its_fixed_profile_exactly":sp.expand(density(0)["combined_local_quadratic_density"])}
    data=top()
    out["actual_finite_local_fourth_coefficient_is_rank_one"]=data["actual_fourth_derivative_coefficient"]-data["expected_rank_one_coefficient"]
    out["no_invertible_two_source_fourth_local_block"]=sp.factor(data["actual_fourth_derivative_coefficient"].det())
    out["actual_rank_one_fourth_order_null_direction"]=sp.simplify(
        data["actual_fourth_derivative_coefficient"]*data["kernel_vector"])
    for order in range(3):
        difference=local.clean(direct_metric_density(order)
            -density(order)["metric_linear_pullback_and_second_map_contact"],extra=(vhat,))
        for output,row in (("N",n),("V",vhat)):
            out[output+f"_independent_nonlinear_metric_map_Euler_{order}"]=local.euler(difference,row,extra=(vhat,))
    return out
