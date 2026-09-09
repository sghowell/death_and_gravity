"""All actual lapse derivatives, retaining the original I and Q definitions."""
from functools import cache

import sympy as sp
from p8_auxiliary_neighborhood import center, model

from . import taylor


@cache
def boundary_jets():
    ctx,e=taylor.context(model.u)
    N=1+e
    c=model.coefficients()
    primitive=-model.N**-2*c["boundary_primitive_s"].subs(model.s,1/model.N)
    # The complete derivative simplifies to QQ(u), although isolated
    # prefactors may contain powers of (1+u^2)^(1/4). Differentiate
    # before evaluating those cancelling algebraic factors.
    primitive_derivatives=tuple(sp.factor(sp.diff(primitive,model.N,n).subs(model.N,1))
                                for n in range(5))
    I_N=taylor.from_derivatives(ctx,e,primitive_derivatives)
    I=taylor.integrate_zero(I_N,e)
    Iphi=taylor.from_derivatives(ctx,e,tuple(sp.factor(sp.diff(value,model.u))
                                           for value in taylor.derivatives(I)))
    x=model.lower.d.x
    mapping={x:-N.power(-2)}
    forcing=taylor.evaluate(c["Q_ODE_forcing"],mapping,ctx)
    coefficient=taylor.evaluate(c["Q_ODE_coefficient"],mapping,ctx)
    dx=2*N.power(-3)
    Q=ctx.jet()
    for degree in range(ctx.n):
        rhs=dx*(forcing-coefficient*Q)
        derivative=rhs.coefficient((1<<degree)-1)
        Q+=ctx.jet(derivative)*e**(degree+1)/sp.factorial(degree+1)
    return {"context":ctx,"epsilon":e,"N":N,"I":I,"Iphi":Iphi,"Q":Q,
            "I_N":I_N,"Q_N_rhs":dx*(forcing-coefficient*Q)}


@cache
def jets():
    data=boundary_jets()
    g=model.generic()
    polynomial=sp.Poly(sp.expand(g["Hamiltonian_over_hat_volume"]),*model.VARIABLES)
    actual=model.specialization()
    boundary={symbol:sum(value*(model.N-1)**n/sp.factorial(n)
                         for n,value in enumerate(taylor.derivatives(data[name])))
              for symbol,name in ((model.I,"I"),(model.Iphi,"Iphi"),(model.Q,"Q"))}
    rows={}
    for powers,coefficient in polynomial.terms():
        expression=coefficient.subs(actual,simultaneous=True).subs(boundary,simultaneous=True)
        rows[powers]=tuple(sp.factor(sp.diff(expression,model.N,n).subs(model.N,1)) for n in range(5))
    return rows


def scalar_coefficients():
    return {powers:row for powers,row in jets().items() if sum(powers)==0}


@cache
def checks():
    data=boundary_jets()
    I,Q=(taylor.derivatives(data[name]) for name in ("I","Q"))
    old=model.adm.clock_lapse()
    out={"original_I_fixed_basepoint":I[0],
         "original_Q_fixed_basepoint":Q[0],
         "original_I_first_lapse_derivative":sp.factor(I[1]-old["I_N"]),
         "original_I_second_lapse_derivative":sp.factor(I[2]-old["I_NN"])}
    bg=model.coefficients()["background"]
    out["original_Q_second_lapse_derivative"]=sp.factor(Q[2]+3*sp.diff(bg["h"],model.u)/bg["h"]**3)
    for degree in range(4):
        out[f"original_Q_ODE_lapse_jet_{degree}"]=sp.factor(
            Q[degree+1]-data["Q_N_rhs"].coefficient((1<<degree)-1))
        out[f"original_I_primitive_lapse_jet_{degree}"]=sp.factor(
            I[degree+1]-data["I_N"].coefficient((1<<degree)-1))
    for degree in range(5):
        out[f"original_boundary_time_derivative_{degree}"]=sp.factor(
            taylor.derivatives(data["Iphi"])[degree]-sp.diff(I[degree],model.u))
    rows=jets()
    zero=(0,)*len(model.VARIABLES)
    out["actual_stationary_background_force"]=rows[zero][1]
    out["actual_background_margin_lapse_Hessian"]=sp.factor(
        rows[zero][2]+2*bg["J"]+8*model.MARGIN/bg["h"]**2)
    point=sp.Poly(sp.expand(center.data()["full_closed_bounce_Hamiltonian"]),*model.VARIABLES)
    reference=dict(point.terms())
    for powers in set(rows)|set(reference):
        expected=reference.get(powers,sp.Integer(0))
        row=rows.get(powers,(sp.Integer(0),)*5)
        for degree in range(5):
            out["actual_closed_bounce_coefficient_"+"_".join(map(str,powers))+f"_N{degree}"]=sp.factor(
                row[degree].subs(model.u,0)-sp.diff(expected,model.N,degree).subs(model.N,1))
    return out
