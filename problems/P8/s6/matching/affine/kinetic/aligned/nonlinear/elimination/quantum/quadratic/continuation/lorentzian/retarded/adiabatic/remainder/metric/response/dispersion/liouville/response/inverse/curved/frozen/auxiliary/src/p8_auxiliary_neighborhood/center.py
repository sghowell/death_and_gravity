"""Closed bounce Hamiltonian retaining the actual time-boundary primitive."""
from functools import cache

import sympy as sp
from p8_s5 import lapse_series

from . import model

N=model.N
r=sp.Symbol("inverse_lapse",positive=True)


@cache
def data():
    eps=model.MARGIN
    Iphi=-6*r**sp.Rational(3,2)-18*r**(-sp.Rational(1,2))+24
    gt=9*(1-4*N**3)/(11+6*N**2-44*N**3)
    gs=9*N*(4+5*N)/(-18+44*N+55*N**2)
    delta=N**-2-1
    background=(-(sp.Rational(23,8)+eps)*N**(-sp.Rational(3,2))
                -(sp.Rational(2849,200)-2*eps)*sp.sqrt(N)
                +(sp.Rational(9,8)-eps)*N**sp.Rational(5,2)+24
                +sp.Rational(1,200)/sp.sqrt(N))
    linear=(2*N**sp.Rational(3,2)*model.shear+sp.sqrt(N)*model.electric/2
            +sp.sqrt(N)*model.magnetic/4+N**sp.Rational(3,2)*model.vector/(2*gs)
            +sp.Rational(1,10)*model.dc/sp.sqrt(N)
            +N**sp.Rational(3,2)*model.chi_gradient/2-model.curvature/(2*sp.sqrt(N)))
    quadratic=(-3*N**sp.Rational(3,2)*(model.dp-delta*model.j)**2/4
               +gt*model.j**2/(2*sp.sqrt(N))+model.dc**2/(2*sp.sqrt(N)))
    return {"actual_boundary_primitive_phi_at_bounce":Iphi,
            "actual_gamma_t_at_bounce":gt,"actual_gamma_s_at_bounce":gs,
            "background":background,"linear_invariant_part":linear,"quadratic_invariant_part":quadratic,
            "full_closed_bounce_Hamiltonian":background+linear+quadratic,
            "lapse_background_derivatives":tuple(sp.factor(sp.diff(background,N,j).subs(N,1)) for j in range(5))}


@cache
def series():
    d=data()
    A=d["lapse_background_derivatives"]
    L=tuple(sp.diff(d["linear_invariant_part"],N,j).subs(N,1) for j in range(4))
    Q=tuple(sp.diff(d["quadratic_invariant_part"],N,j).subs(N,1) for j in range(3))
    return lapse_series.stationary_series(A,L,Q)


@cache
def checks():
    d=data()
    primitive=model.coefficients()["boundary_primitive_s"]
    expected=sp.diff(d["actual_boundary_primitive_phi_at_bounce"],r)
    out={"actual_bounce_boundary_primitive_phi_derivative":sp.factor(
        sp.diff(primitive,model.u).subs({model.u:0,model.s:r})-expected),
         "actual_bounce_boundary_primitive_basepoint":d["actual_boundary_primitive_phi_at_bounce"].subs(r,1)}
    g=model.generic()
    mapping=model.specialization()
    point={model.u:0,model.I:0,model.Q:0,
           model.Iphi:d["actual_boundary_primitive_phi_at_bounce"].subs(r,1/N)}
    # Expand only in invariant variables, never a giant rational polynomial
    # with every time coefficient and canonical jet simultaneously.
    delta=g["Hamiltonian_over_hat_volume"].subs(mapping,simultaneous=True).subs(point)-d["full_closed_bounce_Hamiltonian"]
    out["full_actual_Hamiltonian_matches_closed_bounce_expression"]=sum(
        sp.factor(value)*sp.prod(x**p for x,p in zip(model.VARIABLES,powers))
        for powers,value in sp.Poly(sp.expand(delta),*model.VARIABLES).terms())
    out["bounce_zero_lapse_force"]=d["lapse_background_derivatives"][1]
    out["bounce_actual_margin_lapse_Hessian"]=sp.factor(
        d["lapse_background_derivatives"][2]+sp.Rational(1199,400)+8*model.MARGIN)
    mass=model.gamma()
    out["bounce_temporal_mass_rational_function"]=sp.factor(
        mass["temporal"].subs(mass["p"],1/(2*N))-d["actual_gamma_t_at_bounce"])
    out["bounce_spatial_mass_rational_function"]=sp.factor(
        mass["spatial"].subs(mass["p"],1/(2*N))-d["actual_gamma_s_at_bounce"])
    # The already-derived generic stationary formulas apply to A+L+Q
    # in invariant degree, including the second-order lapse correction.
    out.update({"generic_stationary_"+name:value for name,value in lapse_series.generic_checks().items()})
    return out
