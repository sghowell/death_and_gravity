"""Actual new lapse jets, retaining every unchanged I, I_phi and Q jet."""
from functools import cache

import sympy as sp
from p8_auxiliary_neighborhood import model as old
from p8_coupled_vertices import coefficients as previous

from . import model


def monomial(*indices):
    powers=[0]*len(old.VARIABLES)
    for index in indices:
        powers[index]+=1
    return tuple(powers)


@cache
def jets():
    rows=dict(previous.jets())
    j=old.VARIABLES.index(old.j)
    p=old.VARIABLES.index(old.dp)
    rows[monomial(j)]=(sp.Integer(0),)*5
    rows[monomial(p,j)]=(sp.Integer(0),)*5
    c=old.coefficients()
    coefficient=old.N*c["gamma_t"]/(2*c["U"])
    rows[monomial(j,j)]=tuple(sp.factor(sp.diff(coefficient,old.N,n).subs(old.N,1))
                              for n in range(5))
    return rows


@cache
def checks():
    before=previous.jets()
    after=jets()
    out={}
    for powers in before:
        key="_".join(map(str,powers))
        # All quadratic-in-invariant coefficients at N=1 agree.
        # The invariant map and both canonical boundary generators agree.
        out["all_actual_invariant_coefficients_at_clock_"+key]=sp.factor(after[powers][0]-before[powers][0])
        if sum(powers)==1:
            out["all_linear_invariant_lapse_forces_at_clock_"+key]=sp.factor(after[powers][1]-before[powers][1])
    zero=monomial()
    for n in range(5):
        out[f"full_background_lapse_derivative_unchanged_{n}"]=sp.factor(after[zero][n]-before[zero][n])
    j=old.VARIABLES.index(old.j)
    g=old.generic()
    delta=-old.N*g["d"]**2/(4*g["a"])
    actual=sp.factor(delta.subs(old.specialization(),simultaneous=True))
    for n in range(5):
        # Factor inside positive real polynomial power bases too. This
        # resolves (1+u^2)^(1/4)*((1+u^2)^3)^(1/4) without force=True
        # or changing the fixed analytic branch on the clock tube.
        out[f"actual_new_Gauss_square_lapse_derivative_{n}"]=sp.factor(sp.factor(
            after[monomial(j,j)][n]-before[monomial(j,j)][n]
            -sp.diff(actual,old.N,n).subs(old.N,1),deep=True))
    newpoly=sp.Poly(sp.expand(model.data()["new_Hamiltonian"]),*old.VARIABLES)
    out["new_invariant_polynomial_has_twelve_nonzero_monomials"]=len(newpoly.terms())-12
    out["new_invariant_polynomial_is_at_most_quadratic"]=newpoly.total_degree()-2
    out["actual_changed_Gauss_coefficient_starts_at_second_lapse_degree"]=sp.factor(actual.subs(old.N,1))
    out["actual_changed_Gauss_first_lapse_degree_vanishes"]=sp.factor(sp.diff(actual,old.N).subs(old.N,1))
    return out
