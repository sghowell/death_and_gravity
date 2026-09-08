"""General bundle-rank local heat input and its independent finite action."""
from functools import cache

import sympy as sp
from p8_vector_curvature import heat as frozen

from . import local


@cache
def coefficients():
    D, R, Ric2, Riem2, boxR = local.dimension, frozen.R, frozen.Ric2, frozen.Riem2, frozen.boxR
    values = {0: D, 1: (D/6-1)*R,
              2: (D-12)*R**2/72+(90-D)*Ric2/180+(D-15)*Riem2/180+(D-5)*boxR/30}
    finite = {}
    for n, coefficient in values.items():
        r = sp.Rational((-1)**(2-n), sp.factorial(2-n))
        finite[n] = sp.factor(-2*r*coefficient.subs(D, 3)*(sp.harmonic(2-n)-local.ell)
                              +4*r*sp.diff(coefficient, D).subs(D, 3))
    return {"dimension_dependent_coefficients": values, "finite_Euclidean_action_over_64pi_squared": finite}


@cache
def checks():
    data = coefficients()
    prior = frozen.coefficients()
    out = {"rank_retained_heat_coefficient_reduces_to_frozen_"+str(n):
           sp.factor(data["dimension_dependent_coefficients"][n].subs(local.dimension, 3)-prior["proca"][2*n])
           for n in range(3)}
    out["bundle_rank_derivative_is_minimal_scalar_heat_coefficient"] = sp.factor(
        sp.diff(data["dimension_dependent_coefficients"][2], local.dimension)-prior["scalar"][4])
    expected = {0: 3*local.ell-sp.Rational(5, 2),
                1: (local.ell-sp.Rational(5, 3))*frozen.R,
                2: 2*local.ell*prior["proca"][4]+4*prior["scalar"][4]}
    out.update({"independent_finite_heat_action_Gamma_expansion_"+str(n):
                sp.factor(data["finite_Euclidean_action_over_64pi_squared"][n]-expected[n]) for n in range(3)})
    return out
