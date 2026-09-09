"""The matter boundary and the original regular scalar momenta."""
from functools import cache

import sympy as sp
from p8_auxiliary_neighborhood import model

P0=sp.Symbol("constant_comoving_matter_momentum",real=True)
v,chi,p,pm,vd,sd=sp.symbols("hat_curvature matter scalar_momentum matter_momentum curvature_rate matter_rate",real=True)


@cache
def data():
    old_p=p-3*P0*chi
    old_pm=pm-3*P0*v
    boundary=3*P0*v*chi
    return {"boundary":boundary,"old_curvature_momentum":old_p,
            "old_matter_momentum":old_pm,"spatial_York_scalar_source":p-3*P0*chi,
            "one_form_difference":p*vd+pm*sd-old_p*vd-old_pm*sd}


@cache
def checks():
    d=data()
    bg=model.coefficients()["background"]
    bracket=lambda left,right:sp.diff(left,v)*sp.diff(right,p)-sp.diff(left,p)*sp.diff(right,v)+sp.diff(left,chi)*sp.diff(right,pm)-sp.diff(left,pm)*sp.diff(right,chi)
    out={"matter_boundary_one_form_identity":sp.expand(d["one_form_difference"]-sp.diff(d["boundary"],v)*vd-sp.diff(d["boundary"],chi)*sd),
         "actual_comoving_matter_momentum_constant":sp.factor(sp.diff(bg["a"]**3*bg["ell"],model.u)),
         "actual_comoving_matter_momentum_value":sp.factor(bg["a"]**3*bg["ell"]-sp.Rational(1,10)),
         "old_scalar_canonical_pair":bracket(v,d["old_curvature_momentum"])-1,
         "old_matter_canonical_pair":bracket(chi,d["old_matter_momentum"])-1,
         "cross_momentum_bracket_zero":sp.expand(bracket(d["old_curvature_momentum"],d["old_matter_momentum"]))}
    return out
