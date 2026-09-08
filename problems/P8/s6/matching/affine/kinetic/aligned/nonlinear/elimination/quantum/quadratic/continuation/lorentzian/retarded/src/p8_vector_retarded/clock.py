"""The unchanged clock profiles, canonical preparation, and volume factors."""
from functools import cache

import sympy as sp
from p8_affine_kinetic import scalar
from p8_vector_quadratic import response as prior

from . import canonical, flat

u = scalar.u


@cache
def profiles():
    background = scalar.background()
    alpha, beta = prior.mass_coefficients()
    return {"a": background["a"], "H": background["H"],
            "alpha": alpha, "beta": beta, "initial_time": -sp.Rational(1, 2)}


@cache
def checks():
    out = dict(prior.actual_mass_checks())
    data = profiles()
    out["physical_clock_scale"] = sp.factor(data["a"]-(1+u**2)**2)
    out["physical_clock_Hubble"] = sp.factor(data["H"]-sp.diff(data["a"], u)/data["a"])
    a0, a1 = sp.symbols("scale_at_output scale_at_source", positive=True)
    out["two_density_volumes_cancel_four_field_normalizations"] = (
        a0**3*a1**3*(a0**sp.Rational(-3, 2)*a1**sp.Rational(-3, 2))**2-1)
    frequency, half_rate, d = sp.symbols("Borel_frequency Borel_half_rate canonical_rate", real=True)
    frequency = sp.Symbol("Borel_frequency", positive=True)
    v0 = 1/sp.sqrt(2*frequency)
    p0 = (-half_rate-sp.I*frequency-d)*v0
    Sigma = sp.Matrix([[v0**2, sp.re(v0*sp.conjugate(p0))],
                       [sp.re(p0*v0), sp.expand_complex(p0*sp.conjugate(p0))]])
    out["selected_initial_Wronskian"] = sp.simplify(v0*sp.conjugate(p0)-p0*v0-sp.I)
    out["selected_initial_covariance_determinant"] = sp.factor(Sigma.det()-sp.Rational(1, 4))
    for sector in ("T", "L"):
        expression = canonical.insertion(sector)["J"]
        mapped = expression.subs({canonical.alpha: data["alpha"], canonical.beta: data["beta"]})
        target = (-canonical.m**2*data["beta"]*canonical.v**2/2 if sector == "T" else
                  (data["alpha"]*canonical.q*canonical.p**2/(canonical.m**2+canonical.q)
                   -data["beta"]*(canonical.m**2+canonical.q)*canonical.v**2)/2)
        out[sector+"_actual_mass_only_clock_insertion"] = sp.factor(mapped-target)
    controls = flat.controls()
    replaced = {flat.alpha: data["alpha"], flat.beta: data["beta"]}
    targets = (sp.Rational(113, 8748)*flat.m**4/(1+u**2)**6,
               -sp.Rational(52, 6561)*flat.m**2/(1+u**2)**6,
               sp.Rational(8, 6561)/(1+u**2)**6)
    for j, target in enumerate(targets):
        out["actual_flat_time_pole_"+str(2*j)] = sp.factor(flat.integrated_residue(j).subs(replaced)-target)
    out["actual_nonzero_omitted_contact_UV_residue"] = sp.factor(
        controls["omitted_contact_UV_residue_times_pi_squared"].subs(replaced)
        -flat.m**4/(108*(1+u**2)**6))
    return out
