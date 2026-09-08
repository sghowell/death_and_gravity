"""Original-clock scalar principal potential and continuous momentum correction."""
from functools import cache

import sympy as sp
from p8_vector_state import wkb

from . import reduction

u = wkb.u
r = 1+u**2
a = r**2
mass2 = sp.Symbol("positive_mass_squared", positive=True)


def derivative(value):
    return sp.expand(a*sp.diff(value, u))


@cache
def data():
    U = mass2*a**2
    U1, U2 = derivative(U), derivative(derivative(U))
    b = derivative(a)/a
    curvature = derivative(derivative(a))/a
    numerator = b*U1+U2/2
    return {"U": U, "U1": sp.factor(U1), "U2": sp.factor(U2), "b": sp.factor(b),
            "pump_curvature": sp.factor(curvature), "remainder_first_numerator": sp.factor(numerator),
            "scalar_principal_potential": sp.factor(U-curvature)}


def exact_rational(value, allow_zero=False):
    if type(allow_zero) is not bool:
        raise ValueError("Require a native boolean zero-domain flag")
    if value is None or isinstance(value, (bool, str, float)) or value in (sp.true, sp.false):
        raise ValueError("Require an exact rational scale")
    value = sp.sympify(value)
    if value.is_Rational is not True or value.has(sp.Float) or (value < 0 if allow_zero else value <= 0):
        raise ValueError("Require a finite exact positive mass and nonnegative momentum")
    return value


def bounds(momentum, mass):
    momentum, mass = exact_rational(momentum, True), exact_rational(mass)
    fraction = mass**2/(momentum**2+mass**2)
    return {"momentum": momentum, "mass": mass,
            "clock_nonnegative_remainder_lower": sp.Integer(0),
            "clock_remainder_upper": sp.Rational(359375, 4096)*fraction,
            "coarse_clock_remainder_upper": 88*fraction,
            "clock_scalar_principal_potential_lower": mass**2-sp.Rational(275, 16)}


@cache
def checks():
    z = sp.Symbol("momentum_fraction", nonnegative=True)
    H = 4*u/r
    H1 = sp.diff(H, u)
    item = data()
    substitute = {reduction.U: item["U"], reduction.U1: item["U1"],
                  reduction.U2: item["U2"], reduction.b: item["b"], reduction.Acurv: item["pump_curvature"]}
    actual = reduction.longitudinal_potential()["exact"].subs(substitute).subs(
        reduction.k**2, mass2*a**2*z/(1-z))
    target = a**2*(mass2-z*H1+z*(1-3*z)*H**2)
    out = {"clock_exact_acoustic_potential": sp.factor(actual-target),
           "scalar_principal_positive_curvature_formula": sp.factor(item["pump_curvature"]-4*r**2*(1+7*u**2)),
           "remainder_first_numerator_formula": sp.factor(item["remainder_first_numerator"]-4*mass2*r**6*(1+19*u**2)),
           "positive_remainder_factorization": sp.factor(
               (target-item["scalar_principal_potential"])
               -a**2*(1-z)*(H1+(2+3*z)*H**2)),
           "isotropic_low_momentum_potential": sp.factor(target.subs(z, 0)-item["U"]),
           "high_momentum_principal_potential": sp.factor(target.subs(z, 1)-item["scalar_principal_potential"])}
    # Independent cosmic-time Liouville conversion using the frozen S6.50 rate.
    old = wkb.frequency("longitudinal")["U"].subs(wkb.z, z)
    out["clock_cosmic_time_rate_conversion"] = sp.factor(
        a**2*(mass2-old+H1/2+H**2/4)-target)
    out["clock_transverse_exact_scalar_potential"] = sp.factor(
        a**2*(mass2-wkb.frequency("transverse")["U"]+H1/2+H**2/4)-item["U"])
    return out
