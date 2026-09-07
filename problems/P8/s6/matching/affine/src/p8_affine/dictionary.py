"""Principal metric-affine inverse and lower-order matching algebra.

The source uses (-+++), x=-X_repo.  The sign conversion is a convention,
not a physical metric or matter redefinition.  The source's printed Q2
formula is kept as a negative control; its proposed corrected form must
separately pass the literal connection-action audit before certification.
"""
import json
from fractions import Fraction
from functools import cache
from pathlib import Path

import sympy as sp
from p8.matter import ia_completion

u = sp.Symbol("u", real=True)
x = sp.Symbol("x", negative=True)
P8 = Path(__file__).resolve().parents[5]
WITNESS = P8/"certificates"/"witness-CD_matter.json"


def principal(p, px, c, cx, quartic, variable):
    """Aoki-Shimada 4.5,4.9-4.12, with every X product distinguished.

    Required denominators: p != 0 and Delta != 0. These formulas alone
    do not certify the rank of all unsourced connection components.
    """
    delta = 2*p-c*variable+2*quartic*variable**2
    f = p-c*variable/2
    a1 = -c/2-p*(c-2*quartic*variable)/delta
    a3 = 2*cx+(4*p*quartic+(4*px-c)*(c-2*quartic*variable))/delta
    a4 = (-2*cx+2*px*(3*px-c)/p
          +px*variable*(px*c-4*p*cx)/p**2
          +(c**2-4*p*quartic-2*c*quartic*variable)/delta)
    a5 = (-px*(px*c-4*p*cx)/p**2
          +2*px*(4*p*quartic+(3*px-c)*(c-2*quartic*variable))/(p*delta))
    return {"f": f, "Delta": delta, "alpha1": a1, "alpha2": -a1,
            "alpha3": a3, "alpha4": a4, "alpha5": a5}


@cache
def target():
    """Read the original source-pinned witness, including its scalar F."""
    raw = json.loads(WITNESS.read_text())["covariant_functions"]
    local = {"t": u, "X": -x}
    values = {key: sp.sympify(value, locals=local) for key, value in raw.items()}
    repo_x = sp.Symbol("X", positive=True)
    f2 = values["F2"].subs(x, -repo_x)
    a3 = values["A3"].subs(x, -repo_x)
    a2, a4, a5 = ia_completion(f2, sp.diff(f2, repo_x), sp.Integer(0), a3, repo_x)
    return {"f": -values["F2"], "alpha1": values["A1"],
            "alpha2": a2.subs(repo_x, -x), "alpha3": -values["A3"],
            "alpha4": -a4.subs(repo_x, -x), "alpha5": a5.subs(repo_x, -x),
            "scalar_F": values["F"], "braiding": -values["K"]}


@cache
def lift():
    h = (1+u**2)**3
    w = h-1-x
    f = w/(2*h)
    p = sp.sqrt(w/h)/2
    c = 2*(p-f)/x
    quartic = (sp.Rational(1, 2)-f)/x**2
    return {"h": h, "w": w, "f": f, "p": p, "c": c, "quartic": quartic,
            "px": sp.diff(p, x), "cx": sp.diff(c, x),
            "p_phi": sp.diff(p, u), "f_phi": sp.diff(f, u),
            "Delta": sp.Integer(1)}


@cache
def principal_identities():
    values = lift()
    reduced = principal(*(values[key] for key in ("p", "px", "c", "cx", "quartic")), x)
    expected = target()
    residuals = {name: sp.simplify(reduced[name]-expected[name])
                 for name in ("f", "alpha1", "alpha2", "alpha3", "alpha4", "alpha5")}
    residuals.update({
        "Delta_one": sp.simplify(reduced["Delta"]-1),
        "p_square": sp.simplify(values["p"]**2-values["f"]/2),
        "inverse_ODE": sp.simplify(values["px"]/values["p"]
                                   -sp.diff(values["f"], x)/values["f"]
                                   -x*expected["alpha3"]/(4*values["f"])),
        "clock_normalization": sp.simplify(values["p"].subs(x, -1)-sp.Rational(1, 2)),
    })
    return residuals


def generic_inverse_identities():
    """Derive the inverse for arbitrary nonzero target f and A1=0."""
    z, f, fx, a3, p = sp.symbols("z f fx a3 p", nonzero=True)
    px = p*(fx/f+z*a3/(4*f))
    c = 2*(p-f)/z
    cx = 2*(px-fx)/z-2*(p-f)/z**2
    quartic = (p**2-f**2)/(f*z**2)
    reduced = principal(p, px, c, cx, quartic, z)
    return {"target_f": sp.factor(reduced["f"]-f),
            "target_A1": sp.factor(reduced["alpha1"]),
            "target_A3": sp.factor(reduced["alpha3"]-a3),
            "two_displayed_denominators": sp.factor(reduced["Delta"]-2*p**2/f)}


def rational(value):
    if isinstance(value, (bool, float, str)) or value in (sp.true, sp.false):
        raise TypeError("An exact rational physical input is required")
    if isinstance(value, Fraction):
        return sp.Rational(value.numerator, value.denominator)
    if not isinstance(value, (int, sp.Rational)):
        raise TypeError("An exact rational physical input is required")
    return sp.Rational(value)


def tube_point(time=0, physical_X=1, *, mass_squared=1, time_scale=1):
    """A guarded exact calibration, not a replacement for whole-tube proof."""
    time, physical_X, mass_squared, time_scale = map(
        rational, (time, physical_X, mass_squared, time_scale))
    if not sp.Rational(9, 10) <= physical_X <= sp.Rational(11, 10):
        raise ValueError("The closed original CD clock tube is required")
    if mass_squared <= 0 or time_scale <= 0:
        raise ValueError("Positive physical units are required")
    values = lift()
    substitutions = {u: time, x: -physical_X}
    f = sp.simplify(values["f"].subs(substitutions))
    p = sp.simplify(values["p"].subs(substitutions))
    return {"u": time, "X_repo": physical_X, "x_source": -physical_X,
            "physical_time": time_scale*time, "M_squared": mass_squared,
            "tau": time_scale, "f_physical": mass_squared*f,
            "p_physical": mass_squared*p, "Delta_physical": mass_squared,
            "dimensionless_p_squared": sp.simplify(p**2),
            "additional_quotient_factor": sp.simplify(8*p**2-1),
            "curvature_and_quartic_coefficient_scale": mass_squared,
            "cubic_and_effective_braiding_scale": mass_squared/time_scale,
            "parent_scalar_coefficient_scale": mass_squared/time_scale**2,
            "four_volume_normalized_action_scale": mass_squared*time_scale**2,
            "source_signature_change_not_a_new_matter_metric": True,
            "connection_has_no_added_kinetic_term": True}
