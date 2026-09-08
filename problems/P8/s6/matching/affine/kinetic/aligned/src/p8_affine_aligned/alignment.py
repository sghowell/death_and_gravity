"""Exact source shift in a separately named kinetic operator.

S6.41's mass-centering term is unchanged; only F(T) is replaced by
F(T-B*dphi). The shifted source is quadratic, not zero nonlinearly.
"""
from functools import cache

import sympy as sp
from p8_affine import connection as affine
from p8_affine_retuned import degeneracy as primary
from p8_affine_retuned import dynamics as parent
from p8_affine_retuned import geometry
from p8_affine_ricci import source as original

u = parent.u
X = sp.Symbol("x_source", negative=True)
LOWER = sp.Symbol("Q_lower", real=True)


@cache
def coefficient():
    return sp.factor(-parent.source()["d"]*(1+X)/2)


@cache
def normal_source():
    p, s = affine.P, affine.S
    bg = parent.old.background()
    h, hubble = bg["h"], bg["H"]
    x = -s**2
    fphi = sp.diff(h, u)*(1+x)/(2*h**2)
    cubic = -(LOWER+fphi)/(4*p*x)
    actual = sp.factor(primary.chart()["expected"].subs(affine.F3, cubic)
                       -s*coefficient().subs(X, x))
    expected = (4*p**2-1)*(primary.KHAT-3*hubble*s)+sp.Rational(3, 2)*s*LOWER
    residual = sp.factor(sp.expand(actual-expected).subs(p**2, (h-1+s**2)/(4*h)))
    return {"actual": actual, "expected": expected, "residual": residual,
            "delta": (s**2-1)/h, "physical_K_deviation": primary.KHAT-3*hubble*s,
            "B": coefficient().subs(X, x), "J3": cubic,
            "nonlinear_source_is_not_identically_zero": True}


@cache
def rolling():
    bg = parent.old.background()
    s, n = affine.S, original.n
    actual = parent.source()["actual"].subs({original.H: bg["H"], original.h: bg["h"]})
    shifted = s*coefficient().subs(X, -s**2)
    first_B = -n*sp.diff(shifted, s).subs(s, 1)
    source = geometry.clean(actual-sp.Matrix([first_B, 0, 0, 0]))
    eps, lapse, trace = sp.symbols("eps n delta_K_hat", real=True)
    Q2 = sp.Symbol("Q_second_order", real=True)
    complete = normal_source()["expected"].subs(
        {affine.P**2: (bg["h"]-1+s**2)/(4*bg["h"])})
    expanded = sp.expand(complete.subs({s: 1-eps*lapse,
                                      primary.KHAT: 3*bg["H"]+eps*trace,
                                      LOWER: eps**2*Q2}, simultaneous=True))
    return {"first_variation": source, "first_B_normal": sp.factor(first_B),
            "background_B": coefficient().subs(X, -1), "background_trace": parent.source()["background"],
            "exact_expansion": expanded,
            "zero_order": sp.factor(expanded.coeff(eps, 0)),
            "first_order": sp.factor(expanded.coeff(eps, 1)),
            "second_order": sp.factor(expanded.coeff(eps, 2))}


@cache
def quadratic_mass():
    """Expand a generic coefficient variation; retain the first nonlinear source."""
    eps, volume = sp.symbols("eps volume_first", real=True)
    field = sp.Matrix(sp.symbols("W0:4", real=True))
    source = sp.Matrix(sp.symbols("S2_0:4", real=True))
    entries = sp.symbols("D0:10", real=True)
    variation = sp.zeros(4)
    k = 0
    for i in range(4):
        for j in range(i, 4):
            variation[i, j] = variation[j, i] = entries[k]
            k += 1
    difference = eps*field-eps**2*source
    action = sp.expand((1+eps*volume)*(difference.T*(geometry.ETA+eps*variation)*difference)[0]/2)
    expected = (field.T*geometry.ETA*field)[0]/2
    return {"second_order": action.coeff(eps, 2), "expected": expected,
            "third_order_source_coupling": -(field.T*geometry.ETA*source)[0],
            "residual": sp.expand(action.coeff(eps, 2)-expected)}


@cache
def checks():
    source, roll, quadratic = normal_source(), rolling(), quadratic_mass()
    return {"exact_nonlinear_shifted_source": source["residual"],
            "full_all_component_first_variation": roll["first_variation"],
            "background_shift_zero": roll["background_B"],
            "background_original_trace_zero": roll["background_trace"],
            "source_expansion_zero_order": roll["zero_order"],
            "source_expansion_first_order": roll["first_order"],
            "generic_quadratic_mass_has_no_light_source": quadratic["residual"]}
