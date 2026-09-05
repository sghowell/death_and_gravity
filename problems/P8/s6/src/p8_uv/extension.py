"""Exact D-only smooth splicing algebra; smoothness proof is in the notes."""

from functools import cache

import sympy as sp
from p8.matter import ia_completion
from p8_s5 import nonlinear_d as d_model

X = sp.Symbol("X", real=True)
u = d_model.u
Z, mass = sp.symbols("Z m", positive=True)
u0, coupling = sp.symbols("u0 lambda", real=True)
switch = sp.Function("B")(X)
LEFT, RIGHT = sp.Rational(1, 4), sp.Rational(3, 4)


def flat_step(argument):
    """C-infinity, not real analytic at zero. Do not differentiate a sample."""
    return sp.Piecewise((0, argument <= 0), (sp.exp(-1/argument), True))


def clock_switch(argument):
    argument = sp.sympify(argument)
    lower, upper = flat_step(argument-LEFT), flat_step(RIGHT-argument)
    return lower/(lower+upper)


def vacuum_function():
    return Z*X/2-Z*mass**2*(u-u0)**2/2+coupling*Z**2*X**2


@cache
def family():
    old = d_model.functions()
    old_f = old["F"].subs(d_model.X, X)
    a3 = switch*old["a3"]
    a2, a4, a5 = ia_completion(-sp.Rational(1, 2), sp.S.Zero, sp.S.Zero, a3, X)
    return {"F": switch*old_f+(1-switch)*vacuum_function(),
            "F2": -sp.Rational(1, 2), "K": sp.Integer(0), "A1": sp.Integer(0),
            "A2": a2, "A3": a3, "A4": a4, "A5": a5}


def branch(value):
    if value not in (0, 1):
        raise ValueError("Only the proven constant switch branches are allowed")
    return {key: expr.subs(switch, value) for key, expr in family().items()}


def flat_polynomials(count):
    """r^(n)(y)=exp(-1/y)*P_n(1/y) on y>0."""
    if not isinstance(count, int) or count < 0:
        raise ValueError("Derivative count must be a nonnegative integer")
    z = sp.Symbol("z", positive=True)
    polys = [sp.Integer(1)]
    for _ in range(count):
        polys.append(sp.expand(z**2*(polys[-1]-sp.diff(polys[-1], z))))
    return z, polys


@cache
def checks():
    f, tube, vac = family(), branch(1), branch(0)
    old = d_model.functions()
    targets = {"F": old["F"], "A3": old["a3"], "A4": old["a4"], "A5": old["a5"]}
    residuals = {f"tube_{key}": sp.cancel(tube[key]-value.subs(d_model.X, X))
                 for key, value in targets.items()}
    residuals.update({"vacuum_F": sp.cancel(vac["F"]-vacuum_function()),
                      "Ia_A4": sp.cancel(f["A4"]+f["A3"]-X**2*f["A3"]**2/4),
                      "Ia_A5": sp.cancel(f["A5"]+X*f["A3"]**2)})
    residuals.update({f"vacuum_{key}": vac[key] for key in ("K", "A1", "A2", "A3", "A4", "A5")})
    point = {u: u0, X: 0}
    fv = vac["F"]
    residuals.update({"vacuum_metric_equation": fv.subs(point),
                      "vacuum_scalar_equation": sp.diff(fv, u).subs(point),
                      "vacuum_kinetic": sp.diff(fv, X).subs(point)-Z/2,
                      "vacuum_mass": sp.diff(fv, u, 2).subs(point)+Z*mass**2})
    psi, Y = sp.symbols("psi Y", real=True)
    canonical = fv.subs({u: u0+psi/sp.sqrt(Z), X: Y/Z})
    residuals["canonical_vacuum"] = sp.expand(canonical-(Y/2-mass**2*psi**2/2+coupling*Y**2))
    z, polys = flat_polynomials(8)
    y = sp.Symbol("y", positive=True)
    for n, polynomial in enumerate(polys[:-1]):
        current = sp.exp(-1/y)*polynomial.subs(z, 1/y)
        residuals[f"flat_derivative_{n+1}"] = sp.simplify(
            sp.diff(current, y)-sp.exp(-1/y)*polys[n+1].subs(z, 1/y))
    # Physical Ia scaling is independently obtained from the general formula.
    M, a3 = sp.symbols("M a3", positive=True)
    a2, a4, a5 = ia_completion(-M**2/2, 0, 0, a3, X)
    residuals.update({"physical_A2": a2,
                      "physical_A4": sp.cancel(a4+a3-X**2*a3**2/(4*M**2)),
                      "physical_A5": sp.cancel(a5+X*a3**2/M**2)})
    return {key: sp.cancel(value) for key, value in residuals.items()}


def derive():
    residuals = checks()
    if any(value != 0 for value in residuals.values()):
        raise ValueError("Nonzero vacuum extension residual")
    z, polys = flat_polynomials(8)
    return {"units": "M=tau=1 for algebra; physical restoration in notes/vacuum-extension.md",
            "switch_thresholds": [str(LEFT), str(RIGHT)],
            "clock_tube": ["9/10", "11/10"],
            "tube_margin": str(sp.Rational(9, 10)-RIGHT),
            "switch_formula": "r(X-1/4)/(r(X-1/4)+r(3/4-X)); r(y)=0 if y<=0 else exp(-1/y)",
            "functions": {key: str(value) for key, value in family().items()},
            "vacuum_function": str(vacuum_function()),
            "Ia_denominator_F2_minus_XA1": "-1/2",
            "flat_derivative_polynomials": [str(poly) for poly in polys],
            "flat_polynomial_variable": str(z),
            "exact_residuals": {key: "0" for key in residuals},
            "analytic_branch": "NO_FINITE_CONSTANT_CLOCK_VACUUM_IF_EXACT_TUBE_EQUALITY_AND_CONNECTED_GLOBAL_REAL_ANALYTICITY",
            "scope": "Smooth off-tube algebraic family, not a controlled vacuum-to-bounce or UV matching construction"}
