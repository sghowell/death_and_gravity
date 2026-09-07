"""Regular primary canonical map and secondary scalar constraint.

This chart never divides H. Its instantaneous Hamiltonian momentum block
is not, by itself, a physical high-frequency ghost test: the map mixes
coordinates and momenta and has a time-dependent canonical boundary.
"""

from functools import cache
from math import comb

import sympy as sp
from p8_variable_beta import background

from . import action

Q, R, PR, PE, P0 = sp.symbols("Q R P_R P_E p_0", real=True)
STATE = (R, action.E, action.Z, PR, PE, action.PX)
EXTENDED = (*STATE, P0)
J = sp.zeros(6)
J[:3, 3:] = sp.eye(3)
J[3:, :3] = -sp.eye(3)


@cache
def derive():
    d = action.derive()
    u, c, a, b, h, P, y, kin = (
        d[key] for key in ("u", "c", "a", "b", "h", "P", "y", "kbar"))
    r = 1+2*b*action.K/(3*a**3*P)
    alpha = h/(3*c*a**3*P)
    old = sp.Matrix([Q+r*R+action.E/3+alpha*PR, R+alpha*P0,
                     action.E, action.Z, P0, PR-r*P0, PE-P0/3, action.PX])
    generator = (-sp.diff(r, u)*P0*R-sp.diff(alpha, u)*P0*PR
                 +(r*sp.diff(alpha, u)-alpha*sp.diff(r, u))*P0**2/2)
    boundary = alpha*P0*PR-r*alpha*P0**2/2
    substitution = dict(zip(d["state"], old, strict=True))
    transformed_constraint = sp.factor(d["C_f"].subs(substitution, simultaneous=True))
    expected_D = (-(1+c/y**3)/(6*a**3)+h**2*(1+y/c)**2/(a**3*kin)
                  -2*h**2*y/(3*c**2*a**3*P)+sp.diff(alpha, u))
    return {"u": u, "c": c, "K": action.K, "r": r, "alpha": alpha,
            "old_state": old, "H_time_correction": generator,
            "one_form_boundary": boundary, "transformed_C_f": transformed_constraint,
            "expected_C_f": 3*c*a**3*P*Q, "expected_secondary_D": expected_D,
            "state": STATE, "extended_state": EXTENDED}


@cache
def checks():
    d, r = action.derive(), derive()
    old = r["old_state"]
    new_q = (Q, R, action.E, action.Z)
    new_p = (P0, PR, PE, action.PX)
    new = (*new_q, *new_p)
    boundary, u = r["one_form_boundary"], d["u"]
    residuals = {"regular_primary_constraint": sp.factor(r["transformed_C_f"]-r["expected_C_f"])}
    for variable in new:
        old_coefficient = sum(old[i+4]*sp.diff(old[i], variable) for i in range(4))
        new_coefficient = sum(new_p[i]*sp.diff(new_q[i], variable) for i in range(4))
        residuals[f"one_form_d_{variable}"] = sp.simplify(
            old_coefficient-new_coefficient-sp.diff(boundary, variable))
    residuals["full_time_boundary"] = sp.simplify(
        sum(old[i+4]*sp.diff(old[i], u) for i in range(4))-sp.diff(boundary, u)
        +r["H_time_correction"])
    direction = old.subs(dict.fromkeys((*new_q, PR, PE, action.PX), 0)).subs(P0, 1)
    hpp = 2*d["H"].subs(dict(zip(d["state"], direction, strict=True)), simultaneous=True)
    actual_D = hpp+sp.diff(r["H_time_correction"], P0, 2)
    residuals["secondary_D_from_full_Hessian"] = sp.factor(actual_D-r["expected_secondary_D"])
    residuals["secondary_D_has_no_K"] = sp.diff(r["expected_secondary_D"], action.K)
    residuals["secondary_D_at_center"] = sp.factor(r["expected_secondary_D"].subs(u, 0)+sp.Rational(5, 24))
    return residuals


def jets(time, lapse, order=1):
    """Exact time jets, taken BEFORE evaluation, through order 0, 1 or 2.

Return extended quadratic Hessians M_j, reduced Hessians H_j, and the
eliminated secondary momentum rows L_j with p0=L0 Z. Only H0,H1 are
needed for physical second-order equations; L2 is needed away from H=0.
Generic symbolic lapse is permitted for derivation, not a domain verdict.
"""
    if isinstance(order, bool) or not isinstance(order, (int, sp.Integer)) or order not in (0, 1, 2):
        raise ValueError("Only exact jet orders zero, one or two are supported")
    order = int(order)
    if any(isinstance(value, (bool, float)) for value in (time, lapse)):
        raise TypeError("Only exact real symbolic or rational inputs are supported")
    d = action.derive()
    at = {d["u"]: sp.sympify(time), d["c"]: sp.sympify(lapse)}
    if any(not isinstance(value, sp.Expr) or value.has(sp.Float, sp.oo, -sp.oo, sp.zoo, sp.nan)
           or value.is_real is False for value in at.values()):
        raise ValueError("Only exact finite time/lapse inputs are supported")
    if all(not value.free_symbols for value in at.values()):
        background.evaluate(at[d["u"]], at[d["c"]])
    return _jets_cached(at[d["u"]], at[d["c"]], order)


@cache
def _jets_cached(time, lapse, order):
    """Cache only after validation; bool/float must not alias exact keys."""
    d, r = action.derive(), derive()
    at = {d["u"]: time, d["c"]: lapse}
    u = d["u"]
    T = r["old_state"].subs(Q, 0).jacobian(EXTENDED)
    original = sp.hessian(d["H"], d["state"])
    G = sp.hessian(r["H_time_correction"], EXTENDED)
    tj = tuple(sp.diff(T, u, j).subs(at).applyfunc(sp.cancel) for j in range(order+1))
    hj = tuple(sp.diff(original, u, j).subs(at).applyfunc(sp.cancel) for j in range(order+1))
    gj = tuple(sp.diff(G, u, j).subs(at).applyfunc(sp.cancel) for j in range(order+1))
    extended = []
    for j in range(order+1):
        value = gj[j].copy()
        for left in range(j+1):
            for middle in range(j-left+1):
                right = j-left-middle
                value += comb(j, left)*comb(j-left, middle)*tj[left].T*hj[middle]*tj[right]
        extended.append(value.applyfunc(sp.cancel))
    ds = [m[6, 6] for m in extended]
    if ds[0] == 0:
        raise ValueError("Secondary scalar constraint is singular")
    ls = [-extended[0][6:7, :6]/ds[0]]
    for j in range(1, order+1):
        value = -extended[j][6:7, :6]
        for i in range(1, j+1):
            value -= comb(j, i)*ds[i]*ls[j-i]
        ls.append((value/ds[0]).applyfunc(sp.cancel))
    reduced = []
    for j in range(order+1):
        value = extended[j][:6, :6].copy()
        for i in range(j+1):
            value += comb(j, i)*extended[i][:6, 6:7]*ls[j-i]
        reduced.append(value.applyfunc(sp.cancel))
    return {"M": tuple(extended), "H": tuple(reduced),
            "p0": tuple(row.applyfunc(sp.cancel) for row in ls), "D": tuple(ds)}
