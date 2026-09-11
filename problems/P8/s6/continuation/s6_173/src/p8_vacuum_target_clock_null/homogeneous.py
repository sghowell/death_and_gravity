"""Covariant homogeneous action variation before setting the unit lapse."""

from functools import cache

import sympy as s


@cache
def data():
    t = s.Symbol("t", real=True)
    b, N = s.symbols("log_a lapse", real=True)
    v, vd, vdd, Nd, Ndd, w = s.symbols("v vdot vddot Ndot Nddot chidot", real=True)
    h = s.Function("h")(t)
    f0 = s.Function("F0")(t)
    f1 = s.Function("FX")(t)
    q3, q4, q5 = (s.Function(name)(t) for name in ("A3clock", "A4clock", "A5clock"))
    volume = s.exp(3 * b)
    X = N**-2
    R = -6 * (vd / N**2 - v * Nd / N**3 + 2 * v * v / N**2)
    box = -Nd / N**3 + 3 * v / N**2
    Z = -Nd / N**5
    L3 = Z * box
    L4 = Nd**2 / N**8
    L5 = Z * Z
    # For first metric variation these are exact clock-jet representatives.
    pieces = {
        "Einstein": -R / 2,
        "scalar_F": f0 + f1 * (X - 1),
        "curvature_nonEinstein": -(X - 1) * R / (2 * h),
        "A3": q3 * L3,
        "A4": q4 * L4,
        "A5": q5 * L5,
        "original_M1": w * w / (2 * N * N),
    }
    D = lambda f: (
        s.diff(f, t)
        + s.diff(f, b) * v
        + s.diff(f, v) * vd
        + s.diff(f, vd) * vdd
        + s.diff(f, N) * Nd
        + s.diff(f, Nd) * Ndd
    )
    gauge = {N: 1, Nd: 0, Ndd: 0}
    results = {}
    for name, kernel in pieces.items():
        L = N * volume * kernel
        rho = s.simplify(-(s.diff(L, N) - D(s.diff(L, Nd))).subs(gauge) / volume)
        P = s.simplify(
            (s.diff(L, b) - D(s.diff(L, v)) + D(D(s.diff(L, vd)))).subs(gauge)
            / (3 * volume)
        )
        results[name] = {
            "rho_Euler": rho,
            "P_Euler": P,
            "null_Euler": s.simplify(rho + P),
        }
    expected = {
        "Einstein": (-3 * v * v, 2 * vd + 3 * v * v),
        "scalar_F": (2 * f1 - f0, f0),
        "curvature_nonEinstein": (6 * (vd + 2 * v * v) / h, 0),
        "A3": (-3 * (q3 * vd + s.diff(q3, t) * v + 3 * q3 * v * v), 0),
        "A4": (0, 0),
        "A5": (0, 0),
        "original_M1": (w * w / 2, w * w / 2),
    }
    checks = {}
    for name, (rho, P) in expected.items():
        checks[name + "_energy"] = s.simplify(results[name]["rho_Euler"] - rho)
        checks[name + "_pressure"] = s.simplify(results[name]["P_Euler"] - P)
    checks.update(
        {
            "L3_complete_lapse_expression": s.simplify(
                L3 - Nd * Nd / N**8 + 3 * v * Nd / N**7
            ),
            "L4_complete_lapse_expression": L4 - Nd * Nd / N**8,
            "L5_complete_lapse_expression": L5 - Nd * Nd / N**10,
            "curvature_nonidentity_value_zero": pieces["curvature_nonEinstein"].subs(
                gauge
            ),
            "A3_value_zero_not_first_variation": pieces["A3"].subs(gauge),
        }
    )
    return {
        "symbols": {
            "t": t,
            "h": h,
            "H": v,
            "Hdot": vd,
            "F0": f0,
            "FX": f1,
            "A3": q3,
            "chidot": w,
        },
        "complete_homogeneous_invariants": {
            "X": X,
            "R": R,
            "box_u": box,
            "Z": Z,
            "L3": L3,
            "L4": L4,
            "L5": L5,
        },
        "complete_piecewise_first_metric_variations": results,
        "representative_scope": "F and curvature use exact first clock jets; the A_i first variations need their clock values only. These local representatives are not proposed global replacements of the frozen S6.109 analytic target. The physical action is kappa times this dimensionless kernel.",
        "variation_warning": "The lapse and scale factor are varied before imposing N=1,Ndot=0. The curvature nonidentity and A3 action values vanish on the clock but their metric variations do not. A4,A5 have zero first background variation, not zero constraint/second-variation effects.",
        "checks": checks,
    }
