"""Lapse-metric Ricci contraction and classical scalar-curvature variation."""

from functools import cache

import sympy as sp

t = sp.Symbol("t", real=True)
a, N, Q, kinetic, potential = [sp.Function(name)(t) for name in ("a", "N", "Q", "Tkin", "V")]


@cache
def ricci_from_metric():
    metric = sp.diag(N**2, -a**2, -a**2, -a**2)
    inverse = sp.diag(N**-2, -a**-2, -a**-2, -a**-2)

    def derivative(expr, coordinate):
        return sp.diff(expr, t) if coordinate == 0 else sp.Integer(0)

    christoffel = [[[sp.simplify(sum(inverse[r, s]*(derivative(metric[s, n], m)
                              +derivative(metric[s, m], n)-derivative(metric[m, n], s))
                              for s in range(4))/2)
                    for n in range(4)] for m in range(4)] for r in range(4)]
    ricci = sp.Matrix(4, 4, lambda m, n: sum(
        derivative(christoffel[r][m][n], r)-derivative(christoffel[r][m][r], n)
        +sum(christoffel[r][r][s]*christoffel[s][m][n]
             -christoffel[r][n][s]*christoffel[s][m][r] for s in range(4)) for r in range(4)))
    return sp.simplify(sum(inverse[m, n]*ricci[m, n] for m in range(4) for n in range(4)))


@cache
def equations():
    adot, Ndot = sp.diff(a, t), sp.diff(N, t)
    Qdot, Qddot = sp.diff(Q, t), sp.diff(Q, t, 2)
    raw = -N*a**3*Q*ricci_from_metric()/2
    boundary = 3*Q*a**2*adot/N
    gravity = -3*Q*a*adot**2/N-3*a**2*adot*Qdot/N
    lagrangian = gravity+a**3*kinetic/(2*N)-N*a**3*potential

    def EL(field):
        return sp.diff(lagrangian, field)-sp.diff(sp.diff(lagrangian, sp.diff(field, t)), t)

    cosmic = {N: 1, Ndot: 0, sp.diff(N, t, 2): 0}
    en = sp.simplify((EL(N)/a**3).subs(cosmic))
    ea = sp.simplify((EL(a)/(3*a**2)).subs(cosmic))
    H = sp.diff(a, t)/a
    Hd = sp.diff(H, t)
    raychaudhuri = sp.simplify(ea-en)
    expected_en = 3*Q*H**2+3*H*Qdot-kinetic/2-potential
    expected_ea = 2*Q*Hd+3*Q*H**2+2*H*Qdot+Qddot+kinetic/2-potential
    expected_ray = 2*Q*Hd+kinetic+Qddot-H*Qdot
    M, tau = sp.symbols("M tau", positive=True)
    scale_factor = 1+(t/tau)**2
    target_Hd = sp.factor(sp.diff(sp.diff(scale_factor, t)/scale_factor, t))
    target_stress = sp.factor(-2*M**2*target_Hd)
    residuals = {"Einstein_boundary": sp.simplify(raw-sp.diff(boundary, t)-gravity),
                 "lapse_equation": sp.simplify(en-expected_en),
                 "scale_equation": sp.simplify(ea-expected_ea),
                 "Raychaudhuri_combination": sp.simplify(raychaudhuri-expected_ray),
                 "bounce_null_stress": sp.simplify(target_stress.subs(t, 0)+4*M**2/tau**2)}
    return {"ricci_scalar": ricci_from_metric(), "gravity_boundary": boundary,
            "reduced_lagrangian": lagrangian, "lapse_equation": en, "scale_equation": ea,
            "Raychaudhuri_equation": raychaudhuri, "target_H_dot": target_Hd,
            "target_null_stress": target_stress, "bounce_stress_remainder_minimum": 4*M**2/tau**2,
            "residuals": residuals}


def field_metric_control():
    # A generic three-field LDL^T metric: no diagonal-kinetic assumption.
    l21, l31, l32 = sp.symbols("l21 l31 l32", real=True)
    d1, d2, d3 = sp.symbols("d1 d2 d3", positive=True)
    v = sp.Matrix(sp.symbols("v1:4", real=True))
    L = sp.Matrix([[1, 0, 0], [l21, 1, 0], [l31, l32, 1]])
    metric = L*sp.diag(d1, d2, d3)*L.T
    null_contraction = (v.T*metric*v)[0]
    squares = sum(d*value**2 for d, value in zip((d1, d2, d3), L.T*v))
    return sp.expand(null_contraction-squares), squares


def derive():
    result = equations()
    residual, squares = field_metric_control()
    checks = {**result["residuals"], "mixed_field_metric_squares": residual}
    if any(sp.simplify(value) != 0 for value in checks.values()):
        raise ValueError("Parent geometry or background identity failed")
    return {**{key: str(value) for key, value in result.items() if key != "residuals"},
            "mixed_field_metric_null_squares": str(squares),
            "exact_residuals": {key: "0" for key in checks},
            "minimal_parent_verdict": "NO_EXACT_CLASSICAL_FLAT_BOUNCE_WITH_POSITIVE_FIELD_METRIC_AND_CONSTANT_M",
            "scope": "Physical metric, homogeneous classical parent; no quantum NEC assumption"}
