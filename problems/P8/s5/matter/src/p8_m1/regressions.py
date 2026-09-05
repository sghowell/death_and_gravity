"""Independent literal bounce formulas and intentionally broken reductions."""

from functools import cache

import sympy as sp

from . import nonlinear as model


@cache
def bounce_checks():
    N = model.N
    # Integrate Psi_u,N at u=0 BEFORE forgetting the u derivative:
    # Psi_u,N=9(N^(-5/2)-N^(-1/2)), with Psi_u(1)=0.
    psi_u = 24-6*N**sp.Rational(-3, 2)-18*sp.sqrt(N)
    direct = {
        "A": (sp.Rational(9, 8)*N**sp.Rational(5, 2)-sp.Rational(2849, 200)*sp.sqrt(N)
              -sp.Rational(23, 8)*N**sp.Rational(-3, 2)+24+N**sp.Rational(-1, 2)/200),
        "B": sp.Integer(0), "C": -N**sp.Rational(-1, 2)/2,
        "L": N**sp.Rational(-1, 2)/10,
        "D": -N**sp.Rational(3, 2)/3, "E": 2*N**sp.Rational(3, 2),
        "M": N**sp.Rational(-1, 2)/2, "Z": N**sp.Rational(3, 2)/2,
    }
    jets = model.lapse_jets()
    checks = {f"literal_bounce_{key}_{k}": sp.simplify(sp.diff(expr, N, k).subs(N, 1)-jets[key][k].subs(model.u, 0))
              for key, expr in direct.items() for k in range(5)}
    checks["literal_boundary_u_derivative"] = sp.simplify(sp.diff(psi_u, N)
                                                          -9*(N**sp.Rational(-5, 2)-N**sp.Rational(-1, 2)))
    checks["literal_boundary_u_constant"] = psi_u.subs(N, 1)
    return checks


@cache
def boundary_negative_control():
    omitted = model.lapse_jets(omit_boundary=True)["A"][2]
    residual = sp.factor(omitted+2*model.functions()["J"])
    expected = 18*(1-model.u**2)/model.d**8
    if sp.cancel(residual-expected) != 0 or residual == 0:
        raise ValueError("Dropped-boundary negative control failed")
    return {"omitted_boundary_J_bridge_residual": str(residual), "at_bounce": str(residual.subs(model.u, 0)),
            "accidental_zeros": ["-1", "1"],
            "psi_at_bounce_is_zero_but_its_time_derivative_is_not": True}
