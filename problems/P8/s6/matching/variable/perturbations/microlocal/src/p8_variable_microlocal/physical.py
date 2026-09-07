"""Exact finite-K physical Cauchy map in (pi,xi,chi_B).

This differs from the frozen parent's (xi,chi_B,Psi_f_B) chart. The
relative spatial scalar pi is gauge invariant after the common spatial
gauge; xi and chi_B retain the prescribed physical g metric. A rational
symbol proof, not a point fixture, is needed for uniform principal claims.
"""

from functools import cache

import sympy as sp
from p8_variable_constraints import action, reduction


@cache
def rows():
    d = action.derive()
    a, w, K = d["a"], d["chi_speed"], action.K
    pi = sp.Matrix([[0, -1/K, 0, 0, 0, 0]])
    xi = sp.Matrix([[0, 0, 3*a**3*w/(2*a*K), 0, -3/(2*a*K), 0]])
    chi = sp.Matrix([[0, 0, 1, 0, 0, 0]])+w*xi
    return pi.col_join(xi).col_join(chi)


def equation(time, lapse):
    # Parent validation occurs before caching, including bool/float guards.
    data = reduction.jets(time, lapse, 1)
    return _equation(sp.sympify(time), sp.sympify(lapse), data)


def _equation(time, lapse, data):
    d = action.derive()
    u = d["u"]
    at = {u: time, d["c"]: lapse}
    O = tuple(sp.diff(rows(), u, j).subs(at) for j in range(3))
    A = reduction.J*data["H"][0]
    Ap = reduction.J*data["H"][1]
    C = O[0].col_join(O[1]+O[0]*A).applyfunc(sp.cancel)
    determinant = sp.factor(C.det(method="domain-ge"))
    if determinant == 0:
        raise ValueError("The selected finite-K physical Cauchy map is singular")
    inverse = C.inv().applyfunc(sp.cancel)
    acceleration = (O[2]+2*O[1]*A+O[0]*(Ap+A*A)).applyfunc(sp.cancel)
    equation = (acceleration*inverse).applyfunc(sp.cancel)
    symplectic = (inverse.T*reduction.J*inverse).applyfunc(sp.cancel)
    return {"time": time, "lapse": lapse, "O": O, "C": C,
            "det_C": determinant, "inverse_C": inverse,
            "Mq": equation[:, :3], "Mv": equation[:, 3:],
            "symplectic": symplectic,
            "q_bracket": (O[0]*reduction.J*O[0].T).applyfunc(sp.factor)}


def degree(expression):
    if expression == 0:
        return -sp.oo
    numerator, denominator = sp.fraction(sp.cancel(expression))
    return sp.degree(numerator, action.K)-sp.degree(denominator, action.K)
