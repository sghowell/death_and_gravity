"""Independent full quadratic targets and coupled scalar velocity responses."""

from functools import cache

import sympy as sp

from . import background

q = sp.Symbol("q", positive=True)
Q = sp.Matrix(sp.symbols("v s", real=True))
P = sp.Matrix(sp.symbols("p P", real=True))
H, l, theta, lam, w, J = sp.symbols("H l theta lam w J", real=True)


@cache
def symbolic(chart):
    v, s = Q
    p, pm = P
    density = (-q*v**2+pm**2/2+(q/2-3*l**2/4)*s**2-l*p*s/2
               +(-theta*p/2+lam*q*v+w*pm/2-3*l*theta*s/2)**2/J)
    if chart == "gamma":
        density = density.subs({p: -2*q*v, v: p/(2*q)}, simultaneous=True)-H*v*p
    elif chart != "unitary":
        raise ValueError("Unknown canonical chart")
    A = sp.hessian(density, P).applyfunc(sp.factor)
    B = sp.Matrix(2, 2, lambda i, j: sp.diff(density, P[i], Q[j])).applyfunc(sp.factor)
    inverse = A.inv().applyfunc(sp.factor)
    return {"density": sp.factor(density), "A": A, "B": B,
            "alpha": inverse, "beta": (-inverse*B).applyfunc(sp.factor),
            "kinetic": inverse/2}


@cache
def target(time_point=None, chart="unitary", q_value=q):
    bg = background.functions(time_point)
    mapping = dict(zip((H, l, theta, lam, w, J), (bg[k] for k in ("H", "l", "theta", "lam", "w", "J"))))
    mapping[q] = q_value
    return sp.factor(symbolic(chart)["density"].subs(mapping, simultaneous=True))


@cache
def response(time_point, chart, q_value):
    if q_value <= 0:
        raise ValueError("A velocity response requires q>0")
    bg = background.functions(time_point)
    mapping = dict(zip((H, l, theta, lam, w, J), (bg[k] for k in ("H", "l", "theta", "lam", "w", "J"))))
    mapping[q] = q_value
    result = {key: symbolic(chart)[key].subs(mapping, simultaneous=True).applyfunc(sp.factor)
              for key in ("alpha", "beta", "kinetic")}
    if any(value.has(sp.zoo, sp.nan, sp.oo, -sp.oo) for matrix in result.values() for value in matrix):
        raise ValueError("Singular scalar velocity chart; use the regular phase Hamiltonian")
    K = result["kinetic"]
    if time_point is not None and (K[0, 0] <= 0 or K.det() <= 0):
        raise ValueError("Scalar velocity Hessian is not positive in this chart")
    return result


@cache
def identities():
    data_u, data_g = symbolic("unitary"), symbolic("gamma")
    J0 = J+w**2/2
    D = q*lam**2-J0
    Kg = sp.Matrix([[q*J0/D, -q*lam*w/(2*D)], [-q*lam*w/(2*D), (q*lam**2-J)/(2*D)]])
    expected_principal = sp.Matrix([[J0/lam**2, -w/(2*lam)], [-w/(2*lam), sp.Rational(1, 2)]])
    high_q = data_g["kinetic"].applyfunc(lambda value: sp.limit(value, q, sp.oo))
    checks = {
        "unitary_momentum_Hessian_determinant": sp.cancel(data_u["A"].det()-theta**2/(2*J)),
        "gamma_momentum_Hessian_determinant": sp.cancel(data_g["A"].det()-D/(2*J*q)),
        **{f"gamma_velocity_K_{i}{j}": sp.cancel((data_g["kinetic"]-Kg)[i, j]) for i in range(2) for j in range(2)},
        **{f"gamma_principal_K_{i}{j}": sp.cancel((high_q-expected_principal)[i, j]) for i in range(2) for j in range(2)},
    }
    for label, data in (("unitary", data_u), ("gamma", data_g)):
        checks.update({f"{label}_response_alpha_{i}{j}": sp.cancel((data["A"]*data["alpha"]-sp.eye(2))[i, j])
                       for i in range(2) for j in range(2)})
        checks.update({f"{label}_response_beta_{i}{j}": sp.cancel((data["A"]*data["beta"]+data["B"])[i, j])
                       for i in range(2) for j in range(2)})
    return checks
