"""Algebraic lapse elimination through quartic invariant degree.

The result is NOT the fully spatially constrained scalar/tensor Hamiltonian.
It retains exact curvature and momentum invariants and their spatial geometry.
"""

from functools import cache

import sympy as sp

from . import nonlinear_d as model


def stationary_series(A, L, Q):
    """h=A(n)+L(n)+Q(n); weights L=1, Q=2; inputs are derivatives, not Taylor coefficients."""
    n1 = -L[1]/A[2]
    force2 = A[3]*n1**2/2 + L[2]*n1 + Q[1]
    n2 = -force2/A[2]
    n3 = -(A[3]*n1*n2 + A[4]*n1**3/6 + L[2]*n2
           + L[3]*n1**2/2 + Q[2]*n1)/A[2]
    orders = (A[0], L[0], Q[0]-L[1]**2/(2*A[2]),
              A[3]*n1**3/6+L[2]*n1**2/2+Q[1]*n1,
              A[4]*n1**4/24+L[3]*n1**3/6+Q[2]*n1**2/2-force2**2/(2*A[2]))
    return {"lapse": (n1, n2, n3), "hamiltonian": orders}


@cache
def generic_checks():
    """Expand before substituting the model; check every order, not sample times."""
    n, eps = sp.symbols("n eps")
    A = sp.symbols("A0:5")
    A = (A[0], sp.Integer(0), *A[2:])
    L, Q = sp.symbols("L0:4"), sp.symbols("Q0:3")
    result = stationary_series(A, L, Q)
    h = sum(A[k]*n**k/sp.factorial(k) for k in range(5))
    h += eps*sum(L[k]*n**k/sp.factorial(k) for k in range(4))
    h += eps**2*sum(Q[k]*n**k/sp.factorial(k) for k in range(3))
    # Polynomial substitution is exact; no single Fourier mode or spatial averaging.
    ns = sum(eps**(i+1)*v for i, v in enumerate(result["lapse"]))
    force = sp.Poly(sp.diff(h, n).subs(n, ns), eps)
    reduced = sp.Poly(h.subs(n, ns), eps)
    return {**{f"constraint_order_{k}": sp.cancel(force.nth(k)) for k in range(4)},
            **{f"Hamiltonian_order_{k}": sp.cancel(reduced.nth(k)-result["hamiltonian"][k])
               for k in range(5)}}


@cache
def derive():
    jets = model.lapse_jets()
    # Expand only the invariant polynomial first. Expanding the rational time
    # functions at the same time causes needless denominator/expression blowup.
    symbolic = {key: sp.symbols(f"{key}0:5") for key in jets}
    mapping = {symbolic[key][k]: jets[key][k] for key in jets for k in range(5)}
    L = tuple(symbolic["B"][k]*model.sigma+symbolic["C"][k]*model.rho for k in range(4))
    Q = tuple(symbolic["D"][k]*model.sigma**2+symbolic["E"][k]*model.shear2 for k in range(3))
    raw = stationary_series(symbolic["A"], L, Q)
    variables = (model.sigma, model.rho, model.shear2)
    def specialize(expr):
        poly = sp.Poly(expr, *variables)
        return sp.Add(*(sp.factor(value.subs(mapping))*sp.prod(v**p for v, p in zip(variables, powers))
                        for powers, value in poly.terms()))
    return {key: tuple(specialize(expr) for expr in values) for key, values in raw.items()}


def coefficients(expr):
    """Canonical monomial representation; shear2 has weight two."""
    return {",".join(map(str, powers)): sp.factor(value)
            for powers, value in sp.Poly(expr, model.sigma, model.rho, model.shear2).terms()}
