"""Actual regular canonical light system and physical source readouts."""
from functools import cache

import sympy as sp
from p8_affine_aligned import dynamics

old = dynamics.old
u = old.u
PB, PS = sp.symbols("P_b P_s", real=True)
K_IN2 = sp.Rational(1, 16)
LEFT, RIGHT = -sp.Rational(1, 2), sp.Rational(1, 2)
STATES = sp.Matrix([old.shift, old.matter, PB, PS])


def clean(value):
    return value.applyfunc(sp.cancel) if isinstance(value, sp.MatrixBase) else sp.cancel(value)


@cache
def system():
    bg = old.background()
    a, q = bg["a"], K_IN2/bg["a"]**2
    H = dynamics.first_order()["H"].subs({old.pi: 0, old.sigma: 0})
    canonical_map = {old.v: PB/(2*old.a**3*old.q), old.pm: PS/old.a**3}
    values = {symbol: bg[name] for symbol, name in
              ((old.a, "a"), (old.H, "H"), (old.theta, "theta"),
               (old.lam, "lam"), (old.J, "J"), (old.w, "w"), (old.ell, "ell"))}
    # The time-dependent canonical swap adds -H*b*P_b.
    density = old.a**3*H.subs(canonical_map, simultaneous=True)-old.H*old.shift*PB
    density = density.subs(values, simultaneous=True).subs(old.q, q)
    symplectic = sp.zeros(4)
    symplectic[:2, 2:] = sp.eye(2)
    symplectic[2:, :2] = -sp.eye(2)
    matrix = clean(symplectic*sp.hessian(density, tuple(STATES)))
    lapse = dynamics.first_order()["lapse"].subs(canonical_map, simultaneous=True)
    lapse = lapse.subs(values, simultaneous=True).subs(old.q, q)
    # zeta_metric=v+n/(2h) and omega_N=1/(2h), hence zeta_hat=v.
    # K_hat-3H=3*v'-3H*n+q*b; use the original shift constraint.
    trace = 3*(bg["theta"]-bg["H"])*lapse-3*bg["ell"]*old.matter/2+q*old.shift
    lapse_row = clean(sp.Matrix([[sp.diff(lapse, item) for item in STATES]]))
    trace_row = clean(sp.Matrix([[sp.diff(trace, item) for item in STATES]]))
    return {"a": a, "q": q, "H": density, "A": matrix, "symplectic": symplectic,
            "lapse": lapse, "lapse_row": lapse_row, "trace": trace, "trace_row": trace_row,
            "physical_v": PB/(2*a**3*q)}


@cache
def initial_jets():
    data = system()
    A = [clean(data["A"].diff(u, j).subs(u, LEFT)) for j in range(3)]
    ell = [clean(data["lapse_row"].diff(u, j).subs(u, LEFT)) for j in range(4)]
    rows = [ell[0], ell[1]+ell[0]*A[0],
            ell[2]+2*ell[1]*A[0]+ell[0]*(A[1]+A[0]**2),
            ell[3]+3*ell[2]*A[0]+3*ell[1]*(A[1]+A[0]**2)
            +ell[0]*(A[2]+2*A[1]*A[0]+A[0]*A[1]+A[0]**3)]
    rows = [clean(row) for row in rows]
    prepared = sp.Matrix.vstack(*rows[:3])
    kernel = prepared.nullspace()
    if len(kernel) != 1:
        raise ValueError("The declared preparation does not have the expected one-dimensional kernel")
    vector = kernel[0]/max(abs(value) for value in kernel[0])
    trace0 = (data["trace_row"].subs(u, LEFT)*vector)[0]
    n3 = (rows[3]*vector)[0]
    source3 = -2*n3*trace0/old.background()["h"].subs(u, LEFT)
    return {"A_jets": A, "lapse_observer_jets": rows, "prepared_matrix": prepared,
            "initial_vector": vector, "lapse_third": n3, "trace_initial": trace0,
            "source_third_initial": source3,
            "preparation_residual": prepared*vector,
            "observability_determinant": sp.Matrix.vstack(*rows).det()}
