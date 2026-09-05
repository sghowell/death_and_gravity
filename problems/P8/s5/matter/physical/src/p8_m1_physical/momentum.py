"""Matter-sourced nonlinear York solution; the momentum is a density."""

from itertools import product

import sympy as sp
from p8_physical import jets as j
from p8_physical.momentum import residual as gravity_residual
from p8_physical.momentum import solve_vector, york


def residual(momentum, geo, matter_density, matter_field):
    geometric = gravity_residual(momentum, geo["christoffel"])
    return [geometric[i]-matter_density*sum(geo["inverse"][i][k]*matter_field.derivative(k)
                                          for k in range(3))/2 for i in range(3)]


def derive(context, geo, zeta, tensor, scalar_p, tensor_p, H, matter_density, matter_field, order=None):
    """Solve all 3 constraints, retaining the canonical density pi_chi.

    scalar_p here is the metric momentum BEFORE the mixed matter boundary
    shift, while the external phase variables in vertices.py are after it.
    The frozen linear metric gauge makes York momenta exactly symplectically
    orthogonal to the scalar trace and both TT polarizations after spatial IBP.
    """
    order = context.n-1 if order is None else order
    if not 0 <= order <= context.n:
        raise ValueError("Constraint order is outside this Fourier jet")
    free = j.madd(j.mscale(j.identity(context), -H*(1+zeta)+scalar_p/6),
                  j.madd(tensor_p, j.mscale(tensor, H)))
    total, vector_orders = free, []
    for degree in range(1, order+1):
        source = [value.homogeneous(degree) for value in residual(total, geo, matter_density, matter_field)]
        vector = solve_vector(context, source)
        total = j.madd(total, york(vector))
        vector_orders.append(vector)
    remaining = residual(total, geo, matter_density, matter_field)
    checks = {f"constraint_{degree}_{i}": remaining[i].homogeneous(degree).is_zero()
              for degree in range(order+1) for i in range(3)}
    if not all(checks.values()):
        raise ValueError(f"Matter momentum reduction failed: {checks}")
    return {"momentum": total, "free": free, "vector_orders": vector_orders, "checks": checks}


def symplectic_projector_checks():
    k = sp.Matrix(sp.symbols("kx ky kz", real=True))
    v = sp.Matrix(sp.symbols("vx vy vz", real=True))
    tensor = sp.Matrix(3, 3, lambda i, a: sp.I*(k[i]*v[a]+k[a]*v[i])
                       - (sp.Rational(2, 3)*sp.I*k.dot(v) if i == a else 0))
    # An arbitrary symmetric TT test tensor obeys tr E=0 and E*k=0.
    e = sp.Matrix([[sp.Symbol("e00"), sp.Symbol("e01"), sp.Symbol("e02")],
                   [sp.Symbol("e01"), sp.Symbol("e11"), sp.Symbol("e12")],
                   [sp.Symbol("e02"), sp.Symbol("e12"), sp.Symbol("e22")]])
    contraction = sum(tensor[i, a]*e[i, a] for i, a in product(range(3), repeat=2))
    return {"York_trace": sp.expand(sp.trace(tensor)),
            "York_TT_pairing_identity": sp.expand(contraction-2*sp.I*(v.T*e*k)[0]
                                                   +sp.Rational(2, 3)*sp.I*k.dot(v)*sp.trace(e))}
