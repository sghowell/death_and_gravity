"""Recursive solution of all spatial momentum constraints in flat York form."""

from itertools import product

import sympy as sp

from . import jets as j


class ZeroMomentumConstraint(ValueError):
    """A nonzero source at zero total wavevector needs a separate homogeneous treatment."""


def york(vector):
    divergence = sum(vector[i].derivative(i) for i in range(3))
    return [[vector[i].derivative(k)+vector[k].derivative(i)
             -(sp.Rational(2, 3)*divergence if i == k else 0)
             for k in range(3)] for i in range(3)]


def solve_vector(context, source):
    """Solve (Delta I+grad div/3) W=-source at each nonzero wavevector."""
    result = [{} for _ in range(3)]
    masks = set().union(*(v.data.keys() for v in source))
    zero = context.zero_coefficient
    for mask in masks:
        wave = context.wave[mask]
        k2 = sum(k*k for k in wave)
        values = [v.data.get(mask, zero) for v in source]
        if k2 == 0:
            if any(v != zero for v in values):
                raise ZeroMomentumConstraint(f"Nonzero zero-wavevector momentum source at mask {mask}")
            continue
        longitudinal = sum(context.convert(wave[i])*values[i] for i in range(3))
        for i in range(3):
            result[i][mask] = (values[i]-context.convert(wave[i]/(4*k2))*longitudinal)/context.convert(k2)
    return [j.Jet(context, data) for data in result]


def residual(momentum, connection):
    return [sum(momentum[i][k].derivative(k) for k in range(3))
            +sum(connection[i][a][b]*momentum[a][b] for a, b in product(range(3), repeat=2))
            for i in range(3)]


def derive(context, geo, zeta, tensor, scalar_p, tensor_p, H, order=None):
    """Canonical shifts included; pi_bar=a^2 pi with local a=1 at evaluation."""
    order = context.n-1 if order is None else order
    unit = j.identity(context)
    free = j.madd(j.mscale(unit, -H*(1+zeta)+scalar_p/6),
                  j.madd(tensor_p, j.mscale(tensor, H)))
    total = free
    vector_orders = []
    for degree in range(1, order+1):
        source = [v.homogeneous(degree) for v in residual(total, geo["christoffel"])]
        vector = solve_vector(context, source)
        total = j.madd(total, york(vector))
        vector_orders.append(vector)
    remainder = residual(total, geo["christoffel"])
    checks = {f"constraint_{degree}_{i}": remainder[i].homogeneous(degree).is_zero()
              for degree in range(1, order+1) for i in range(3)}
    if not all(checks.values()):
        raise ValueError(f"Momentum reduction failed: {checks}")
    return {"momentum": total, "vector_orders": vector_orders, "checks": checks}


def projector_checks():
    k = sp.Matrix(sp.symbols("kx ky kz", real=True))
    k2 = (k.T*k)[0]
    identity = sp.eye(3)
    operator = -k2*identity-k*k.T/3
    inverse = -(identity-k*k.T/(4*k2))/k2
    return {"inverse": (operator*inverse-identity).applyfunc(sp.cancel),
            "determinant": sp.factor(operator.det()), "nonzero_domain": "k^2>0"}
