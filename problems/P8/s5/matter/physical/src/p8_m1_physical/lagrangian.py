"""Unnormalized physical velocity vertices including matrix Legendre contact.

The full coupled scalar inverse Hessian, including off-diagonal terms, and
both tensor polarizations are contracted at every quartic internal momentum.
This contact correction is not the cubic-exchange scattering amplitude.
"""

from functools import cache

import sympy as sp

from . import quadratic
from .vertices import Leg, hamiltonian_kernel, tensor_basis


@cache
def stationary_checks():
    # Independent two-component arbitrary symmetric quadratic completion.
    a, b, c = sp.symbols("a b c", real=True)
    A = sp.Matrix([[a, b], [b, c]])
    g = sp.Matrix(sp.symbols("g0 g1", real=True))
    dp = -A.inv()*g
    actual = -(dp.T*A*dp)[0]/2-(g.T*dp)[0]
    expected = (g.T*A.inv()*g)[0]/2
    return {"two_scalar_quartic_Legendre": sp.cancel(actual-expected),
            "off_diagonal_contact": sp.cancel(sp.diff(expected, g[0], g[1])+b/(a*c-b**2))}


@cache
def kernel(legs, time_point=None, chart="unitary"):
    if len(legs) not in (2, 3, 4) or any(leg.kind not in ("s", "s_dot", "m", "m_dot", "t", "t_dot") for leg in legs):
        raise ValueError("Velocity kernels need 2--4 scalar/matter/tensor coordinate or velocity legs")
    base = hamiltonian_kernel(legs, time_point, chart, "velocity")
    raw = -base["kernel"]
    if len(legs) == 2:
        raw += base["symplectic"]
    corrections = []
    if len(legs) == 4:
        partitions = (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2)))
        for left, right in partitions:
            wave = tuple(sum(sp.Rational(legs[i].wave[j]) for i in left) for j in range(3))
            opposite = tuple(-value for value in wave)
            q = sum(value**2 for value in wave)
            alpha = quadratic.response(time_point, chart, q)["alpha"]

            def source(indices, momentum_wave, kind, E=None):
                inputs = tuple(legs[i] for i in indices)+(Leg(momentum_wave, kind, E),)
                return hamiltonian_kernel(inputs, time_point, chart, "velocity")["kernel"]

            scalar_left = sp.Matrix([source(left, opposite, kind) for kind in ("p", "P")])
            scalar_right = sp.Matrix([source(right, wave, kind) for kind in ("p", "P")])
            scalar = sp.factor((scalar_left.T*alpha*scalar_right)[0])
            diagonal_only = sum(scalar_left[i]*alpha[i, i]*scalar_right[i] for i in range(2))
            tensor = sp.Integer(0)
            tensor_pieces = []
            for E in tensor_basis(wave):
                norm = sp.trace(sp.Matrix(E)**2)
                piece = sp.factor(norm*source(left, opposite, "pi", E)*source(right, wave, "pi", E)/4)
                tensor_pieces.append(piece)
                tensor += piece
            raw += scalar+tensor
            corrections.append({"partition": (left, right), "scalar": scalar,
                                "scalar_off_diagonal": sp.factor(scalar-diagonal_only),
                                "tensor": sp.factor(tensor), "tensor_polarizations": tuple(tensor_pieces)})
    return {"kernel": sp.factor(raw), "minus_H": -base["kernel"],
            "Legendre_corrections": corrections, "chart": chart, "normalized": False,
            "not_a_scattering_amplitude": True}
