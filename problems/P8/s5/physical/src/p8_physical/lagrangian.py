"""Velocity kernels including the quartic correction from nonlinear momenta.

L3=-H3(P0); L4=-H4(P0)+(1/2) H3_P (H2_PP)^-1 H3_P. The second term is a
Legendre/contact correction, NOT the cubic-exchange scattering diagram.
Exact canonical kinetic normalization includes the time derivative of Z.
"""

from functools import cache

import sympy as sp
from p8_s5 import compact

from .vertices import Leg, chart_functions, hamiltonian_kernel, tensor_basis


@cache
def stationary_checks():
    eps, p, v, a, b, c, g0, g1, g2, g3, h = sp.symbols("eps p v a b c g0 g1 g2 g3 h")
    H2 = a*p**2/2+b*p+c
    H3 = g0+g1*p+g2*p**2+g3*p**3
    H4 = h*p**4
    p0 = (v-b)/a
    p1 = -sp.diff(H3, p).subs(p, p0)/a
    H = H2+eps*H3+eps**2*H4
    actual = sp.Poly(sp.expand((p*v-H).subs(p, p0+eps*p1)), eps).nth(2)
    expected = -H4.subs(p, p0)+sp.diff(H3, p).subs(p, p0)**2/(2*a)
    return {"quartic_Legendre": sp.cancel(actual-expected)}


@cache
def kernel(legs, time_point=None, chart="unitary", normalize=True):
    if len(legs) not in (2, 3, 4) or any(leg.kind not in ("s", "s_dot", "t", "t_dot") for leg in legs):
        raise ValueError("Velocity kernels need 2--4 coordinate/velocity external legs")
    base = hamiltonian_kernel(legs, time_point, chart, "velocity", normalize)
    raw = -base["rational_kernel"]
    if len(legs) == 2:
        raw += base["rational_symplectic"]
    corrections = []
    if len(legs) == 4:
        # Three unordered 2+2 partitions: ordered pairs occur twice and cancel 1/2.
        partitions = (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2)))
        q_symbol = sp.Symbol("q", positive=True)
        for left, right in partitions:
            wave = tuple(sum(sp.Rational(legs[i].wave[j]) for i in left) for j in range(3))
            opposite = tuple(-v for v in wave)
            q = sum(v*v for v in wave)
            inverse_scalar = chart_functions(chart, q_symbol)["alpha"].subs(q_symbol, q)
            if time_point is not None:
                inverse_scalar = inverse_scalar.subs(compact.x, time_point)
                if inverse_scalar.has(sp.zoo, sp.nan, sp.oo, -sp.oo) or inverse_scalar <= 0:
                    raise ValueError("Internal momentum is outside the positive-kinetic velocity chart")
            # Source legs carry raw canonical momentum, so their normalization is one.
            def source(indices, momentum_wave, kind, E=None):
                inputs = tuple(legs[i] for i in indices)+(Leg(momentum_wave, kind, E),)
                return hamiltonian_kernel(inputs, time_point, chart, "velocity", normalize)["rational_kernel"]
            scalar = inverse_scalar*source(left, opposite, "p")*source(right, wave, "p")
            tensor = sp.Integer(0)
            for E in tensor_basis(wave):
                norm = sp.trace(sp.Matrix(E)**2)
                tensor += norm*source(left, opposite, "pi", E)*source(right, wave, "pi", E)/4
            value = sp.factor(scalar+tensor)
            raw += value
            corrections.append({"partition": (left, right), "rational_scalar": sp.factor(scalar),
                                "rational_tensor": sp.factor(tensor), "rational_total": value})
    raw = sp.factor(raw)
    return {"kernel": sp.factor(base["normalization"]*raw), "rational_kernel": raw,
            "normalization": base["normalization"], "rational_minus_H": -base["rational_kernel"],
            "Legendre_corrections": corrections, "chart": chart, "normalized": normalize}
