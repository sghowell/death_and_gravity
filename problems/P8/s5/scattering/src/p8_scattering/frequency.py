"""Direct frequency-seeded jets, retaining the exact quartic Legendre term.

Y=a^(3/2) Qc is the oscillator coordinate. An incoming frequency w inserts
dot Qc=(-i*w-3H/2)*Y at the evaluation point a=1. The Qc normalization
derivative is also retained. This constructs local symbols, not an S matrix.
"""

from dataclasses import dataclass, replace
from functools import cache

import sympy as sp
from p8_physical import geometry, momentum
from p8_physical import jets as j
from p8_physical.vertices import (
    chart_functions,
    invariant_coefficients,
    nonexceptional,
    polarization,
    tensor_basis,
)
from p8_s5 import compact
from sympy.polys.domains import QQ_I

PARTITIONS = (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2)))


@dataclass(frozen=True)
class Leg:
    wave: tuple
    kind: str
    frequency: object = 0
    polarization: tuple | None = None


class Context(j.Context):
    def __init__(self, waves, parameters):
        super().__init__(waves)
        if parameters:
            self.domain = QQ_I.frac_field(*parameters)
            self.zero_coefficient = self.domain.zero
            self.one_coefficient = self.domain.one


def response(chart, wave, time_point):
    q = sp.Symbol("q", positive=True)
    wave2 = sum(sp.Rational(k)**2 for k in wave)
    data = {key: value.subs(q, wave2) for key, value in chart_functions(chart, q).items()}
    if time_point is not None:
        data = {key: sp.cancel(value.subs(compact.x, time_point)) for key, value in data.items()}
        if any(v.has(sp.zoo, sp.nan, sp.oo, -sp.oo) for v in data.values()) or data["K"] <= 0:
            raise ValueError("Nonpositive or singular scalar velocity chart")
    return data


@cache
def hamiltonian(legs, time_point, chart="gamma", volume=True):
    if len(legs) not in (2, 3, 4):
        raise ValueError("Use two through four frequency legs")
    nonexceptional(legs)
    if time_point is not None:
        time_point = sp.Rational(time_point)
        if abs(time_point) > 1:
            raise ValueError("Compact time must be in [-1,1]")
    x = compact.x if time_point is None else time_point
    H, Adot = 2*x, 2+8*x*x
    parameters = set().union(*(sp.sympify(leg.frequency).free_symbols for leg in legs))
    if time_point is None:
        parameters.add(compact.x)
    c = Context([leg.wave for leg in legs], sorted(parameters, key=str))
    zeta, p = c.jet(), c.jet()
    tensor, tensor_p = j.zeros(c), j.zeros(c)
    scalar_coordinate, scalar_momentum = c.jet(), c.jet()
    normalization = sp.Integer(1)
    for i, leg in enumerate(legs):
        amplitude = c.leg(i)
        q = sum(sp.Rational(k)**2 for k in leg.wave)
        nu = -sp.I*sp.sympify(leg.frequency)-(3*H/2 if volume else 0)
        if leg.kind in ("s", "p"):
            if leg.kind == "s":
                data = response(chart, leg.wave, time_point)
                S, P = amplitude, (data["alpha"]*(nu-data["logZdot"])+data["beta"])*amplitude
                normalization /= sp.sqrt(2*data["K"])
            else:
                S, P = c.jet(), amplitude
            scalar_coordinate += S
            scalar_momentum += P
            if chart == "unitary":
                zeta += S
                p += P
            elif chart == "gamma":
                zeta += P/(2*q)
                p += -2*q*S
            else:
                raise ValueError("Use unitary or gamma chart")
        elif leg.kind in ("t", "pi"):
            E, norm = polarization(leg)
            if leg.kind == "t":
                normalization *= 2/sp.sqrt(norm)
                coordinate, conjugate = amplitude, nu*amplitude/4
            else:
                coordinate, conjugate = c.jet(), amplitude/norm
            tensor = j.madd(tensor, [[E[a, b]*coordinate for b in range(3)] for a in range(3)])
            tensor_p = j.madd(tensor_p, [[E[a, b]*conjugate for b in range(3)] for a in range(3)])
        else:
            raise ValueError(f"Unknown frequency leg kind: {leg.kind}")
    geo = geometry.derive(c, zeta, tensor)
    mom = momentum.derive(c, geo, zeta, tensor, p, tensor_p, H)
    mixed = j.mscale(j.mmul(mom["momentum"], geo["metric"]), geo["volume"].power(-1))
    tr = j.trace(mixed)
    sigma, shear2 = tr+3*H, j.trace(j.mmul(mixed, mixed))-tr*tr/3
    density = c.jet()
    for order in invariant_coefficients()[:len(legs)+1]:
        for (a, b, d), coefficient in order:
            if time_point is not None:
                coefficient = coefficient.subs(compact.x, time_point)
            density += coefficient*sigma**a*geo["curvature"]**b*shear2**d
    h = geo["volume"]*density-2*H*j.contract(mom["momentum"], geo["metric"])
    h -= Adot*(6*zeta+3*zeta*zeta-j.contract(tensor, tensor)/2)
    if chart == "gamma":
        h -= H*scalar_coordinate*scalar_momentum
    return {"rational_kernel": sp.factor(h.coefficient(c.full)), "normalization": sp.simplify(normalization),
            "constraints": mom["checks"]}


@cache
def vertex(legs, time_point, chart="gamma", volume=True):
    if len(legs) not in (3, 4) or any(leg.kind not in ("s", "t") for leg in legs):
        raise ValueError("Use three/four scalar or tensor frequency legs")
    base = hamiltonian(legs, time_point, chart, volume)
    contact = -base["rational_kernel"]
    corrections = []
    if len(legs) == 4:
        for left, right in PARTITIONS:
            wave = tuple(sum(sp.Rational(legs[i].wave[j]) for i in left) for j in range(3))
            opposite = tuple(-k for k in wave)
            alpha = response(chart, wave, time_point)["alpha"]
            def source(indices, internal_wave, kind, E=None):
                inputs = tuple(legs[i] for i in indices)+(Leg(internal_wave, kind, 0, E),)
                return hamiltonian(inputs, time_point, chart, volume)["rational_kernel"]
            scalar = sp.factor(alpha*source(left, opposite, "p")*source(right, wave, "p"))
            tensor = sp.Integer(0)
            for E in tensor_basis(wave):
                norm = sp.trace(sp.Matrix(E)**2)
                tensor += norm*source(left, opposite, "pi", E)*source(right, wave, "pi", E)/4
            corrections.append({"partition": (left, right), "scalar": scalar, "tensor": sp.factor(tensor)})
    rational = sp.factor(contact+sum(e["scalar"]+e["tensor"] for e in corrections))
    return {"kernel": sp.factor(base["normalization"]*rational), "rational_kernel": rational,
            "normalization": base["normalization"], "rational_minus_H4": contact,
            "Legendre_corrections": corrections}


def reverse(leg):
    return replace(leg, wave=tuple(-k for k in leg.wave), frequency=-leg.frequency)
