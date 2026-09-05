"""Canonical physical phase-space kernels for nonexceptional spatial momenta.

Metric gauge: hat_h=a^2[(1+2 zeta)I+gamma_TT]. Momentum constraints are
solved, not set to their linear values. Background and quadratic canonical
shifts, plus the time-dependent a^2 coordinate change, are retained exactly.
The local background unit ell is held CONSTANT within each evaluation patch.
"""

from dataclasses import dataclass
from functools import cache
from itertools import combinations

import sympy as sp
from p8_s5 import compact, lapse_series

from . import geometry, momentum
from . import jets as j


@dataclass(frozen=True)
class Leg:
    wave: tuple
    kind: str  # s, p, t, pi; s_dot/t_dot additionally in velocity representation
    polarization: tuple | None = None


def tensor_basis(wave):
    """Two rational, transverse, tracefree, orthogonal polarizations; norms explicit."""
    k = sp.Matrix(wave)
    k2 = k.dot(k)
    if k2 == 0:
        raise ValueError("Tensor polarization at zero momentum needs a separate convention")
    e = sp.Matrix([-k[1], k[0], 0]) if k[0]**2+k[1]**2 else sp.Matrix([1, 0, 0])
    f = k.cross(e)
    e2 = e.dot(e)
    plus = (e*e.T-f*f.T/k2)/e2
    cross = (e*f.T+f*e.T)/e2
    return tuple(tuple(tuple(v for v in row) for row in matrix.tolist()) for matrix in (plus, cross))


def polarization(leg):
    if leg.polarization is None:
        raise ValueError("Supply a tensor polarization")
    matrix = sp.Matrix(leg.polarization)
    if matrix.shape != (3, 3) or matrix != matrix.T or sp.trace(matrix) != 0:
        raise ValueError("Polarization must be symmetric and tracefree")
    if matrix*sp.Matrix(leg.wave) != sp.zeros(3, 1):
        raise ValueError("Polarization must be transverse")
    norm = sp.trace(matrix*matrix)
    if norm <= 0:
        raise ValueError("Tensor polarization must be nonzero and real")
    return matrix, norm


def nonexceptional(legs):
    waves = tuple(tuple(map(sp.Rational, leg.wave)) for leg in legs)
    if any(sum(k[i] for k in waves) != 0 for i in range(3)):
        raise ValueError("Integrated kernels require total spatial momentum zero")
    # A homogeneous scalar one-leg entry is allowed solely as a tadpole check.
    if len(legs) == 1:
        return
    for size in range(1, len(legs)):
        for subset in combinations(range(len(legs)), size):
            if all(sum(waves[k][i] for k in subset) == 0 for i in range(3)):
                raise momentum.ZeroMomentumConstraint("A proper subset has zero momentum; do not drop its constraint")


@cache
def invariant_coefficients():
    result = []
    for order in lapse_series.derive()["hamiltonian"]:
        coefficients = []
        for key, value in lapse_series.coefficients(order).items():
            p, r, w = map(int, key.split(","))
            exponent = 1-sp.Rational(p, 2)-r-w
            coefficients.append(((p, r, w), compact.normalized(value, exponent, p % 2)))
        result.append(tuple(coefficients))
    return tuple(result)


@cache
def chart_functions(chart, q):
    x, P = compact.x, compact.P
    J, theta, lam, H = 2*P, 2*x, 1-2*(1-x**2)**3, 2*x
    if chart == "unitary":
        K = J/theta**2
        alpha = 2*K
        beta = 2*lam*q/theta
        weight = 0
    elif chart == "gamma":
        K = q*J/(q*lam**2-J)
        alpha = 2*K
        beta = -alpha*(theta*lam*q/J-H)
        weight = 1
    else:
        raise ValueError("Use unitary or gamma chart")
    # Local constant-unit time derivative: x'=1-x^2; the compact q parameter
    # has moving-unit derivative -2*x*q; remove normalization weight 2*w*x.
    logZdot = sp.cancel(((1-x**2)*sp.diff(K, x)-2*x*q*sp.diff(K, q)-2*weight*x*K)/(2*K))
    return {"K": sp.factor(K), "alpha": sp.factor(alpha), "beta": sp.factor(beta),
            "logZdot": sp.factor(logZdot)}


def fields(context, legs, time_point, chart, representation, normalize):
    if representation not in ("phase", "velocity"):
        raise ValueError("Use phase or velocity representation")
    if normalize and representation != "velocity":
        raise ValueError("Configuration normalization belongs to velocity kernels")
    s, p, s_dot = context.jet(), context.jet(), context.jet()
    tensor, tensor_p, tensor_dot = j.zeros(context), j.zeros(context), j.zeros(context)
    normalization = sp.Integer(1)
    q_symbol = sp.Symbol("q", positive=True)
    for index, leg in enumerate(legs):
        amplitude = context.leg(index)
        q = sum(sp.Rational(k)**2 for k in leg.wave)
        def at_time(expr, q=q):
            expr = expr.subs(q_symbol, q)
            return sp.cancel(expr if time_point is None else expr.subs(compact.x, time_point))
        if leg.kind in ("s", "s_dot", "p"):
            if representation == "velocity" and leg.kind != "p":
                response = {key: at_time(value) for key, value in chart_functions(chart, q_symbol).items()}
                if any(v.has(sp.zoo, sp.nan, sp.oo, -sp.oo) for v in response.values()):
                    raise ValueError("Singular scalar velocity chart")
                if time_point is not None and response["K"] <= 0:
                    raise ValueError("Scalar kinetic coefficient is not positive in this velocity chart")
                if normalize:
                    normalization /= sp.sqrt(2*response["K"])
                velocity_mix = -response["logZdot"] if normalize else 0
                if leg.kind == "s":
                    s += amplitude
                    s_dot += velocity_mix*amplitude
                    p += (response["beta"]+response["alpha"]*velocity_mix)*amplitude
                else:
                    s_dot += amplitude
                    p += response["alpha"]*amplitude
            elif leg.kind == "s":
                s += amplitude
            elif leg.kind == "p":
                p += amplitude
            else:
                raise ValueError("Velocity leg requested in phase representation")
        elif leg.kind in ("t", "t_dot", "pi"):
            matrix, norm = polarization(leg)
            if normalize and leg.kind != "pi":
                normalization *= 2/sp.sqrt(norm)
            if leg.kind == "t":
                tensor = j.madd(tensor, [[matrix[i, k]*amplitude for k in range(3)] for i in range(3)])
            elif leg.kind == "pi":
                tensor_p = j.madd(tensor_p, [[matrix[i, k]*amplitude/norm for k in range(3)] for i in range(3)])
            elif representation == "velocity":
                tensor_p = j.madd(tensor_p, [[matrix[i, k]*amplitude/4 for k in range(3)] for i in range(3)])
                tensor_dot = j.madd(tensor_dot, [[matrix[i, k]*amplitude for k in range(3)] for i in range(3)])
            else:
                raise ValueError("Velocity leg requested in phase representation")
        else:
            raise ValueError(f"Unknown external leg kind: {leg.kind}")
    if chart == "unitary":
        zeta, scalar_p = s, p
        generator = context.jet()
    elif chart == "gamma":
        def rescale_modes(jet, inverse=False):
            values = {}
            for mask, coefficient in jet.data.items():
                q = sum(k*k for k in context.wave[mask])
                if q == 0:
                    raise momentum.ZeroMomentumConstraint("Gamma canonical swap requires nonzero momentum")
                factor = 1/(2*q) if inverse else -2*q
                values[mask] = coefficient*context.convert(factor)
            return j.Jet(context, values)
        zeta, scalar_p = rescale_modes(p, True), rescale_modes(s)
        H = 2*(compact.x if time_point is None else time_point)
        generator = -H*s*p  # -Cdot b zeta, C=2*a^3*q, Cdot=H*C
    else:
        raise ValueError("Use unitary or gamma chart")
    symplectic = p*s_dot+j.contract(tensor_p, tensor_dot)
    return zeta, tensor, scalar_p, tensor_p, generator, normalization, symplectic


@cache
def hamiltonian_kernel(legs, time_point=None, chart="unitary", representation="phase", normalize=False):
    """Return the labelled integrated coefficient; no 1/n! is inserted."""
    if not 1 <= len(legs) <= 4:
        raise ValueError("Only first through fourth order are supported")
    nonexceptional(legs)
    if time_point is not None:
        time_point = sp.Rational(time_point)
        if abs(time_point) > 1:
            raise ValueError("Compact time is in [-1,1]")
    context = j.Context([leg.wave for leg in legs], parameter=compact.x if time_point is None else None)
    x = compact.x if time_point is None else time_point
    H, Adot = 2*x, 2+8*x**2  # Adot/a^3=Hdot+3H^2 in the fixed local units
    zeta, tensor, scalar_p, tensor_p, generator, normalization, symplectic = fields(
        context, legs, time_point, chart, representation, normalize)
    geo = geometry.derive(context, zeta, tensor)
    mom = momentum.derive(context, geo, zeta, tensor, scalar_p, tensor_p, H)
    mixed = j.mscale(j.mmul(mom["momentum"], geo["metric"]), geo["volume"].power(-1))
    tr = j.trace(mixed)
    sigma = tr+3*H
    shear2 = j.trace(j.mmul(mixed, mixed))-tr**2/3
    if not sigma.homogeneous(0).is_zero() or not shear2.homogeneous(1).is_zero():
        raise ValueError("Momentum invariant background/order check failed")
    density = context.jet()
    for order in invariant_coefficients()[:len(legs)+1]:
        for (p, r, w), coefficient in order:
            if time_point is not None:
                coefficient = coefficient.subs(compact.x, time_point)
            density += coefficient*sigma**p*geo["curvature"]**r*shear2**w
    h = geo["volume"]*density-2*H*j.contract(mom["momentum"], geo["metric"])
    h -= Adot*(6*zeta+3*zeta**2-j.contract(tensor, tensor)/2)
    h += generator
    return {"rational_kernel": sp.factor(h.coefficient(context.full)),
            "rational_symplectic": sp.factor(symplectic.coefficient(context.full)),
            "normalization": normalization,
            "kernel": sp.factor(normalization*h.coefficient(context.full)),
            "constraint_checks": mom["checks"],
            "chart": chart, "representation": representation, "normalized": normalize}
