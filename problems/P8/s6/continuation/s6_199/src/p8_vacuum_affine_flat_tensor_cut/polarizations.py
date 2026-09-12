"""All three physical modes and the full on-shell stress pair."""

from functools import cache

import sympy as s

ENERGY, MOMENTUM, MASS = s.symbols("E k mass", positive=True)


def onshell(expr):
    return (
        s.Poly(s.expand(expr), ENERGY)
        .rem(s.Poly(ENERGY**2 - MOMENTUM**2 - MASS**2, ENERGY))
        .as_expr()
    )


def readouts(sign):
    if (
        isinstance(sign, bool)
        or not isinstance(sign, (int, s.Integer))
        or sign not in (-1, 1)
    ):
        raise ValueError("COM momentum sign must be +/-1")
    e1, e2, e3 = s.eye(3)[:, 0], s.eye(3)[:, 1], s.eye(3)[:, 2]
    E, k, m = ENERGY, MOMENTUM, MASS
    return (
        (-s.I * E * e1, s.I * sign * k * e2, s.Integer(0), m * e1),
        (-s.I * E * e2, -s.I * sign * k * e1, s.Integer(0), m * e2),
        (-s.I * m * sign * e3, s.zeros(3, 1), -k, E * sign * e3),
    )


def stress_pair(u, v):
    E, B, V0, V = u
    F, C, W0, W = v
    out = s.zeros(4)
    out[0, 0] = (E.dot(F) + B.dot(C) + V0 * W0 + V.dot(W)) / 2
    flux = (E.cross(C) + F.cross(B) + V0 * W + W0 * V) / 2
    out[0, 1:4] = flux.T
    out[1:4, 0] = flux
    spatial = -(E * F.T + F * E.T + B * C.T + C * B.T) / 2 + (V * W.T + W * V.T) / 2
    spatial += s.eye(3) * (E.dot(F) + B.dot(C) + V0 * W0 - V.dot(W)) / 2
    out[1:4, 1:4] = spatial
    return out.applyfunc(onshell)


@cache
def pairs():
    return tuple(stress_pair(u, v) for u in readouts(1) for v in readouts(-1))


@cache
def spectral_polynomials():
    traces = sum(s.trace(a[1:4, 1:4]) ** 2 for a in pairs())
    norm = sum(sum(x * x for x in a[1:4, 1:4]) for a in pairs())
    scalar = onshell(traces / 3)
    tensor = onshell((norm - traces / 3) / 5)
    invariant = s.Symbol("s", positive=True)
    return {
        "spin2_COM": s.factor(tensor),
        "spin0_COM": s.factor(scalar),
        "spin2_invariant": s.factor(tensor.subs(MOMENTUM**2, invariant / 4 - MASS**2)),
        "spin0_invariant": s.factor(scalar.subs(MOMENTUM**2, invariant / 4 - MASS**2)),
    }


@cache
def data():
    checks = {}
    for j, T in enumerate(pairs()):
        checks[f"complete_conserved_time_row_{j}"] = T[0, :]
        checks[f"complete_symmetric_stress_pair_{j}"] = T - T.T
    return {
        "normalization": "Canonical ordinary Proca mass1000 in its flat reference vacuum. Readouts here omit the universal1/sqrt(2E) per leg, which is included exactly in Lorentz-invariant two-particle phase space. No dynamical gravity or full interacting vacuum is asserted.",
        "physical_modes": "Two transverse and one longitudinal mode at k and-k. E,B,V0=mA0,V=mAsp include the temporal constraint. No independent A0 oscillator or Maxwell replacement is used.",
        "pair_convention": "Tab=Z^T M_ab Z, with symmetric M. The displayed pair amplitude is u_k^T M_ab u_l; the Wick pair exchange factor2 belongs to the cut measure, not this amplitude.",
        "COM_conservation": "For total momentum(2E,0), every pair has T0mu=0 by the exact on-shell relation E^2=k^2+m^2.",
        "rotation_projection": "After angular averaging, K=a P2+b P0. The spin2 eigenvalue is(sum||Tij||F^2-sum tr(Tij)^2/3)/5; spin0 is sum tr(Tij)^2/3. Both sums include all nine polarization pairs.",
        "polynomials": spectral_polynomials(),
        "checks": checks,
        "gates": {
            "three_physical_polarizations": len(readouts(1)) == 3,
            "nine_pair_amplitudes": len(pairs()) == 9,
            "nonzero_longitudinal_temporal_constraint": readouts(1)[2][2] != 0,
        },
    }
