"""The angular Einstein equation is implied by the (tt, tr, rr) equations and div T = 0.

Lemma.  Let E^{mn} = G^{mn} - 8 pi T^{mn} for the spherically symmetric metric of
``einstein_euler`` with a perfect fluid.  If nabla_m T^{mn} = 0 and E^{tt} = E^{tr} = E^{rr}
= 0 on an open set (so their derivatives vanish too), then E^{theta theta} = E^{phi phi} = 0.
Proof: nabla_m G^{mn} = 0 identically (``contracted_bianchi``), hence nabla_m E^{mn} = 0;
for a spherically symmetric symmetric tensor with vanishing (t,r)-block this divergence is
c * E^{theta theta} with c = -2 r / a^2 != 0 (``angular_equation_coefficient``); and
E^{phi phi} = E^{theta theta}/sin^2 theta (``angular_structure``).
"""
from __future__ import annotations

import sympy as sp

from .einstein_euler import (
    a,
    christoffel,
    divergence,
    einstein,
    metric,
    r,
    ricci,
    stress_tensor,
    t,
    th,
)


def contracted_bianchi():
    """nabla_m G^{mn} for the general metric (third derivatives of alpha, a): must be 0."""
    g, ginv = metric()
    Gam = christoffel(g, ginv)
    G, _ = einstein(g, ginv, ricci(Gam))
    G_up = (ginv * G * ginv).applyfunc(sp.cancel)
    return [sp.cancel(v) for v in divergence(G_up, Gam)]


def angular_equation_coefficient():
    """For E^{mn} = diag(0, 0, Eth, Eth/sin^2), returns (nabla_m E^{mr} / Eth, nabla_m E^{mt})."""
    Eth = sp.Function("Eth")(t, r)
    E = sp.zeros(4)
    E[2, 2] = Eth
    E[3, 3] = Eth / sp.sin(th) ** 2
    g, ginv = metric()
    D = divergence(E, christoffel(g, ginv))
    return sp.cancel(D[1] / Eth), sp.cancel(D[0])


def angular_structure():
    """E_{phi phi} - sin^2(theta) E_{theta theta} for the actual E = G - 8 pi T: must be 0."""
    g, ginv = metric()
    Gam = christoffel(g, ginv)
    G, _ = einstein(g, ginv, ricci(Gam))
    _, T_dn, _ = stress_tensor(g, ginv)
    E = G - 8 * sp.pi * T_dn
    return sp.cancel(E[3, 3] - sp.sin(th) ** 2 * E[2, 2])


EXPECTED_COEFFICIENT = -2 * r / a**2
