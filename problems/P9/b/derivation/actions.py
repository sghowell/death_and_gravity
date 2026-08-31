"""The frozen action and its GPV-form equivalent, as exact eps-jets.

L86  : transcription of GLV14 eq. (86), the FROZEN five-alpha quadratic action
       (alpha_B in GLV convention; alpha_H defined by their eq. (82); alpha_M
       enters only through dM^2/dt = alpha_M H M^2 downstream).  Fetched
       2026-08-31 from arXiv:1411.3712 (ar5iv), transcribed term by term.
L87  : GLV14 eq. (87) unitary-gauge operator action (GPV operator basis;
       1210.0201's own names: m33 = mbar_1^3, m42 = Mbar_2^2/2 = -Mbar_3^2/2,
       mt42 = 2 mu_1^2).  Used to DERIVE background equations + dictionary.
Lmat : minimally coupled k-essence matter P(Y), modelling dust as c_m^2 -> 0.

Every jet is exact through O(eps^2); no background equation is imposed here.
"""
from __future__ import annotations

from functools import cache

import sympy as sp

from . import geometry as geo
from .tools import (
    M2,
    M24,
    P0,
    P1,
    P2,
    H,
    Lam,
    Ms2,
    aBg,
    aHh,
    aK,
    aT,
    bgsubs,
    cc,
    cut,
    dsf,
    eps,
    f_,
    k,
    m33,
    m42,
    mt42,
    mul,
    sgb,
    xavg,
)


def _e1(e):
    return sp.expand(e).coeff(eps, 1)


def _e2(e):
    return sp.expand(e).coeff(eps, 2)


@cache
def L86():
    """Frozen GLV14 (86) integrand (eps^2 coefficient, x-averaged)."""
    d = geo.build()
    a = geo.a
    dK1ud = [[_e1(d["Kud"][i][j]) for j in range(3)] for i in range(3)]
    # eps^1 of K^i_j already has no background piece (K̄ = H delta at eps^0)
    dK1 = sum(dK1ud[i][i] for i in range(3))
    dKdK = sum(dK1ud[i][j] * dK1ud[j][i] for i in range(3) for j in range(3))
    R1 = _e1(d["R3"])
    R2 = _e2(d["R3"])
    z1 = _e1(d["zeta"])            # first-order zeta (with its cos factor)
    dN1 = _e1(d["dN"])
    dsqh1_over_a3 = 3 * z1        # GLV (74) delta sqrt(h) = 3 a^3 zeta
    integrand = a**3 * M2 / 2 * (
        dKdK - dK1**2
        + (1 + aT) * (R1 * dsqh1_over_a3 + R2)
        + aK * H**2 * dN1**2
        + 4 * aBg * H * dK1 * dN1
        + (1 + aHh) * R1 * dN1)
    return xavg(sp.expand(integrand))


@cache
def L87_jet():
    """GLV14 (87) integrand as a full eps-jet (eps^0, eps^1, eps^2 kept)."""
    d = geo.build()
    R4 = geo.ricci4()
    g00 = -d["invN2"]              # upper 00 component of the ADM inverse
    dg00 = cut(g00 + 1)
    dK = cut(bgsubs(d["trK"]) - 3 * H)          # eps^0 part vanishes exactly
    dKud = [[cut(bgsubs(d["Kud"][i][j]) - (H if i == j else 0)) for j in range(3)]
            for i in range(3)]
    dK2m = cut(mul(dK, dK) - sum(mul(dKud[i][j], dKud[j][i])
                                 for i in range(3) for j in range(3)))
    sqg = d["sqg"]
    L = cut(mul(sqg, cut(Ms2 / 2 * f_ * R4 - Lam - cc * g00))
            + mul(sqg, cut(M24 / 2 * mul(dg00, dg00) - m33 / 2 * mul(dK, dg00)
                           - m42 * dK2m + mt42 / 2 * mul(d["R3"], dg00))))
    return L


@cache
def Lmat_jet():
    """Matter integrand sqrt(-g) P(Y) as a full eps-jet (k-essence P)."""
    d = geo.build()
    Y = geo.matter_Y(sgb, dsf)
    Y0 = sp.expand(Y).coeff(eps, 0)
    dY = cut(Y - Y0)
    Pjet = cut(P0 + P1 * dY + P2 / 2 * mul(dY, dY))
    return cut(mul(d["sqg"], Pjet))


def tadpole(Ljet):
    """eps^1 part at the homogeneous (k = 0) mode."""
    return sp.expand(Ljet).coeff(eps, 1).subs(k, 0)


def quad(Ljet):
    """eps^2 part, x-averaged (k kept general)."""
    return xavg(sp.expand(Ljet).coeff(eps, 2))
