"""Exact ADM geometry of the perturbed flat-FLRW scalar sector (single mode).

Everything is derived from the metric, no transcription: 3d Christoffels and
Ricci by plain index loops (P4 house style), extrinsic curvature from
K_ij = (h_ij' - D_i N_j - D_j N_i)/(2N), 4d curvature from the full ADM
metric by 4d index loops.  GLV14 (73)-(75) are verified downstream as
known-answer tests of these jets.

All returned quantities are eps-jets (exact Taylor polynomials, O(eps^2)).
"""
from __future__ import annotations

from functools import cache

import sympy as sp

from .tools import (
    COORDS3,
    COORDS4,
    a,
    cut,
    cx,
    dnf,
    eps,
    mul,
    psf,
    t,
    zf,
)


def _zeta():
    return eps * zf * cx


def expjet(c):
    """e^{c zeta} as an eps-jet (exact Taylor)."""
    z = _zeta()
    return cut(1 + c * z + c**2 * z**2 / 2)


@cache
def build():
    """All geometric jets for the scalar sector.  Cached (expensive)."""
    z = _zeta()
    dN = eps * dnf * cx
    N = 1 + dN
    invN = cut(1 - dN + dN**2)          # 1/N jet
    invN2 = cut(1 - 2 * dN + 3 * dN**2)  # 1/N^2 jet
    psi = eps * psf * cx

    # 3-metric h_ij = a^2 e^{2 zeta} delta_ij and inverse
    E2, Em2 = expjet(2), expjet(-2)
    h = [[a**2 * E2 if i == j else sp.Integer(0) for j in range(3)] for i in range(3)]
    hinv = [[Em2 / a**2 if i == j else sp.Integer(0) for j in range(3)] for i in range(3)]
    sqh = a**3 * expjet(3)

    # 3d Christoffels, Ricci (index loops; fields depend on (t, x) only)
    G3 = [[[cut(sp.Rational(1, 2) * sum(
        hinv[l][s] * (sp.diff(h[s][j], COORDS3[i]) + sp.diff(h[s][i], COORDS3[j])
                      - sp.diff(h[i][j], COORDS3[s])) for s in range(3)))
        for j in range(3)] for i in range(3)] for l in range(3)]
    Ric3 = [[cut(sum(sp.diff(G3[l][i][j], COORDS3[l]) - sp.diff(G3[l][i][l], COORDS3[j])
                     + sum(G3[l][l][s] * G3[s][i][j] - G3[l][j][s] * G3[s][i][l]
                           for s in range(3)) for l in range(3)))
             for j in range(3)] for i in range(3)]
    R3 = cut(sum(hinv[i][j] * Ric3[i][j] for i in range(3) for j in range(3)))

    # shift: N^i = delta^{ij} d_j psi (GLV14 eq. 73), N_i = h_ij N^j
    Nup = [sp.diff(psi, c) for c in COORDS3]
    Ndown = [cut(sum(h[i][j] * Nup[j] for j in range(3))) for i in range(3)]

    # extrinsic curvature
    DN = [[cut(sp.diff(Ndown[j], COORDS3[i])
               - sum(G3[l][i][j] * Ndown[l] for l in range(3)))
           for j in range(3)] for i in range(3)]
    Kdd = [[mul(invN, sp.Rational(1, 2) * cut(sp.diff(h[i][j], t) - DN[i][j] - DN[j][i]))
            for j in range(3)] for i in range(3)]
    Kud = [[cut(sum(mul(hinv[i][l], Kdd[l][j]) for l in range(3)))
            for j in range(3)] for i in range(3)]
    trK = cut(sum(Kud[i][i] for i in range(3)))

    # 4d ADM metric and closed-form inverse
    g4 = [[sp.Integer(0)] * 4 for _ in range(4)]
    g4[0][0] = cut(-N**2 + sum(Ndown[i] * Nup[i] for i in range(3)))
    for i in range(3):
        g4[0][i + 1] = g4[i + 1][0] = Ndown[i]
        for j in range(3):
            g4[i + 1][j + 1] = h[i][j]
    g4i = [[sp.Integer(0)] * 4 for _ in range(4)]
    g4i[0][0] = -invN2
    for i in range(3):
        g4i[0][i + 1] = g4i[i + 1][0] = mul(invN2, Nup[i])
        for j in range(3):
            g4i[i + 1][j + 1] = cut(hinv[i][j] - mul(invN2, Nup[i], Nup[j]))
    sqg = mul(N, sqh)

    return dict(N=N, invN=invN, invN2=invN2, sqh=sqh, sqg=sqg, h=h, hinv=hinv,
                R3=R3, Ric3=Ric3, Kdd=Kdd, Kud=Kud, trK=trK, g4=g4, g4i=g4i,
                Nup=Nup, Ndown=Ndown, zeta=z, dN=dN, psi=psi)


@cache
def ricci4():
    """4d Ricci scalar of the ADM metric as an eps-jet (index loops)."""
    d = build()
    g4, g4i = d["g4"], d["g4i"]
    G4 = [[[cut(sp.Rational(1, 2) * sum(
        g4i[l][s] * (sp.diff(g4[s][n], COORDS4[m]) + sp.diff(g4[s][m], COORDS4[n])
                     - sp.diff(g4[m][n], COORDS4[s])) for s in range(4)))
        for n in range(4)] for m in range(4)] for l in range(4)]
    Ric = [[cut(sum(sp.diff(G4[l][m][n], COORDS4[l]) - sp.diff(G4[l][m][l], COORDS4[n])
                    + sum(G4[l][l][s] * G4[s][m][n] - G4[l][n][s] * G4[s][m][l]
                          for s in range(4)) for l in range(4)))
            for n in range(4)] for m in range(4)]
    return cut(sum(g4i[m][n] * Ric[m][n] for m in range(4) for n in range(4)))


def matter_Y(field_bg, dfield):
    """Y = -(1/2) g^{mu nu} d_mu s d_nu s for s = field_bg(t) + eps*dfield(t)*cos(kx)."""
    d = build()
    s = field_bg + eps * dfield * cx
    ds = [sp.diff(s, c) for c in COORDS3]
    u = cut(sp.diff(s, t) - sum(d["Nup"][i] * ds[i] for i in range(3)))
    return cut(sp.Rational(1, 2) * (mul(d["invN2"], u, u)
               - sum(mul(d["hinv"][i][j], ds[i], ds[j]) for i in range(3) for j in range(3))))
