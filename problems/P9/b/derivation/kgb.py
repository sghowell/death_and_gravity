"""Covariant kinetic-gravity-braiding in unitary gauge: independent route.

L = sqrt(-g) [ MPl^2/2 R4 + K(X) + G(X) box(phi) ],  X = -(1/2) (d phi)^2,
phi = phi_b(t) (unitary gauge), signature (-,+,+,+).  NOTE: DPSV
(arXiv:1008.0048) use signature (+,-,-,-), where box(phi) flips sign; their
model K = -X, G = mu X (their eq. 59) is L_ours = -X - mu X box(phi), i.e.
G_ours = -mu X.  The alpha known-answer (BS14 A.7 restricted) fixes all signs.

Everything here is derived by the same jet machinery: box(phi) from
(1/sqrt(-g)) d_mu (sqrt(-g) g^{mu nu} d_nu phi), tadpoles give the background
(Friedmann + shift-charge conservation, DPSV (40)-(41) as known-answer),
the eps^2 part is matched against the frozen L86 to DERIVE alpha_K, alpha_B.
"""
from __future__ import annotations

from functools import cache

import sympy as sp

from . import actions as A
from . import geometry as geo
from .tools import (
    G0,
    G1,
    G2,
    K0,
    K1,
    K2,
    H,
    MPl2,
    bgsubs,
    cut,
    jetpow,
    mul,
    phb,
    t,
)

Xb = phb.diff(t) ** 2 / 2


@cache
def box_phi():
    """box(phi) as an eps-jet for phi = phi_b(t) on the perturbed ADM metric."""
    d = geo.build()
    sqg = d["sqg"]
    invsqg = jetpow(geo.a**3, cut(sqg - geo.a**3), -1)
    pd = phb.diff(t)
    inner_t = mul(sqg, d["g4i"][0][0], pd)
    inner_x = mul(sqg, d["g4i"][1][0], pd)
    from .tools import x
    return cut(mul(invsqg, cut(sp.diff(inner_t, t) + sp.diff(inner_x, x))))


@cache
def L_kgb_jet():
    """sqrt(-g) [MPl^2/2 R4 + K(X) + G(X) box(phi)] as an eps-jet."""
    d = geo.build()
    R4 = geo.ricci4()
    Xjet = sp.expand(Xb * d["invN2"])       # X = Xbar/N^2 exactly (unitary gauge)
    dX = cut(Xjet - Xb)
    Kjet = cut(K0 + K1 * dX + K2 / 2 * mul(dX, dX))
    Gjet = cut(G0 + G1 * dX + G2 / 2 * mul(dX, dX))
    return cut(mul(d["sqg"], cut(MPl2 / 2 * R4 + Kjet + mul(Gjet, box_phi()))))


@cache
def background():
    """Tadpole equations of L_kgb (+ dust handled separately in tests).

    Returns (E_dN, E_zeta, rho_kgb, p_kgb, J):
      E_dN = 0   <=>  3 MPl^2 H^2 = rho_kgb          (Friedmann, vacuum part)
      E_z  = 0   <=>  MPl^2(2 Hdot + 3 H^2) = -p_kgb
      J = phidot (K_X - 3 H phidot G_X): Bianchi gives the exact identity
      rhodot + 3H(rho+p) = +phidot (Jdot + 3HJ)  [sign fixed by the test],
      i.e. DPSV (40)-(41) under the signature map G_DPSV = -G_ours.
    """
    from .background import euler_lagrange
    from .tools import dnf, zf
    L1 = A.tadpole(L_kgb_jet())
    E_dn = euler_lagrange(L1, dnf)
    E_z = euler_lagrange(L1, zf)
    a = geo.a
    rho = sp.expand(sp.cancel(3 * MPl2 * H**2 - E_dn / a**3))
    p_kgb = sp.expand(sp.cancel(E_z / (3 * a**3)
                                - MPl2 * (2 * sp.Derivative(H, t) + 3 * H**2)))
    J = phb.diff(t) * (K1 - 3 * H * phb.diff(t) * G1)   # DPSV (41) under G -> -G
    return E_dn, E_z, rho, p_kgb, J


def quad():
    """eps^2 x-averaged quadratic Lagrangian of L_kgb (gravity+scalar only)."""
    return sp.expand(bgsubs(A.quad(L_kgb_jet())))


# expected alphas: BS14 (A.7) restricted to K(X), G3 = -G(X) (KYY convention),
# G4 = MPl^2/2 — transcribed as KNOWN-ANSWER expected values, tested not input.
def expected_alphas():
    pd = phb.diff(t)
    M2e = MPl2
    aK_e = (2 * Xb * (K1 + 2 * Xb * K2) - 12 * pd * Xb * H * (G1 + Xb * G2)) / (M2e * H**2)
    aB_BS_e = -2 * pd * Xb * G1 / (H * M2e)
    return M2e, sp.expand(aK_e), sp.expand(aB_BS_e)
