"""GPV operator dictionary: L87 == L86 under the alpha-map (deliverable v).

The pinned table (GLV14, fetched 2026-08-31, table between their (89), (90)):

    M^2        = Ms^2 f + 2 m4^2
    M^2 H a_M  = d/dt M^2                      (definition, GLV (68) = BS (3.4))
    M^2 H^2 a_K = 2 c + 4 M_2^4
    M^2 H a_B^GLV = (Ms^2 fdot - m_3^3)/2
    M^2 a_T    = -2 m4^2
    M^2 a_H    = 2 (mtilde_4^2 - m4^2)

Verification: substitute the dictionary + the DERIVED background solutions for
c(t), Lambda(t) into the eps^2 parts, and show L87_quad - L86_quad is a total
time derivative (IBP canonicaliser with explicit certificate).  1210.0201's
own operator names: m_3^3 = mbar_1^3, m4^2 = Mbar_2^2/2 = -Mbar_3^2/2,
mtilde_4^2 = 2 mu_1^2 (naming follows GLV14 (87) / GLPV 1304.4840).
"""
from __future__ import annotations

from functools import cache

import sympy as sp

from . import actions as A
from . import ibp
from .tools import (
    M2,
    M24,
    H,
    Lam,
    Ms2,
    aBg,
    aHh,
    aK,
    aT,
    bgsubs,
    cc,
    f_,
    m33,
    m42,
    mt42,
    pm,
    rhom,
    subs_fun,
    t,
)

# the pinned dictionary (expected values; the identity test proves them)
M2_dict = Ms2 * f_ + 2 * m42
DICT = {
    M2: M2_dict,
    aK: (2 * cc + 4 * M24) / (M2_dict * H**2),
    aBg: (Ms2 * sp.Derivative(f_, t) - m33) / (2 * M2_dict * H),
    aT: -2 * m42 / M2_dict,
    aHh: 2 * (mt42 - m42) / M2_dict,
}

# background solutions for c, Lambda from the derived (88)-(89)
_fd = sp.Derivative(f_, t)
_fdd = sp.Derivative(f_, (t, 2))
_Hd = sp.Derivative(H, t)
_sum88 = 3 * Ms2 * (f_ * H**2 + _fd * H) - rhom            # c + Lambda
_dif89 = Ms2 * (2 * f_ * _Hd + 3 * f_ * H**2 + 2 * _fd * H + _fdd) + pm  # Lambda - c
C_SOL = sp.together((_sum88 - _dif89) / 2)
LAM_SOL = sp.together((_sum88 + _dif89) / 2)


def csub(e):
    """Substitute the solved background c(t), Lambda(t) (and their derivatives)."""
    for fn, val in ((cc, C_SOL), (Lam, LAM_SOL)):
        for o in (2, 1):
            e = e.subs(sp.Derivative(fn, (t, o)), bgsubs(sp.diff(val, (t, o))))
        e = e.subs(fn, val)
    return sp.expand(bgsubs(e))


@cache
def match():
    """(difference, Lc, F): L87_quad - L86_quad|dict with (88)-(89) imposed,
    canonicalised mod d/dt.  Lc equals exactly the matter-coupling remainder
    a^3 (3/2 rho_m zeta dN - 9/4 p_m zeta^2) (x-averaged) — see the test."""
    from .tools import dnf, psf, zf
    L87q = sp.expand(bgsubs(A.quad(A.L87_jet())))
    L86d = subs_fun(sp.expand(bgsubs(A.L86())), DICT)
    diff = csub(L87q - L86d)
    Lc, F = ibp.canon(diff, [zf, dnf, psf])
    return diff, Lc, F


def frow_limit():
    """f(R) limit (deliverable vi): c = M24 = m33 = m42 = mt42 = 0.

    Returns (alphaK, alphaB_BS + alphaM, alphaT, alphaH) — all must be 0.
    """
    lim = {cc: sp.Integer(0), M24: sp.Integer(0), m33: sp.Integer(0),
           m42: sp.Integer(0), mt42: sp.Integer(0)}
    aK_l = DICT[aK].subs(lim)
    aB_BS = -2 * DICT[aBg].subs(lim)               # BS14 convention map
    M2_l = M2_dict.subs(lim)
    aM_l = sp.Derivative(M2_l, t).doit() / (H * M2_l)   # alpha_M = M2dot/(H M2)
    aT_l = DICT[aT].subs(lim)
    aH_l = DICT[aHh].subs(lim)
    return (sp.cancel(aK_l), sp.cancel(sp.together(aB_BS + aM_l)),
            sp.cancel(aT_l), sp.cancel(aH_l))
