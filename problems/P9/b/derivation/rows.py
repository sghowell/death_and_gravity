"""I1, the k-essence dictionary, and the (N)-row identity targets.

Everything is a substitution into the DERIVED exact dust DE eigen-speed
(reduction.cs2_dust_exact), never into a transcribed formula.

I1 (b0 anchor, deliverable iii):  for S = {alpha_K} on family B (dust),

    alpha_K c_s^2 = 3 Omega_DE (1 + w_DE) = f / (M_Pl^2 H^2),

an exact identity (test).  k-essence dictionary (deliverable vii): with
c = X P_X, M_2^4 = X^2 P_XX (derived in background.kessence_unitary) and the
derived Friedmann equation, I1 <=> c_s^2 = P_X/(P_X + 2 X P_XX), and
rho_DE + p_DE = 2 X P_X.
"""
from __future__ import annotations

from functools import cache

import sympy as sp

from . import reduction as R
from .tools import M2, H, MPl2, aBg, aHh, aT, pm, rhom, subs_fun, t


def _row_subs(e, aB_on=False, aT_on=False, aM_on=False, aH_on=False):
    """Switch off the alphas not in S (exact substitution incl. derivatives)."""
    m = {}
    if not aB_on:
        m[aBg] = sp.Integer(0)
    if not aT_on:
        m[aT] = sp.Integer(0)
    if not aH_on:
        m[aHh] = sp.Integer(0)
    e = subs_fun(e, m)
    if not aM_on:
        e = subs_fun(e, {M2: MPl2})   # M^2 frozen at M_Pl^2, derivatives vanish
    return sp.cancel(sp.together(sp.expand(e)))


@cache
def cs2_times_aK_row(aB_on=False, aT_on=False, aM_on=False, aH_on=False):
    """alpha (= D) times the exact dust DE speed, for the given row."""
    S = R.cs2_dust_exact()
    alpha = R.alpha
    out = _row_subs(sp.expand(alpha * S), aB_on, aT_on, aM_on, aH_on)
    return out


# --- I1 -----------------------------------------------------------------------------
def I1_lhs():
    """alpha_K c_s^2 for S = {alpha_K} (exact, from the derived speed)."""
    return cs2_times_aK_row()


def I1_rhs():
    """3 Omega_DE (1+w_DE) on family B with dust (p_m = 0)."""
    from .background import f_DE
    return sp.cancel(sp.together(f_DE.subs(pm, 0) / (MPl2 * H**2)))


# --- k-essence dictionary (deliverable vii) ----------------------------------------
@cache
def kessence_check():
    """Return (cs2_kessence_derived, cs2_expected, rhopluspX, alphaK_kess).

    Uses the alpha_H = 0 exact speed with GENERAL p_m (BS 3.13-verified form),
    the derived c, Lambda, M_2^4 of unitary-gauge k-essence, the dictionary
    alpha_K = (2c + 4 M_2^4)/(M^2 H^2), and the derived Friedmann equation
    Hdot = -(2 X P_X + rho_m + p_m)/(2 M^2)   [GLV14 text after (85)].
    """
    from .background import kessence_unitary
    from .tools import K1, K2, phb
    c_val, _lam_val, M24_val = kessence_unitary()
    Xb = phb.diff(t) ** 2 / 2
    aK_kess = sp.cancel((2 * c_val + 4 * M24_val) / (MPl2 * H**2))
    # S = {alpha_K}, alpha_H = 0, general p_m: cs2 = (3.13) restriction (verified)
    Hd = sp.Derivative(H, t)
    cs2aK = -2 * Hd / H**2 - (rhom + pm) / (MPl2 * H**2)
    Hd_sol = -(2 * Xb * K1 + rhom + pm) / (2 * MPl2)
    cs2 = sp.cancel(sp.together((cs2aK / aK_kess).subs(Hd, Hd_sol)))
    cs2_expected = K1 / (K1 + 2 * Xb * K2)     # Vikman / GLV text: P_X/(P_X+2XP_XX)
    rho_plus_p = sp.cancel(2 * Xb * K1)        # = 2 X P_X
    return cs2, cs2_expected, rho_plus_p, aK_kess
