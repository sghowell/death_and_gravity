"""Constraint elimination -> single-variable action, D and c_s^2 (deliverable ii).

Total quadratic action of the frozen theory with matter (DERIVED, not assumed):
the GPV-route total [L87 + Lmat]_eps^2 with the derived (88)-(89) substituted
equals, modulo an exact d/dt (IBP certificate in tests),

    L86 (frozen five-alpha action) + [sqrt(-g) P(Y)]_eps^2
      + a^3 [3 rho_m zeta dN - (9/2) p_m zeta^2]   (pre-average normalisation),

i.e. the x-averaged extra term is (3/2) a^3 rho_m z dN - (9/4) a^3 p_m z^2.
The remainder is the background-equation-proportional piece GLV absorbed when
deriving (86); it vanishes in vacuum, where (86) alone was verified exactly
against (77), (79)-(80), (83).  Reduction:
solve the dN and psi constraints exactly, canonicalise by IBP (explicit
certificate), extract kinetic/gradient matrices, UV speeds at k -> oo.
"""
from __future__ import annotations

from functools import cache

import sympy as sp

from . import actions as A
from . import ibp
from .background import fluid
from .tools import (
    M2,
    H,
    aBg,
    aHh,
    aK,
    aT,
    bgsubs,
    cm2,
    dnf,
    dsf,
    k,
    pm,
    psf,
    rhom,
    t,
    zf,
)

c2 = sp.Symbol("c2", real=True)     # c_s^2 placeholder in dispersion relations
alpha = aK + 6 * aBg**2             # GLV (80); equals BS14 D under the map


@cache
def total_quad():
    """L86 + Lmat_eps2 + matter-coupling remainder (x-averaged, derived)."""
    from . import geometry as geo
    from .background import pm_of_P, rhom_of_P
    Lm2 = sp.expand(bgsubs(A.quad(A.Lmat_jet())))
    extra = (sp.Rational(3, 2) * geo.a**3 * rhom_of_P * zf * dnf
             - sp.Rational(9, 4) * geo.a**3 * pm_of_P * zf**2)
    return sp.expand(bgsubs(A.L86()) + Lm2 + extra)


def solve_constraints(L2):
    """Exact elimination of dN and psi (both algebraic)."""
    E_dn = sp.expand(sp.diff(L2, dnf))
    E_ps = sp.expand(sp.diff(L2, psf))
    A11, A12 = E_dn.coeff(dnf), E_dn.coeff(psf)
    A21, A22 = E_ps.coeff(dnf), E_ps.coeff(psf)
    b1 = -sp.expand(E_dn - A11 * dnf - A12 * psf)
    b2 = -sp.expand(E_ps - A21 * dnf - A22 * psf)
    det = sp.cancel(A11 * A22 - A12 * A21)
    dn_sol = sp.cancel(sp.together((b1 * A22 - b2 * A12) / det))
    ps_sol = sp.cancel(sp.together((A11 * b2 - A21 * b1) / det))
    return dn_sol, ps_sol


@cache
def reduced():
    """(L_canonical, certificate_ok, dn_sol, ps_sol) for the two-field system."""
    L2 = total_quad()
    dn_sol, ps_sol = solve_constraints(L2)
    Lred = sp.expand(L2.subs([(dnf, dn_sol), (psf, ps_sol)]))
    Lc, F = ibp.canon(Lred, [zf, dsf])
    ok = ibp.check_canon(Lred, Lc, F)
    return Lc, ok, dn_sol, ps_sol


@cache
def uv_matrices(mode="gen"):
    """Leading k->oo kinetic matrix and k^2 potential matrix (fluid variables).

    mode="ah0" sets alpha_H = 0 (and its derivative) before extraction.
    """
    Lc, ok, _, _ = reduced()
    assert ok
    if mode == "ah0":
        Lc = Lc.subs(sp.Derivative(aHh, t), 0).subs(aHh, 0)
    Tk, U, _ = ibp.quad_matrices(Lc, [zf, dsf])
    Tinf = Tk.applyfunc(lambda e: ibp.kleading(fluid(e), k))
    G2 = U.applyfunc(lambda e: ibp.kleading(fluid(sp.cancel(sp.together(e) / k**2)), k))
    return Tinf, G2


def dispersion(mode="gen"):
    """det(c2 * Tinf + a^2 * G2): with L = T vd vd + U v v the frequency-domain
    EOM is (omega^2 T + U) v = 0, so c^2 = omega^2 a^2/k^2 are the roots of
    det(c2 T + a^2 U/k^2) (the antisymmetric W term is O(k^0), subleading)."""
    Tinf, G2 = uv_matrices(mode)
    from . import geometry as geo
    return sp.expand(sp.together((c2 * Tinf + G2 * geo.a**2).det()))


@cache
def sum_product(mode="gen"):
    """(sum, product) of the two UV speeds-squared (exact, fluid variables)."""
    p = sp.Poly(dispersion(mode), c2)
    A2, A1, A0 = [sp.cancel(sp.together(q)) for q in p.all_coeffs()]
    return sp.cancel(sp.together(-A1 / A2)), sp.cancel(sp.together(A0 / A2))


@cache
def cs2_dust_exact():
    """Exact DE eigen-speed for dust (c_m^2 -> 0, p_m = 0), alpha_H free.

    Dust: the product of speeds -> 0 (one exactly-soft dust mode); the DE
    speed is the surviving root = sum|_{c_m^2=0, p_m=0}.  For alpha_H = 0 it
    equals BS14 (3.13) = GLV14 (85); for alpha_H != 0 it differs by an exact
    kinetic-matter-mixing term (see notes/s1-derivation.md).
    """
    S, _ = sum_product("gen")
    return sp.cancel(sp.together(sp.expand(S.subs({cm2: 0, pm: 0}))))


# ---- published expected values (KNOWN-ANSWER, transcribed from fetches) -----------
_Hd = sp.Derivative(H, t)
_aM = sp.Derivative(M2, t) / (H * M2)               # alpha_M (GLV 68 / BS 3.4)

# GLV14 (80) gradient coefficient and (85) sound speed (fetched 2026-08-31)
Lgrad_glv80 = 2 * M2 * (1 + aT - (1 + aHh) / (1 + aBg) * (1 + _aM - _Hd / H**2)
                        - sp.diff((1 + aHh) / (1 + aBg), t) / H)
cs2_glv85 = (-(1 + aBg)**2 / alpha * Lgrad_glv80 / M2
             - (1 + aHh)**2 / alpha * (rhom + pm) / (M2 * H**2))

# BS14 (3.13) with alpha_B^BS = -2 alpha_B^GLV (the convention map, deliverable iv)
aB_BS = -2 * aBg
D_BS = aK + sp.Rational(3, 2) * aB_BS**2
QS_BS = 2 * M2 * D_BS / (2 - aB_BS)**2
cs2_bs313 = -((2 - aB_BS) * (_Hd - H**2 * aB_BS * (1 + aT) / 2 - H**2 * (_aM - aT))
              - H * sp.diff(aB_BS, t) + (rhom + pm) / M2) / (H**2 * D_BS)
