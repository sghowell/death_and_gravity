"""Background (tadpole) equations and family-B definitions.

Derived facts (each an identity test in tests/):
  * eps^1 (homogeneous) variation of L87 + Lmat gives GLV14 (88)-(89) exactly
    [known-answer], plus the matter equation of motion d/dt(a^3 sigma' P_Y)=0.
  * matter eom  =>  sigma'' = -3 H sigma' c_m^2  and energy conservation
    rho_m' + 3H(rho_m + p_m) = 0.
  * family B (FORMULATION 1.3): rho_DE = 3 M_Pl^2 H^2 - rho_m,
    p_DE = -M_Pl^2 (2 Hdot + 3 H^2) - p_m, f = rho_DE + p_DE; DE conservation
    follows identically.
  * unitary-gauge k-essence P(phi, X): c = X P_X, Lambda = X P_X - P,
    M_2^4 = X^2 P_XX  (derived by jet expansion, feeds the (vii) dictionary).
"""
from __future__ import annotations

from functools import cache

import sympy as sp

from . import actions as A
from .tools import P0, P1, P2, H, MPl2, bgsubs, cm2, ddt, eps, pm, rhom, sgb, t

Yb = sgb.diff(t) ** 2 / 2
rhom_of_P = 2 * Yb * P1 - P0          # rho_m = 2 Y P_Y - P
pm_of_P = P0                          # p_m = P
cm2_of_P = P1 / (P1 + 2 * Yb * P2)    # matter sound speed
sig_dd_sol = -3 * H * sgb.diff(t) * cm2   # from the derived matter eom

FLUID_SUBS = {P0: pm, P1: (rhom + pm) / (2 * Yb),
              P2: (rhom + pm) * (1 - cm2) / (4 * Yb**2 * cm2)}


def euler_lagrange(L, v):
    """EL operator for L(v, v', v'') with background rules applied."""
    e = sp.expand(L)
    r = sp.diff(e, v) - ddt(e.coeff(sp.Derivative(v, t)))
    c2 = e.coeff(sp.Derivative(v, (t, 2)))
    if c2 != 0:
        r += ddt(ddt(c2))
    return sp.expand(bgsubs(r))


@cache
def tadpole_equations():
    """(E_dN, E_zeta, E_dsigma) of the homogeneous eps^1 action L87 + Lmat."""
    from .tools import dnf, dsf, zf
    L1 = sp.expand(A.tadpole(A.L87_jet()) + A.tadpole(A.Lmat_jet()))
    return (euler_lagrange(L1, dnf), euler_lagrange(L1, zf),
            euler_lagrange(L1, dsf))


def fluid(e):
    """Translate P-Taylor background data to (rho_m, p_m, c_m^2)."""
    e = e.subs(sp.Derivative(sgb, (t, 2)), sig_dd_sol)
    e = e.subs(FLUID_SUBS)
    return sp.cancel(sp.together(sp.expand(e)))


# --- family B (frozen definitions, FORMULATION 1.3; dust has p_m = 0) --------------
rho_DE = 3 * MPl2 * H**2 - rhom
p_DE = -MPl2 * (2 * sp.Derivative(H, t) + 3 * H**2) - pm
f_DE = rho_DE + p_DE                  # = -2 MPl2 Hdot - (rho_m + p_m)
Omega_DE = rho_DE / (3 * MPl2 * H**2)
one_plus_w_DE = f_DE / rho_DE


@cache
def kessence_unitary():
    """c(t), Lambda(t), M_2^4(t) of unitary-gauge k-essence P(phi(t), X).

    Derived by expanding P(X), X = -(1/2) g^{munu} d_mu phi d_nu phi with
    phi = phi_b(t), in delta g^00 = 1 - 1/N^2 (exact jets); matched against
    the -Lambda - c g^00 + (M_2^4/2)(dg00)^2 operators of L87.
    """
    from . import geometry as geo
    from .tools import K0, K1, K2, phb
    d = geo.build()
    Xb = phb.diff(t) ** 2 / 2
    Xjet = sp.expand(Xb * d["invN2"])         # X = Xbar / N^2 exactly
    dX = sp.expand(Xjet - Xb)
    Pjet = sp.expand(K0 + K1 * dX + K2 / 2 * dX**2)   # P(X) Taylor (exact, eps^2)
    dg00 = sp.expand(1 - d["invN2"])          # delta g^00
    # match Pjet == -Lambda - c*g00 + M24/2 dg00^2 with g00 = -1 + dg00
    c_val = sp.together(-sp.expand(Pjet).coeff(eps, 1) / sp.expand(dg00).coeff(eps, 1))
    lam_val = sp.together(c_val - K0)
    resid2 = sp.expand(Pjet.coeff(eps, 2) - (-c_val) * sp.expand(dg00).coeff(eps, 2))
    M24_val = sp.together(2 * resid2 / sp.expand(dg00).coeff(eps, 1) ** 2)
    return c_val, lam_val, M24_val
