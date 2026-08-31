"""DPSV kinetic-gravity-braiding stable phantom crossing (known-answer check).

Model (DPSV arXiv:1008.0048 eq. (59), their signature (+,-,-,-)): K = -X,
G = mu X.  In our frozen conventions (signature (-,+,+,+), L = K + G box phi)
this is K(X) = -X, G(X) = -mu X (map derived in derivation/kgb.py; DPSV
(40)-(41) shift-charge conservation verified there as an exact identity).

Everything numerical here uses the DERIVED quantities:
  rho, p, J        from derivation.kgb.background()
  alpha_K, alpha_B from the exact quadratic-action match (test_derivation)
  D, c_s^2         from the derived reduction (BS14 (3.13) form, alpha_H=0,
                   itself a verified known-answer of the chain).
Units: M_Pl^2 = 1, mu = 1.  This is a known-answer test of the machinery
(digest D2), NOT a witness for the ladder.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

MU = 1.0


# --- derived background quantities for K = -X, G = -mu X (MPl2 = 1) ----------------
def rho_phi(pd, h):
    x = pd * pd / 2
    return -x + 3 * MU * h * pd**3


def p_phi(pd, h, pdd):
    x = pd * pd / 2
    return -x - MU * pd * pd * pdd


def shift_charge(pd, h):
    return pd * (-1.0 + 3 * MU * h * pd)


def solve_hdot_pdd(pd, h, rhom):
    """Linear solve of {pressure equation, shift-charge conservation}.

    2 Hdot = -(rho_phi + p_phi + rho_m);  d/dt J + 3 H J = 0 with
    J = J(pd, h):  J_pd * pdd + J_h * Hdot + 3 H J = 0.
    """
    x = pd * pd / 2
    # eq1: 2 Hdot - mu pd^2 pdd = 2x - 3 mu h pd^3 - rhom
    a11, a12 = 2.0, -MU * pd * pd
    b1 = 2 * x - 3 * MU * h * pd**3 - rhom
    # eq2: (3 mu pd^2) Hdot + (-1 + 6 mu h pd) pdd = -3 h J
    a21, a22 = 3 * MU * pd * pd, (-1.0 + 6 * MU * h * pd)
    b2 = -3 * h * shift_charge(pd, h)
    det = a11 * a22 - a12 * a21
    hdot = (b1 * a22 - b2 * a12) / det
    pdd = (a11 * b2 - a21 * b1) / det
    return hdot, pdd


def alphas(pd, h):
    """alpha_K, alpha_B^BS from the derived match (M^2 = 1)."""
    x = pd * pd / 2
    aK = (-2 * x + 12 * MU * pd * x * h) / (h * h)
    aB = 2 * MU * pd * x / h
    return aK, aB


def stability(pd, h, rhom):
    """(D, c_s^2) from the derived reduction (dust, alpha_T = alpha_M = alpha_H = 0)."""
    hdot, pdd = solve_hdot_pdd(pd, h, rhom)
    aK, aB = alphas(pd, h)
    D = aK + 1.5 * aB * aB
    x = pd * pd / 2
    # d(alpha_B)/dt = 2 mu (1.5 pd^2 pdd / h - x pd hdot / h^2)  [chain rule]
    aBdot = 2 * MU * (1.5 * pd * pd * pdd / h - pd * x * hdot / (h * h))
    cs2 = -((2 - aB) * (hdot - 0.5 * h * h * aB) - h * aBdot + rhom) / (h * h * D)
    return D, cs2


def w_de(pd, h, rhom):
    hdot, _ = solve_hdot_pdd(pd, h, rhom)
    rho_de = 3 * h * h - rhom
    p_de = -(2 * hdot + 3 * h * h)
    return p_de / rho_de, rho_de


def run(pd0=0.50, om0=0.90, h0=1.0, n_efolds=8.0, n=4000):
    """Integrate d/dN (pd, rhom, h) in N = ln a; return trajectory arrays."""
    rhom0 = 3 * h0 * h0 * om0
    # consistency: rho_phi(pd0, h0) must equal 3 h0^2 (1 - om0): solve h0 instead
    from scipy.optimize import brentq
    f = lambda h: rho_phi(pd0, h) + rhom0 - 3 * h * h
    h0 = brentq(f, 1e-3, 1e3)

    def rhs(N, y):
        pd, rhom, h = y
        hdot, pdd = solve_hdot_pdd(pd, h, rhom)
        return [pdd / h, -3 * rhom / 1.0, hdot / h]

    sol = solve_ivp(rhs, (0, n_efolds), [pd0, rhom0, h0], rtol=1e-10, atol=1e-12,
                    dense_output=True, max_step=0.01)
    N = np.linspace(0, n_efolds, n)
    pd, rhom, h = sol.sol(N)
    W = np.array([w_de(*v)[0] for v in zip(pd, h, rhom)])
    DD, CS2 = np.transpose([stability(*v) for v in zip(pd, h, rhom)])
    return dict(N=N, pd=pd, rhom=rhom, h=h, w=W, D=DD, cs2=CS2, sol=sol)


def crossing_report(res):
    """First w = -1 crossing with D > 0, c_s^2 > 0 in a bracket around it."""
    N, W, D, CS2 = res["N"], res["w"], res["D"], res["cs2"]
    rho_de = np.array([w_de(*v)[1] for v in zip(res["pd"], res["h"], res["rhom"])])
    s = np.sign(1 + W)
    idx = np.where((s[:-1] > 0) & (s[1:] < 0) & (rho_de[:-1] > 0))[0]
    if len(idx) == 0:
        return None
    i = idx[0]
    lo, hi = max(0, i - 200), min(len(N) - 1, i + 200)
    return dict(N_c=float(N[i]), w_before=float(W[i - 5]), w_after=float(W[i + 5]),
                D_min=float(D[lo:hi].min()), cs2_min=float(CS2[lo:hi].min()),
                rho_de_c=float(rho_de[i]), window=(float(N[lo]), float(N[hi])))


if __name__ == "__main__":
    res = run()
    rep = crossing_report(res)
    print("crossing:", rep)
