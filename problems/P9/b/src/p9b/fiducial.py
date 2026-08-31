"""F* fiducial background and S2 near-witnesses for the no-alpha_B ladder rows.

Fiducial F* (FORMULATION 1.3): CPL w_DE = w0 + wa(1-a), (w0, wa) = (-4/5, -4/5),
Omega_m0 = 3/10, dust + DE only, crossing at a_c = 3/4.  Units: H0 = 1,
M_Pl^2 = 1.  Then rho_DE(a)/rho_DE0 = a^(9/5) exp((12/5)(1-a)) exactly.

The single symbolic expression MASTER_TARGET below encodes alpha_K * c_s^2
(c_s^2 = the exact dust DE eigen-speed pinned in FORMULATION rev. v1.1) for
EVERY no-alpha_B row, in the variables
    s_aT = alpha_T, s_aM = alpha_M, s_aH = alpha_H, s_dH = d(alpha_H)/dx,
    s_hd = Hdot/H^2, s_rm = rho_m/(M^2 H^2);
rows switch off absent alphas (and set M^2 = M_Pl^2, i.e. m2 = 1, when alpha_M
is off).  It is matched EXACTLY against the derived targets
derivation.rows.cs2_times_aK_row(...) for all 8 rows in
tests/test_s2_rows.py (test_master_target_matches_derived_rows); all numeric
evaluation goes through its lambdified form TARGET_FN, so the float
near-witnesses cannot drift from the derivation chain.

Near-witness conventions (S2, float level — S3 owns Arb certification):
alpha_K = 1 (so D = alpha_K = 1 > 0 and c_s^2 = target); alpha_T, alpha_M
constant; alpha_H constant or linear in x = ln a; M^2(a) = a^alpha_M
(normalisation M^2 = M_Pl^2 at t0: a = 1, the right endpoint of every J).
"""
from __future__ import annotations

import numpy as np
import sympy as sp

# --- fiducial constants (exact rationals as floats) --------------------------------
W0, WA, OM0 = -0.8, -0.8, 0.3
A_C = 0.75                       # crossing scale factor (exact: 3/4)

# --- master row target (see module docstring) --------------------------------------
s_aT, s_aM, s_aH, s_dH, s_hd, s_rm = sp.symbols(
    "s_aT s_aM s_aH s_dH s_hd s_rm", real=True)
MASTER_TARGET = (-2 * s_aT + 2 * s_aH + 2 * (1 + s_aH) * s_aM
                 - 2 * (1 + s_aH) * s_hd + 2 * s_dH - (1 + 2 * s_aH) * s_rm)
TARGET_FN = sp.lambdify((s_aT, s_aM, s_aH, s_dH, s_hd, s_rm),
                        MASTER_TARGET, "numpy")


# --- F* background -----------------------------------------------------------------
def wde(a):
    return W0 + WA * (1.0 - a)


def g_de(a):
    """rho_DE(a)/rho_DE0 for CPL(-4/5, -4/5): a^(9/5) e^((12/5)(1-a))."""
    return a ** 1.8 * np.exp(2.4 * (1.0 - a))


def E2(a):
    """H^2/H0^2 = Omega_m0 a^-3 + Omega_DE0 g_de(a)."""
    return OM0 * a ** -3 + (1.0 - OM0) * g_de(a)


def Om(a):
    return OM0 * a ** -3 / E2(a)


def Ode(a):
    return (1.0 - OM0) * g_de(a) / E2(a)


def fhat(a):
    """f/(M_Pl^2 H^2) = 3 Omega_DE (1 + w_DE); sign = sign(4a - 3) on F*."""
    return 3.0 * Ode(a) * (1.0 + wde(a))


def hdot_over_H2(a):
    """Hdot/H^2 = -(3/2)[Omega_m + (1 + w_DE) Omega_DE] (dust + DE)."""
    return -1.5 * (Om(a) + (1.0 + wde(a)) * Ode(a))


# --- near-witnesses per row --------------------------------------------------------
# Each entry: J = (a1, a2); constants aT, aM; alpha_H = aH0 + aHs * (x - ln a2)
# (aHs = d alpha_H/dx); alpha_K = 1 everywhere.  Rows absent -> zeros.
WITNESSES = {
    "KT":   dict(J=(0.60, 1.0), aT=-0.10),
    "KM":   dict(J=(0.60, 1.0), aM=0.50),
    "KH":   dict(J=(0.70, 1.0), aH0=1.0),                 # constant alpha_H
    "KH2":  dict(J=(0.60, 1.0), aH0=np.log(1 / 0.60), aHs=1.0),  # slope variant
    "KTM":  dict(J=(0.60, 1.0), aT=-0.10, aM=0.50),
    "KTH":  dict(J=(0.70, 1.0), aT=-0.10, aH0=1.0),
    "KMH":  dict(J=(0.70, 1.0), aM=0.50, aH0=1.0),
    "KTMH": dict(J=(0.70, 1.0), aT=-0.10, aM=0.50, aH0=1.0),
}


def evaluate(row, n=4001):
    """Evaluate the row's near-witness on its J-grid; return a report dict.

    target = alpha_K c_s^2 (alpha_K = 1), through TARGET_FN only.
    """
    w = WITNESSES[row]
    a1, a2 = w["J"]
    a = np.linspace(a1, a2, n)
    x = np.log(a)
    aT = w.get("aT", 0.0) + 0 * a
    aM = w.get("aM", 0.0)
    aH = w.get("aH0", 0.0) + w.get("aHs", 0.0) * (x - np.log(a2))
    daHdx = w.get("aHs", 0.0) + 0 * a
    m2 = a ** aM                          # M^2/M_Pl^2 (== 1 when alpha_M off)
    rm = 3.0 * Om(a) / m2                 # rho_m/(M^2 H^2)
    target = TARGET_FN(aT, aM + 0 * a, aH, daHdx, hdot_over_H2(a), rm)
    # transversality: d fhat/dx at a_c (analytic in a, finite grid step in x)
    ac = np.array([A_C - 1e-6, A_C + 1e-6])
    dfhat_dx = np.diff(fhat(ac))[0] / np.diff(np.log(ac))[0]
    return dict(
        row=row, J=(a1, a2), n=n,
        min_target=float(target.min()), max_target=float(target.max()),
        argmin_a=float(a[target.argmin()]),
        min_cs2=float(target.min()),      # alpha_K = 1
        D=1.0,                            # D = alpha_K (alpha_B = 0)
        alpha_K=1.0,
        max_abs_aT=float(np.abs(aT).max()), max_abs_aM=abs(aM),
        max_abs_aH=float(np.abs(aH).max()),
        min_cT2=float(1.0 + aT.min()),
        min_m2=float(m2.min()), max_m2=float(m2.max()),
        fhat_left=float(fhat(np.array([0.5 * (a1 + A_C)]))[0]),
        fhat_right=float(fhat(np.array([0.5 * (A_C + a2)]))[0]),
        dfhat_dx_at_crossing=float(dfhat_dx),
        min_Ode=float(Ode(a).min()),
        crossing_interior=bool(a1 < A_C < a2),
    )


def report_all():
    return {row: evaluate(row) for row in WITNESSES}


if __name__ == "__main__":
    for row, r in report_all().items():
        print(f"{row:5s} J={r['J']}  min c_s^2 = {r['min_cs2']:+.4f} "
              f"(at a = {r['argmin_a']:.3f})  fhat: {r['fhat_left']:+.4f} -> "
              f"{r['fhat_right']:+.4f}  cT2_min = {r['min_cT2']:.2f}  "
              f"m2 in [{r['min_m2']:.3f}, {r['max_m2']:.3f}]")
