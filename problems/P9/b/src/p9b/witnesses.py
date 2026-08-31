"""S3 witness registry and the braiding-extended master target (all 15 rows).

Master target (alpha * c_s^2, alpha = D = alpha_K + (3/2)(alpha_B^BS)^2) for
EVERY ladder row S != {alpha_K}, extending fiducial.MASTER_TARGET with the
braiding.  Variables (all scalars; ' = d/dx, x = ln a; BS convention):

    s_aB = alpha_B^BS,  s_dB = (alpha_B^BS)',  s_aT = alpha_T,
    s_aM = alpha_M,     s_aH = alpha_H,        s_dH = alpha_H',
    s_hd = Hdot/H^2,    s_rm = rho_m/(M^2 H^2).

With g = 1 + alpha_B^GLV = 1 - s_aB/2 (map alpha_B^BS = -2 alpha_B^GLV):

    MASTER_TARGET_B = -2 g^2 (1+s_aT) + 2 g (1+s_aH)(1+s_aM-s_hd)
                      + 2 g s_dH + (1+s_aH) s_dB - (1+2 s_aH) s_rm.

This is matched EXACTLY (sympy) against the derived exact dust DE eigen-speed
rows.cs2_times_aK_row(aB_on=..., ...) for all 8 flag combinations in
tests/test_s3_witnesses.py (the FORMULATION rev-v1.1 eigen-speed pin, with
kinetic matter mixing, is inherited from the derivation chain); at s_aB =
s_dB = 0 it reduces exactly to fiducial.MASTER_TARGET.  The lambdified
TARGET_POLY_FN is a pure polynomial (dyadic coefficients only), so it
evaluates exactly both on floats and on python-flint arb balls; all
transcendentals live in the background functions (certify.py).

Witnesses (S3, frozen here; all D3 bounds respected):
alpha_K = 1; alpha_T = -1/10, alpha_M = 1/2, alpha_H = 1 (constants) where
the row switches them on; alpha_B^BS(a) = a/2 (rational in a; magnitude
guided by the S1 KGB check, aB_BS ~ 0.5, and its slope term (aB_BS)' = a/2);
M^2/M_Pl^2 = a^alpha_M (normalisation M^2(a=1) = M_Pl^2).  J = [3/5, 1]
(rows without alpha_H) or [7/10, 1] (rows with alpha_H), both containing the
F* crossing a_c = 3/4 in the interior.
"""
from __future__ import annotations

from fractions import Fraction as F

import sympy as sp

# --- scalar symbols and the master target ------------------------------------------
s_aB, s_dB, s_aT, s_aM, s_aH, s_dH, s_hd, s_rm = sp.symbols(
    "s_aB s_dB s_aT s_aM s_aH s_dH s_hd s_rm", real=True)
SCALARS = (s_aB, s_dB, s_aT, s_aM, s_aH, s_dH, s_hd, s_rm)

_g = 1 - s_aB / 2
MASTER_TARGET_B = (-2 * _g**2 * (1 + s_aT) + 2 * _g * (1 + s_aH) * (1 + s_aM - s_hd)
                   + 2 * _g * s_dH + (1 + s_aH) * s_dB - (1 + 2 * s_aH) * s_rm)

# pure dyadic polynomial: works on floats and on arb balls alike
TARGET_POLY_FN = sp.lambdify(SCALARS, MASTER_TARGET_B, "math")

A_C = F(3, 4)          # F* crossing scale factor (exact)
GAP = F(1, 50)         # IVT gap half-width around a_c for the fhat sign sweeps
ALPHA_K = F(1)         # alpha_K = 1 in every witness


def _row(J1, b1=None, aT=None, aM=None, aH=None, c_min=None):
    return dict(J=(J1, F(1)), b1=b1, aT=aT, aM=aM, aH=aH,
                c_min=c_min, D_min=F(11, 10) if b1 is not None else F(1))


# c_min: certified strict lower bound for the exact eigen-speed c_s^2 on J
# (rational, frozen below the float minimum found by scan(); Arb-certified in
# certify.py).  D_min: 1 exact for alpha_B = 0 rows (D = alpha_K = 1);
# 11/10 Arb-certified for braiding rows (D = 1 + (3/2)(a/2)^2 >= 1.135 on J).
REG = {
    "KB":    _row(F(3, 5),  b1=F(1, 2),                                  c_min=F(1, 10)),
    "KT":    _row(F(3, 5),              aT=F(-1, 10),                    c_min=F(1, 20)),
    "KM":    _row(F(3, 5),                            aM=F(1, 2),        c_min=F(1, 4)),
    "KH":    _row(F(7, 10),                                      aH=F(1), c_min=F(1, 4)),
    "KBT":   _row(F(3, 5),  b1=F(1, 2), aT=F(-1, 10),                    c_min=F(1, 5)),
    "KBM":   _row(F(3, 5),  b1=F(1, 2),               aM=F(1, 2),        c_min=F(7, 20)),
    "KBH":   _row(F(7, 10), b1=F(1, 2),                          aH=F(1), c_min=F(3, 10)),
    "KTM":   _row(F(3, 5),              aT=F(-1, 10), aM=F(1, 2),        c_min=F(2, 5)),
    "KTH":   _row(F(7, 10),             aT=F(-1, 10),            aH=F(1), c_min=F(2, 5)),
    "KMH":   _row(F(7, 10),                           aM=F(1, 2), aH=F(1), c_min=F(5, 4)),
    "KBTM":  _row(F(3, 5),  b1=F(1, 2), aT=F(-1, 10), aM=F(1, 2),        c_min=F(2, 5)),
    "KBTH":  _row(F(7, 10), b1=F(1, 2), aT=F(-1, 10),            aH=F(1), c_min=F(2, 5)),
    "KBMH":  _row(F(7, 10), b1=F(1, 2),               aM=F(1, 2), aH=F(1), c_min=F(4, 5)),
    "KTMH":  _row(F(7, 10),             aT=F(-1, 10), aM=F(1, 2), aH=F(1), c_min=F(7, 5)),
    "KBTMH": _row(F(7, 10), b1=F(1, 2), aT=F(-1, 10), aM=F(1, 2), aH=F(1), c_min=F(1)),
}

MINIMAL_ROWS = ("KB", "KT", "KM", "KH")     # the two-operator rows


def flags(row):
    """(B_on, T_on, M_on, H_on) from the row name (canonical K[B][T][M][H])."""
    return tuple(ch in row for ch in "BTMH")


def witness_desc(row):
    """JSON-ready exact description of the frozen witness."""
    w = REG[row]
    return {
        "alpha_K": "1",
        "alpha_B_BS": "a/2" if w["b1"] is not None else "0",
        "alpha_T": str(w["aT"]) if w["aT"] is not None else "0",
        "alpha_M": str(w["aM"]) if w["aM"] is not None else "0",
        "alpha_H": str(w["aH"]) if w["aH"] is not None else "0",
        "M2_over_MPl2": f"a^({w['aM']})" if w["aM"] is not None else "1",
        "J": [str(w["J"][0]), str(w["J"][1])],
    }


# --- float scan (proposer / regression; NOT the certificate) -----------------------
def scan(row, n=4001):
    """Float grid evaluation of c_s^2 = target/D on J (numpy; mirrors S2's
    fiducial.evaluate but covers the braiding rows).  Returns a report dict."""
    import numpy as np

    from . import fiducial as fid

    w = REG[row]
    a1, a2 = float(w["J"][0]), float(w["J"][1])
    a = np.linspace(a1, a2, n)
    b1 = float(w["b1"]) if w["b1"] is not None else 0.0
    aT = float(w["aT"]) if w["aT"] is not None else 0.0
    aM = float(w["aM"]) if w["aM"] is not None else 0.0
    aH = float(w["aH"]) if w["aH"] is not None else 0.0
    b = b1 * a
    db = b1 * a                    # d(b1*a)/dx = b1 * a
    m2 = a**aM
    rm = 3.0 * fid.Om(a) / m2
    hd = fid.hdot_over_H2(a)
    tgt = np.array([TARGET_POLY_FN(b[i], db[i], aT, aM, aH, 0.0, hd[i], rm[i])
                    for i in range(len(a))])
    D = 1.0 + 1.5 * b**2
    cs2 = tgt / D
    return dict(
        row=row, J=(a1, a2), n=n,
        min_cs2=float(cs2.min()), argmin_a=float(a[cs2.argmin()]),
        min_D=float(D.min()), max_D=float(D.max()),
        max_abs_aB=float(np.abs(b).max()),
        min_cT2=1.0 + aT, min_m2=float(m2.min()),
        fhat_left=float(fid.fhat(np.array([0.5 * (a1 + float(A_C))]))[0]),
        fhat_right=float(fid.fhat(np.array([0.5 * (float(A_C) + a2)]))[0]),
        min_Ode=float(fid.Ode(a).min()),
        crossing_interior=bool(a1 < float(A_C) < a2),
        c_min=float(w["c_min"]), D_min=float(w["D_min"]),
    )


if __name__ == "__main__":
    for r in REG:
        s = scan(r)
        print(f"{r:6s} J=({s['J'][0]:.2f},1)  min cs2 = {s['min_cs2']:+.4f}  "
              f"c_min = {s['c_min']:.3f}  D in [{s['min_D']:.3f},{s['max_D']:.3f}]")
