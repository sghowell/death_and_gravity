"""Float cross-checks of S2 (reviewer script, 2026-10-17; evidence, not proof).

A. Complex scaling: eigenvalues lambda of H = -d^2/ds^2 + i V_l(r(s)) on the rotated line
   x = x_P + e^{i pi/4} s (FD, Dirichlet on [-60, 30]); lambda = i w^2, so w = sqrt(-i lambda).
   The PS QNMs (S1 float solver values recorded in the certificates) must appear as eigenvalues,
   and no eigenvalue may lie in the window W_l.  This tests Lemma 3's direction and sign conventions
   on the actual contour, independently of the coercivity machinery.
B. Ground-state margin lambda_min(w) of -(c u')' + q u (the quantity the certificate bounds below)
   at the four corners of W_l, along the window edges, at the PS mode, and on the column
   Re w = w_l, Im w = -t gamma0 (t = 0 .. 1).  Compared with mu_cert (A3 generalised: mu_cert must
   not exceed min over W_l of lambda_min).
C. Where the certificate loses: float min over s of E(s; w^2) with the stored ground state, at the
   corners of the w^2 hull (what the certificate must cover) and at the true window corners.
Own contour (scipy DOP853), own potential, own discretisation; parameters from params.json.
Run (repo root):  nice -n 10 .venv/bin/python problems/P2/certificates/review_s2/indep_float_spectrum.py
"""
import json
import time
from fractions import Fraction
from pathlib import Path

import numpy as np
import scipy.sparse as sps
from scipy.integrate import solve_ivp
from scipy.linalg import eigh_tridiagonal
from scipy.sparse.linalg import eigs

HERE = Path(__file__).resolve().parent
P2 = HERE.parents[1]
J = json.loads((P2 / "certificates" / "params.json").read_text())["points"]["PB"]["quantities"]
g = lambda k: float(J[k]["approx"])
rn, rm, rp, rc = g("r_n"), g("r_minus"), g("r_plus"), g("r_c")
Q2, rP, OmP, gam0, km = g("Q2"), g("r_P"), g("Omega_P"), g("gamma0"), g("kappa_minus")
beta = float(J["strip_depth_beta"]["approx"])
E1 = np.exp(1j * np.pi / 4)
T0 = time.time()


def F(r):
    return -(r - rn) * (r - rm) * (r - rp) * (r - rc) / (50 * r * r)


def V(r, l):
    L = l * (l + 1)
    D = -r ** 4 / 50 + r * r - 2 * r + Q2
    Dp = -4 * r ** 3 / 50 + 2 * r - 2
    Fr = D / r ** 2
    Fp = (r * Dp - 2 * D) / r ** 3
    return Fr * (L / r ** 2 + Fp / r)


def contour(s):
    """r(s) on a sorted grid s (both signs)."""
    out = np.empty(s.shape, complex)
    f = lambda t, y: [(E1 * F(y[0] + 1j * y[1])).real, (E1 * F(y[0] + 1j * y[1])).imag]
    for side, S in ((-1, s.min()), (1, s.max())):
        sol = solve_ivp(f, (0.0, S), [rP, 0.0], method="DOP853", rtol=1e-12, atol=1e-14, dense_output=True)
        mask = s <= 0 if side < 0 else s >= 0
        y = sol.sol(s[mask])
        out[mask] = y[0] + 1j * y[1]
    return out


S_LO, S_HI = -60.0, 30.0
H = 0.004
s = np.arange(S_LO, S_HI + H / 2, H)
r = contour(s)
res = {"grid": {"s": [S_LO, S_HI], "h": H, "n": int(s.size)}}

# ------------------------------------------------------------------------------ A. complex scaling
res["A"] = {}
for l, ps in ((20, (100.823703512, -1.68357264593)), (50, (248.396138233, -1.6832715886)),
              (100, (494.340406988, -1.68322053463))):
    Vs = V(r, l)
    n = s.size
    main = 2 / H ** 2 + 1j * Vs
    off = -np.ones(n - 1) / H ** 2
    Hm = sps.diags([off, main, off], [-1, 0, 1], format="csc")
    wl = (l + 0.5) * OmP
    a1, a2, b1 = wl - 8 * gam0, wl + 6 * gam0, -beta
    found = []
    # shifts: the PS fundamental and two overtones, and a sweep across the window
    targets = [complex(ps[0], ps[1] * (2 * k + 1)) * km for k in range(3)]
    targets += [complex(a, b) for a in np.linspace(a1, a2, 7) for b in (b1, b1 / 2, 0.0)]
    for wt in targets:
        lam, _ = eigs(Hm, k=6, sigma=1j * wt * wt, which="LM", tol=1e-12, maxiter=5000)
        for la in lam:
            w = np.sqrt(-1j * la)
            w = w if w.real >= 0 else -w
            found.append(w)
    found = np.unique(np.round(np.array(found), 9))
    in_window = [w for w in found if a1 <= w.real <= a2 and b1 <= w.imag <= 0]
    near = sorted(found, key=lambda w: abs(w - complex(ps[0], ps[1]) * km))[:1]
    ov1 = sorted(found, key=lambda w: abs(w - complex(ps[0], 3 * ps[1]) * km))[:1]
    res["A"][f"l{l}"] = {
        "PS_fundamental_S1_over_kappa_minus": list(ps),
        "nearest_eigen_w_over_kappa_minus": [near[0].real / km, near[0].imag / km],
        "abs_difference_over_kappa_minus": abs(near[0] - complex(ps[0], ps[1]) * km) / km,
        "first_overtone_candidate_over_kappa_minus": [ov1[0].real / km, ov1[0].imag / km],
        "eigenvalues_found": int(len(found)),
        "eigen_w_in_window": [[w.real / km, w.imag / km] for w in in_window],
        "min_dist_of_any_eigen_w_to_window_over_kappa_minus": float(min(
            np.hypot(np.clip(w.real, a1, a2) - w.real, np.clip(w.imag, b1, 0) - w.imag) for w in found) / km)}
    print("A", l, res["A"][f"l{l}"], flush=True)


# ------------------------------------------------------------------------------ B. ground-state margins
def lam_min(l, w, phinf=-1.4375, s1=5.0, eta=0.5):
    Vs = V(r, l)
    ph = phinf * s ** 4 / (s1 ** 4 + s ** 4)
    dph = phinf * 4 * s ** 3 * s1 ** 4 / (s1 ** 4 + s ** 4) ** 2
    c = np.cos(ph) - eta * np.abs(dph)
    q = np.real(np.exp(-1j * ph) * 1j * (Vs - w * w)) - np.abs(dph) / (4 * eta)
    ch = 0.5 * (c[1:] + c[:-1])
    d = np.empty(s.size)
    d[1:-1] = (ch[1:] + ch[:-1]) / H ** 2
    d[0], d[-1] = (ch[0] + c[0]) / H ** 2, (ch[-1] + c[-1]) / H ** 2
    d += q
    e = -ch / H ** 2
    return float(eigh_tridiagonal(d, e, select="i", select_range=(0, 0), eigvals_only=True)[0])


res["B"] = {}
for l, ps in ((20, (100.823703512, -1.68357264593)), (50, (248.396138233, -1.6832715886)),
              (100, (494.340406988, -1.68322053463))):
    st = json.loads((P2 / "certificates" / "m3" / f"m3_PB_l{l}.json").read_text())
    muHO = 2 * (l + 0.5) * OmP * gam0
    wl = (l + 0.5) * OmP
    a1, a2, b1 = wl - 8 * gam0, wl + 6 * gam0, -beta
    corners = {"low_left (a1, -beta)": complex(a1, b1), "low_right (a2, -beta)": complex(a2, b1),
               "top_left (a1, 0)": complex(a1, 0), "top_right (a2, 0)": complex(a2, 0)}
    cm = {k: lam_min(l, w) / muHO for k, w in corners.items()}
    edge = [lam_min(l, complex(a, b1)) / muHO for a in np.linspace(a1, a2, 15)]
    col = {f"{t:.2f}": lam_min(l, complex(wl, -t * gam0)) / muHO for t in (0, 0.25, 0.5, 0.75, 1.0)}
    psm = lam_min(l, complex(ps[0], ps[1]) * km) / muHO
    mu_cert = float(st["mu_cert"].split()[0].strip("[")) / muHO
    res["B"][f"l{l}"] = {"lambda_over_muHO_corners": cm, "min_over_bottom_edge": min(edge),
                         "argmin_bottom_edge_Re_over_kappa_minus": float(np.linspace(a1, a2, 15)[int(np.argmin(edge))] / km),
                         "author_float_low_corner": st["float_lambda_min"] / muHO,
                         "column_Re=w_l_Im=-t_gamma0": col, "at_PS_mode": psm, "mu_cert_over_muHO": mu_cert,
                         "mu_cert_le_min_window": bool(mu_cert <= min(min(cm.values()), min(edge)))}
    print("B", l, res["B"][f"l{l}"], flush=True)


# ------------------------------------------------------------------------------ C. where the loss is
def Gp_of(nodes, x):
    xs = np.array([float(a) for a, _ in nodes])
    ys = np.array([float(b) for _, b in nodes])
    gp = np.interp(x, xs, ys, left=0.0, right=0.0)
    gpp = np.zeros_like(x)
    k = np.searchsorted(xs, x) - 1
    ok = (k >= 0) & (k < len(xs) - 1)
    gpp[ok] = (ys[k[ok] + 1] - ys[k[ok]]) / (xs[k[ok] + 1] - xs[k[ok]])
    return gp, gpp


res["C"] = {}
for l in (20, 50, 100):
    st = json.loads((P2 / "certificates" / "m3" / f"m3_PB_l{l}.json").read_text())
    nodes = [(Fraction(a), Fraction(b)) for a, b in st["ground_state_nodes"]]
    muHO = 2 * (l + 0.5) * OmP * gam0
    wl = (l + 0.5) * OmP
    a1, a2, b1 = wl - 8 * gam0, wl + 6 * gam0, -beta
    m = (s > -8) & (s < 8)
    x = s[m]
    Vs = V(r[m], l)
    phinf, s1, eta = -1.4375, 5.0, 0.5
    ph = phinf * x ** 4 / (s1 ** 4 + x ** 4)
    d1 = phinf * 4 * x ** 3 * s1 ** 4 / (s1 ** 4 + x ** 4) ** 2
    d2 = 4 * phinf * s1 ** 4 * x ** 2 * (3 * s1 ** 4 - 5 * x ** 4) / (s1 ** 4 + x ** 4) ** 3
    sg = np.where(x > 0, -1.0, 1.0)
    c = np.cos(ph) - eta * sg * d1
    cp = -np.sin(ph) * d1 - eta * sg * d2
    gp, gpp = Gp_of(nodes, x)
    base = cp * gp - c * gp ** 2 + c * gpp - sg * d1 / (4 * eta)

    def Emin(w2):
        return float(np.min(np.real(np.exp(-1j * ph) * 1j * (Vs - w2)) + base))
    hull = {"A_lo,B_lo": complex(a1 ** 2 - beta ** 2, -2 * a2 * beta), "A_lo,B_hi": complex(a1 ** 2 - beta ** 2, 0),
            "A_hi,B_lo": complex(a2 ** 2, -2 * a2 * beta), "A_hi,B_hi": complex(a2 ** 2, 0)}
    win = {"(a1,-beta)": complex(a1, b1) ** 2, "(a2,-beta)": complex(a2, b1) ** 2,
           "(a1,0)": complex(a1, 0) ** 2, "(a2,0)": complex(a2, 0) ** 2}
    eh = {k: Emin(v) / muHO for k, v in hull.items()}
    # E is affine in w^2, so min over s of min over the hull = min over s of the pointwise min of 4 corners
    pm = np.min(np.stack([np.real(np.exp(-1j * ph) * 1j * (Vs - v)) + base for v in hull.values()]), axis=0)
    pw = np.min(np.stack([np.real(np.exp(-1j * ph) * 1j * (Vs - v)) + base for v in win.values()]), axis=0)
    res["C"][f"l{l}"] = {"float_min_E_over_hull_over_muHO": float(pm.min() / muHO),
                         "float_min_E_over_true_window_corners_over_muHO": float(pw.min() / muHO),
                         "per_hull_corner": eh,
                         "mu_cert_over_muHO": float(st["mu_cert"].split()[0].strip("[")) / muHO}
    print("C", l, res["C"][f"l{l}"], flush=True)

res["seconds"] = round(time.time() - T0, 1)
(HERE / "indep_float_spectrum.json").write_text(json.dumps(res, indent=1) + "\n")
print("done", res["seconds"], "s")
