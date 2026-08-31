"""Theorem C, S4-1 route (a): FLOAT leading term of E(kappa) at large |kappa|.

Law fitted and validated here (see notes/s4-reduction.md section 4):
    E(kappa) = c_s D_inf e^{-J} kappa^{-q} e^{Phi+ kappa} (1 + d/kappa + O(1/kappa^2)),
with Phi+ = int_{-inf}^0 mu+ dx, q the sonic-layer exponent, D_inf the t -> 0 limit of
det[r1, r2, e+] on the rows (Ahat, nhat, vhat) (kappa-free), and J = the zone integral
int_{-inf}^{-0.05} (B~_44 + 3) dx in
the Levinson frame (T columns: unit-norm scaled fluid eigenvectors, v-component > 0;
exponent integrals based at x = 0).  So c_s = C e^J / D_inf with C the fitted constant
of E ~ C kappa^{-q} e^{Phi+ kappa}.  Nothing here is validated arithmetic.

Run from the repo root (PYTHONPATH=problems/P4/src):
    uv run python -m p4.validated.tc_leading scan     # E on 3 rays, |kappa| 30..300 (~min)
    uv run python -m p4.validated.tc_leading fit      # joint fit of (Phi+, q, C, d per ray)
    uv run python -m p4.validated.tc_leading factor   # D_inf, J, c_s decomposition
Results in problems/P4/results/theorem_c/tc_leading.json (scan appends E values).
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

from p4 import css, perturb, shoot, taylor

V0 = 0.112439401388092
OUT = os.path.join(os.path.dirname(__file__), "..", "..", "..", "results", "theorem_c")
RAYS = {"real": 0.0, "diag": np.pi / 4, "imag": np.pi / 2}
RADII = np.geomspace(30.0, 300.0, 13)


def E_of(kappa):
    """E(kappa) = A_p(x_end) e^{x_end}; sonic start delta, K adapted to |kappa|."""
    r = abs(kappa)
    delta = min(0.04, 8.0 / r)
    prob = perturb.Problem(V0, x_end=-8.0, delta=delta, K=44)
    return complex(prob.E(np.array([kappa]))[0])


def run_scan(path):
    data = {"rays": {}}
    for tag, th in RAYS.items():
        ks, Es = [], []
        for r in RADII:
            kap = r * np.exp(1j * th)
            E = E_of(kap)
            ks.append([kap.real, kap.imag])
            Es.append([E.real, E.imag])
            print(f"{tag} |k|={r:7.2f}  E={E:.6e}", flush=True)
        data["rays"][tag] = {"kappa": ks, "E": Es}
    os.makedirs(path, exist_ok=True)
    with open(os.path.join(path, "tc_leading.json"), "w") as fh:
        json.dump(data, fh, indent=1)


def load(path):
    with open(os.path.join(path, "tc_leading.json")) as fh:
        return json.load(fh)


def fit(data, use_d=True):
    """ln|E| = ln|C| + Phi+ Re(kappa) - q ln|kappa| + Re(d_ray/kappa), real least
    squares jointly over the rays (shared |C|, Phi+, q; complex d per ray).  The
    complex constant C is then estimated per sample as E kappa^q e^{-Phi+ kappa}
    (principal branch) -- its constancy per ray validates the one-term law."""
    rows, rhs, tags = [], [], list(data["rays"])
    for j, tag in enumerate(tags):
        d = data["rays"][tag]
        ks = np.array([complex(a, b) for a, b in d["kappa"]])
        Es = np.array([complex(a, b) for a, b in d["E"]])
        for k, E in zip(ks, Es):
            row = np.zeros(3 + 2 * len(tags))
            row[0], row[1], row[2] = 1.0, k.real, -np.log(abs(k))
            if use_d:
                row[3 + 2 * j], row[4 + 2 * j] = (1 / k).real, -(1 / k).imag
            rows.append(row)
            rhs.append(np.log(abs(E)))
    M, v = np.array(rows), np.array(rhs)
    sol, *_ = np.linalg.lstsq(M, v, rcond=None)
    out = {"lnC": sol[0], "Phi": sol[1], "q": sol[2],
           "d": {t: complex(sol[3 + 2 * j], sol[4 + 2 * j]) for j, t in enumerate(tags)},
           "max_resid": float(np.abs(M @ sol - v).max()), "C": {}}
    for tag in tags:
        d = data["rays"][tag]
        ks = np.array([complex(a, b) for a, b in d["kappa"]])
        Es = np.array([complex(a, b) for a, b in d["E"]])
        Cs = Es * ks ** complex(out["q"]) * np.exp(-out["Phi"] * ks)
        out["C"][tag] = (Cs[-1], float(np.abs(np.diff(Cs[-5:])).max() / abs(Cs[-1])))
    return out


def scaled_BC(x, u, du):
    """B~(x), C~(x) of the centre-scaled 4D system q~' = (B~ - kappa C~) q~."""
    A, N, W, V = u
    q = css.coeffs(A, N, W, V)
    Pfl = np.array([[q["a"], q["b"]], [q["c"], q["d"]]])
    Psfl = np.array([[q["sa"], q["sb"]], [q["sc"], q["sd"]]])
    ell = np.exp(np.array([2 * x, -x, 2 * x, x]))
    from p4.tps import Dual
    B = np.zeros((4, 4))
    for l in range(4):
        y = [Dual(u[i], 1.0 if i == l else 0.0) for i in range(4)]
        qq = css.coeffs(*y, Am1=Dual(A - 1, 1.0 if l == 0 else 0.0))
        r3 = qq["a"] * Dual(du[2], 0.0) + qq["b"] * Dual(du[3], 0.0) + qq["e"]
        r4 = qq["c"] * Dual(du[2], 0.0) + qq["d"] * Dual(du[3], 0.0) + qq["f"]
        B[0, l] = (y[0] * qq["FA"]).b
        B[1, l] = (y[1] * qq["FN"]).b
        rr = np.linalg.solve(Pfl, -np.array([r3.b, r4.b]))
        B[2, l], B[3, l] = rr
    C = np.zeros((4, 4))
    C[2:, 2:] = np.linalg.solve(Pfl, Psfl)
    Bs = (B * ell[None, :] / ell[:, None]) - np.diag([2.0, -1.0, 2.0, 1.0])
    Cs = C * ell[None, :] / ell[:, None]
    return Bs, Cs


def eplus(Cs):
    """Unit eigenvector of the scaled fluid block for mu+ (largest eigenvalue)."""
    w, vec = np.linalg.eig(Cs[2:, 2:])
    i = int(np.argmax(w.real))
    v = vec[:, i].real
    v /= np.linalg.norm(v)
    if v[1] < 0:
        v = -v
    out = np.zeros(4)
    out[2:] = v
    return w[i].real, out


def run_factor():
    """Phi+, J = int(B~44 + 3), D_inf, and c_s = C/(D_inf e^J)."""
    sh = shoot.shoot(V0, 1, x_end=-16.0, keep_sol=True)
    bgs = taylor.background_series(V0, 1, K=36)

    def u_du(x):
        if x > -0.05:
            u = taylor.eval_series(bgs, x)
        else:
            y = sh.sol.sol(x)
            N, W, V = y[0] * np.exp(-x), y[1] * np.exp(2 * x), y[2] * np.exp(x)
            u = np.array([1 + css.Am1_constraint(N, W, V), N, W, V])
        return u, css.rhs_plain(x, u)

    xs = np.concatenate([np.linspace(-16, -0.1, 3000), -np.geomspace(0.1, 1e-6, 800)[1:]])
    mu, B44 = np.empty(len(xs)), np.empty(len(xs))
    T4 = np.empty((len(xs), 4))
    for i, x in enumerate(xs):
        u, du = u_du(x)
        Bs, Cs = scaled_BC(x, u, du)
        mu[i], T4[i] = eplus(Cs)
        B44[i] = np.nan                       # filled below with T'
    dT = np.gradient(T4, xs, axis=0)
    for i, x in enumerate(xs):
        u, du = u_du(x)
        Bs, Cs = scaled_BC(x, u, du)
        B44[i] = T4[i] @ (Bs @ T4[i] - dT[i])      # T orthonormal col: T^{-1} row = T^T
    Phi = np.trapezoid(mu, xs) + mu[0]             # tail int_{-inf}^{-16} ~ mu(-16)
    zone = xs <= -0.05                             # layer boundary delta* = 0.05:
    J = np.trapezoid(B44[zone] + 3.0, xs[zone]) + (B44[0] + 3.0)   # J = int_(-inf,-0.05](B44+3)
    res = [B44[i] * xs[i] for i in (-1, -200, -400)]
    print(f"Phi+ = {Phi:.6f}   J(0.05) = {J:.6f}   B44(-16) = {B44[0]:.6f}")
    print(f"sonic tail: x*B44 at x = {xs[-1]:.1e}, {xs[-200]:.1e}, {xs[-400]:.1e}: "
          f"{res[0]:.4f}, {res[1]:.4f}, {res[2]:.4f}  (~ +0.03/x, log absorbed in the layer)")
    # D_inf := lim_{t->0} det[r1, r2, e+] rows (Ahat, nhat, vhat); r_i(0) from the
    # certified regular centre family (kappa = 0), e+ evaluated deep in the centre zone.
    from flint import arb

    from p4.validated import centre, lincentre, linscaled
    from p4.validated.arbseries import precision
    with precision(256):
        a = arb('-0.2123656467659762832750918714540807905889') + arb('7.031e-12') + arb(0, 8.87e-16)
        m_ = arb('8.901323275379966931515526907200000000000') + arb('2.3757e-9') + arb(0, 4.61e-14)
        ce = centre.centre_expansion(m_ * (2 * a).exp(), nhat=(-a).exp(), K=30)
        ce.certify()
        rf = lincentre.RegularFamily(linscaled.full_system(), ce, 0.0, K=30)
        r0 = [[complex(v) for v in (b[0][i, 0] for i in range(4))] for b in rf.basis]
    M = np.array([[r0[0][i], r0[1][i], T4[0][i]] for i in (0, 1, 3)], dtype=complex)
    D_inf = complex(np.linalg.det(M))
    print(f"r1(0) = {np.round(np.array(r0[0]).real, 6)}   r2(0) = {np.round(np.array(r0[1]).real, 6)}")
    print(f"e+(-16) = {np.round(T4[0], 6)}   D_inf = {D_inf.real:.6f}")
    try:
        r = fit(load(OUT))
        C = np.exp(r["lnC"])
        cs = C * np.exp(J) / D_inf.real
        print(f"C = {C:.5f}  =>  c_s = C e^J / D_inf = {cs:.5f}   (q = {r['q']:.4f})")
    except FileNotFoundError:
        pass
    return Phi, J, D_inf, xs, B44, mu, T4


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "fit"
    if which == "scan":
        run_scan(OUT)
    elif which == "fit":
        data = load(OUT)
        r = fit(data)
        print(f"joint |E| fit: Phi+ = {r['Phi']:.6f}   q = {r['q']:.4f}   |C| = {np.exp(r['lnC']):.5f}")
        for t, dv in r["d"].items():
            Ct, drift = r["C"][t]
            print(f"  ray {t}: d = {dv:.3f};  C_est(|k|=300) = {Ct:.5f} (|C|={abs(Ct):.5f}, "
                  f"drift over last 5 samples {drift:.1e})")
        print(f"  max |resid(ln|E|)| = {r['max_resid']:.2e}")
        r2 = fit(data, use_d=False)
        print(f"without 1/kappa term: Phi+ = {r2['Phi']:.6f}, q = {r2['q']:.4f}, resid {r2['max_resid']:.2e}")
    elif which == "factor":
        run_factor()
