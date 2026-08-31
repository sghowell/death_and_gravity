"""Theorem C scoping (S4 pre-stage): FLOAT numerical sizing of E(kappa) at large |kappa|
and of the WKB slopes mu_{1,2}(x).  Nothing here is validated; see notes/theorem-c-plan.md.

Run from the repo root:  PYTHONPATH=problems/P4/src uv run python -m p4.validated.tc_sizing
(~10 min; pass "mu" or "rays" or "arcs" to run one part).
"""
from __future__ import annotations

import sys

import numpy as np

from p4 import css, perturb, shoot, taylor

V0 = 0.112439401388092


def mu_slopes(x, NWV):
    """Eigenvalues of P_fl^{-1} P_s,fl (the kappa-matrix of the fluid block), sorted."""
    N, W, V = NWV
    q = css.coeffs(1.0 + css.Am1_constraint(N, W, V), N, W, V)
    P = np.array([[q["a"], q["b"]], [q["c"], q["d"]]])
    Ps = np.array([[q["sa"], q["sb"]], [q["sc"], q["sd"]]])
    return np.sort(np.linalg.eigvals(np.linalg.solve(P, Ps)).real)


def run_mu():
    sh = shoot.shoot(V0, 1, x_end=-8.0, keep_sol=True)
    bg = taylor.background_series(V0, 1, K=36)
    plain = lambda x: (lambda y: (y[0] * np.exp(-x), y[1] * np.exp(2 * x), y[2] * np.exp(x)))(sh.sol.sol(x))
    print("# x, mu-, mu+ near 0 (Taylor bg): x*mu- -> 0.90989 (= -d sigma/d kappa), mu+ -> 0.43425")
    for x in (-1e-1, -1e-2, -1e-3, -1e-4):
        m = mu_slopes(x, taylor.eval_series(bg, x)[1:])
        print(f"{x:9.1e}  {m[0]:+12.6f}  {m[1]:+.6f}   x*mu- = {x*m[0]:+.6f}")
    print("# interior/centre: mu+- ~ +-1.4011 e^x as x -> -inf")
    for x in (-0.5, -1, -2, -4, -8):
        m = mu_slopes(x, plain(x))
        print(f"{x:5.1f}  {m[0]:+.6f}  {m[1]:+.6e}   mu+*e^-x = {m[1]*np.exp(-x):+.6f}")
    xs = np.linspace(-8.0, -0.1, 400)
    I = np.trapezoid([mu_slopes(x, plain(x))[1] for x in xs], xs)
    xs2 = -np.logspace(-1, -6, 300)                  # increasing: -0.1 -> -1e-6
    I2 = np.trapezoid([mu_slopes(x, taylor.eval_series(bg, x)[1:])[1] for x in xs2], xs2)
    tail = mu_slopes(-8.0, plain(-8.0))[1]           # mu+ ~ c e^x: int_{-inf}^{-8} = mu+(-8)
    print(f"Phi+ = int mu+ over (-inf, 0] ~ {I + I2 + tail:.4f}  (pieces {I:.4f} + {I2:.4f} + {tail:.5f})")


def scan(prob, tag, kappas, chunk=12):
    kappas = np.asarray(kappas, dtype=complex)
    E = np.concatenate([prob.E(kappas[i:i + chunk]) for i in range(0, len(kappas), chunk)])
    a = np.abs(E)
    j = int(np.argmin(a))
    dips = [f"{kappas[i]:.2f}:{a[i]:.2e}" for i in range(1, len(a) - 1) if a[i] < a[i - 1] and a[i] < a[i + 1]]
    print(f"## {tag}: min|E| = {a[j]:.3e} at {kappas[j]:.3f} (|kappa|*min = {abs(kappas[j])*a[j]:.3f}), "
          f"max = {a.max():.3e}, interior local minima: {dips or 'none'}")
    return kappas, E


def run_rays(prob):
    scan(prob, "ray A: 15 + t, t = 0..45", 15.0 + np.arange(0, 46, 1.5))
    scan(prob, "ray B: s + 14i, s = 0..15", np.arange(0, 15.01, 0.75) + 14j)
    k, E = scan(prob, "ray C: i tau, tau = 14..60", 1j * np.arange(14, 60.01, 1.0))
    d = np.diff(np.unwrap(np.angle(E)))
    print(f"   ray C: d(arg E)/d tau = {d[0]:.3f} (tau=14) -> {d[-1]:.3f} (tau=60); "
          f"tau*|E| = {abs(k[0])*abs(E[0]):.3f} -> {abs(k[-1])*abs(E[-1]):.3f}")
    for s in (0.5, 1.3, 2.0):
        scan(prob, f"strip: {s} + i tau, tau = 14..60", s + 1j * np.arange(14, 60.01, 1.0))


def run_arcs(prob):
    for r in (20.0, 30.0, 50.0):
        scan(prob, f"arc |kappa| = {r:.0f}, first quadrant", r * np.exp(1j * np.linspace(0, np.pi / 2, 61)))


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("mu", "all"):
        run_mu()
    if which in ("rays", "arcs", "all"):
        prob = perturb.Problem(V0, x_end=-8.0)
        if which in ("rays", "all"):
            run_rays(prob)
        if which in ("arcs", "all"):
            run_arcs(prob)
