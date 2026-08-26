"""Stage D: monotonicity diagnostics of the EC profile on the closed sound cone
[centre, sonic point] = {x <= 0}, and plots.

Quantities (all functions of x = ln(-r/t)):
  V        fluid velocity relative to the static (r = const) observer (KHA's V;
           = HM01's V_R up to sign);
  v_rel    fluid velocity relative to the x = const (z = const) line,
           v_rel = (1+NV)/(N+V)  (HM01's |V_z|; = 1/sqrt3 at the sonic point);
  rho_hat  4 pi rho t^2 = omega e^{-2x} / A   (density in similarity units);
  w        1/3 - v_rel^2  (signed distance to the sonic point, GHJS's w = 1 - v^2);
  Omega    4 pi r^2 rho = omega / A.
"""
from __future__ import annotations

import numpy as np

from . import css, shoot

V0_EC = 0.112439401388092


def ec_profile(V0=V0_EC, branch=1, x_min=-8.0, x_max=6.0, n=2001):
    """Dense EC profile on [x_min, x_max] (inside from the scaled shot, outside from the
    plain continuation). Returns dict of arrays."""
    inner = shoot.shoot(V0, branch, x_end=x_min, keep_sol=True)
    outer = shoot.continue_outward(V0, branch, x_max=x_max)
    x = np.linspace(x_min, x_max, n)
    y = np.empty((4, n))
    for i, xi in enumerate(x):
        if xi <= -inner.delta:
            y[:, i] = css.from_scaled(xi, css.scaled3_to_full(xi, inner.sol.sol(xi)))
        elif xi >= outer.t[0]:
            N, W, V = outer.sol(xi)
            y[:, i] = [1 + css.Am1_constraint(N, W, V), N, W, V]
        else:                                               # |x| < delta: sonic series
            _, _, coef = shoot.start_data(V0, branch)
            y[:, i] = shoot.taylor.eval_series(coef, xi)
    A, N, W, V = y
    return dict(x=x, A=A, N=N, W=W, V=V, v_rel=css.v_rel(N, V), rho_hat=W * np.exp(-2 * x) / A,
                Omega=W / A, w=1.0 / 3.0 - css.v_rel(N, V) ** 2, mu=1 - 1 / A)


def monotonicity(prof, x_lo=-8.0, x_hi=0.0):
    """Sign statistics of d/dx of the diagnostic quantities on [x_lo, x_hi]."""
    m = (prof["x"] >= x_lo) & (prof["x"] <= x_hi)
    x = prof["x"][m]
    out = {}
    for key in ("V", "v_rel", "rho_hat", "W", "w", "Omega", "A", "mu"):
        f = prof[key][m]
        df = np.gradient(f, x)
        out[key] = dict(min_d=float(df.min()), max_d=float(df.max()),
                        frac_pos=float(np.mean(df > 0)), frac_neg=float(np.mean(df < 0)),
                        f_lo=float(f[0]), f_hi=float(f[-1]),
                        sign_changes_of_f=int(np.sum(np.sign(f[1:]) != np.sign(f[:-1]))))
    return out


def zero_of_V(prof):
    x, V = prof["x"], prof["V"]
    i = np.where(np.sign(V[1:]) != np.sign(V[:-1]))[0]
    return x[i] - V[i] * (x[i + 1] - x[i]) / (V[i + 1] - V[i])


def save_plots(prof, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    x = prof["x"]
    fig, ax = plt.subplots(2, 2, figsize=(10, 7))
    ax[0, 0].plot(x, prof["A"], label="A = a²")
    ax[0, 0].plot(x, prof["mu"], label="2m/r = 1 − 1/A")
    ax[0, 0].plot(x, prof["Omega"], label="Ω = 4πr²ρ")
    ax[0, 0].set_ylim(-0.1, 2.5)
    ax[0, 1].plot(x, prof["V"], label="V (rel. static observer)")
    ax[0, 1].plot(x, prof["v_rel"], label="v_rel (rel. x = const line)")
    ax[0, 1].axhline(1 / np.sqrt(3), ls=":", c="k", lw=0.8)
    ax[0, 1].axhline(0, ls=":", c="k", lw=0.8)
    ax[0, 1].set_ylim(-1.2, 1.3)
    ax[1, 0].semilogy(x, prof["N"], label="N = α/(a e^x)")
    ax[1, 0].semilogy(x, prof["W"], label="ω = 4πr²a²ρ")
    ax[1, 1].semilogy(x, prof["rho_hat"], label="ρ̂ = 4πρt² = ω e^{−2x}/A")
    for a in ax.ravel():
        a.axvline(0, c="gray", lw=0.6)
        a.set_xlabel("x = ln(−r/t)")
        a.legend(fontsize=8)
    fig.suptitle("Evans–Coleman CSS profile, k = 1/3 (sonic point at x = 0)")
    fig.tight_layout()
    fig.savefig(f"{outdir}/ec_profile.png", dpi=130)
    plt.close(fig)

    m = x <= 0
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.8))
    ax[0].plot(x[m], np.abs(prof["V"][m]), label="|V|")
    ax[0].plot(x[m], np.abs(prof["v_rel"][m]), label="|v_rel|")
    ax[0].axhline(1 / np.sqrt(3), ls=":", c="k", lw=0.8)
    ax[1].semilogy(x[m], prof["rho_hat"][m], label="ρ̂ = ωe^{−2x}/A")
    ax[1].semilogy(x[m], prof["W"][m] * np.exp(-2 * x[m]), label="ωe^{−2x}")
    ax[2].plot(x[m], prof["w"][m], label="1/3 − v_rel²")
    ax[2].plot(x[m], np.gradient(prof["w"][m], x[m]), label="d/dx (1/3 − v_rel²)")
    ax[2].axhline(0, ls=":", c="k", lw=0.8)
    for a in ax:
        a.set_xlabel("x = ln(−r/t)")
        a.legend(fontsize=8)
    fig.suptitle("Monotonicity diagnostics on the sound cone x ≤ 0")
    fig.tight_layout()
    fig.savefig(f"{outdir}/ec_monotonicity.png", dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    import json
    import os
    import sys

    outdir = sys.argv[1] if len(sys.argv) > 1 else "problems/P4/results"
    os.makedirs(outdir, exist_ok=True)
    prof = ec_profile()
    mono = monotonicity(prof)
    print(json.dumps(mono, indent=1))
    print("zeros of V:", zero_of_V(prof))
    save_plots(prof, outdir)
    np.savez(f"{outdir}/ec_profile.npz", **prof)
    with open(f"{outdir}/ec_monotonicity.json", "w") as f:
        json.dump(mono, f, indent=1)
