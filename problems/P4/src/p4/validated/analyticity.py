"""Theorem B, Stage 3 (b): the hypotheses under which E is analytic on the closed rectangle R, so
that the winding number of ``modecount`` counts its zeros (argument principle).

E(kappa) = e^{3 x_d} det[r_1(kappa), r_2(kappa), Phi(x_d, x_0; kappa) p(x_0; kappa)]_{(A^, n~, v~)} with
 (S) p(x_0; kappa): the 4D sonic series (``linsonic4``).  On each box of a cover of R the box run
     certifies: every level matrix M_n(kappa) = n D + E(kappa) is invertible for 1 <= n <= K (ball
     solves) and for n > K (rank-1 bound, Re sigma(kappa) < K + 1) -- non-resonance on the whole box --
     and the tail bound Z < 1 holds uniformly on the box with nu > |x_0|, so p(x_0; .) is a uniform
     limit of the analytic partial sums sum_{n<=N} p_n(kappa) x_0^n (p_n rational in kappa without
     poles on the box): analytic on the box.
 (P) Phi: the scaled 4D system has kappa-free P, invertible on every sub-box of every tube step
     (``tube_regular``: ball inverses), and G affine in kappa: Phi(x_d, x_0; .) is entire.
 (C) r_1, r_2: on each box the ``lincentre.RegularFamily`` certificate over the box kappa proves
     P_0 invertible, (n P_0 - G_0(kappa)) invertible for 1 <= n <= K (ball solves) and n > K
     (||P_0^{-1} G_0|| < K + 1: no positive-integer exponent), rank G_0 = 2 (kappa-independent minor)
     and Z < 1: the exponent-0 family is 2-dimensional and its basis is analytic in kappa on the box.
Hence E is analytic on an open neighbourhood of R (each closed tile is certified on a 5% larger box,
so the open certified boxes cover R); zeros of E <=> p~(x_d) in the regular plane.
The covers are adaptive bisections of R (``cover``); ``certify_analyticity`` runs everything.
"""
from __future__ import annotations

import json
import multiprocessing as mp
import time

from flint import acb, acb_mat, arb, fmpq

from . import lincentre, linsonic4, linscaled
from .arbseries import precision
from .linsonic4 import _upper
from .linstep import norm_inf
from .linsys import kappa_box

_G = {}
ENLARGE = 1.05          # tiles of half-width w are certified on half-width 1.05 w, so the open boxes cover R


def tube_regular(tube, S=None, nsub=8, prec=256):
    """Ball-invert P (kappa-free) on every sub-box of every step: (ok, n_steps, max ||P^-1||)."""
    S = S or linscaled.full_system()
    worst, n = arb(0), 0
    with precision(prec):
        for sd in tube.steps:
            data = tube.system_data(sd, S, nsub)
            for b in range(nsub):
                Pi = acb_mat([[acb(e) for e in row] for row in data["Pbox"][0][b]]).inv()     # raises if not invertible
                worst = worst.max(norm_inf(Pi.tolist()))
            n += 1
    return True, n, float(worst)


def sonic_box(bg, cx, cy, w, K=40, x0=-0.05, prec=256):
    """(ok, details) of the 4D sonic certificate on the box cx + i cy +/- w (+/- i w)."""
    with precision(prec):
        try:
            ex = linsonic4.linear_sonic_expansion4(bg, acb(cx, cy), width=w, m=0, K=K)
            cert = ex.certify()
        except ZeroDivisionError:
            return False, dict(reason="singular level matrix over the box")
        if not cert.ok:
            return False, dict(reason=cert.details.get("reason"), Z=float(cert.Z1 + cert.Z2))
        sg = cert.details["sigma"]
        ok = bool(cert.nu > abs(x0))
        return ok, dict(nu=float(cert.nu), Z=float(cert.Z1 + cert.Z2), eps=float(cert.eps),
                        re_sigma_max=float(_upper(sg.real)), g=float(cert.details["g"]))


def centre_box(S, ce, cx, cy, w, K=50, nu=0.06, prec=256):
    with precision(prec):
        try:
            rf = lincentre.RegularFamily(S, ce, kappa_box(acb(cx, cy), w), K=K, prec=prec)
            ok, eps, det = rf.certify(nu, prec=prec)
        except ZeroDivisionError:
            return False, dict(reason="singular level matrix over the box")
        if not ok:
            return False, dict(reason="Z >= 1 or g >= K+1", **{k: v for k, v in det.items() if k != "eps"})
        return True, dict(Z=det["Z1"] + det["Z2"], g=det["g"], eps=max(det["eps"]))


def _laurent(poly, k, z0):
    """t-order-k coefficient of t^{d_N} poly(N = n/t, W = w t^2, V = v t) at (n, w, v) = z0 (arb balls);
    the t-order of a monomial N^a W^b V^e is d_N - a + 2b + e >= 0 (exact bookkeeping)."""
    dN, tot = poly.degrees()[2], arb(0)
    for exps, c in poly.terms():
        et, eA, a, b, e = map(int, exps)
        assert et == 0 and eA == 0
        if dN - a + 2 * b + e == k:
            tot += arb(fmpq(c)) * z0[0] ** a * z0[1] ** b * z0[2] ** e
    return tot


def constraint_exponents(bg, ce, x_d=-3.0, prec=256):
    """Ball certificate of the conditions under which the exact identity  16 S D c' = Lambda_lin c
    (``linsys.linear_constraint_propagation``) forces the linearised constraint c to vanish identically:
    (a) sonic point: D(u_0) = 0 exactly (A1 closed forms), D_1 != 0 and gamma := Lambda_lin(u_0)/(16 S_0 D_1)
        not a positive integer  =>  every solution analytic at x = 0 with c(0) = 0 (the order-0 constraint of
        the 4D recursion) has c == 0 on its disc of convergence  (p~ in Sigma);
    (b) centre: with N = n/t, W = w t^2, V = v t, c = t^2 c~ and  theta c~ = (lambda(t) - 2) c~ with lambda =
        Lambda_lin/(16 S D) analytic at t = 0 (Lambda_lin has no N^3 monomial: exact), lambda(0) = [N^2]/(16 (-1))
        = -1, so rho := lambda(0) - 2 = -3 is no nonnegative integer  =>  every solution analytic in t at the
        centre (the exponent-0 family r_1, r_2) has c~ == 0  (r_i in Sigma);
    (c) e_W not in Sigma(x_d):  l_Sigma . e_W = 2 A T~ != 0 at u(x_d).
    Lambda_lin and D are kappa-free, so (a)-(c) hold for every kappa at once."""
    from flint import fmpq
    from .arbseries import Series
    from .linsys import LinSystem, linear_constraint_propagation, plain_from_scaled
    from .polysys import _PolyEvaluator
    with precision(prec):
        L = LinSystem()
        lam, D = linear_constraint_propagation(L)
        u = [s.with_cap(2) for s in bg.series()]
        ev = _PolyEvaluator([Series.var(2)] + u)
        D0, D1, S0, W0 = ev(D)[0], ev(D)[1], ev(L.S)[0], u[2][0]
        gamma = ev(lam)[0] / (16 * S0 * D1)
        ok_a = bool(D0.contains(arb(0)) and S0 != 0 and W0 != 0 and D1 != 0) and bool(gamma < 1 or not gamma.contains_integer())
        dN = lam.degrees()[2]
        ok_reg = not any(dN - int(e[2]) + 2 * int(e[3]) + int(e[4]) < dN - 2 for e, _ in lam.terms())
        z0 = ce.balls[0]
        Dc = _laurent(D, 0, z0)
        lam0 = _laurent(lam, dN - 2, z0) / (16 * Dc)
        rho = lam0 - 2
        ok_b = ok_reg and bool(Dc != 0) and bool(rho < 0 or not rho.contains_integer())
        A, N, W, V = plain_from_scaled(x_d, *ce.eval(x_d))
        lw = 2 * A * (1 + V * V / 3 + N * V * arb(fmpq(4, 3)))
        ok_c = bool(lw != 0)
        return dict(ok=ok_a and ok_b and ok_c, sonic_ok=ok_a, centre_ok=ok_b, eW_ok=ok_c, gamma_sonic=gamma, D1=D1,
                    D0=D0, rho_centre=rho, lambda_centre=lam0, D_centre=Dc, lSigma_eW=lw, n_terms=len(lam.coeffs()))


def _init():
    from . import linmatch
    from .modecount import V0_EC, W_V0, certified_centre
    with precision(256):
        _G["bg"] = linmatch.box_background(V0_EC, W_V0, K=41)
        _G["ce"] = certified_centre()
        _G["S4"] = linscaled.full_system()


def _check(box):
    """Bisect until both certificates hold; returns the list of certified sub-boxes."""
    cx, cy, w, depth = box
    we = w * ENLARGE                       # certified on the enlarged box: every point of the closed tile is interior
    oks, dets = sonic_box(_G["bg"], cx, cy, we), centre_box(_G["S4"], _G["ce"], cx, cy, we)
    if oks[0] and dets[0]:
        return [dict(cx=cx, cy=cy, w=w, w_cert=we, sonic=oks[1], centre=dets[1])]
    if depth >= 8:
        raise RuntimeError(f"box {box} not certifiable: {oks[1]} {dets[1]}")
    out = []
    for sx in (-1, 1):
        for sy in (-1, 1):
            out.extend(_check((cx + sx * w / 2, cy + sy * w / 2, w / 2, depth + 1)))
    return out


def cover(rect, w0=0.5, workers=2, verbose=True):
    """Certified cover of the closed rectangle by squares (half-width <= w0) for (S) and (C)."""
    a, b, c, d = rect
    boxes = []
    x = a
    while x < b - 1e-12:
        y = c
        while y < d - 1e-12:
            boxes.append((x + w0, y + w0, w0, 0))
            y += 2 * w0
        x += 2 * w0
    t0, out = time.time(), []
    with mp.get_context("spawn").Pool(workers, initializer=_init) as pool:
        for i, r in enumerate(pool.imap_unordered(_check, boxes)):
            out.extend(r)
            if verbose and (i + 1) % 100 == 0:
                print(f"[{i + 1}/{len(boxes)}] {len(out)} boxes, {time.time() - t0:.0f}s", flush=True)
    return out


def certify_analyticity(tube_path, rect=(0.0, 15.0, -14.0, 14.0), w0=0.5, workers=2, out=None):
    """(P), (S), (C) on R; summary dict (saved to ``out`` as JSON if given)."""
    from . import lintube
    t0 = time.time()
    from . import linmatch
    from .modecount import V0_EC, W_V0, certified_centre
    tube = lintube.Tube.load(tube_path)
    ok_t, nsteps, pinv = tube_regular(tube)
    last = tube.steps[-1]
    x_d = float(last.x - arb(last.h))
    with precision(256):
        cx = constraint_exponents(linmatch.box_background(V0_EC, W_V0, K=41), certified_centre(), x_d)
    cx = {k: (str(v) if isinstance(v, arb) else v) for k, v in cx.items()}
    boxes = cover(rect, w0, workers)
    summ = dict(rect=rect, w0=w0, enlarge=ENLARGE, x_d=x_d, tube_info=tube.info, constants=dict(V0=V0_EC, w_V0=W_V0, K_sonic=40, K_centre=50, nu_centre=0.06, x0=-0.05),
                constraint=cx, tube_regular=ok_t, n_steps=nsteps, max_Pinv=pinv, n_boxes=len(boxes),
                min_w=min(b["w"] for b in boxes), min_nu=min(b["sonic"]["nu"] for b in boxes),
                max_re_sigma=max(b["sonic"]["re_sigma_max"] for b in boxes), max_Z_sonic=max(b["sonic"]["Z"] for b in boxes),
                max_Z_centre=max(b["centre"]["Z"] for b in boxes), max_g_centre=max(b["centre"]["g"] for b in boxes),
                time=time.time() - t0, boxes=boxes)
    if out:
        json.dump(summ, open(out, "w"))
    return summ


if __name__ == "__main__":                       # PYTHONPATH=problems/P4/src uv run python -m p4.validated.analyticity TUBE.json
    import argparse
    import os
    from .modecount import RESULTS_DIR
    ap = argparse.ArgumentParser(description="analyticity cover of R + constraint-exponent certificate (results JSON)")
    ap.add_argument("tube")
    ap.add_argument("--rect", nargs=4, type=float, default=(0.0, 15.0, -14.0, 14.0), metavar=("A", "B", "C", "D"))
    ap.add_argument("--w0", type=float, default=0.5)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--out", default=os.path.join(RESULTS_DIR, "analyticity_R.json"))
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    s = certify_analyticity(a.tube, rect=tuple(a.rect), w0=a.w0, workers=a.workers, out=a.out)
    print({k: v for k, v in s.items() if k not in ("boxes", "tube_info")}, "->", a.out)
