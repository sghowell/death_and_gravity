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

from flint import acb, acb_mat, arb

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
    tube = lintube.Tube.load(tube_path)
    ok_t, nsteps, pinv = tube_regular(tube)
    boxes = cover(rect, w0, workers)
    summ = dict(rect=rect, tube_regular=ok_t, n_steps=nsteps, max_Pinv=pinv, n_boxes=len(boxes),
                min_w=min(b["w"] for b in boxes), min_nu=min(b["sonic"]["nu"] for b in boxes),
                max_re_sigma=max(b["sonic"]["re_sigma_max"] for b in boxes), max_Z_sonic=max(b["sonic"]["Z"] for b in boxes),
                max_Z_centre=max(b["centre"]["Z"] for b in boxes), max_g_centre=max(b["centre"]["g"] for b in boxes),
                time=time.time() - t0, boxes=boxes)
    if out:
        json.dump(summ, open(out, "w"))
    return summ
