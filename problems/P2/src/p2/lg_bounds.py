"""P2-2: certified high-frequency thresholds W_l (Liouville-Green with Olver-type error control).

Statement certified per (point X, l):  W_l(omega) != 0 for all omega with |Re omega| >= W_l and
-(1/2 + delta) kappa_- <= Im omega <= 0   (Im omega >= 0 is P2-0; Re omega < 0 by the exact
symmetry W_l(-conj w) = conj W_l(w)).  Proof: notes/s1-p2-2-lg.md.  Summary (x = r*, ' = d/dx):

  Phi'' + p^2 Phi = 0,  p = sqrt(w^2 - V_l) (principal branch; Re(w^2 - V) >= a^2 - b^2 - V > 0),
  w = a + i b, a >= a0, -beta <= b <= 0, beta = (1/2 + delta) kappa_-.
  Jost solutions (= nonzero multiples of r R_+ and r R_c, because beta < kappa_+, kappa_c):
    Phi_+ = p^{-1/2} e^{-i S}(1 + h_+),  Phi_c = p^{-1/2} e^{+i S}(1 + h_c),  S' = p,
  with Volterra equations whose kernels are bounded by (1 + e^{2 int |Im p|})/2 and e^{2 int |Im p|}
  (the reference points are at the horizons, where these solutions are dominant when b < 0;
  the weights are integrable against theta/p = V''/(4p^3) + 5 V'^2/(16 p^5) ~ e^{2 kappa |x|}).
  For a matching point x_X:
    |h_+(x_X)| <= e^{Psi_+} - 1,   |h_+'(x_X)/p| <= e^{Psi_+} Psit_+,
    Psi_+  = int_{-inf}^{x_X} (1 + e^{2 int_t^{x_X} |Im p|})/2 |theta/p| dt,
    Psit_+ = int_{-inf}^{x_X} e^{2 int_t^{x_X} |Im p|} |theta/p| dt,      (c side likewise)
  and the Wronskian is 2i (const)[(1+h_+)(1+h_c) + ((1+h_+) h_c' - h_+'(1+h_c))/(2 i p)], which
  cannot vanish if   e^{Psi_+ + Psi_c} (1 + (Psit_+ + Psit_c)/2) < 2.            (*)
Uniform majorants (valid for all a >= a0, b in [-beta, 0]):
    |p|^2 >= m(r) := a0^2 - beta^2 - V(r) > 0,
    |theta/p| dx = J(r) dr,  J <= |(F V_r)_r|/(4 m^{3/2}) + 5 F V_r^2/(16 m^{5/2}),
    |Im p| <= beta G(r),  G = max(1, a0/sqrt(m))  (G = 1 where beta^2 + V <= 0),
    int dx = int dr/F = r*(r_hi) - r*(r_lo),  r* = sum_i A_i ln|r - r_i|  (exact antiderivative).
Quadrature: upper Riemann sums on r-cells (Arb interval evaluation, outward rounding); the two
end cells [r_+, r_+ + e] and [r_c - e, r_c] are bounded analytically:
    weight <= E e^{2 beta G B} ((r_1 - r_+)/(r - r_+))^s,  s = beta G/kappa_+ < 1,
    int (1 + weight)/2 J <= J_max e (1 + E e^{2 beta G B}/(1 - s))/2.

Replay:  .venv/bin/python -m p2.lg_bounds verify [PB|PA] [lmax]    (cwd problems/P2/src)
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

from flint import arb, arb_poly, fmpq

from .arbutil import git_commit, precision
from .params import Point

CERT_DIR = Path(__file__).resolve().parents[2] / "certificates" / "lg"
PREC = 128


class LG:
    def __init__(self, pt: Point, l: int, prec=PREC):
        self.pt, self.l, self.prec = pt, l, prec
        with precision(prec):
            L = l * (l + 1)
            D = arb_poly([pt.Q2, arb(-2), arb(1), arb(0), -pt.Lam / 3])
            Dp = D.derivative()
            r = arb_poly([0, 1])
            N = D * (L * r * r + r * Dp - 2 * D)          # V = N / r^6
            self.D, self.Dp, self.N = D, Dp, N
            self.Np, self.Npp = N.derivative(), N.derivative().derivative()
            self.beta = pt.beta
            self.roots = pt.roots
            self.A = pt.Acoef

    # -- pointwise (interval) quantities --------------------------------------------------------
    def F(self, r):
        return self.D(r) / (r * r)

    def V(self, r):
        return self.N(r) / r ** 6

    def parts(self, r):
        """(F, V, (F V_r)_r, F V_r^2) on the arb interval r."""
        D, Dp, N, Np, Npp = self.D(r), self.Dp(r), self.N(r), self.Np(r), self.Npp(r)
        F = D / (r * r)
        Fr = (r * Dp - 2 * D) / r ** 3
        V = N / r ** 6
        Vr = (r * Np - 6 * N) / r ** 7
        Vrr = (r * r * Npp - 12 * r * Np + 42 * N) / r ** 8
        return F, V, Fr * Vr + F * Vrr, abs(F) * Vr * Vr

    def rstar(self, r):
        """r* = sum A_i ln|r - r_i| (exact antiderivative of 1/F up to a constant)."""
        return sum((Ai * abs(r - ri).log() for Ai, ri in zip(self.A, self.roots)), arb(0))

    def cell(self, lo, hi, a0):
        """Upper bounds on the cell [lo, hi]: J_max, G_max (or None if m is not > 0)."""
        X = (lo + hi) / 2 + arb(0, 1) * ((hi - lo) / 2)
        _F, V, FVr_r, FVr2 = self.parts(X)
        m = a0 * a0 - self.beta ** 2 - V
        if not (m > 0):
            return None
        J = abs(FVr_r) / (4 * m ** arb(1.5)) + 5 * FVr2 / (16 * m ** arb(2.5))
        if (self.beta ** 2 + V) < 0:
            G = arb(1)
        else:
            G = arb((a0 / m.sqrt()).upper()).max(arb(1))
        return arb(J.upper()), G

    # -- the two sides ----------------------------------------------------------------------------
    def grid(self, side, rX, n_geo, n_uni, eps_rel, qgeo):
        """Cell endpoints (exact rationals) from the horizon end towards rX."""
        pt = self.pt
        rh = pt.rp if side == "+" else pt.rc
        span = abs(float((arb(rX) - rh).mid()))
        e0 = eps_rel * span
        # geometric part: distances e0, e0/q, ... up to span * 0.05, then uniform to rX
        ds = [e0]
        while ds[-1] < 0.05 * span and len(ds) < n_geo:
            ds.append(ds[-1] / qgeo)
        d_end = ds[-1]
        uni = [d_end + (span - d_end) * k / n_uni for k in range(1, n_uni + 1)]
        dists = ds + uni[:-1]
        sgn = 1 if side == "+" else -1
        rhf = float(rh.mid())
        pts = [fmpq(round((rhf + sgn * d) * 2 ** 40), 2 ** 40) for d in dists]
        pts.append(fmpq(rX))
        return pts

    def side_bounds(self, side, a0, rX, n_geo=400, n_uni=2000, eps_rel=1e-7, qgeo=0.93):
        """Upper bounds (Psi, Psit) for one side, uniform over a >= a0, b in [-beta, 0]."""
        pt = self.pt
        beta = self.beta
        pts = self.grid(side, rX, n_geo, n_uni, eps_rel, qgeo)
        # ordered from horizon outward: cells [pts[k], pts[k+1]] (as r-intervals)
        cells = []
        for k in range(len(pts) - 1):
            lo_, hi_ = (arb(pts[k]), arb(pts[k + 1])) if side == "+" else (arb(pts[k + 1]), arb(pts[k]))
            c = self.cell(lo_, hi_, a0)
            if c is None:
                return None
            dx = abs(self.rstar(arb(pts[k + 1])) - self.rstar(arb(pts[k])))
            cells.append((c[0], c[1], dx, abs(hi_ - lo_)))
        # cumulative weights: for cell k, I_k = sum_{j>=k} beta G_j dx_j (from cell k to rX)
        n = len(cells)
        I = [arb(0)] * (n + 1)
        for k in reversed(range(n)):
            J, G, dx, w = cells[k]
            I[k] = I[k + 1] + beta * G * dx
        Psi, Psit = arb(0), arb(0)
        for k in range(n):
            J, G, dx, w = cells[k]
            wt = (2 * I[k]).exp()
            Psi += (1 + wt) / 2 * J * w
            Psit += wt * J * w
        # end cell [r_h, r_1] analytic
        rh = pt.rp if side == "+" else pt.rc
        kap = pt.kp if side == "+" else pt.kc
        r1 = arb(pts[0])
        lo_, hi_ = (rh, r1) if side == "+" else (r1, rh)
        X = (lo_.lower() + hi_.upper()) / 2 + arb(0, 1) * ((hi_.upper() - lo_.lower()) / 2)
        _F, V, FVr_r, FVr2 = self.parts(X)
        m = a0 * a0 - beta ** 2 - V
        if not (m > 0):
            return None
        Jmax = arb((abs(FVr_r) / (4 * m ** arb(1.5)) + 5 * FVr2 / (16 * m ** arb(2.5))).upper())
        Gmax = arb((a0 / m.sqrt()).upper()).max(arb(1))
        # r* - A_h ln|r - r_h| is smooth on the end cell: its oscillation B bounds the remainder
        others = [(Ai, ri) for Ai, ri in zip(self.A, self.roots) if not ri.overlaps(rh)]
        reg = lambda r: sum((Ai * abs(r - ri).log() for Ai, ri in others), arb(0))
        B = arb(abs(reg(X) - reg(r1)).upper())
        s = beta * Gmax / kap
        if not (s < 1):
            return None
        E1 = (2 * I[0]).exp()
        e = abs(r1 - rh)
        amp = E1 * (2 * beta * Gmax * B).exp() / (1 - s)
        Psi += Jmax * e * (1 + amp) / 2
        Psit += Jmax * e * amp
        return {"Psi": arb(Psi.upper()), "Psit": arb(Psit.upper()), "cells": n,
                "weight_exponent_total": float((2 * I[0]).upper()), "s_end": float(s.upper())}

    def criterion(self, a0, rX, **kw):
        with precision(self.prec):
            a0 = arb(a0)
            bp = self.side_bounds("+", a0, rX, **kw)
            bc = self.side_bounds("c", a0, rX, **kw)
            if bp is None or bc is None:
                return None
            val = (bp["Psi"] + bc["Psi"]).exp() * (1 + (bp["Psit"] + bc["Psit"]) / 2)
            return {"value": val, "ok": bool(val < 2), "plus": bp, "c": bc}


def vmax_float(lg: LG):
    """Float location/value of max V_l on (r_+, r_c) (design only)."""
    import mpmath as mp
    pt = lg.pt
    lo_, hi_ = float(pt.rp.mid()), float(pt.rc.mid())
    best = (None, -1e9)
    for k in range(1, 4000):
        r = lo_ + (hi_ - lo_) * k / 4000
        v = float(lg.V(arb(r)).mid())
        if v > best[1]:
            best = (r, v)
    del mp
    return best


def find_threshold(pt: Point, l: int, tol_rel=2e-3, verbose=False):
    """Smallest a0 (exact dyadic, bisection on the monotone criterion) with (*) certified."""
    lg = LG(pt, l)
    with precision(PREC):
        rVmax, Vmax = vmax_float(lg)
        beta = float(pt.beta.mid())
        base = (Vmax + beta * beta) ** 0.5
        cands = sorted({fmpq(round(rVmax * 2 ** 20), 2 ** 20),
                        fmpq(round(float(pt.rP.mid()) * 2 ** 20), 2 ** 20)})
        best = None
        for rX in cands:
            lo_, hi_ = base * 1.0001, base * 1.5 + 10 * float(pt.km.mid())
            while True:
                res = lg.criterion(fmpq(int(hi_ * 2 ** 24) + 1, 2 ** 24), rX)
                if res is not None and res["ok"]:
                    break
                lo_, hi_ = hi_, hi_ * 1.5
            for _ in range(40):
                if hi_ - lo_ < tol_rel * hi_:
                    break
                mid = (lo_ + hi_) / 2
                res = lg.criterion(fmpq(int(mid * 2 ** 24) + 1, 2 ** 24), rX)
                if res is not None and res["ok"]:
                    hi_ = mid
                else:
                    lo_ = mid
            if best is None or hi_ < best[0]:
                best = (hi_, rX)
            if verbose:
                print(f"  l={l} rX={float(rX):.5f} threshold ~ {hi_:.6g}", flush=True)
        a0 = fmpq(int(best[0] * 2 ** 24) + 1, 2 ** 24)
        return a0, best[1], lg, (rVmax, Vmax)


def certify(name: str, l: int, verbose=False):
    t0, c0 = time.time(), time.process_time()
    pt = Point(name)
    a0, rX, lg, (rVmax, Vmax) = find_threshold(pt, l, verbose=verbose)
    t_search = time.time() - t0
    t1, c1 = time.time(), time.process_time()
    res = lg.criterion(a0, rX)
    t_cert, c_cert = time.time() - t1, time.process_time() - c1
    with precision(PREC):
        km = pt.km
        a0b = arb(a0)
        out = {
            "claim": f"P2-2 at {name}, l = {l}: W_l(omega) != 0 for |Re omega| >= W_l, "
                     f"-(1/2+delta) kappa_- <= Im omega <= 0",
            "point": name, "qfrac": str(pt.qfrac), "delta": str(pt.delta), "l": l,
            "W_l_exact": str(a0), "W_l_float": float(a0b.mid()),
            "W_l_over_kappa_minus": (a0b / km).str(12),
            "W_l_over_l_half_OmegaP": (a0b / ((l + arb(1) / 2) * pt.OmegaP)).str(12),
            "sqrt_Vmax_plus_beta2_over_kappa_minus_float": (Vmax + float(pt.beta.mid()) ** 2) ** 0.5 / float(km.mid()),
            "planning_value_over_kappa_minus": float(((arb(5) / 4) * (l + arb(1) / 2) * pt.OmegaP / km + 2).mid()),
            "matching_point_exact": str(rX), "r_Vmax_float": rVmax,
            "strip_depth_beta": pt.beta.str(20),
            "criterion": {"statement": "exp(Psi_+ + Psi_c) (1 + (Psit_+ + Psit_c)/2) < 2",
                          "value_upper": res["value"].str(12), "certified": res["ok"],
                          "Psi_plus": res["plus"]["Psi"].str(8), "Psit_plus": res["plus"]["Psit"].str(8),
                          "Psi_c": res["c"]["Psi"].str(8), "Psit_c": res["c"]["Psit"].str(8),
                          "cells": [res["plus"]["cells"], res["c"]["cells"]],
                          "end_cell_exponents_s": [res["plus"]["s_end"], res["c"]["s_end"]],
                          "weight_exponents": [res["plus"]["weight_exponent_total"], res["c"]["weight_exponent_total"]]},
            "grid": {"n_geo": 400, "n_uni": 2000, "eps_rel": 1e-7, "qgeo": 0.93, "dyadic_bits": 40},
            "prec_bits": PREC, "commit": git_commit(),
            "timing_s": {"search_wall": round(t_search, 2), "certify_wall": round(t_cert, 2),
                         "certify_cpu": round(c_cert, 2), "total_cpu": round(time.process_time() - c0, 2)}}
    return out


def run(name, lmin, lmax, verbose=False):
    CERT_DIR.mkdir(parents=True, exist_ok=True)
    for l in range(lmin, lmax + 1):
        out = certify(name, l, verbose)
        (CERT_DIR / f"lg_{name}_l{l}.json").write_text(json.dumps(out, indent=1) + "\n")
        print(f"{name} l={l}: W_l/k_- = {out['W_l_over_kappa_minus']}  (/(l+1/2)Omega_P = "
              f"{out['W_l_over_l_half_OmegaP']})  crit {out['criterion']['value_upper']} "
              f"cpu {out['timing_s']['total_cpu']}s", flush=True)


def verify(name, lmax):
    """Replay: recompute the criterion at the stored (W_l, matching point) for each l."""
    bad = []
    pt = Point(name)
    for l in range(lmax + 1):
        f = CERT_DIR / f"lg_{name}_l{l}.json"
        if not f.exists():
            bad.append(f"missing {f.name}")
            continue
        st = json.loads(f.read_text())
        lg = LG(pt, l)
        a0 = fmpq(*[int(x) for x in st["W_l_exact"].split("/")]) if "/" in st["W_l_exact"] else fmpq(int(st["W_l_exact"]))
        rX = fmpq(*[int(x) for x in st["matching_point_exact"].split("/")])
        res = lg.criterion(a0, rX)
        if res is None or not res["ok"]:
            bad.append(f"l={l}: criterion not certified on replay")
    return bad


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "run"
    name = sys.argv[2] if len(sys.argv) > 2 else "PB"
    if cmd == "run":
        lmin, lmax = int(sys.argv[3]), int(sys.argv[4])
        run(name, lmin, lmax, verbose=True)
    else:
        lmax = int(sys.argv[3]) if len(sys.argv) > 3 else 50
        b = verify(name, lmax)
        print("lg replay", name, "OK" if not b else b)
        sys.exit(1 if b else 0)
