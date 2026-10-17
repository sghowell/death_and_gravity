"""Independent Arb re-check of the S2 window certificates (reviewer script, 2026-10-17).

For l in {20, 50, 100} (or the l's given on the command line) and the ground-state nodes stored in
certificates/m3/m3_PB_l{l}.json, prove E(s) >= mu on s in [-64, 32] and on both tails, for every
w^2 in the reviewer's own interval hull of the window, with code that shares nothing with m3.py:
  * parameters: the S1 certified balls in certificates/params.json (parsed, not recomputed);
  * r(s): Krawczyk enclosure of the root of the branch-tracked equation r*(r) = x_P + e^{i pi/4}s
    (point s for the midpoint value, interval s for the whole cell), not the author's ODE Taylor model;
  * V_l: expanded polynomial Delta(r) = -r^4/50 + r^2 - 2r + Q^2, V = F (L/r^2 + F'/r);
  * derivatives: forward-mode dual numbers (value, d/ds), not the hand-derived E' of m3.py;
  * cells: the reviewer's own cover (node intervals; G' = 0 regions in pieces of <= 1/16), bisected
    until the lower bound reaches the author's mu_cert or depth 18;
  * tails: own triangle-inequality disc bound and own lower bound of q there.
Mean-value form per cell: E(s) >= E(s_m) - sup_cell |E'| w/2  (E is C^1 inside a cell; G' linear).
Run (repo root):  nice -n 10 .venv/bin/python problems/P2/certificates/review_s2/indep_cells_arb.py [l ...]
"""
import json
import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np
from flint import acb, arb, ctx, fmpq
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
P2 = HERE.parents[1]
ctx.prec = 128

J = json.loads((P2 / "certificates" / "params.json").read_text())["points"]["PB"]["quantities"]
ball = lambda k: arb(J[k]["ball"])
ROOTS = [ball("r_n"), ball("r_minus"), ball("r_plus"), ball("r_c")]
AB = [arb(J["A_coefficients"][k]["ball"]) for k in ("r_n", "r_minus", "r_plus", "r_c")]
Q2 = ball("Q2")
OMP, G0, BETA, KM, RP = ball("Omega_P"), ball("gamma0"), arb(J["strip_depth_beta"]["ball"]), ball("kappa_minus"), ball("r_P")
XP = sum((a * abs(acb(RP) - ri).log() for a, ri in zip(AB, ROOTS)), arb(0))
EI4 = acb(1, 1) / arb(2).sqrt()
I_ = acb(0, 1)


def q_(x):
    x = Fraction(x)
    return arb(fmpq(x.numerator, x.denominator))


# ------------------------------------------------------------------------------------- dual numbers
class D:
    __slots__ = ("v", "d")

    def __init__(self, v, d=0):
        self.v, self.d = v, d

    def __add__(self, o):
        o = o if isinstance(o, D) else D(o)
        return D(self.v + o.v, self.d + o.d)
    __radd__ = __add__

    def __sub__(self, o):
        o = o if isinstance(o, D) else D(o)
        return D(self.v - o.v, self.d - o.d)

    def __rsub__(self, o):
        return D(o) - self

    def __mul__(self, o):
        o = o if isinstance(o, D) else D(o)
        return D(self.v * o.v, self.v * o.d + self.d * o.v)
    __rmul__ = __mul__

    def __truediv__(self, o):
        o = o if isinstance(o, D) else D(o)
        return D(self.v / o.v, (self.d * o.v - self.v * o.d) / (o.v * o.v))

    def __rtruediv__(self, o):
        return D(o) / self

    def __neg__(self):
        return D(-self.v, -self.d)

    def __pow__(self, n):
        assert isinstance(n, int) and n >= 1
        out = self
        for _ in range(n - 1):
            out = out * self
        return out


def dcos(x):
    return D(x.v.cos(), -x.v.sin() * x.d)


def dsin(x):
    return D(x.v.sin(), x.v.cos() * x.d)


def dre(x):
    return D(x.v.real, x.d.real if isinstance(x.d, acb) else arb(0))


def dim(x):
    return D(x.v.imag, x.d.imag if isinstance(x.d, acb) else arb(0))


# ------------------------------------------------------------------------------------- model
def Delta(r):
    return -(r ** 4) / 50 + r * r - 2 * r + Q2


def dDelta(r):
    return -4 * r ** 3 / 50 + 2 * r - 2


def F_val(r):
    """F = Delta/r^2 in factored form with the S1 root balls (tight on boxes; used for dr/ds and
    in the Krawczyk derivative f' = 1/F).  V_l itself uses the expanded polynomial below."""
    p = acb(-1) / 50
    for ri in ROOTS:
        p = p * (r - ri)
    return p / (r * r)


def V_dual(r, L):
    """V_l(r) = F (L/r^2 + F'/r) with r a dual number (dual arithmetic gives dV/ds)."""
    Dl = Delta(r)
    F = Dl / (r * r)
    Fp = (r * dDelta(r) - 2 * Dl) / (r * r * r)
    return F * (L / (r * r) + Fp / r)


class Phase:
    def __init__(self, phinf, s1):
        self.pi, self.s14 = arb(phinf), arb(s1) ** 4

    def all(self, s):
        """(phi, phi', phi'') as duals of s (closed forms; duals add one more derivative)."""
        s2 = s * s
        s4 = s2 * s2
        den = self.s14 + s4
        phi = self.pi * s4 / den
        p1 = 4 * self.pi * self.s14 * s2 * s / (den * den)
        p2 = 4 * self.pi * self.s14 * s2 * (3 * self.s14 - 5 * s4) / (den * den * den)
        return phi, p1, p2


# ------------------------------------------------------------------------------------- r(s) enclosure
class Contour:
    """Float path (design only: guesses and branch args) + Krawczyk enclosures (rigorous)."""

    def __init__(self):
        rs = [float(x.mid()) for x in ROOTS]
        e = np.exp(1j * np.pi / 4)
        F = lambda r: -(r - rs[0]) * (r - rs[1]) * (r - rs[2]) * (r - rs[3]) / (50 * r * r)
        f = lambda s, y: [(e * F(y[0] + 1j * y[1])).real, (e * F(y[0] + 1j * y[1])).imag]
        self.sol = {}
        self.grid = {}
        for side, S in ((-1, -64.0), (1, 32.0)):
            sol = solve_ivp(f, (0.0, S), [float(RP.mid()), 0.0], method="DOP853", rtol=1e-13, atol=1e-15,
                            dense_output=True)
            ss = np.linspace(0.0, S, int(abs(S) * 256) + 1)
            y = sol.sol(ss)
            r = y[0] + 1j * y[1]
            z = np.stack([r - rs[0], r - rs[1], r - rs[2], rs[3] - r])
            self.sol[side] = sol
            self.grid[side] = (ss, np.unwrap(np.angle(z), axis=1))
        self.rs = rs

    def guess(self, s):
        side = 1 if s >= 0 else -1
        y = self.sol[side].sol(s)
        r = y[0] + 1j * y[1]
        ss, args = self.grid[side]
        k = int(round(abs(s) * 256))
        z = np.array([r - self.rs[0], r - self.rs[1], r - self.rs[2], self.rs[3] - r])
        a0 = np.angle(z)
        al = a0 + 2 * np.pi * np.round((args[:, k] - a0) / (2 * np.pi))
        return r, al

    def enclose(self, s_lo, s_hi, depth=0):
        """Rigorous enclosure of r(s) for all s in [s_lo, s_hi] (Fractions; equal = point).
        Krawczyk with epsilon-inflation; if it does not contract, split the s-interval and take the
        hull of the two enclosures."""
        sm = (s_lo + s_hi) / 2
        r0, al = self.guess(float(sm))
        S = q_(sm) if s_lo == s_hi else (q_(s_lo) + q_(s_hi)) / 2 + arb(0, 1) * ((q_(s_hi) - q_(s_lo)) / 2)
        X = acb(XP) + EI4 * S
        sig = (1, 1, 1, -1)
        rot = [acb(arb(0), arb(-float(a))).exp() for a in al]
        m = acb(r0.real, r0.imag)

        def f(r):
            return sum((AB[i] * ((sig[i] * (r - ROOTS[i])) * rot[i]).log() + I_ * AB[i] * arb(float(al[i]))
                        for i in range(4)), acb(0)) - X

        def fp(r):                                   # f' = 1/F (factored: no cancellation)
            return 1 / F_val(r)
        Y = 1 / acb(complex(fp(m).mid()))
        fm = f(m)
        Yf = Y * fm
        rho = 1.2 * max(float(Yf.real.rad()), float(Yf.imag.rad())) + 1e-10   # float guess error ~1e-12
        for _ in range(30):
            R = acb(m.real + arb(0, 1) * rho, m.imag + arb(0, 1) * rho)
            if not all(bool(((sig[i] * (R - ROOTS[i])) * rot[i]).real > 0) for i in range(4)):
                break
            K = m - Yf + (1 - Y * fp(R)) * (R - m)
            if not (K.real.is_finite() and K.imag.is_finite()):
                break
            if (K.real.lower() > R.real.lower() and K.real.upper() < R.real.upper()
                    and K.imag.lower() > R.imag.lower() and K.imag.upper() < R.imag.upper()):
                return K
            ext = max(abs(float((K.real - m.real).lower())), abs(float((K.real - m.real).upper())),
                      abs(float((K.imag - m.imag).lower())), abs(float((K.imag - m.imag).upper())))
            new = 1.1 * ext
            if new > 1e3 * rho or new != new:
                break
            rho = max(new, 1.05 * rho)
        if s_lo == s_hi or depth >= 12:
            raise ArithmeticError(f"Krawczyk failed on [{float(s_lo)}, {float(s_hi)}]")
        a = self.enclose(s_lo, sm, depth + 1)
        b = self.enclose(sm, s_hi, depth + 1)
        return a.union(b)


# ------------------------------------------------------------------------------------- E on a cell
class Check:
    def __init__(self, l, design, nodes, contour):
        self.l, self.L = l, l * (l + 1)
        self.eta = arb(design["eta"])
        self.ph = Phase(design["phi_inf"], design["s1"])
        self.sgn_neg_phinf = design["phi_inf"] < 0
        self.nodes = nodes
        self.ct = contour
        wl = (l + arb(1) / 2) * OMP
        a1, a2, b1 = wl - 8 * G0, wl + 6 * G0, -BETA
        A1, A2, B1 = arb(a1.lower()), arb(a2.upper()), arb(b1.lower())
        Alo, Ahi = arb((A1 * A1 - B1 * B1).lower()), arb((A2 * A2).upper())
        Blo, Bhi = arb((2 * A2 * B1).lower()), arb(0)
        self.A = (Alo + Ahi) / 2 + arb(0, 1) * ((Ahi - Alo) / 2)
        self.B = (Blo + Bhi) / 2 + arb(0, 1) * ((Bhi - Blo) / 2)
        self.window = (a1, a2, b1)
        self.decay_ok = bool(a1 + b1 > 0)

    def E(self, s, r, rp, Ga, Gpp, sa, sgn):
        """Dual E at s (dual), r = dual (value, dr/ds); G' = Ga + Gpp (s - sa)."""
        phi, p1, p2 = self.ph.all(s)
        ap1 = p1 * sgn
        c = dcos(phi) - self.eta * ap1
        cp = -(dsin(phi) * p1) - self.eta * (p2 * sgn)
        V = V_dual(D(r, rp), self.L)
        V1, V2 = dre(V), dim(V)
        q = (V1 - self.A) * dsin(phi) - (V2 - self.B) * dcos(phi) - ap1 / (4 * self.eta)
        Gp = Ga + Gpp * (s - sa)
        E = q + cp * Gp - c * Gp * Gp + c * Gpp
        return E, c

    def cell(self, sa, sb, Ga, Gb):
        """(lower bound of E on [sa, sb], lower bound of c there)."""
        w = sb - sa
        sm = (sa + sb) / 2
        sgn = (-1 if sm > 0 else 1) * (1 if self.sgn_neg_phinf else -1)
        Gpp = (q_(Gb) - q_(Ga)) / q_(w)
        rm = self.ct.enclose(sm, sm)
        Em, _ = self.E(D(q_(sm), arb(1)), rm, EI4 * F_val(rm), q_(Ga), Gpp, q_(sa), sgn)
        S = (q_(sa) + q_(sb)) / 2 + arb(0, 1) * (q_(w) / 2)
        Rc = self.ct.enclose(sa, sb)
        Ec, cc = self.E(D(S, arb(1)), Rc, EI4 * F_val(Rc), q_(Ga), Gpp, q_(sa), sgn)
        low = Em.v.lower() - abs(Ec.d).upper() * q_(w) / 2
        return arb(low), arb(cc.v.lower())

    def cover(self):
        xs = [x for x, _ in self.nodes]
        cells = []

        def fill(a, b, h=Fraction(1, 16)):
            n = max(1, int(np.ceil(float((b - a) / h))))
            pts = [a + (b - a) * Fraction(k, n) for k in range(n + 1)]
            if a < 0 < b and 0 not in pts:
                pts = sorted(pts + [Fraction(0)])
            for u, v in zip(pts, pts[1:]):
                cells.append((u, v, Fraction(0), Fraction(0)))
        fill(Fraction(-64), xs[0])
        for (x0, y0), (x1, y1) in zip(self.nodes, self.nodes[1:]):
            if x0 < 0 < x1:                                   # never straddle s = 0 (sign of phi')
                y = y0 + (y1 - y0) * (0 - x0) / (x1 - x0)
                cells += [(x0, Fraction(0), y0, y), (Fraction(0), x1, y, y1)]
            else:
                cells.append((x0, x1, y0, y1))
        fill(xs[-1], Fraction(32))
        return cells

    def run(self, target, max_depth=18):
        mu, cmin, n, worst, below = None, None, 0, None, 0
        for (a, b, ya, yb) in self.cover():
            stack = [(a, b, ya, yb, 0)]
            while stack:
                a_, b_, ya_, yb_, d = stack.pop()
                low, cl = self.cell(a_, b_, ya_, yb_)
                n += 1
                if low >= target or d >= max_depth:
                    if mu is None or low < mu:
                        mu, worst = low, (float(a_), float(b_), d)
                    cmin = cl if cmin is None or cl < cmin else cmin
                    below += int(not (low >= target))
                    continue
                m, ym = (a_ + b_) / 2, (ya_ + yb_) / 2
                stack += [(m, b_, ym, yb_, d + 1), (a_, m, ya_, ym, d + 1)]
        return mu, cmin, n, worst, below

    def tails(self, dL, dR, design):
        """E = q on the tails (G' = 0): q >= -|V|_disc + min(-A sin phi + B cos phi) - |phi'(S)|/(4 eta)."""
        out = {}
        for name, S, rh, d in (("left", -64, ROOTS[2], dL), ("right", 32, ROOTS[3], dR)):
            # sup |V| on the square [r_h - d, r_h + d] x [-d, d] (contains the disc): the expanded
            # polynomial has large dependency on a wide box, so subdivide it into n x n sub-boxes
            n = 48
            Vmax = arb(0)
            for i in range(n):
                xr = rh - d + (2 * d) * (2 * i + 1) / (2 * n) + arb(0, 1) * (d / n)
                for j in range(n):
                    yi = -d + (2 * d) * (2 * j + 1) / (2 * n) + arb(0, 1) * (d / n)
                    Vb = V_dual(D(acb(xr, yi), acb(0)), self.L).v
                    Vmax = arb(Vmax.upper()).max(arb(abs(Vb).upper()))
            V = acb(arb(0, 1) * Vmax)                         # |V| <= Vmax on the whole square
            phS, p1S, _ = self.ph.all(D(arb(S), arb(1)))
            pinf = arb(design["phi_inf"])
            lo, hi = (pinf, phS.v) if pinf < phS.v else (phS.v, pinf)
            PH = (lo.lower() + hi.upper()) / 2 + arb(0, 1) * ((hi.upper() - lo.lower()) / 2)
            term = -(self.A * PH.sin()) + self.B * PH.cos()
            qlow = term.lower() - abs(V).upper() - abs(p1S.v).upper() / (4 * self.eta)
            c_low = PH.cos().lower() - self.eta * abs(p1S.v).upper()
            out[name] = {"q_lower": arb(qlow).str(8), "c_lower": arb(c_low).str(6), "q": arb(qlow), "c": arb(c_low)}
        return out


def tail_disc(ih, d0, S):
    rh, Ah = ROOTS[ih], AB[ih]
    smin, smax = arb(0), arb(0)
    for i in range(4):
        if i == ih:
            continue
        dist = abs(rh - ROOTS[i])
        t1, t2 = AB[i] * (dist - d0).log(), AB[i] * (dist + d0).log()
        smin += arb(min(t1.lower(), t2.lower()))
        smax += arb(max(t1.upper(), t2.upper()))
    xS = XP + arb(S) * arb(2).sqrt() / 2
    d = ((xS - smin) / Ah).exp() if ih == 2 else ((xS - smax) / Ah).exp()
    return arb(d.upper())


def main(ls):
    t0 = time.time()
    ct = Contour()
    dL, dR = tail_disc(2, arb("0.06"), -64), tail_disc(3, arb("0.6"), 32)
    # starting points of the first-exit argument lie in the discs of radius d0
    sL = abs(ct.enclose(Fraction(-64), Fraction(-64)) - acb(ROOTS[2]))
    sR = abs(ct.enclose(Fraction(32), Fraction(32)) - acb(ROOTS[3]))
    res = {"tail_discs": {"d_left": dL.str(8), "d_right": dR.str(8), "start_left_in_0.06": bool(sL < arb("0.06")),
                          "start_right_in_0.6": bool(sR < arb("0.6"))}}
    print(res, flush=True)
    for l in ls:
        st = json.loads((P2 / "certificates" / "m3" / f"m3_PB_l{l}.json").read_text())
        nodes = [(Fraction(x), Fraction(y)) for x, y in st["ground_state_nodes"]]
        mu_auth = arb(st["mu_cert"])
        chk = Check(l, st["design"], nodes, ct)
        t1 = time.time()
        mu, cmin, n, worst, below = chk.run(arb(mu_auth.lower()))
        tl = chk.tails(dL, dR, st["design"])
        mu_all = min([mu, tl["left"]["q"], tl["right"]["q"]], key=lambda z: z.lower())
        muHO = 2 * (l + arb(1) / 2) * OMP * G0
        row = {"l": l, "cells_evaluated": n, "cells_below_author_mu_at_max_depth": below,
               "mu_reviewer_cells_lower": mu.str(10), "worst_cell": worst,
               "mu_author_cert": st["mu_cert"], "reviewer_ge_author": bool(mu >= arb(mu_auth.lower())),
               "tails": {k: {kk: vv for kk, vv in v.items() if kk in ("q_lower", "c_lower")} for k, v in tl.items()},
               "mu_reviewer_all": arb(mu_all.lower()).str(10), "m_ratio_reviewer": (arb(mu_all.lower()) / muHO).str(8),
               "c_min_cells_lower": cmin.str(8), "c_nonneg": bool(cmin > 0 and tl["left"]["c"] > 0 and tl["right"]["c"] > 0),
               "decay_condition_a1_plus_b1_gt_0": chk.decay_ok,
               "window_Re_over_kappa_minus": [(chk.window[0] / KM).str(8), (chk.window[1] / KM).str(8)],
               "hull_Re_w2": chk.A.str(8), "hull_Im_w2": chk.B.str(8),
               "certified_positive": bool(mu_all > 0 and cmin > 0 and chk.decay_ok),
               "seconds": round(time.time() - t1, 1)}
        res[f"l{l}"] = row
        print(json.dumps(row), flush=True)
    res["seconds_total"] = round(time.time() - t0, 1)
    name = "indep_cells_arb" + ("" if ls == [20, 50, 100] else "_" + "_".join(map(str, ls)))
    (HERE / f"{name}.json").write_text(json.dumps(res, indent=1) + "\n")


if __name__ == "__main__":
    main([int(x) for x in sys.argv[1:]] or [20, 50, 100])
