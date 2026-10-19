# ruff: noqa: N999  (module name fixed by the S4 pre-registration, section 6)
"""S4 regimes T (barrier top) and A (above the barrier): global M3 on S2's contour, uniformly in h.

Pre-registration: notes/s4-preregistration.md §3.3, §3.4, §5 and Amendment 4.  Proof: notes/s4-proof.md.

On S2's validated contour x = x_P + e^{i pi/4} s, with h = 1/(l + 1/2), omega = zeta/h + i nu and the
exact splitting V_l = U/h^2 + Wt (U = F/r^2, Wt = F F_r/r - U/4), a phase phi(s), eta = 1/2,
c = cos(phi) - eta |phi'| and the semiclassical ground state G' = Gamma/h, S2's
E = q + c'G' - cG'^2 + cG'' satisfies exactly

    h^2 E = E0(s; zeta) + h E1(s; zeta, nu) + h^2 E2(s; nu),
    E0 = Re(i e^{-i phi}(U - zeta^2)) - c Gamma^2,
    E1 = 2 zeta nu cos(phi) + (c Gamma)',
    E2 = Re(i e^{-i phi}(Wt + nu^2)) - |phi'|/(4 eta).

If E0 >= 0 and F := E0/h0 + E1 + h0 min(0, E2) >= mu > 0 on the whole contour, then E >= mu/h for
every h in (0, h0] (E0/h >= E0/h0 and h E2 >= h0 min(0, E2)).  On a box zeta in [zl, zh],
nu in [nl, nh] (nh <= 0): E0 is nondecreasing in zeta^2 (sin phi <= 0), E1 is nonincreasing in zeta
and nondecreasing in nu (cos phi > 0), E2 is nonincreasing in nu^2.  So the box bound is exact:
E0 at zl, E1 at (zh, nl), E2 at nl.

Per s-cell [a, b] (a, b exact dyadics inside one contour step), lower bounds by the mean-value form
f(s_m) - sup_cell |f'| w/2 and, for E0, also the second-order form
f(s_m) - |f'(s_m)| w/2 - sup_cell |f''| w^2/8.  E0(0) = 0 exactly; on cells with an endpoint at
s = 0, E0 = s^2 E0~ with  U - U_P in s^2/2 conv{U_ss([0, s])}  (U_s(0) = 0), so E0 >= 0 follows from
E0~ >= 0 on the cell.  Tails |s| >= |S| (Gamma = 0) by S2's horizon discs.

T (Amendment 4): phi = phi_inf s^4/(s1^4 + s^4), Gamma = (1 - eps_g) a_d s theta(s).
A: phi = phi0 + (phi_inf - phi0) s^4/(s1^4 + s^4), phi0 = -1/5, Gamma = 0, zeta = zeta_A only; then
monotone in Re omega (Lemma W7) because phi <= phi0 < arg omega.
"""
from __future__ import annotations

import json
import os
import sys
import time
from fractions import Fraction
from pathlib import Path

from flint import acb, arb, arb_series

from .arbutil import Timer, git_commit, precision
from .m3 import S_LEFT, S_RIGHT, fa, geometry_certificate, theta_unit
from .params import Point

PREC = 128
H0 = Fraction(2, 41)                       # h0 = 1/20.5 (l0 = 20)
STEP = Fraction(1, 16)                     # S2 contour step = base cell
ROOT = Path(__file__).resolve().parents[2]
CERT_DIR = ROOT / "certificates" / "part1"
RES_DIR = ROOT / "results" / "part1"

DESIGN_T = {
    "name": "T", "s1": "7", "p": 4, "phi_inf": "-23/16", "phi0": "0", "eta": "1/2",
    "eps_g": "1/20", "a_d": "32178552/1073741824", "taper": ["2", "4"],
    "zeta_over_OmegaP": ["4/5", "2"], "zeta_box_over_OmegaP": "1/20",
    "max_depth": 12, "target_over_a_star": "1/10",
}
DESIGN_A = {
    "name": "A", "s1": "7", "p": 4, "phi_inf": "-23/16", "phi0": "-1/5", "eta": "1/2",
    "eps_g": None, "a_d": "0", "taper": ["2", "4"],
    "zeta_over_OmegaP": ["3/2", "3/2"], "zeta_box_over_OmegaP": "0",
    "max_depth": 12, "target_over_a_star": "1/10",
}
DEPTH_CHARS = "0123456789abcdefghij"


def Fr(x) -> Fraction:
    return Fraction(x)


def pot_parts(pt: Point, r):
    """U, U_r, U_xx, F, Wt, Wt_r at an acb r (ball or box); Delta in factored form."""
    lam3 = pt.Lam / 3
    D = acb(-lam3)
    for ri in pt.roots:
        D = D * (r - ri)
    Dp = -4 * lam3 * r ** 3 + 2 * r - 2
    Dpp = -12 * lam3 * r * r + 2
    r2 = r * r
    r3 = r2 * r
    r4 = r2 * r2
    F = D / r2
    Fr_ = (Dp * r - 2 * D) / r3
    Frr = (Dpp * r2 - 4 * Dp * r + 6 * D) / r4
    U = F / r2
    Ur = (Dp * r - 4 * D) / (r4 * r)
    Urr = (Dpp * r2 - 8 * Dp * r + 20 * D) / (r4 * r2)
    W = F * Fr_ / r - U / 4
    Wr = (Fr_ * Fr_ + F * Frr) / r - F * Fr_ / r2 - Ur / 4
    Uxx = F * (Fr_ * Ur + F * Urr)
    return {"U": U, "Ur": Ur, "Uxx": Uxx, "F": F, "W": W, "Wr": Wr}


class Design:
    """Exact design constants as arbs."""

    def __init__(self, pt: Point, d: dict):
        self.d = d
        self.s1 = fa(Fr(d["s1"]))
        self.s14 = self.s1 ** 4
        self.phinf = fa(Fr(d["phi_inf"]))
        self.phi0 = fa(Fr(d["phi0"]))
        self.eta = fa(Fr(d["eta"]))
        self.gs = d["eps_g"] is not None
        self.k = fa((1 - Fr(d["eps_g"])) * Fr(d["a_d"])) if self.gs else arb(0)
        self.tc, self.td = (fa(Fr(x)) for x in d["taper"])
        self.tc_f, self.td_f = (Fr(x) for x in d["taper"])
        self.a_d = fa(Fr(d["a_d"])) if self.gs else None

    def phase(self, s):
        """(phi, phi', phi'', phi''') at s (arb point or interval)."""
        x = arb_series([s, 1], prec=4)
        x4 = x ** 4
        f = self.phi0 + (self.phinf - self.phi0) * x4 / (self.s14 + x4)
        c = f.coeffs()
        c = [c[k] if k < len(c) else arb(0) for k in range(4)]
        return c[0], c[1], 2 * c[2], 6 * c[3]

    def zone(self, a: Fraction, b: Fraction) -> int:
        """0: |s| <= tc, 1: tc <= |s| <= td, 2: |s| >= td (cells never straddle the breakpoints)."""
        lo_, hi_ = (a, b) if a >= 0 else (-b, -a)
        if hi_ <= self.tc_f:
            return 0
        if lo_ >= self.td_f:
            return 2
        if lo_ >= self.tc_f and hi_ <= self.td_f:
            return 1
        raise ValueError(f"cell [{a}, {b}] straddles a taper breakpoint")

    def gamma(self, s, zone: int, side: int):
        """(Gamma, Gamma', Gamma'') on a cell of the given zone and side (sign of s)."""
        if not self.gs or zone == 2:
            z = arb(0)
            return z, z, z
        k = self.k
        if zone == 0:
            return k * s, k, arb(0)
        pi = arb.pi()
        u = pi * (side * s - self.tc) / (self.td - self.tc)
        w = self.td - self.tc
        su, cu = u.sin(), u.cos()
        th = (1 + cu) / 2
        thp = -(pi / (2 * w)) * su * side
        thpp = -(pi * pi / (2 * w * w)) * cu
        return k * s * th, k * (th + s * thp), k * (2 * thp + s * thpp)


class TChecker:
    def __init__(self, pt: Point, rays, tails, design: dict):
        self.pt, self.rays, self.tails = pt, rays, tails
        self.D = Design(pt, design)
        self.e4 = theta_unit()                 # e^{i pi/4}
        self.I = acb(0, 1)
        self.UP = pot_parts(pt, acb(pt.rP))["U"].real

    # -- E0, E1, E2 and derivatives ----------------------------------------------------------
    def parts(self, s, r, zl, zh, nl, sgn, zone, side):
        D, I = self.D, self.I
        phi, d1, d2, d3 = D.phase(s)
        eta = D.eta
        ad1 = d1 * sgn
        cphi, sphi = phi.cos(), phi.sin()
        c = cphi - eta * ad1
        cp = -sphi * d1 - eta * sgn * d2
        cpp = -cphi * d1 * d1 - sphi * d2 - eta * sgn * d3
        G, Gp, Gpp = D.gamma(s, zone, side)
        P = pot_parts(self.pt, r)
        rot = acb(sphi, cphi)                  # i e^{-i phi}
        Us = self.e4 * P["F"] * P["Ur"]
        Uss = I * P["Uxx"]
        A0 = P["U"] - zl * zl
        E0 = (rot * A0).real - c * G * G
        E0p = (rot * (Us - I * d1 * A0)).real - cp * G * G - 2 * c * G * Gp
        E0pp = ((rot * (Uss - 2 * I * d1 * Us - I * d2 * A0 - d1 * d1 * A0)).real
                - cpp * G * G - 4 * cp * G * Gp - 2 * c * Gp * Gp - 2 * c * G * Gpp)
        E1 = 2 * zh * nl * cphi + cp * G + c * Gp
        E1p = -2 * zh * nl * sphi * d1 + cpp * G + 2 * cp * Gp + c * Gpp
        A2 = P["W"] + nl * nl
        Ws = self.e4 * P["F"] * P["Wr"]
        E2 = (rot * A2).real - ad1 / (4 * eta)
        E2p = (rot * (Ws - I * d1 * A2)).real - sgn * d2 / (4 * eta)
        return {"E0": E0, "E0p": E0p, "E0pp": E0pp, "E1": E1, "E1p": E1p, "E2": E2, "E2p": E2p, "c": c}

    def E0_tilde(self, cell, R_hull, zl):
        """Lower bound of E0/s^2 on a cell with an endpoint at 0 (T design: phi0 = 0, theta = 1)."""
        D = self.D
        phi, d1, _, _ = D.phase(cell)
        sgn = -1 if cell.upper() <= 0 else 1
        sgn = -sgn                              # sign(phi') = -sign(s) for phi_inf < phi0
        c = phi.cos() - D.eta * d1 * sgn
        M = (self.I * pot_parts(self.pt, R_hull)["Uxx"]) / 2     # (U - U_P)/s^2 in 1/2 conv U_ss
        rot = acb(phi.sin(), phi.cos())
        s2 = cell * cell
        sin_over_s2 = phi.sinc() * (D.phinf - D.phi0) * s2 / (D.s14 + s2 * s2)
        Et = (rot * M).real + (self.UP - zl * zl) * sin_over_s2 - c * D.k * D.k
        return Et, c

    # -- one cell ---------------------------------------------------------------------------------
    def check_cell(self, a: Fraction, b: Fraction, box, h0):
        """Certified lower bounds on [a, b]: (F_low, E0_low, c_low)."""
        zl, zh, nl = box
        D = self.D
        side = 1 if (a + b) > 0 else -1
        ray = self.rays[side]
        sgn = -side if D.phinf < D.phi0 else side
        zone = D.zone(a, b)
        a_, b_ = fa(a), fa(b)
        w = b_ - a_
        sm = (a_ + b_) / 2
        cell = sm + arb(0, 1) * (w / 2)
        pm = self.parts(sm, ray.at(sm), zl, zh, nl, sgn, zone, side)
        pc = self.parts(cell, ray.enclose(a_, b_), zl, zh, nl, sgn, zone, side)
        hw = w / 2
        E1_low = pm["E1"].lower() - abs(pc["E1p"]).upper() * hw
        E2_low = pm["E2"].lower() - abs(pc["E2p"]).upper() * hw
        c_low = pc["c"].lower()
        if (a == 0 or b == 0) and D.gs:
            R_hull = ray.enclose(a_, b_)        # the cell is [0, b] or [a, 0]: hull of [0, s] inside
            Et, c_ = self.E0_tilde(cell, R_hull, zl)
            E0_low = arb(0) if Et.lower() >= 0 else arb(-1)
            c_low = arb(c_low).min(c_.lower())
        else:
            e1 = pm["E0"].lower() - abs(pc["E0p"]).upper() * hw
            e2 = pm["E0"].lower() - abs(pm["E0p"]).upper() * hw - abs(pc["E0pp"]).upper() * (w * w / 8)
            E0_low = arb(e1).max(arb(e2))
        h0_ = fa(h0)
        F_low = E0_low / h0_ + E1_low + h0_ * arb(E2_low).min(arb(0))
        return arb(F_low.lower()), arb(E0_low.lower()), arb(c_low)

    # -- tails -----------------------------------------------------------------------------------
    def tail(self, side, box, h0):
        """Lower bounds on |s| >= |S| (Gamma = 0): (F_low, E0_low, c_low)."""
        zl, zh, nl = box
        pt, D = self.pt, self.D
        S = S_LEFT if side < 0 else S_RIGHT
        d = self.tails[0] if side < 0 else self.tails[1]
        rh = pt.rp if side < 0 else pt.rc
        disc = acb(rh + arb(0, 1) * d, arb(0, 1) * d)
        P = pot_parts(pt, disc)
        phS, d1S, _, _ = D.phase(arb(S))
        # beyond S, |phi| increases from |phi(S)| to |phi_inf| < pi/2 and |phi'| decreases
        sin_abs_min = abs(phS.sin()).lower()
        cos_max = phS.cos().upper()
        cos_min = D.phinf.cos().lower()
        E0_low = zl * zl * arb(sin_abs_min) - abs(P["U"]).upper()
        E1_low = 2 * zh * nl * arb(cos_max)
        E2_low = -abs(P["W"]).upper() - nl * nl - abs(d1S).upper() / (4 * D.eta)
        c_low = arb(cos_min) - D.eta * abs(d1S).upper()
        h0_ = fa(h0)
        F_low = E0_low / h0_ + E1_low + h0_ * arb(E2_low.lower()).min(arb(0))
        return arb(F_low.lower()), arb(E0_low.lower()), arb(c_low.lower())

    # -- a (zeta, nu) box --------------------------------------------------------------------------
    def base_cells(self):
        out = []
        n_left, n_right = int(-Fraction(S_LEFT) / STEP), int(Fraction(S_RIGHT) / STEP)
        for k in range(n_left, 0, -1):
            out.append((-k * STEP, -(k - 1) * STEP))
        for k in range(n_right):
            out.append((k * STEP, (k + 1) * STEP))
        return out

    def cells_from_leaves(self, leaves: str):
        """The partition encoded by a depth string (leaves in left-first order within base cells)."""
        ds = [DEPTH_CHARS.index(ch) for ch in leaves]
        pos = 0
        out = []

        def build(a, b, d):
            nonlocal pos
            if ds[pos] == d:
                out.append((a, b, d))
                pos += 1
                return
            m = (a + b) / 2
            build(a, m, d + 1)
            build(m, b, d + 1)

        for (a, b) in self.base_cells():
            build(a, b, 0)
        if pos != len(ds):
            raise ValueError("leaf string does not match the base cells")
        return out

    def run_box(self, box, h0, target, max_depth, leaves=None, sample=None):
        """Certify one box.  Production: bisect cells that miss the target (or E0 >= 0, c >= 0) up to
        max_depth.  Replay: evaluate the stored partition `leaves` (only cells with index in `sample`
        if given)."""
        mu, e0min, cmin = None, None, None
        n_eval = 0
        worst = None
        fails = []
        depth_str = []

        def record(a_, b_, dpt, F_low, E0_low, c_low):
            nonlocal mu, e0min, cmin, worst
            depth_str.append(DEPTH_CHARS[dpt])
            if not (F_low > 0 and E0_low >= 0 and c_low >= 0):
                fails.append([str(a_), str(b_), F_low.str(5), E0_low.str(5), c_low.str(5)])
            if mu is None or F_low < mu:
                mu, worst = F_low, (float(a_), float(b_))
            if e0min is None or E0_low < e0min:
                e0min = E0_low
            if cmin is None or c_low < cmin:
                cmin = c_low

        if leaves is not None:
            for i, (a_, b_, dpt) in enumerate(self.cells_from_leaves(leaves)):
                if sample is not None and i not in sample:
                    continue
                F_low, E0_low, c_low = self.check_cell(a_, b_, box, h0)
                n_eval += 1
                record(a_, b_, dpt, F_low, E0_low, c_low)
        else:
            for (a, b) in self.base_cells():
                stack = [(a, b, 0)]
                while stack:
                    a_, b_, dpt = stack.pop()
                    F_low, E0_low, c_low = self.check_cell(a_, b_, box, h0)
                    n_eval += 1
                    ok = bool(F_low >= target and E0_low >= 0 and c_low >= 0)
                    if not ok and dpt < max_depth:
                        m = (a_ + b_) / 2
                        stack.append((m, b_, dpt + 1))
                        stack.append((a_, m, dpt + 1))
                        continue
                    record(a_, b_, dpt, F_low, E0_low, c_low)
        tl = self.tail(-1, box, h0)
        tr = self.tail(1, box, h0)
        return {"mu_cells": mu, "E0_low_min": e0min, "c_min": cmin, "cells": n_eval, "worst_cell": worst,
                "fails": fails, "leaves": "".join(depth_str), "tail_left": tl, "tail_right": tr}


# ---------------------------------------------------------------------------------------------
# boxes
# ---------------------------------------------------------------------------------------------
def zeta_boxes(design: dict):
    lo_, hi_ = (Fr(x) for x in design["zeta_over_OmegaP"])
    w = Fr(design["zeta_box_over_OmegaP"])
    if w == 0:
        return [(lo_, hi_)]
    n = int((hi_ - lo_) / w)
    assert lo_ + n * w == hi_
    return [(lo_ + k * w, lo_ + (k + 1) * w) for k in range(n)]


def box_arbs(pt: Point, zb, nl=None):
    """(zl, zh, nl) arbs for zeta in [zb0, zb1] Omega_P and nu in [nl, 0] (default nl = -d_delta)."""
    zl = fa(zb[0]) * pt.OmegaP
    zh = fa(zb[1]) * pt.OmegaP
    if nl is None:
        nl = -pt.beta
    # endpoint balls: zl may be taken smaller and zh, |nl| larger (bounds are monotone)
    return arb(zl.lower()), arb(zh.upper()), arb(nl.lower())


def setup(pt: Point):
    geo, rays, tails = geometry_certificate(pt)
    return geo, rays, tails


def a_star(pt: Point):
    """a* = sqrt(|U''|/2) = 2 Omega_P gamma0 (the oscillator constant)."""
    return 2 * pt.OmegaP * pt.gamma0


def certify_box(ch: TChecker, pt: Point, zb, design: dict, h0=H0, nl=None, leaves=None, sample=None):
    box = box_arbs(pt, zb, nl)
    target = fa(Fr(design["target_over_a_star"])) * a_star(pt)
    with Timer() as tm:
        res = ch.run_box(box, h0, arb(target.lower()), design["max_depth"], leaves=leaves, sample=sample)
    mu_all = arb(res["mu_cells"]).min(res["tail_left"][0]).min(res["tail_right"][0])
    e0_all = arb(res["E0_low_min"]).min(res["tail_left"][1]).min(res["tail_right"][1])
    c_all = arb(res["c_min"]).min(res["tail_left"][2]).min(res["tail_right"][2])
    zl, _zh, nl_ = box
    # Lemma 3 decay on the contour tails: Re omega + Im omega > 0 for all h <= h0
    decay = bool(zl / fa(h0) + nl_ > 0)
    ok = bool(mu_all > 0 and e0_all >= 0 and c_all >= 0 and not res["fails"] and decay)
    return {"zeta_over_OmegaP": [str(zb[0]), str(zb[1])], "nu_low": nl_.str(10), "h0": str(h0),
            "mu_cert": mu_all.str(8), "mu_over_a_star": (mu_all / a_star(pt)).str(5),
            "E0_low_min": e0_all.str(5), "c_min": c_all.str(6),
            "tail_left": [x.str(5) for x in res["tail_left"]], "tail_right": [x.str(5) for x in res["tail_right"]],
            "cells": res["cells"], "worst_cell": res["worst_cell"], "fails": res["fails"][:20],
            "n_fails": len(res["fails"]), "jost_decay": decay, "pass": ok, "leaves": res["leaves"],
            "timing_s": {"wall": round(tm.wall, 1), "cpu": round(tm.cpu, 1)}}


def certify_box_abs(ch: TChecker, pt: Point, zl, zh, nl, design: dict, h0=H0):
    """The T/A checker on an absolute box zeta in [zl, zh], nu in [nl, 0] (arbs), h in (0, h0] (controls)."""
    box = (arb(zl.lower()), arb(zh.upper()), arb(nl.lower()))
    target = fa(Fr(design["target_over_a_star"])) * a_star(pt)
    res = ch.run_box(box, h0, arb(target.lower()), design["max_depth"])
    mu_all = arb(res["mu_cells"]).min(res["tail_left"][0]).min(res["tail_right"][0])
    e0_all = arb(res["E0_low_min"]).min(res["tail_left"][1]).min(res["tail_right"][1])
    ok = bool(mu_all > 0 and e0_all >= 0 and not res["fails"])
    return {"zeta": [box[0].str(10), box[1].str(10)], "nu_low": box[2].str(10), "h0": str(h0),
            "mu_cert_lower": float(arb(mu_all.lower())), "mu_over_a_star": float((mu_all / a_star(pt)).mid()),
            "worst_cell": res["worst_cell"], "n_fails": len(res["fails"]), "cells": res["cells"], "certified": ok}


def monotonicity_A(pt: Point, design: dict, h0=H0):
    """Lemma W7: phi(s) <= phi0 < -arctan(d_delta h0/zeta_A) <= arg omega on A, and phi > -pi."""
    zA = fa(Fr(design["zeta_over_OmegaP"][0])) * pt.OmegaP
    bound = -((pt.beta * fa(h0)) / zA).atan()
    phi0, phinf = fa(Fr(design["phi0"])), fa(Fr(design["phi_inf"]))
    return {"phi0": design["phi0"], "minus_arctan_bound": bound.str(10),
            "phi_le_phi0_by_construction": bool(phinf < phi0),
            "phi0_lt_bound": bool(phi0 < bound), "phi_inf_gt_minus_pi": bool(phinf > -arb.pi()),
            "pass": bool(phinf < phi0 < bound and phinf > -arb.pi())}


# ---------------------------------------------------------------------------------------------
# resumable runs
# ---------------------------------------------------------------------------------------------
def run_regime(name: str, worker: int = 0, nworkers: int = 1):
    design = DESIGN_T if name == "T" else DESIGN_A
    RES_DIR.mkdir(parents=True, exist_ok=True)
    out = RES_DIR / f"{name}_PB.w{worker}.jsonl"
    done = set()
    if out.exists():
        for line in out.read_text().splitlines():
            if line.strip():
                done.add(tuple(json.loads(line)["zeta_over_OmegaP"]))
    pt = Point("PB", 256)
    with precision(PREC):
        _, rays, tails = setup(pt)
        ch = TChecker(pt, rays, tails, design)
        for k, zb in enumerate(zeta_boxes(design)):
            if k % nworkers != worker or (str(zb[0]), str(zb[1])) in done:
                continue
            row = certify_box(ch, pt, zb, design)
            row["commit"] = git_commit()
            row["finished"] = time.strftime("%Y-%m-%dT%H:%M:%S")
            with out.open("a") as f:
                f.write(json.dumps(row) + "\n")
            print(name, zb, row["pass"], row["mu_over_a_star"], row["cells"], row["timing_s"], flush=True)


def assemble(name: str):
    """Collect the worker JSONL files into certificates/part1/{name}_PB.json."""
    design = DESIGN_T if name == "T" else DESIGN_A
    rows = {}
    for f in sorted(RES_DIR.glob(f"{name}_PB.w*.jsonl")):
        for line in f.read_text().splitlines():
            if line.strip():
                r = json.loads(line)
                rows[tuple(r["zeta_over_OmegaP"])] = r
    boxes = [(str(a), str(b)) for a, b in zeta_boxes(design)]
    missing = [b for b in boxes if b not in rows]
    pt = Point("PB", 256)
    with precision(PREC):
        geo, _, _ = setup(pt)
        extra = {}
        if name == "A":
            extra["monotonicity_W7"] = monotonicity_A(pt, design)
    ok = (not missing) and all(rows[b]["pass"] for b in boxes) and geo["A0_pass"] and \
        (extra.get("monotonicity_W7", {"pass": True})["pass"])
    mus = [rows[b]["mu_cert"] for b in boxes if b in rows]
    cert = {"regime": name, "point": "PB", "h0": str(H0), "precision_bits": PREC, "design": design,
            "geometry_A0_pass": geo["A0_pass"], "boxes": [rows[b] for b in boxes if b in rows],
            "missing_boxes": missing, "mu_min": min(mus, key=lambda s: float(arb(s.split("+/-")[0].strip("[ ")).mid())) if mus else None,
            "cells_total": sum(rows[b]["cells"] for b in boxes if b in rows),
            "cpu_s_total": round(sum(rows[b]["timing_s"]["cpu"] for b in boxes if b in rows), 1),
            "commits": sorted({rows[b]["commit"] for b in boxes if b in rows}), "pass": ok, **extra}
    CERT_DIR.mkdir(parents=True, exist_ok=True)
    (CERT_DIR / f"{name}_PB.json").write_text(json.dumps(cert, indent=1) + "\n")
    return cert


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "run":
        run_regime(sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 0,
                   int(sys.argv[4]) if len(sys.argv) > 4 else 1)
    elif cmd == "assemble":
        c = assemble(sys.argv[2])
        print(c["regime"], c["pass"], c["mu_min"], c["missing_boxes"], c["cpu_s_total"])
    else:
        raise SystemExit(f"unknown command {cmd} (pid {os.getpid()})")
