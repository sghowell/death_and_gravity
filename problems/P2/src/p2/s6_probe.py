"""S6 Phase-1 probes for P2-6 (the two-sided beta at P_B).  Pre-registration: notes/s6-preregistration.md.

These runs size the work and test the foundations.  They are probes, not certificates: their
outputs go to results/s6/ and nothing here is cited at a ledger level.

  E   the E_+ point omega0 = -i kappa_+ (inside the deeper strip):
      * sympy: at omega0 the Z_+ recursion's first denominator D_1 = Delta'(r_+) - 2 i w r_+^2
        vanishes, and the first coefficient C_0 = -(2 i w r_+ + L) = -(2 kappa_+ r_+ + L) does not, so
        the exponent-0 solution has a log term for every l;
      * Arb: Wt_l := W(Z2, Z_c) at omega0, Z2 the exponent-1 Frobenius solution (Z2 = t + ...), the
        Jost-normalised Wronskian form of zsol; Res_{omega0} W_l = i (2 kappa_+ r_+ + L)/(2 r_+^2) Wt_l;
      * cross-check: (w - omega0) W_l(w) at w = omega0 + eps from zsol, against that residue;
      * the bound behind the written argument for l >= 1: L + r F' > 0 on [r_+, r_c] for L >= 2.
  M3  S2's barrier-top windows (route M3) at depth e kappa_- for l = 20, 50, 100.
  TA  S4's regimes T (24 zeta-boxes) and A (zeta_A = 3/2) at depth e, uniformly in h <= 2/41;
      and A at fixed h = 1/(l + 1/2) for l <= 19 (the part-(ii) replacement), scanning zeta_A.
  S   S4's regime S (both sides, all zeta-boxes and the closing box) at depth e.
  L   S4's regime L on the band boxes [0, 8] x [-e', -11/20] at depth e, with a pole-free series
      R~ = (1 - i w/kappa_+) R_+ in the boxes that touch omega0.

Run (cwd problems/P2, PYTHONPATH=src):  python -m p2.s6_probe {E|M3|TA|AFIX|S+|Sc|L} [args]
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from pathlib import Path

from flint import acb, arb, ctx

from .arbutil import A, Timer, git_commit, precision
from .params import Point

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "s6"
PREC = 128
EDGE = Fraction(23, 25)                 # e = 0.92 (FORMULATION P2-6)
I = acb(0, 1)


def deep_point(depth: Fraction = EDGE) -> Point:
    """Point P_B with the strip depth d = depth * kappa_- (every S2/S4 regime reads pt.beta)."""
    pt = Point("PB", 256)
    pt.beta = A(depth) * pt.km
    return pt


def _write(name: str, obj) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / name).write_text(json.dumps(obj, indent=1) + "\n")


# ---------------------------------------------------------------------------------------------
# E: the point omega0 = -i kappa_+
# ---------------------------------------------------------------------------------------------
def e_sympy() -> dict:
    """The first Frobenius step of the Z_+ equation at r_+ (zsol's equation, s = +1)."""
    import sympy as sp
    t, w, L, rp, d1, d2, d3, d4 = sp.symbols("t w L r_p d1 d2 d3 d4")
    a1 = sp.symbols("a1")
    Delta = d1 * t + d2 * t ** 2 + d3 * t ** 3 + d4 * t ** 4          # Delta(r_+ + t), Delta(r_+) = 0
    r = rp + t
    Z = 1 + a1 * t
    expr = sp.expand(Delta * sp.diff(Z, t, 2) + (sp.diff(Delta, t) - 2 * sp.I * w * r ** 2) * sp.diff(Z, t)
                     - (2 * sp.I * w * r + L) * Z)
    c0 = sp.simplify(expr.coeff(t, 0))
    sol = sp.solve(sp.Eq(c0, 0), a1)[0]
    D1 = d1 - 2 * sp.I * w * rp ** 2
    C0 = -(2 * sp.I * w * rp + L)
    ok_form = sp.simplify(sol - (-C0 / D1)) == 0
    kp = sp.symbols("kappa_p", positive=True)
    w0 = -sp.I * kp
    sub = {d1: 2 * kp * rp ** 2, w: w0}
    D1_0 = sp.simplify(D1.subs(sub))
    C0_0 = sp.simplify(C0.subs(sub))
    return {"t0_coefficient": str(c0), "a1_equals_minus_C0_over_D1": bool(ok_form),
            "D1_at_omega0": str(D1_0), "C0_at_omega0": str(C0_0),
            "C0_at_omega0_is_minus_(2 kappa r + L)": bool(sp.simplify(C0_0 + 2 * kp * rp + L) == 0),
            "pass": bool(ok_form and D1_0 == 0 and sp.simplify(C0_0 + 2 * kp * rp + L) == 0)}


def frob_z2(pt: Point, w, L, t, min_terms=20):
    """(Z, Z') at r_+ + t (t exact) for the exponent-1 solution Z2 = t + sum_{m>=2} a_m t^m of the
    Z_+ equation at omega0 = -i kappa_+ (zsol's recurrence, started at a_0 = 0, a_1 = 1).  At omega0,
    D_1 = 0 exactly and D_m = Delta'(r_+) m (m - 1) for m >= 2; w is the ball for omega0."""
    from .arbw import _geom_tails, _solve_q
    from .zsol import N_MAX, _absu, zpolys
    rh = pt.rp
    Ap, B, C = zpolys(pt, rh, acb(w), L, 1, singular=True)
    A1 = Ap[1]
    ta = arb(t)
    a = [acb(0), acb(1)]
    S = acb(0) + acb(ta)
    Sd = acb(1)
    tk = arb(ta)
    small = 0
    eps = arb(2) ** (-(ctx.prec - 8))
    for m in range(2, N_MAX):
        acc = acb(0)
        for n in range(max(0, m - 3), m):
            j2, j1, j0 = m + 1 - n, m - n, m - 1 - n
            term = acb(0)
            if j2 < len(Ap):
                term += Ap[j2] * (n * (n - 1))
            if j1 < len(B):
                term += B[j1] * n
            if 0 <= j0 < len(C):
                term += C[j0]
            acc += a[n] * term
        Dm = A1 * (m * (m - 1))                       # exact indicial factor at omega0
        am = -acc / Dm
        a.append(am)
        Sd += am * m * tk
        tk = tk * ta
        tv = am * tk
        S += tv
        if m >= min_terms and _absu(tv) < eps * abs(S).lower():
            small += 1
            if small >= 6:
                break
        else:
            small = 0
    M = len(a)
    za = arb(1)                                       # |second exponent| = 1
    Nn = arb(M)
    A1l = abs(A1).lower()
    g = []
    for k in range(1, 4):
        Ak = _absu(Ap[k + 1]) if k + 1 < len(Ap) else arb(0)
        Bk = _absu(B[k]) if k < len(B) else arb(0)
        Ck = _absu(C[k - 1]) if k - 1 < len(C) else arb(0)
        g.append(arb(((Ak * Nn * Nn + Bk * Nn + Ck) / (A1l * Nn * (Nn - za))).upper()))
    q = arb(_solve_q(g))
    if not (sum((gk * q ** (-(k + 1)) for k, gk in enumerate(g)), arb(0)) <= 1):
        raise ArithmeticError("frob_z2 majorant root not verified")
    K = arb(max((_absu(a[m]) / q ** m for m in range(M - 3, M)), key=lambda x: x.upper()).upper())
    s0, s1 = _geom_tails(K, q, abs(ta), M, arb(0))
    S += acb(arb(0, s0.upper()), arb(0, s0.upper()))
    Sd += acb(arb(0, s1.upper()), arb(0, s1.upper()))
    return S, Sd, M


def wt_at_omega0(pt: Point, l: int) -> dict:
    """Wt_l = [Delta (Z2 Zc' - Z2' Zc) + 2 i w r^2 Z2 Zc]/(r_+ r_c) at omega0 and r = r_m (zsol's form)."""
    from . import zsol
    from .arbw import _rat, majorant_radius
    L = l * (l + 1)
    w0 = acb(0, -1) * pt.kp
    rs = [float(x.mid()) for x in pt.roots]
    xr = majorant_radius(rs[2], [rs[0], rs[1], rs[3]])
    t = _rat(min(0.2 * 2.5 * xr, 2.5 * xr), 30)
    Z, dZ, nterms = frob_z2(pt, w0, L, t)
    rm = _rat(pt.rP.mid(), 20)
    start = pt.rp + arb(t)
    pts = zsol.path(pt, start, rm, float(abs(w0).upper()), 0.5, 1.0)
    c = start
    for nxt in pts:
        h = arb(nxt) - c
        Z, dZ = zsol.taylor_z(pt, c, Z, dZ, w0, L, h, 1)
        c = arb(nxt)
    Zc, dZc, ec, sm = zsol.propagate(pt, "c", w0, None, L, rmatch=rm)
    Zc = Zc + acb(arb(0, ec.upper()), arb(0, ec.upper()))
    dZc = dZc + acb(arb(0, (ec * sm).upper()), arb(0, (ec * sm).upper()))
    rma = arb(rm)
    D = acb(pt.Delta(rma))
    wt = (D * (Z * dZc - dZ * Zc) + 2 * I * w0 * acb(rma * rma) * Z * Zc) / acb(pt.rp * pt.rc)
    res = I * (2 * pt.kp * pt.rp + L) / (2 * pt.rp * pt.rp) * wt
    return {"l": l, "Wt": wt.str(12), "abs_Wt_lower": float(abs(wt).lower()), "Wt_nonzero": bool(not wt.contains(0)),
            "Re_Wt": wt.real.str(10), "Im_Wt_contains_0": bool(wt.imag.contains(0)),
            "residue_of_W_l": res.str(12), "frobenius_terms": nterms, "steps_plus": len(pts)}


def v_positivity(pt: Point) -> dict:
    """L + r F'(r) > 0 on [r_+, r_c] for every L >= 2: r F' = (2/r)(M - Q^2/r) - (2 Lambda/3) r^2 with
    r_+ > Q^2/M (so the first term is >= 0) and (2 Lambda/3) r_c^2 < 2."""
    m1 = pt.rp - pt.Q2                       # M = 1
    m2 = 2 - 2 * pt.Lam / 3 * pt.rc * pt.rc
    return {"r_plus_minus_Q2_over_M": m1.str(10), "two_minus_2Lambda_rc2_over_3": m2.str(10),
            "pass": bool(m1 > 0 and m2 > 0)}


def run_E(lmax=19, lfloat=60):
    pt = Point("PB", 256)
    out = {"omega0": "-i kappa_+", "omega0_over_kappa_minus": (pt.kp / pt.km).str(12), "commit": git_commit()}
    with Timer() as tm:
        out["sympy_first_step"] = e_sympy()
        with precision(PREC):
            out["V_l_positive_for_l_ge_1"] = v_positivity(pt)
            rows = []
            for l in range(0, lmax + 1):
                rows.append(wt_at_omega0(pt, l))
                print("E", l, rows[-1]["Wt"], rows[-1]["Wt_nonzero"], flush=True)
            out["Wt_l_le_19"] = rows
            more = []
            for l in list(range(lmax + 1, lfloat + 1, 5)) + [100]:
                more.append(wt_at_omega0(pt, l))
                print("E", l, more[-1]["Wt"], more[-1]["Wt_nonzero"], flush=True)
            out["Wt_l_ge_20_sample"] = more
            # residue cross-check from the other side: eps W_l(omega0 + eps) -> Res
            chk = []
            from . import zsol
            for l in (0, 1, 5, 19):
                w0 = acb(0, -1) * pt.kp
                for er in (Fraction(1, 10 ** 5), Fraction(1, 10 ** 7)):
                    eps = acb(1, 1) / arb(2).sqrt() * A(er) * pt.km
                    Wv, _ = zsol.wronskian_z(pt, l, w0 + eps)
                    prod = eps * Wv
                    r = next(x for x in rows if x["l"] == l)
                    chk.append({"l": l, "eps_over_kappa_minus": str(er), "eps_times_W": prod.str(10),
                                "residue": r["residue_of_W_l"]})
            out["residue_cross_check"] = chk
    out["timing_s"] = {"wall": round(tm.wall, 1), "cpu": round(tm.cpu, 1)}
    _write("E_PB.json", out)
    return out


# ---------------------------------------------------------------------------------------------
# M3: S2's barrier-top windows at depth e
# ---------------------------------------------------------------------------------------------
def run_M3(ls=(20, 50, 100), depth=EDGE):
    from . import m3
    pt = deep_point(depth)
    out = {"depth_over_kappa_minus": str(depth), "commit": git_commit(), "rows": {}}
    with precision(m3.PREC):
        geo, rays, tails = m3.geometry_certificate(pt)
        out["A0_pass"] = geo["A0_pass"]
        for l in ls:
            c = m3.certify_l(pt, l, rays, tails)
            out["rows"][l] = {k: c[k] for k in ("m_ratio", "mu_cert", "mu_HO", "float_lambda_min", "A1_pass",
                                                "jost_decay_condition", "window_Im_over_kappa_minus", "timing_s")}
            print("M3", l, c["m_ratio"], c["A1_pass"], c["timing_s"], flush=True)
    _write(f"M3_PB_e{depth.numerator}_{depth.denominator}.json", out)
    return out


# ---------------------------------------------------------------------------------------------
# TA: S4's T and A at depth e (uniform in h <= 2/41); A at fixed h for l <= 19
# ---------------------------------------------------------------------------------------------
def run_TA(depth=EDGE):
    from . import part1_T as P
    pt = deep_point(depth)
    out = {"depth_over_kappa_minus": str(depth), "commit": git_commit()}
    with precision(P.PREC):
        _, rays, tails = P.setup(pt)
        chT = P.TChecker(pt, rays, tails, P.DESIGN_T)
        rows = []
        for zb in P.zeta_boxes(P.DESIGN_T):
            r = P.certify_box(chT, pt, zb, P.DESIGN_T)
            rows.append({k: r[k] for k in ("zeta_over_OmegaP", "mu_over_a_star", "E0_low_min", "c_min", "cells",
                                           "n_fails", "pass", "timing_s")})
            print("T", zb, r["pass"], r["mu_over_a_star"], r["cells"], flush=True)
        out["T"] = rows
        chA = P.TChecker(pt, rays, tails, P.DESIGN_A)
        r = P.certify_box(chA, pt, (Fraction(3, 2), Fraction(3, 2)), P.DESIGN_A)
        out["A"] = {k: r[k] for k in ("mu_over_a_star", "E0_low_min", "pass", "timing_s")}
        out["A_monotonicity_W7"] = P.monotonicity_A(pt, P.DESIGN_A)
        print("A", out["A"], out["A_monotonicity_W7"]["pass"], flush=True)
    _write(f"TA_PB_e{depth.numerator}_{depth.denominator}.json", out)
    return out


def run_AFIX(ls=None, depth=EDGE):
    """A at fixed h_l = 1/(l + 1/2) (uniform in h <= h_l): smallest zeta_A on a grid that passes."""
    from . import part1_T as P
    from .m3 import fa
    pt = deep_point(depth)
    ls = list(range(20)) if ls is None else ls
    grid = [Fraction(k, 4) for k in range(4, 49)]          # zeta_A / Omega_P = 1, 1.25, ..., 12
    out = {"depth_over_kappa_minus": str(depth), "commit": git_commit(), "rows": {}}
    with precision(P.PREC):
        _, rays, tails = P.setup(pt)
        chA = P.TChecker(pt, rays, tails, P.DESIGN_A)
        for l in ls:
            h = Fraction(2, 2 * l + 1)
            found = None
            tried = []
            lo_, hi_ = 0, len(grid) - 1
            # the A bound is monotone in zeta_A in practice; bisect on the grid
            ok_hi = None
            while lo_ <= hi_:
                mid = (lo_ + hi_) // 2
                z = grid[mid]
                r = P.certify_box(chA, pt, (z, z), P.DESIGN_A, h0=h)
                zA = fa(z) * pt.OmegaP
                w7 = bool(fa(Fraction(-1, 5)) < -((pt.beta * fa(h)) / zA).atan())
                passed = bool(r["pass"] and w7)
                tried.append({"zeta_A": str(z), "pass": passed, "mu_over_a_star": r["mu_over_a_star"]})
                if passed:
                    found, ok_hi = z, r
                    hi_ = mid - 1
                else:
                    lo_ = mid + 1
            a_l = (fa(found) * pt.OmegaP / fa(h) / pt.km).str(6) if found is not None else None
            out["rows"][l] = {"h": str(h), "zeta_A_min_over_OmegaP": str(found) if found is not None else None,
                              "a_l_over_kappa_minus": a_l, "tried": tried}
            print("AFIX", l, found, a_l, flush=True)
    _write(f"AFIX_PB_e{depth.numerator}_{depth.denominator}.json", out)
    return out


# ---------------------------------------------------------------------------------------------
# S: S4's regime S at depth e
# ---------------------------------------------------------------------------------------------
def run_S(side_name: str, depth=EDGE):
    from . import throat as T
    from .part1 import load_bulk
    side = 1 if side_name == "+" else -1
    W_bound, _ = load_bulk()
    pt = deep_point(depth)
    g = T.FGeo(pt)
    out = {"side": side_name, "depth_over_kappa_minus": str(depth), "commit": git_commit(), "rows": []}
    with precision(T.PREC), Timer() as tm:
        try:
            c = T.certify_closing(pt, g, side, W_bound)
            out["closing"] = {k: c.get(k) for k in ("pass", "mu_cells", "mu_tail", "error", "design", "eps_max")}
        except Exception as ex:                                      # noqa: BLE001 (probe)
            out["closing"] = {"pass": False, "error": repr(ex)}
        print("S", side_name, "closing", out["closing"].get("pass"), flush=True)
        for a, b in T.zboxes():
            with Timer() as tb:
                leaves = T.certify_zbox_refined(pt, g, side, a, b, W_bound)
            row = {"zeta": [str(a), str(b)], "pass": all(x["pass"] for x in leaves), "n_leaves": len(leaves),
                   "float_margins": [x.get("float_margin") for x in leaves],
                   "mu_min": min((e["mu_cells"] for x in leaves for e in x.get("eps_boxes", []) if e.get("mu_cells")),
                                 default=None, key=lambda s: float(s.strip("[]").split(",")[0])),
                   "errors": [x.get("error") for x in leaves if x.get("error")][:3],
                   "cpu": round(tb.cpu, 1)}
            out["rows"].append(row)
            print("S", side_name, (str(a), str(b)), row["pass"], row["n_leaves"], row["mu_min"], row["cpu"], flush=True)
    out["timing_s"] = {"wall": round(tm.wall, 1), "cpu": round(tm.cpu, 1)}
    tag = "plus" if side > 0 else "c"
    _write(f"S_PB_{tag}_e{depth.numerator}_{depth.denominator}.json", out)
    return out


# ---------------------------------------------------------------------------------------------
# L: S4's regime L on the band, with a pole-free series near omega0
# ---------------------------------------------------------------------------------------------
def frob_tau_pf(A_, B_, C_, rho, tau, min_terms=30):
    """lowfreq.frob_tau for R~ = (2 rho + 1) R (side +): a_0 = 2 rho + 1 and the n = 1 step with the
    factor (2 rho + 1) cancelled exactly, a_1 = -(A_3 rho (rho - 1) + B_2 rho + C_1)/A_2."""
    from .affine import AF, up
    from .arbw import _geom_tails, _solve_q
    from .lowfreq import N_MAX
    A2 = A_[2]
    ta = arb(tau)
    a0 = rho * 2 + 1
    term0 = rho.const(0)
    if 3 < len(A_):
        term0 = term0 + A_[3] * (rho * (rho + (-1)))
    if 2 < len(B_):
        term0 = term0 + B_[2] * rho
    if 1 < len(C_):
        term0 = term0 + C_[1]
    a1 = -(term0 / A2)
    a = [a0, a1]
    S = a0 + a1 * ta
    Sd = a0 * rho * (1 / acb(ta)) + a1 * (rho + 1)
    tk = arb(ta)
    small = 0
    eps = arb(2) ** (-(ctx.prec - 12))
    for n in range(2, N_MAX):
        acc = rho.const(0)
        for k in range(max(0, n - 6), n):
            j2, j1, j0 = n - k + 2, n - k + 1, n - k
            term = rho.const(0)
            if j2 < len(A_):
                term = term + A_[j2] * ((rho + k) * (rho + (k - 1)))
            if j1 < len(B_):
                term = term + B_[j1] * (rho + k)
            if j0 < len(C_):
                term = term + C_[j0]
            acc = acc + a[k] * term
        an = -(acc / (A2 * ((rho * 2 + n) * n)))
        a.append(an)
        tk = tk * ta
        S = S + an * tk
        Sd = Sd + an * ((rho + n) * (tk / ta))
        if n >= min_terms and an.absu() * tk < eps * abs(S.c0).lower():
            small += 1
            if small >= 8:
                break
        else:
            small = 0
    M = len(a)
    Rr = rho.absu()
    Nn = arb(M)
    if not (Nn > 2 * Rr + 1):
        raise ArithmeticError("frob_tau_pf: N too small relative to |rho|")
    A2l = abs(A2.c0).lower() - A2.dev()
    g = []
    for j in range(1, 7):
        Aj = A_[j + 2].absu() if j + 2 < len(A_) else arb(0)
        Bj = B_[j + 1].absu() if j + 1 < len(B_) else arb(0)
        Cj = C_[j].absu() if j < len(C_) else arb(0)
        gj = (Aj * (1 + Rr / Nn) ** 2 / (1 - 2 * Rr / Nn) + Bj * (Nn + Rr) / (Nn * (Nn - 2 * Rr))
              + Cj / (Nn * (Nn - 2 * Rr))) / A2l
        g.append(up(gj))
    q = arb(_solve_q(g))
    if not (sum((gj * q ** (-(j + 1)) for j, gj in enumerate(g)), arb(0)) <= 1):
        raise ArithmeticError("frob_tau_pf majorant root not verified")
    K = max((a[m].absu() / q ** m for m in range(M - 6, M)), key=lambda x: x.upper())
    s0, s1 = _geom_tails(up(K), q, abs(ta), M, Rr)
    S = S + AF(0, 0, 0, s0, S.Rw, S.Rh)
    Sd = Sd + AF(0, 0, 0, s1, S.Rw, S.Rh)
    return S, Sd, {"terms": M}


def run_L(depth=EDGE, im_top=Fraction(-11, 20), only_pole=False):
    """Band boxes [0, 8] x [-ceil20(depth), im_top] (side 1/20).  Boxes touching omega0 use frob_tau_pf."""
    from . import lowfreq as LF
    from .part1 import load_bulk
    W_bound, _ = load_bulk()
    pt = deep_point(depth)
    n_lo = -(-depth * 20 // 1)                          # ceil(20 depth)
    im_lo = -Fraction(int(n_lo), 20)
    kpk = float((pt.kp / pt.km).mid())
    out = {"depth_over_kappa_minus": str(depth), "band": [str(im_lo), str(im_top)], "commit": git_commit(),
           "fails": [], "worst_plus": None, "worst_c": None, "pole_boxes": []}
    orig = LF.frob_tau

    def touches_pole(b):
        return b[0] <= 0 <= b[1] and float(b[2]) <= -kpk <= float(b[3])

    with precision(LF.PREC), Timer() as tm:
        ch = LF.LChecker(pt, W_bound)
        boxes = []
        nre = int(Fraction(8) / LF.BOX)
        j = 0
        while True:
            hi_ = im_top - j * LF.BOX
            lo_ = hi_ - LF.BOX
            if hi_ <= im_lo:
                break
            for i in range(nre):
                boxes.append((i * LF.BOX, (i + 1) * LF.BOX, max(lo_, im_lo), hi_))
            j += 1
        n_ok = 0
        for b in boxes:
            pole = touches_pole(b)
            if only_pole and not pole:
                continue
            if pole:
                LF.frob_tau = lambda A_, B_, C_, rho, tau, min_terms=30: frob_tau_pf(A_, B_, C_, rho, tau, min_terms)
            try:
                rec = ch.box(b)
            finally:
                LF.frob_tau = orig
            good = all(rec.get(s, {}).get("y_sign_ok") and rec.get(s, {}).get("exit_rule_ok") for s in ("plus", "c"))
            n_ok += int(good)
            for s, key in (("plus", "worst_plus"), ("c", "worst_c")):
                v = rec.get(s, {}).get("Re_y_bound")
                if v is not None:
                    f = float(arb(v).lower()) if s == "plus" else float(arb(v).upper())
                    cur = out[key]
                    if cur is None or (f < cur[0] if s == "plus" else f > cur[0]):
                        out[key] = [f, [str(x) for x in b]]
            if pole:
                out["pole_boxes"].append({"box": [str(x) for x in b], "rec": rec, "ok": good})
            if not good:
                out["fails"].append({"box": [str(x) for x in b], "rec": rec})
        out["n_boxes"] = len(boxes) if not only_pole else len(out["pole_boxes"])
        out["n_ok"] = n_ok
    out["timing_s"] = {"wall": round(tm.wall, 1), "cpu": round(tm.cpu, 1)}
    _write(f"L_PB_e{depth.numerator}_{depth.denominator}{'_pole' if only_pole else ''}.json", out)
    print("L", out["n_ok"], "/", out["n_boxes"], "worst", out["worst_plus"], out["worst_c"], out["timing_s"], flush=True)
    return out


def main(argv):
    cmd = argv[1]
    t0 = time.time()
    if cmd == "E":
        run_E()
    elif cmd == "M3":
        run_M3(tuple(int(x) for x in argv[2:]) or (20, 50, 100))
    elif cmd == "TA":
        run_TA()
    elif cmd == "AFIX":
        run_AFIX([int(x) for x in argv[2:]] or None)
    elif cmd in ("S+", "Sc"):
        run_S(cmd[1])
    elif cmd == "L":
        run_L(only_pole=(len(argv) > 2 and argv[2] == "pole"))
    else:
        raise SystemExit(f"unknown probe {cmd}")
    print("done", cmd, round(time.time() - t0, 1), "s", flush=True)


if __name__ == "__main__":
    main(sys.argv)


# ---------------------------------------------------------------------------------------------
# float evaluator of the T/A bound (design study only; S4 part1_T formulas)
# ---------------------------------------------------------------------------------------------
class FloatTA:
    """F(s) = E0/h0 + E1 + h0 min(0, E2) on S2's contour, in floats (part1_T.TChecker.parts)."""

    def __init__(self, pt: Point, h=1.0 / 256, s_lo=-64.0, s_hi=32.0):
        import numpy as np
        from .m3 import float_contour
        self.np = np
        s, r = float_contour(pt, s_lo, s_hi, h)
        self.s, self.r = s, r
        lam3 = float((pt.Lam / 3).mid())
        Q2 = float(pt.Q2.mid())
        D = -lam3 * r ** 4 + r ** 2 - 2 * r + Q2
        Dp = -4 * lam3 * r ** 3 + 2 * r - 2
        F = D / r ** 2
        Fr = (Dp * r - 2 * D) / r ** 3
        self.U = F / r ** 2
        self.W = F * Fr / r - self.U / 4
        self.OmP = float(pt.OmegaP.mid())
        self.km = float(pt.km.mid())
        self.astar = float((2 * pt.OmegaP * pt.gamma0).mid())

    def F(self, d: dict, zl, zh, nl, h0):
        np = self.np
        s = self.s
        s1, phinf, phi0, eta = float(Fraction(d["s1"])), float(Fraction(d["phi_inf"])), float(Fraction(d["phi0"])), float(Fraction(d["eta"]))
        phi = phi0 + (phinf - phi0) * s ** 4 / (s1 ** 4 + s ** 4)
        ds = s[1] - s[0]
        d1 = np.gradient(phi, ds)
        c = np.cos(phi) - eta * np.abs(d1)
        cp = np.gradient(c, ds)
        if d["eps_g"] is None:
            G = np.zeros_like(s)
        else:
            k = (1 - float(Fraction(d["eps_g"]))) * float(Fraction(d["a_d"]))
            tc, td = (float(Fraction(x)) for x in d["taper"])
            a = np.abs(s)
            th = np.where(a <= tc, 1.0, np.where(a >= td, 0.0, (1 + np.cos(np.pi * (a - tc) / (td - tc))) / 2))
            G = k * s * th
        Gp = np.gradient(G, ds)
        rot = np.sin(phi) + 1j * np.cos(phi)
        E0 = np.real(rot * (self.U - zl * zl)) - c * G * G
        E1 = 2 * zh * nl * np.cos(phi) + cp * G + c * Gp
        E2 = np.real(rot * (self.W + nl * nl)) - np.abs(d1) / (4 * eta)
        Fv = E0 / h0 + E1 + h0 * np.minimum(0, E2)
        i = int(np.argmin(Fv))
        return {"F_min_over_astar": float(Fv[i] / self.astar), "s_at_min": float(s[i]),
                "E0_min": float(E0.min()), "c_min": float(c.min())}

    def box(self, d: dict, zb, depth, h0):
        zl, zh = float(zb[0]) * self.OmP, float(zb[1]) * self.OmP
        nl = -float(depth) * self.km
        return self.F(d, zl, zh, nl, float(h0))


# ---------------------------------------------------------------------------------------------
# part (iii): band-box segment plan and timing samples; small pole-aware KAT boxes
# ---------------------------------------------------------------------------------------------
def band_plan(pt: Point, a_k: float, b_k=Fraction(23, 25), top_k=Fraction(1, 2)):
    """Segments of the right half of [-a, a] x [-b, -top] (kappa_- units), S3's planner."""
    import itertools
    from flint import fmpq
    from . import part3
    km = float(pt.km.mid())
    a = part3._dy_up(a_k * km)
    b = part3._dy_up(float(b_k) * km)
    t = part3._dy_down(float(top_k) * km)
    corners = [(fmpq(0), -b), (a, -b), (a, -t), (fmpq(0), -t)]
    segs = []
    for k, (z0, z1) in enumerate(itertools.pairwise(corners)):
        for j, s in enumerate(part3.plan_edge(pt, z0, z1, part3.RHO_MAX_K * km)):
            s.update({"edge": k, "idx": j})
            segs.append(s)
    return segs


def run_P3(spec):
    """spec: list of (l, a_k, [indices of segments to time])."""
    from . import part3
    pt = Point("PB", 256)
    out = {"commit": git_commit(), "rows": []}
    km = float(pt.km.mid())
    with precision(part3.PREC):
        for l, a_k, idxs in spec:
            segs = band_plan(pt, a_k)
            row = {"l": l, "a_over_kappa_minus": a_k, "planned_segments": len(segs), "samples": []}
            for i in idxs:
                s = segs[i]
                with Timer() as tm:
                    try:
                        res = part3.segment_with_split(pt, l, s, part3.RHO_MIN_K * km)
                        ok = True
                    except ArithmeticError as ex:
                        res, ok = [], str(ex)[:120]
                row["samples"].append({"i": i, "c_over_kappa_minus": [float(s["c"][0]) / km, float(s["c"][1]) / km],
                                       "rho_over_kappa_minus": float(s["rho"]) / km, "ok": ok, "pieces": len(res),
                                       "increments": [r["increment"][:24] for r in res], "cpu": round(tm.cpu, 1)})
                print("P3", l, i, row["samples"][-1], flush=True)
            out["rows"].append(row)
    _write("P3_timing_PB.json", out)
    return out


def run_KATBOX(boxes):
    """Small pole-aware winding boxes (full boundary): (l, (x0, x1, y0, y1)) in kappa_- units."""
    from . import part3
    out = {"commit": git_commit(), "rows": []}
    for l, bx in boxes:
        with Timer() as tm:
            r = part3.rectangle_winding("PB", l, bx, workers=1)
        r["cpu_total"] = round(tm.cpu, 1)
        out["rows"].append(r)
        print("KATBOX", l, bx, r.get("winding"), r.get("winding_ball"), r.get("segments"), r["cpu_total"], flush=True)
    _write("KATBOX_PB.json", out)
    return out


# ---------------------------------------------------------------------------------------------
# design records (T redesign, A thresholds, the S grid extension, L with the pole-free series)
# ---------------------------------------------------------------------------------------------
T_REDESIGN = {"s1": "12", "phi_inf": "-1", "taper": ["3", "6"]}        # the rest as S4's DESIGN_T
S_GRID_EXT = {"w": ["1/2", "3/2", "5/2"], "y_a": ["1/40", "3/10", "3/5"], "phi_n": [5, 6],
              "phase_centre": ["-3/2", "1/2"], "phase_width": ["3/4", "2"]}


def t_design():
    from . import part1_T as P
    d = dict(P.DESIGN_T)
    d.update({"s1": T_REDESIGN["s1"], "phi_inf": T_REDESIGN["phi_inf"], "taper": list(T_REDESIGN["taper"])})
    return d


def extend_s_grid():
    from . import throat as T
    T.GRID_W = T.GRID_W + [Fraction(x) for x in S_GRID_EXT["w"]]
    T.GRID_YA = T.GRID_YA + [Fraction(x) for x in S_GRID_EXT["y_a"]]
    T.GRID_PHI = T.GRID_PHI + list(S_GRID_EXT["phi_n"])
    T.GRID_PC = T.GRID_PC + [Fraction(x) for x in S_GRID_EXT["phase_centre"]]
    T.GRID_PW = T.GRID_PW + [Fraction(x) for x in S_GRID_EXT["phase_width"]]


def run_DESIGN(depth=EDGE):
    import math
    import random
    from . import lowfreq as LF
    from . import part1_T as P
    from . import throat as T
    from .m3 import fa
    from .part1 import load_bulk
    pt = deep_point(depth)
    out = {"depth_over_kappa_minus": str(depth), "commit": git_commit()}
    f = FloatTA(pt)
    dT = t_design()
    zbs = [(Fraction(3, 5) + Fraction(k, 20), Fraction(3, 5) + Fraction(k + 1, 20)) for k in range(28)]
    out["T_float"] = [{"zeta": [str(a), str(b)],
                       "S4_design": round(f.box(P.DESIGN_T, (a, b), depth, P.H0)["F_min_over_astar"], 4),
                       "redesign": round(f.box(dT, (a, b), depth, P.H0)["F_min_over_astar"], 4)} for a, b in zbs]
    out["T_redesign"] = T_REDESIGN
    rows = []
    for l in range(20):
        h = Fraction(2, 2 * l + 1)
        zA = None
        for k in range(400):
            z = Fraction(100 + 5 * k, 100)
            v = f.box(P.DESIGN_A, (z, z), depth, h)
            w7 = -0.2 < -math.atan(float(depth) * f.km * float(h) / (float(z) * f.OmP))
            if v["F_min_over_astar"] > 0.02 and v["E0_min"] >= 0 and w7:
                zA = z
                break
        rows.append({"l": l, "h": str(h), "zeta_A_float": str(zA),
                     "a_prime_over_kappa_minus": round(float(zA) * f.OmP / float(h) / f.km, 2)})
    out["A_fixed_h_float"] = rows
    W_bound, _ = load_bulk()
    with precision(128):
        _, rays, tails = P.setup(pt)
        chT = P.TChecker(pt, rays, tails, dT)
        out["T_redesign_arb"] = []
        for zb in [(Fraction(4, 5), Fraction(17, 20)), (Fraction(9, 10), Fraction(19, 20)), (Fraction(5, 4), Fraction(13, 10))]:
            r = P.certify_box(chT, pt, zb, dT)
            out["T_redesign_arb"].append({"zeta": [str(zb[0]), str(zb[1])], "pass": r["pass"], "mu_over_a_star": r["mu_over_a_star"],
                                          "cells": r["cells"]})
        dA = dict(P.DESIGN_A)
        dA["zeta_over_OmegaP"] = ["6/5", "6/5"]
        chA = P.TChecker(pt, rays, tails, dA)
        r = P.certify_box(chA, pt, (Fraction(6, 5), Fraction(6, 5)), dA)
        out["A_6_5_arb"] = {"pass": r["pass"], "mu_over_a_star": r["mu_over_a_star"], "W7": P.monotonicity_A(pt, dA)["pass"]}
        chA0 = P.TChecker(pt, rays, tails, P.DESIGN_A)
        out["A_fixed_h_arb"] = []
        for l in (0, 1, 2, 5, 10, 19):
            h = Fraction(2, 2 * l + 1)
            z = Fraction(rows[l]["zeta_A_float"])
            r = P.certify_box(chA0, pt, (z, z), P.DESIGN_A, h0=h)
            w7 = bool(fa(Fraction(-1, 5)) < -((pt.beta * fa(h)) / (fa(z) * pt.OmegaP)).atan())
            out["A_fixed_h_arb"].append({"l": l, "zeta_A": str(z), "pass": bool(r["pass"] and w7),
                                         "mu_over_a_star": r["mu_over_a_star"], "cells": r["cells"]})
            print("A fixed h", l, z, r["pass"], r["mu_over_a_star"], flush=True)
        g = T.FGeo(pt)
        extend_s_grid()
        leaves = T.certify_zbox_refined(pt, g, 1, Fraction(3, 50), Fraction(2, 25), W_bound)
        out["S_plus_0.06_0.08_extended_grid"] = [{"zeta": x["zeta_over_OmegaP"], "pass": x["pass"], "design": x.get("design"),
                                                  "mu": [e["mu_cells"] for e in x.get("eps_boxes", [])]} for x in leaves]
        out["S_grid_extension"] = S_GRID_EXT
        random.seed(20261027)
        boxes = []
        for _ in range(60):
            i, j = random.randrange(160), random.randrange(8)
            boxes.append((Fraction(i, 20), Fraction(i + 1, 20), Fraction(-11, 20) - Fraction(j + 1, 20), Fraction(-11, 20) - Fraction(j, 20)))
        boxes += [(Fraction(0), Fraction(1, 20), Fraction(-17, 20), Fraction(-4, 5)),
                  (Fraction(1, 20), Fraction(1, 10), Fraction(-17, 20), Fraction(-4, 5)),
                  (Fraction(0), Fraction(1, 20), Fraction(-9, 10), Fraction(-17, 20))]
        orig = LF.frob_tau
        LF.frob_tau = lambda A_, B_, C_, rho, tau, min_terms=30: frob_tau_pf(A_, B_, C_, rho, tau, min_terms)
        ch = LF.LChecker(pt, W_bound)
        wp, wc, bad = None, None, 0
        try:
            for b in boxes:
                rec = ch.box(b)
                yp = float(arb(rec["plus"]["Re_y_bound"]).lower())
                yc = float(arb(rec["c"]["Re_y_bound"]).upper())
                ok = all(rec[s]["y_sign_ok"] and rec[s]["exit_rule_ok"] for s in ("plus", "c"))
                bad += int(not ok)
                wp = yp if wp is None else min(wp, yp)
                wc = yc if wc is None else max(wc, yc)
        finally:
            LF.frob_tau = orig
        out["L_pole_free_sample"] = {"boxes": len(boxes), "failed": bad, "min_Re_y_plus": wp, "max_Re_y_c": wc,
                                     "seed": 20261027}
    _write(f"DESIGN_PB_e{depth.numerator}_{depth.denominator}.json", out)
    return out
