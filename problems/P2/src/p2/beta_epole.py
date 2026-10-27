"""S6 Lemma E at P_B: the E_+ point omega0 = -i kappa_+ is not a resonance (pre-registration O1, A1).

Normalisation (A1).  zsol's Z_+ equation  Delta Z'' + (Delta' - 2 i w r^2) Z' - (2 i w r + L) Z = 0,
Z_+ = sum_{m>=0} a_m t^m, a_0 = 1, t = r - r_+ (zsol.frob_z).  The m-th step is D_m a_m + acc_m = 0 with
D_m = m (Delta'(r_+)(m - 1) + B_0), B_0 = Delta'(r_+) - 2 i w r_+^2; at omega0, B_0 = 0 exactly, so
D_1 = 0 while acc_1 = C_0 = -(2 i w r_+ + L) = -(2 kappa_+ r_+ + L) != 0 (a log term), and
D_m = Delta'(r_+) m (m - 1) != 0 for m >= 2.  Z2 := t + sum_{m>=2} b_m t^m is the exponent-1 solution at
omega0.  Near omega0, Z_+ = [i (2 kappa_+ r_+ + L)/(2 r_+^2)] Z2/(w - omega0) + O(1), hence
    Res_{omega0} W_l = i (2 kappa_+ r_+ + L)/(2 r_+^2) Wt_l,   Wt_l := W(Z2, Z_c) in zsol's form
    [Delta (Z Zc' - Z' Zc) + 2 i w r^2 Z Zc]/(r_+ r_c) at r = r_m,
and omega0 is a resonance in sector l iff Wt_l = 0.

Certified here (128 bits): the first step (sympy); Wt_l(omega0) != 0 for l <= 19 (and sample values
beyond); the inequalities behind the written argument for l >= 1 (V_l > 0 on (r_+, r_c));
KAT B5: eps W_l(omega0 + eps) -> Res, with zsol and with the independent S1 code arbw.

Run:  python -m p2.beta_epole            (writes certificates/beta/epole_PB.json)
"""
from __future__ import annotations

import sys
import time
from fractions import Fraction

from flint import acb, arb, ctx

from .arbutil import A, Timer, precision
from .beta_common import CERT, PREC, commit_tag, write_json
from .params import Point

I = acb(0, 1)
SAMPLES = [20, 25, 30, 40, 50, 60, 100]
B5_L = [0, 1, 5, 19]
B5_EPS = [Fraction(1, 10 ** 5), Fraction(1, 10 ** 7)]


def sympy_first_step() -> dict:
    """The t^0 coefficient of zsol's Z_+ equation for Z = 1 + a1 t + O(t^2); its value at omega0."""
    import sympy as sp
    t, w, L, rp, d1, d2, d3, d4, a1 = sp.symbols("t w L r_p d1 d2 d3 d4 a1")
    Delta = d1 * t + d2 * t ** 2 + d3 * t ** 3 + d4 * t ** 4
    r = rp + t
    Z = 1 + a1 * t
    expr = sp.expand(Delta * sp.diff(Z, t, 2) + (sp.diff(Delta, t) - 2 * sp.I * w * r ** 2) * sp.diff(Z, t)
                     - (2 * sp.I * w * r + L) * Z)
    c0 = sp.expand(expr.coeff(t, 0))
    D1 = d1 - 2 * sp.I * w * rp ** 2
    C0 = -(2 * sp.I * w * rp + L)
    ok_form = sp.expand(c0 - (D1 * a1 + C0)) == 0
    kp = sp.symbols("kappa_p", positive=True)
    sub = {d1: 2 * kp * rp ** 2, w: -sp.I * kp}
    D1_0 = sp.simplify(D1.subs(sub))
    C0_0 = sp.simplify(C0.subs(sub))
    ok0 = bool(D1_0 == 0 and sp.simplify(C0_0 + 2 * kp * rp + L) == 0)
    # the general step: for Z = t^m the lowest-order term is t^(m-1) times m (d1 (m - 1) + B0)
    m = sp.symbols("m", positive=True, integer=True)
    B0 = d1 - 2 * sp.I * w * rp ** 2
    Dm = m * (d1 * (m - 1) + B0)
    Dm0 = sp.factor(sp.simplify(Dm.subs(sub)))
    ok_m = sp.simplify(Dm0 - 2 * kp * rp ** 2 * m * (m - 1)) == 0
    return {"t0_coefficient": str(c0), "equals_D1_a1_plus_C0": bool(ok_form),
            "D1_at_omega0": str(D1_0), "C0_at_omega0": str(C0_0), "Dm_at_omega0": str(Dm0),
            "Dm_at_omega0_equals_Delta'_m_(m-1)": bool(ok_m), "pass": bool(ok_form and ok0 and ok_m)}


def frob_z2(pt: Point, w, L, t, min_terms=20):
    """(Z, Z') at r_+ + t (t exact) for Z2 = t + sum_{m>=2} b_m t^m at omega0 (w a ball containing
    omega0; the indicial factors use the exact exponents 0, 1 of omega0)."""
    from .arbw import _geom_tails, _solve_q
    from .zsol import N_MAX, _absu, zpolys
    Ap, B, C = zpolys(pt, pt.rp, acb(w), L, 1, singular=True)
    A1 = Ap[1]
    ta = arb(t)
    a = [acb(0), acb(1)]
    S = acb(ta)
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
        am = -acc / (A1 * (m * (m - 1)))
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
    za = arb(1)
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


def wt_ball(pt: Point, l: int):
    """(Wt_l(omega0) as an acb ball, Res_{omega0} W_l as an acb ball, info)."""
    from . import zsol
    from .arbw import _rat, majorant_radius
    L = l * (l + 1)
    w0 = -I * pt.kp
    rs = [float(x.mid()) for x in pt.roots]
    xr = majorant_radius(rs[2], [rs[0], rs[1], rs[3]])
    t = _rat(0.5 * xr, 30)
    Z, dZ, nterms = frob_z2(pt, w0, L, t)
    rm = _rat(pt.rP.mid(), 20)
    start = pt.rp + arb(t)
    pts = zsol.path(pt, start, rm, float(abs(w0).upper()), 0.5, 1.0)
    c = start
    for nxt in pts:
        Z, dZ = zsol.taylor_z(pt, c, Z, dZ, w0, L, arb(nxt) - c, 1)
        c = arb(nxt)
    Zc, dZc, ec, sm = zsol.propagate(pt, "c", w0, None, L, rmatch=rm)
    Zc = Zc + acb(arb(0, ec.upper()), arb(0, ec.upper()))
    dZc = dZc + acb(arb(0, (ec * sm).upper()), arb(0, (ec * sm).upper()))
    rma = arb(rm)
    D = acb(pt.Delta(rma))
    w_ = (D * (Z * dZc - dZ * Zc) + 2 * I * w0 * acb(rma * rma) * Z * Zc) / acb(pt.rp * pt.rc)
    pref = I * (2 * pt.kp * pt.rp + L) / (2 * pt.rp * pt.rp)
    info = {"frobenius_terms_Z2": nterms, "taylor_steps_plus": len(pts), "r_match": str(rm), "t_start": str(t),
            "residue_prefactor_i(2k r+L)/(2r^2)": pref.str(14)}
    return w_, pref * w_, info


def wt(pt: Point, l: int) -> dict:
    """The record of Wt_l(omega0) and Res_{omega0} W_l."""
    w_, res, info = wt_ball(pt, l)
    return {"l": l, "Wt": w_.str(14), "Wt_nonzero": bool(not w_.contains(0)),
            "Wt_real_ball_excludes_0": bool(not w_.real.contains(0)), "Wt_imag_contains_0": bool(w_.imag.contains(0)),
            "residue_of_W_l": res.str(14), **info}


def v_positivity(pt: Point) -> dict:
    """L + r F' > 0 on [r_+, r_c] for L >= 2:  r F' = (2/r)(M - Q^2/r) - (2 Lambda/3) r^2, with
    r_+ > Q^2/M (the first term is >= 0 on [r_+, r_c]) and (2 Lambda/3) r_c^2 < 2; also kappa_+ < kappa_c
    (the limits of the written argument)."""
    m1 = pt.rp - pt.Q2                      # M = 1
    m2 = 2 - 2 * pt.Lam / 3 * pt.rc * pt.rc
    m3 = pt.kc - pt.kp
    return {"r_plus_minus_Q2_over_M": m1.str(12), "two_minus_2Lambda_rc2_over_3": m2.str(12),
            "kappa_c_minus_kappa_plus": m3.str(12), "pass": bool(m1 > 0 and m2 > 0 and m3 > 0)}


def b5(pt: Point) -> dict:
    """KAT B5: eps W_l(omega0 + eps) against Res_{omega0} W_l, with zsol and with the independent S1 code
    arbw; the relative difference must be <= 100 eps (it is O(eps): the regular part of W_l times eps)."""
    from . import arbw, zsol
    w0 = -I * pt.kp
    rows = []
    for l in B5_L:
        _, res, _ = wt_ball(pt, l)
        for er in B5_EPS:
            eps = acb(1, 1) / arb(2).sqrt() * A(er) * pt.km
            Wz = zsol.wronskian_z(pt, l, w0 + eps)[0]
            Wa = arbw.wronskian(pt, l, w0 + eps)
            dz = abs(eps * Wz - res) / abs(res)
            da = abs(eps * Wa - res) / abs(res)
            bound = 100 * A(er)
            rows.append({"l": l, "eps_over_kappa_minus": str(er), "direction": "(1+i)/sqrt(2)",
                         "eps_W_zsol": (eps * Wz).str(12), "eps_W_arbw": (eps * Wa).str(12), "residue": res.str(12),
                         "rel_diff_zsol": dz.str(4), "rel_diff_arbw": da.str(4),
                         "pass": bool(dz < bound and da < bound)})
    return {"rows": rows, "criterion": "|eps W - Res|/|Res| < 100 eps for both codes",
            "pass": all(r["pass"] for r in rows)}


def certificate() -> dict:
    pt = Point("PB", 256)
    out = {"claim": "Lemma E at P_B (notes/s6-preregistration.md O1, A1): omega0 = -i kappa_+ is not a resonance "
                    "of the forward problem on Omega_23 in any sector l; W_l has a simple pole there with residue "
                    "i(2 kappa_+ r_+ + L)/(2 r_+^2) Wt_l. Certified here: the first Frobenius step (sympy), Wt_l(omega0) "
                    "!= 0 for l <= 19 (Arb), the inequalities of the written argument for l >= 1, and KAT B5",
           "point": "PB", "omega0_over_kappa_minus": None, "normalisation": "Z_+ = sum a_m t^m, a_0 = 1, t = r - r_+ "
           "(zsol.frob_z); Z2 = t + sum_{m>=2} b_m t^m; Wt_l = [Delta (Z2 Zc' - Z2' Zc) + 2 i w r^2 Z2 Zc]/(r_+ r_c) "
           "at r = r_m (zsol's Jost-normalised form)", "prec_bits": PREC}
    with Timer() as tm:
        out["sympy_first_step"] = sympy_first_step()
        with precision(PREC):
            out["omega0_over_kappa_minus"] = (pt.kp / pt.km).str(15)
            out["V_l_positive_for_l_ge_1"] = v_positivity(pt)
            rows = {}
            for l in range(20):
                rows[l] = wt(pt, l)
                print("Wt", l, rows[l]["Wt"][:40], rows[l]["Wt_nonzero"], flush=True)
            out["Wt_l_le_19"] = [rows[l] for l in range(20)]
            out["Wt_samples_l_ge_20"] = [wt(pt, l) for l in SAMPLES]
            out["all_l_le_19_nonzero"] = all(rows[l]["Wt_nonzero"] for l in range(20))
            out["B5"] = b5(pt)
    out["timing_s"] = {"wall": round(tm.wall, 1), "cpu": round(tm.cpu, 1)}
    out["pass"] = bool(out["sympy_first_step"]["pass"] and out["V_l_positive_for_l_ge_1"]["pass"]
                       and out["all_l_le_19_nonzero"] and out["B5"]["pass"])
    out["commit"] = commit_tag()
    return out


def main(argv):
    rec = certificate()
    write_json(CERT / "epole_PB.json", rec)
    print("Lemma E certificate pass:", rec["pass"], "B5:", rec["B5"]["pass"], rec["timing_s"], rec["commit"])


if __name__ == "__main__":
    main(sys.argv)
