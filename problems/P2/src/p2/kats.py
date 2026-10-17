"""Known-answer tests KAT-1 ... KAT-7 (FORMULATION.md section 5).  Floats are evidence only.

Pass criteria were fixed in this docstring BEFORE any KAT was run (S1, 2026-10-16) and are not
to be loosened after the fact:

KAT-1  [CCDHJ18] Table I at P_A, all nine printed values (l = 0 NE; l = 1 PS and dS; l = 2, 3,
       4, 5, 10, 100 PS).  PASS iff for every value, Re and Im of the computed omega/kappa_-
       round to the printed 10 decimals (|component difference| <= 5e-11), and the solver's
       internal consistency check (dps 50, K 6  vs  dps 70, K 3) differs by < 1e-20.
KAT-2  Eikonal consistency at P_A.  (a) [certified, test_params] gamma0/kappa_- rounds to
       0.7643026660 and is < 0.7643094446.  (b) D_l = (l+1/2)^2 (-Im w_l/k_- - gamma0/k_-) > 0
       for l in {10, 20, 50, 100} and |D_l/D_100 - 1| <= 0.05 for l in {10, 20, 50} (O(l^-2)).
       (c) E_l = (l+1/2)(Re w_l/k_- - (l+1/2) Omega_P/k_-) satisfies |E_l/E_100 - 1| <= 0.05
       for l in {10, 20, 50} (Re w/(l+1/2) -> Omega_P with O(l^-2)).  (d) At l = 100:
       (l+1/2) Omega_P/k_- rounds to 222.0614 and Re w/k_- rounds to 222.0603 (4 decimals).
KAT-3  [CCDHJ18] Fig. 3, Lambda M^2 = 0.06, l = 10 PS proxy.  q1: beta_PS10(q1) = 1/2;
       q2: beta_NE0(q2) = beta_PS10(q2).  PASS iff round(q1, 4) = 0.9920 and
       round(q2, 4) = 0.9967, and at q1 the PS l=10 mode is the smallest of
       {PS l=10, NE l=0, dS l=1} in -Im w/k_- (so that it sets beta there).
KAT-4  [Hin25] Thm 1.1.  At q in {0.999, 0.9999, 0.99999} (Lambda = 3/50) the NE modes
       (l, n) in {(0,0), (0,1), (1,0), (2,0)}: e(q) = |w/k_- + i (lambda_l^+ + n)| with varkappa
       at the current r_+, r_c.  PASS iff, for every mode, e strictly decreases along the three
       q values and e(0.99999) <= 0.5 e(0.9999).
KAT-5  [CCDHJ18] supplement, near-Nariai at Lambda M^2 = 0.14 (> 1/9): omega/kappa_+ ->
       sqrt(-1 + 2 l(l+1) Upsilon)/2 - i/2 (n = 0), Upsilon = r_N/(r_N - 3/2) (also checked
       against the printed closed form in Lambda).  For l in {1, 2} and
       eps = (r_c - r_+)/r_N in {0.04, 0.02, 0.01, 0.005}: e(eps) = |w/k_+ - w_N|.
       PASS iff e strictly decreases along the four eps and e(0.005) <= 0.65 e(0.01).
KAT-6  (validated, winding.py)  winding number 1 of W_l around small boxes containing the
       Table I l = 1 PS mode (P_A) and the P_B l = 0 NE mode; Krawczyk enclosures contain the
       float values (and the Table I printed value lies within 5e-11 per component).
KAT-7  (spectral.py)  hyperboloidal Chebyshev collocation agrees with the Frobenius/Taylor
       solver to >= 10 significant digits (|diff| <= 1e-10 |w|) on all eight P_B modes of
       FORMULATION 1.3; the spectral values must also be converged in N (two N differ < 1e-12 |w|).
"""
from __future__ import annotations

import itertools
import json
import sys
import time
from pathlib import Path

import mpmath as mp

from .arbutil import git_commit
from .qnm_float import Geom, Solver, qnm

RESULTS = Path(__file__).resolve().parents[2] / "results"

TABLE_I = [  # CCDHJ18 Table I (supplement), omega_0/kappa_-, Lambda M^2 = 0.06, Q/Q_max = 0.996
    ("NE l=0", 0, "0", "-0.8539013779"),
    ("PS l=1", 1, "3.2426164126", "-0.7958326323"),
    ("dS l=1", 1, "0", "-1.5003853731"),
    ("PS l=2", 2, "5.4796067815", "-0.7754179185"),
    ("PS l=3", 3, "7.7016152057", "-0.7699348317"),
    ("PS l=4", 4, "9.9181960834", "-0.7676996545"),
    ("PS l=5", 5, "12.1322536641", "-0.7665731858"),
    ("PS l=10", 10, "23.1897597770", "-0.7649242108"),
    ("PS l=100", 100, "222.0602900249", "-0.7643094446"),
]

PB_SPECTRUM = [  # FORMULATION 1.3 (S0 floats), omega/kappa_- at P_B
    ("NE l=0 n=0", 0, "0", "-0.909381"),
    ("NE l=0 n=1", 0, "0", "-1.855200"),
    ("NE l=1", 1, "0", "-1.911925"),
    ("dS l=1", 1, "0", "-3.324477"),
    ("PS l=1", 1, "7.214864", "-1.753306"),
    ("NE l=2", 2, "0", "-2.891324"),
    ("PS l=2", 2, "12.196453", "-1.707868"),
    ("PS l=10", 10, "51.62353", "-1.68459"),
]


def _round_ok(x, printed):
    nd = len(printed.split(".")[1]) if "." in printed else 0
    return abs(x - mp.mpf(printed)) <= mp.mpf(5) * mp.mpf(10) ** (-(nd + 1)), abs(x - mp.mpf(printed))


def _c(z):
    return [mp.nstr(mp.re(z), 25), mp.nstr(mp.im(z), 25)]


# ----------------------------------------------------------------------------------------- KAT-1
def kat1():
    g = Geom.frozen("PA")
    rows, ok_all = [], True
    for label, l, re_s, im_s in TABLE_I:
        t0 = time.time()
        ref = mp.mpc(re_s, im_s)
        v, err = qnm(g, l, ref, dps=50, K=6)
        okr, dr = _round_ok(mp.re(v), re_s)
        oki, di = _round_ok(mp.im(v), im_s)
        ok = okr and oki and err < mp.mpf("1e-20")
        ok_all &= ok
        rows.append({"mode": label, "l": l, "printed": [re_s, im_s], "computed": _c(v),
                     "abs_diff_re": mp.nstr(dr, 3), "abs_diff_im": mp.nstr(di, 3),
                     "self_consistency": mp.nstr(err, 3), "pass": bool(ok), "seconds": round(time.time() - t0, 1)})
        print("KAT-1", label, rows[-1]["computed"], rows[-1]["abs_diff_re"], rows[-1]["abs_diff_im"], ok, flush=True)
    return {"pass": bool(ok_all), "rows": rows}


# ----------------------------------------------------------------------------------------- KAT-2
def kat2(kat1_rows=None):
    g = Geom.frozen("PA")
    km = g.km
    gam = g.gamma0 / km
    Om = g.OmegaP / km
    vals = {}
    known = {row["l"]: mp.mpc(*row["computed"]) for row in (kat1_rows or []) if row["mode"].startswith("PS")}
    for l in (10, 20, 50, 100):
        if l in known:
            vals[l] = known[l]
            continue
        seed = mp.mpc((l + mp.mpf(1) / 2) * Om - mp.mpf("0.112") / (l + mp.mpf(1) / 2),
                      -(gam + mp.mpf("0.0685") / (l + mp.mpf(1) / 2) ** 2))
        vals[l], _ = qnm(g, l, seed, dps=50, K=6)
        print("KAT-2 l", l, _c(vals[l]), flush=True)
    D = {l: (l + mp.mpf(1) / 2) ** 2 * (-mp.im(v) - gam) for l, v in vals.items()}
    E = {l: (l + mp.mpf(1) / 2) * (mp.re(v) - (l + mp.mpf(1) / 2) * Om) for l, v in vals.items()}
    ok_b = all(D[l] > 0 for l in D) and all(abs(D[l] / D[100] - 1) <= mp.mpf("0.05") for l in (10, 20, 50))
    ok_c = all(abs(E[l] / E[100] - 1) <= mp.mpf("0.05") for l in (10, 20, 50))
    pred = mp.mpf("100.5") * Om
    ok_d = _round_ok(pred, "222.0614")[0] and _round_ok(mp.re(vals[100]), "222.0603")[0]
    return {"pass": bool(ok_b and ok_c and ok_d),
            "gamma0_over_km": mp.nstr(gam, 20), "Omega_P_over_km": mp.nstr(Om, 20),
            "modes": {str(l): _c(v) for l, v in vals.items()},
            "D_l": {str(l): mp.nstr(x, 8) for l, x in D.items()}, "E_l": {str(l): mp.nstr(x, 8) for l, x in E.items()},
            "eikonal_Re_pred_l100": mp.nstr(pred, 10),
            "checks": {"b_Ol-2_imag": bool(ok_b), "c_Re_over_l": bool(ok_c), "d_l100_digits": bool(ok_d)}}


# ----------------------------------------------------------------------------------------- KAT-3
def _mode_abs(q, l, seed_abs, dps=40, guard="0.03"):
    """Mode at Q/Q_max = q from a seed in ABSOLUTE units (PS modes move little in absolute units;
    NE modes scale with kappa_-, so their seeds are rescaled by the caller).  Mode-hopping guard:
    the result must lie within `guard` (relative) of the seed."""
    guard = mp.mpf(guard)
    with mp.workdps(dps):
        g = Geom.from_qfrac(mp.mpf(3) / 50, q)
        w = Solver(g, K=6).find(l, mp.mpc(seed_abs), tol=mp.mpf("1e-18"))
        if abs(w - seed_abs) > guard * abs(seed_abs):
            raise ValueError(f"mode hop at q={mp.nstr(q, 8)}, l={l}: seed {seed_abs} -> {w}")
        return w, g


def kat3():
    gA = Geom.frozen("PA")
    # P_A values (Table I), absolute units
    st = {"ps": mp.mpc("23.1897597770", "-0.7649242108") * gA.km, "ne": mp.mpc(0, "-0.8539013779") * gA.km,
          "ne_km": None}

    def ps(q):
        w, g = _mode_abs(q, 10, st["ps"])
        st["ps"] = w
        return w / g.km, g

    def ne(q, km_new):
        # NE modes scale with kappa_-: rescale the seed
        seed = st["ne"] / st["ne_km"] * km_new if st["ne_km"] else st["ne"]
        w, g = _mode_abs(q, 0, seed, guard="0.05")
        st["ne"], st["ne_km"] = w, g.km
        return w / g.km

    st["ne_km"] = gA.km

    def f_ps(q):
        v, _ = ps(q)
        return -mp.im(v) - mp.mpf(1) / 2

    def f_ne(q):
        vp, g = ps(q)
        vn = ne(q, g.km)
        return mp.im(vp) - mp.im(vn)          # = beta_NE - beta_PS

    for k in range(1, 9):                    # 0.9955 ... 0.9920 in steps of 0.0005
        f_ps(mp.mpf("0.996") - k * mp.mpf("0.0005"))
    q1 = mp.findroot(f_ps, (mp.mpf("0.992"), mp.mpf("0.9921")), solver="secant", tol=mp.mpf("1e-30"), verify=False)
    if abs(f_ps(q1)) > mp.mpf("1e-12"):
        raise ValueError("KAT-3: q1 not converged")
    v_ps_q1, g1 = ps(q1)
    # NE l=0 and dS l=1 at q1 by continuation from P_A in small steps
    st_ne, st_ne_km = mp.mpc(0, "-0.8539013779") * gA.km, gA.km
    st_ds = mp.mpc(0, "-1.5003853731") * gA.km
    qq = mp.mpf("0.996")
    while qq > q1:
        qq = max(q1, qq - mp.mpf("0.0005"))
        gq = Geom.from_qfrac(mp.mpf(3) / 50, qq)
        st_ne, _ = _mode_abs(qq, 0, st_ne / st_ne_km * gq.km, guard="0.05")
        st_ne_km = gq.km
        st_ds, _ = _mode_abs(qq, 1, st_ds, guard="0.05")
    vn_q1, vd_q1 = st_ne / g1.km, st_ds / g1.km
    dom_ok = (-mp.im(v_ps_q1) < -mp.im(vn_q1)) and (-mp.im(v_ps_q1) < -mp.im(vd_q1))
    # NE takeover: march up from P_A in steps of 0.0001
    st["ps"], st["ne"], st["ne_km"] = mp.mpc("23.1897597770", "-0.7649242108") * gA.km, mp.mpc(0, "-0.8539013779") * gA.km, gA.km
    for k in range(1, 8):
        f_ne(mp.mpf("0.996") + k * mp.mpf("0.0001"))
    q2 = mp.findroot(f_ne, (mp.mpf("0.9966"), mp.mpf("0.9967")), solver="secant", tol=mp.mpf("1e-30"), verify=False)
    if abs(f_ne(q2)) > mp.mpf("1e-12"):
        raise ValueError("KAT-3: q2 not converged")
    v_ps_q2, g2 = ps(q2)
    v_ne_q2 = st["ne"] / g2.km
    # dS l=1 at q2 (continuation upward)
    st_ds = mp.mpc(0, "-1.5003853731") * gA.km
    qq = mp.mpf("0.996")
    while qq < q2:
        qq = min(q2, qq + mp.mpf("0.0002"))
        st_ds, _ = _mode_abs(qq, 1, st_ds, guard="0.05")
    vd_q2 = st_ds / g2.km
    ok1 = round(float(q1), 4) == 0.992
    ok2 = round(float(q2), 4) == 0.9967
    return {"pass": bool(ok1 and ok2 and dom_ok),
            "q1_beta_half": mp.nstr(q1, 12), "q2_NE_takeover": mp.nstr(q2, 12),
            "round4": [round(float(q1), 4), round(float(q2), 4)],
            "at_q1": {"PS10": _c(v_ps_q1), "NE0": _c(vn_q1), "dS1": _c(vd_q1), "PS_dominant": bool(dom_ok)},
            "at_q2": {"PS10": _c(v_ps_q2), "NE0": _c(v_ne_q2), "dS1": _c(vd_q2),
                      "dS_above": bool(-mp.im(vd_q2) > -mp.im(v_ne_q2))}}


# ----------------------------------------------------------------------------------------- KAT-4
def kat4():
    modes = [(0, 0, mp.mpc(0, "-0.909381")), (0, 1, mp.mpc(0, "-1.855200")),
             (1, 0, mp.mpc(0, "-1.911925")), (2, 0, mp.mpc(0, "-2.891324"))]
    qs = ["0.999", "0.9999", "0.99999"]
    path = ["0.999", "0.9995", "0.9998", "0.9999", "0.99995", "0.99998", "0.99999"]   # continuation in q
    table, ok_all = [], True
    for l, n, seed in modes:
        errs, vals, lims = [], [], []
        s = seed
        for qq in path:
            g = Geom.from_qfrac(mp.mpf(3) / 50, mp.mpf(qq))
            vk = (g.rc ** 2 + 2 * g.rp * g.rc - 3 * g.rp ** 2) / (g.rc ** 2 + 2 * g.rp * g.rc + 3 * g.rp ** 2)
            lam = (1 + mp.sqrt(1 + 4 * l * (l + 1) / vk)) / 2
            target = mp.mpc(0, -(lam + n))
            v, err = qnm(g, l, s, dps=40, K=6, check=qq in qs)
            s = v
            if qq not in qs:
                continue
            vals.append(v)
            lims.append(target)
            errs.append(abs(v - target))
            print("KAT-4", l, n, qq, _c(v), mp.nstr(errs[-1], 6), "self", mp.nstr(err, 3),
                  "k+/k-", mp.nstr(g.kp / g.km, 6), flush=True)
        dec = errs[0] > errs[1] > errs[2]
        rate = errs[2] <= errs[1] / 2
        ok = dec and rate
        ok_all &= ok
        table.append({"l": l, "n": n, "q": qs, "omega_over_km": [_c(v) for v in vals],
                      "limit": [_c(t) for t in lims], "error": [mp.nstr(e, 8) for e in errs],
                      "ratio_last": mp.nstr(errs[2] / errs[1], 6), "monotone": bool(dec), "halving": bool(rate),
                      "pass": bool(ok)})
    return {"pass": bool(ok_all), "modes": table}


# ----------------------------------------------------------------------------------------- KAT-5
def _nariai_geom(Lam, eps):
    """Q^2 slightly above Q_min^2 with (r_c - r_+)/r_N = eps (solved by secant)."""
    Lam = mp.mpf(Lam)
    _r_ext, r_N = Geom.extremal_roots(Lam)
    Qmin2 = r_N - Lam * r_N ** 4 / 3

    def gap(t):
        g = Geom(Lam, Qmin2 + t)
        return (g.rc - g.rp) / r_N - eps

    # near Nariai, r_c - r_+ ~ sqrt(t): bisection in log t (robust)
    lo_, hi_ = mp.mpf(10) ** -40, mp.mpf("1e-2")
    while gap(hi_) < 0:
        hi_ *= 2
    for _ in range(200):
        mid = mp.sqrt(lo_ * hi_)
        if gap(mid) < 0:
            lo_ = mid
        else:
            hi_ = mid
        if hi_ / lo_ - 1 < mp.mpf(10) ** -25:
            break
    return Geom(Lam, Qmin2 + hi_), r_N, Qmin2


def kat5():
    Lam = mp.mpf("0.14")
    _r_ext, r_N = Geom.extremal_roots(Lam)
    ups = r_N / (r_N - mp.mpf(3) / 2)
    # printed closed form (CCDHJ18 supplement)
    gam = (-3 * Lam ** 2 + Lam ** mp.mpf(1.5) * mp.sqrt(9 * Lam - 2)) ** (mp.mpf(1) / 3)
    ups_cf = (gam ** 2 + mp.cbrt(2) * Lam) / (gam ** 2 + (mp.cbrt(2) - 3 / mp.cbrt(2) * gam) * Lam)
    out = {"Lambda": "0.14", "r_N": mp.nstr(r_N, 20), "Upsilon": mp.nstr(ups, 20),
           "Upsilon_closed_form": mp.nstr(ups_cf, 20), "Upsilon_forms_agree": bool(abs(ups - ups_cf) < mp.mpf("1e-15")),
           "modes": []}
    ok_all = out["Upsilon_forms_agree"]
    for l in (1, 2):
        wN = mp.sqrt(-1 + 2 * l * (l + 1) * ups) / 2 - mp.mpc(0, mp.mpf(1) / 2)
        errs, vals = [], []
        seed = None
        for eps in ("0.04", "0.02", "0.01", "0.005"):
            g, _, _ = _nariai_geom(Lam, mp.mpf(eps))
            s = wN * g.kp / g.km if seed is None else seed * g.kp / g.km
            v, _err = qnm(g, l, s, dps=40, K=6)
            vk = v * g.km / g.kp
            seed = vk
            vals.append(vk)
            errs.append(abs(vk - wN))
            print("KAT-5", l, eps, _c(vk), mp.nstr(errs[-1], 6), flush=True)
        dec = all(x > y for x, y in itertools.pairwise(errs))
        rate = errs[3] <= mp.mpf("0.65") * errs[2]
        ok = dec and rate
        ok_all &= ok
        out["modes"].append({"l": l, "n": 0, "omega_N_over_kp": _c(wN), "eps": ["0.04", "0.02", "0.01", "0.005"],
                             "omega_over_kp": [_c(v) for v in vals], "error": [mp.nstr(e, 8) for e in errs],
                             "ratio_last": mp.nstr(errs[3] / errs[2], 6), "monotone": bool(dec),
                             "rate": bool(rate), "pass": bool(ok)})
    out["pass"] = bool(ok_all)
    return out


# ----------------------------------------------------------------------------------------- KAT-7
def kat7():
    """Spectral collocation (spectral.py) vs Frobenius/Taylor (qnm_float) on the P_B spectrum."""
    import numpy as np

    from .spectral import Collocation
    g = Geom.frozen("PB")
    rows, ok_all = [], True
    cache = {}
    for label, l, re_s, im_s in PB_SPECTRUM:
        seed = complex(float(re_s), float(im_s))
        Ns = (80, 100, 120) if l >= 10 else (60, 80, 100)
        spec = []
        for N in Ns:
            key = (l, N)
            if key not in cache:
                cache[key] = Collocation(g, l, N, prec=256)
            c = cache[key]
            ev = c.eig_double()
            k = int(np.argmin(np.abs(ev - seed)))
            _, wb, _last = c.refine(ev[k])
            spec.append(mp.mpc(mp.mpf(wb.real.mid().str(40, radius=False)), mp.mpf(wb.imag.mid().str(40, radius=False))))
        conv = abs(spec[-1] - spec[-2]) / abs(spec[-1])
        frob, err = qnm(g, l, mp.mpc(re_s, im_s), dps=50, K=6)
        rel = abs(spec[-1] - frob) / abs(frob)
        ok = conv < mp.mpf("1e-12") and rel <= mp.mpf("1e-10")
        ok_all &= ok
        rows.append({"mode": label, "l": l, "S0_value": [re_s, im_s], "spectral_N": list(Ns),
                     "spectral": [_c(v) for v in spec], "spectral_convergence_rel": mp.nstr(conv, 3),
                     "frobenius": _c(frob), "frobenius_self_consistency": mp.nstr(err, 3),
                     "rel_diff": mp.nstr(rel, 3), "digits_agree": mp.nstr(-mp.log10(rel), 4) if rel > 0 else "inf",
                     "pass": bool(ok)})
        print("KAT-7", label, rows[-1]["spectral"][-1], rows[-1]["frobenius"], rows[-1]["rel_diff"], ok, flush=True)
    # extra evidence: converged spectral eigenvalues inside the P_B strip box (float, not a proof)
    scan = {}
    for l in (0, 1, 2, 10):
        e1 = Collocation(g, l, 80, prec=128).eig_double()
        e2 = Collocation(g, l, 100, prec=128).eig_double()
        keep = []
        for z in e2:
            if z.imag > -0.6 and abs(z) < 100:
                d = np.min(np.abs(e1 - z))
                if d < 1e-6 * max(1.0, abs(z)):
                    keep.append([round(z.real, 8), round(z.imag, 8)])
        scan[str(l)] = keep
    return {"pass": bool(ok_all), "rows": rows,
            "strip_scan_Im_gt_-0.6_kappa_minus": {"note": "double-precision eigenvalues stable between N=80 and N=100",
                                                   "eigenvalues_over_kappa_minus": scan}}


def run(which):
    mp.mp.dps = 50
    res = {"commit": git_commit(), "criteria": "see p2/kats.py docstring (fixed before running)"}
    t0 = time.time()
    k1 = None
    for k in which:
        t = time.time()
        if k == "1":
            res["KAT-1"] = k1 = kat1()
        elif k == "2":
            res["KAT-2"] = kat2(k1["rows"] if k1 else None)
        elif k == "3":
            res["KAT-3"] = kat3()
        elif k == "4":
            res["KAT-4"] = kat4()
        elif k == "5":
            res["KAT-5"] = kat5()
        elif k == "7":
            res["KAT-7"] = kat7()
        res[f"KAT-{k}"]["seconds"] = round(time.time() - t, 1)
        print(f"KAT-{k} pass={res[f'KAT-{k}']['pass']} ({res[f'KAT-{k}']['seconds']} s)", flush=True)
    res["total_seconds"] = round(time.time() - t0, 1)
    return res


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "12345"
    out = run(which)
    RESULTS.mkdir(exist_ok=True)
    (RESULTS / f"kat_float_{which}.json").write_text(json.dumps(out, indent=1) + "\n")
