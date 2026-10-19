"""Bulk facts and the affine arithmetic (S4 review).

(1) Independent bulk facts (mpmath, no p2 code): min of W~ on [r_+, r_c] by golden-section search from a dense grid,
    against W_bound = 99/12500; U' sign on (r_+, r_P) and (r_P, r_c) on a dense grid; the identity
    Delta' r - 4 Delta = -2 (r^2 - 3 r + 2 Q^2) and the splitting V_l = (l+1/2)^2 U + W~ in sympy (re-derived here).
(2) A randomized enclosure test of p2.affine.AF (the author's class): random expression trees (products, quotients,
    sums, scalar multiples, powers by repeated products, nested inverses) over random affine inputs, evaluated at
    random admissible points (dw in the disc |dw| <= Rw, dh real |dh| <= Rh, each input's error term drawn in its disc)
    in mpmath at 50 digits; every value must lie in the output form's disc |c0 - value| <= dev (with c0's ball radius).
"""
from __future__ import annotations

import json
import random
import sys
from pathlib import Path

import mpmath as mp
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from indep_params import PB  # noqa: E402

W = Path(__file__).resolve().parents[2]


def bulk():
    p = PB(40)
    mp.mp.dps = 40
    n = 20000
    best = None
    for k in range(n + 1):
        r = p.rp + (p.rc - p.rp) * k / n
        v = p.Wt(r)
        if best is None or v < best[0]:
            best = (v, r)
    r0 = best[1]
    f = lambda r: p.Wt(r)
    rmin = mp.findroot(lambda r: mp.diff(f, r), r0)
    wmin = f(rmin)
    # U monotone: sign of U' on a grid
    up_ok = all(mp.diff(p.U, p.rp + (p.rP - p.rp) * k / 2000) > 0 for k in range(1, 2000))
    dn_ok = all(mp.diff(p.U, p.rP + (p.rc - p.rP) * k / 2000) < 0 for k in range(1, 2000))
    r, Lam, Q2, l = sp.symbols("r Lambda Q2 l")
    D = -Lam / 3 * r ** 4 + r ** 2 - 2 * r + Q2
    ident = sp.expand(sp.diff(D, r) * r - 4 * D + 2 * (r ** 2 - 3 * r + 2 * Q2)) == 0
    F = D / r ** 2
    V = F * (l * (l + 1) / r ** 2 + sp.diff(F, r) / r)
    Wt = F * sp.diff(F, r) / r - F / r ** 2 / 4
    split = sp.simplify(V - ((l + sp.Rational(1, 2)) ** 2 * F / r ** 2 + Wt)) == 0
    rPp = 3 - p.rP
    return {"Wt_min": mp.nstr(wmin, 12), "at_r": mp.nstr(rmin, 8), "W_bound": "99/12500",
            "W_bound_ok": bool(-wmin <= mp.mpf(99) / 12500), "margin": mp.nstr(mp.mpf(99) / 12500 + wmin, 4),
            "U_increasing_grid": up_ok, "U_decreasing_grid": dn_ok, "inner_root_3_minus_rP": mp.nstr(rPp, 12),
            "inner_root_lt_rplus": bool(rPp < p.rp), "sympy_identity": bool(ident), "sympy_splitting": bool(split)}


def affine_test(n_trees=300, n_pts=60, seed=11):
    sys.path.insert(0, str(W / "src"))
    from flint import acb, arb

    from p2.affine import AF
    from p2.arbutil import precision
    rng = random.Random(seed)
    mp.mp.dps = 50
    worst = 0.0
    fails = 0
    n_eval = 0
    with precision(128):
        for t in range(n_trees):
            Rw = rng.choice([1e-4, 1e-3, 1e-2, 0.05])
            Rh = rng.choice([1e-6, 1e-4, 1e-3, 1e-2])
            ins = []
            for _ in range(3):
                c = [complex(rng.uniform(-2, 2), rng.uniform(-2, 2)) for _ in range(3)]
                e = rng.choice([0.0, 1e-8, 1e-5])
                ins.append((c, e))
            Rwa, Rha = arb(repr(Rw)), arb(repr(Rh))
            X = [AF(acb(c[0].real, c[0].imag), acb(c[1].real, c[1].imag), acb(c[2].real, c[2].imag), e, Rwa, Rha) for c, e in ins]
            # random expression program over the three inputs
            prog = []
            for _ in range(rng.randint(3, 9)):
                op = rng.choice(["add", "mul", "div", "inv", "smul", "sub", "sq"])
                prog.append((op, rng.randrange(10 ** 6), rng.randrange(10 ** 6), complex(rng.uniform(-3, 3), rng.uniform(-3, 3))))

            def run_af(vals):
                st = list(vals)
                for op, i, j, s in prog:
                    a, b = st[i % len(st)], st[j % len(st)]
                    try:
                        if op == "add":
                            st.append(a + b)
                        elif op == "sub":
                            st.append(a - b)
                        elif op == "mul":
                            st.append(a * b)
                        elif op == "div":
                            st.append(a / b)
                        elif op == "inv":
                            st.append(a.inv())
                        elif op == "smul":
                            st.append(a * acb(s.real, s.imag))
                        elif op == "sq":
                            st.append(a * a)
                    except ZeroDivisionError:
                        st.append(a)
                return st

            def run_pt(vals, af_states):
                st = list(vals)
                k = 3
                for op, i, j, s in prog:
                    a, b = st[i % len(st)], st[j % len(st)]
                    # mirror the AF path's fallbacks: if the AF op raised, the AF appended a
                    afo = af_states[k]
                    if op == "add":
                        v = a + b
                    elif op == "sub":
                        v = a - b
                    elif op == "mul":
                        v = a * b
                    elif op == "div":
                        v = a / b if afo is not None else a
                    elif op == "inv":
                        v = 1 / a if afo is not None else a
                    elif op == "smul":
                        v = a * mp.mpc(s.real, s.imag)
                    else:
                        v = a * a
                    st.append(v)
                    k += 1
                return st

            # record which div/inv raised
            st_af = list(X)
            raised = [None, None, None]
            for op, i, j, s in prog:
                a, b = st_af[i % len(st_af)], st_af[j % len(st_af)]
                try:
                    if op == "add":
                        r = a + b
                    elif op == "sub":
                        r = a - b
                    elif op == "mul":
                        r = a * b
                    elif op == "div":
                        r = a / b
                    elif op == "inv":
                        r = a.inv()
                    elif op == "smul":
                        r = a * acb(s.real, s.imag)
                    else:
                        r = a * a
                    raised.append(True)
                except ZeroDivisionError:
                    r = a
                    raised.append(None)
                st_af.append(r)
            for _ in range(n_pts):
                ang, rad = rng.uniform(0, 2 * mp.pi), Rw * mp.sqrt(rng.random())
                dw = mp.mpc(rad * mp.cos(ang), rad * mp.sin(ang))
                dh = mp.mpf(rng.uniform(-Rh, Rh))
                vals = []
                for (c, e) in ins:
                    a2, r2 = rng.uniform(0, 2 * mp.pi), e * mp.sqrt(rng.random())
                    E = mp.mpc(r2 * mp.cos(a2), r2 * mp.sin(a2))
                    vals.append(mp.mpc(c[0]) + mp.mpc(c[1]) * dw + mp.mpc(c[2]) * dh + E)
                st = run_pt(vals, raised)
                for af, v in zip(st_af, st):
                    c0 = mp.mpc(float(af.c0.real.mid()), float(af.c0.imag.mid()))
                    r0 = float(abs(af.c0).rad()) if False else float(af.c0.real.rad()) + float(af.c0.imag.rad())
                    dev = float(af.dev().upper()) + r0
                    dist = float(abs(v - c0))
                    n_eval += 1
                    if dist > dev * (1 + 1e-12) + 1e-30:
                        fails += 1
                    if dev > 0:
                        worst = max(worst, dist / dev)
    return {"trees": n_trees, "points_per_tree": n_pts, "evaluations": n_eval, "violations": fails,
            "max_distance_over_dev": worst}


if __name__ == "__main__":
    out = {"bulk": bulk(), "affine": affine_test()}
    print(json.dumps(out, indent=1))
    Path(sys.argv[1] if len(sys.argv) > 1 else "bulk_and_affine.json").write_text(json.dumps(out, indent=1) + "\n")
