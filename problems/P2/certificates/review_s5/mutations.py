#!/usr/bin/env python
"""S5 review: extra negative controls (mutations) for the author's assembly check p2.assembly.check.

N1-N8 are the author's controls.  M1-M22 probe the checker on items N1-N8 do not isolate: a gap from
the P2-2 side, a missing lg certificate, a broken or gapped path, a wrong increment, a certificate for
the wrong l, an uncertified lg criterion, a narrower lg strip, uncertified winding flags, cover tiling
and design flags, a missing regime file, the Im range of L and T, an oversized disc, beta and params.
B1-B2 document what the check does NOT catch by design (it trusts each ingredient certificate's own
replay).  Each mutation acts on an in-memory deep copy; nothing on disk is changed.

Run:  python mutations.py   (writes mutations_out.json next to this file)
"""
from __future__ import annotations

import copy
import json
import sys
import time
from fractions import Fraction as Fr
from pathlib import Path

HERE = Path(__file__).resolve().parent
P2 = HERE.parents[1]
sys.path.insert(0, str(P2 / "src"))

from p2 import assembly  # noqa: E402  (the code under test)
from p2.params import Point  # noqa: E402

OUT = HERE / "mutations_out.json"


def shift_coord(wd, old, new):
    for s in wd["pieces"]:
        for k in ("p0", "p1", "c"):
            s[k] = [new if v == old else v for v in s[k]]


def m1(J):   # P2-2 threshold raised: a gap from the P2-2 side at l = 11
    lg = J["lg:11"]
    lg["W_l_exact"] = str(Fr(lg["W_l_exact"]) + Fr(1, 2 ** 24))


def m2(J):   # one piece removed from the middle of the l = 4 path
    p = J["wind:4"]["pieces"]
    del p[len(p) // 2]


def m3(J):   # one increment off by pi at l = 9
    s = J["wind:9"]["pieces"][5]
    m = s["increment"].strip("[]").split("+/-")
    s["increment"] = f"[{float(m[0]) + 3.14159265358979} +/- {m[1].strip()}]"


def m4(J):   # the l = 6 certificate stored in l = 7's slot
    J["wind:7"] = copy.deepcopy(J["wind:6"])


def m5(J):   # lg certificate for l = 12 missing
    del J["lg:12"]


def m6(J):   # lg criterion not below 2 at l = 3
    J["lg:3"]["criterion"]["value_upper"] = "[2.00010000000 +/- 1.00e-12]"


def m7(J):   # lg certificate for a narrower strip (0.5 kappa_-) at l = 5
    km = Fr(J["params"]["points"]["PB"]["quantities"]["kappa_minus"]["approx"])
    J["lg:5"]["strip_depth_beta"] = f"[{float(km / 2):.20f} +/- 1e-22]"


def m8(J):   # winding certificate flagged uncertified at l = 0
    J["wind:0"]["certified"] = False


def m9(J):   # an S tiling flag false in the cover
    t = J["cover"]["tiling"]
    k = [x for x in t if x.startswith("S plus")][0]
    t[k] = False


def m10(J):  # regime certificate L_PB_plus missing
    del J["regime:L_PB_plus"]


def m11(J):  # L's Im range narrower: [-1/2, 0]
    J["regime:L_PB_c"]["design"]["Im_range_over_km"] = ["-1/2", "0"]


def m12(J):  # a T box whose lower Im edge is -0.5 kappa_- (does not reach -d_delta)
    km = Fr(J["params"]["points"]["PB"]["quantities"]["kappa_minus"]["approx"])
    J["regime:T_PB"]["boxes"][3]["nu_low"] = f"[{-float(km / 2):.20f} +/- 1e-22]"


def m13(J):  # an oversized disc at l = 2 (endpoints still within rho/2, but the box reaches -i kappa_+)
    J["wind:2"]["pieces"][0]["rho"] = "1"


def m14(J):  # a 2^-40 gap between two pieces at l = 10
    p = J["wind:10"]["pieces"]
    x, y = p[3]["p1"]
    p[3]["p1"] = [str(Fr(x) - Fr(1, 2 ** 40)) if Fr(x) != 0 else x, y] if p[3]["edge"] != 1 else [x, str(Fr(y) - Fr(1, 2 ** 40))]


def m15(J):  # l = 13 region b raised (still >= d_delta) but the path not moved
    J["wind:13"]["region"]["b"] = str(Fr(J["wind:13"]["region"]["b"]) + Fr(1, 2 ** 20))


def m16(J):  # krawczyk enclosure moved to -0.95 i kappa_-
    J["krawczyk"]["enclosure_over_kappa_minus"]["im"] = "[-0.95 +/- 1e-30]"


def m17(J):  # params.json kappa_+ inconsistent with the recomputed point (0.6 kappa_-)
    q = J["params"]["points"]["PB"]["quantities"]
    km = Fr(q["kappa_minus"]["approx"])
    q["kappa_plus"] = {"ball": f"[{float(km * Fr(6, 10)):.25f} +/- 1e-30]"}


def m18(J):  # stored winding 0 at l = 0 (the zero at 0 not counted)
    J["wind:0"]["winding"] = 0


def m19(J):  # cover: zeta_SU = 1/10 (S no longer meets T)
    J["cover"]["constants"]["zeta_SU"] = "1/10"


def m20(J):  # cover: h0 = 2/39 (l0 = 19; overlapping, so harmless, but not the registered value)
    J["cover"]["constants"]["h0"] = "2/39"


def m21(J):  # a P2-0 identity step recorded as failed
    J["p2_0"]["steps"][0]["ok"] = False


def m22(J):  # the symmetry identity Z4 recorded as failed
    J["zeq"]["Z4"] = False


def b1(J):   # blind spot by design: a_l raised consistently in lg, the box and the path at l = 8
    lg, wd = J["lg:8"], J["wind:8"]
    old = wd["region"]["a"]
    new = str(Fr(old) + 1)
    lg["W_l_exact"] = new
    wd["region"]["a"] = new
    shift_coord(wd, old, new)


def b3(J):   # blind spot by design: the lg matching point (analytic content) changed at l = 6
    J["lg:6"]["matching_point_exact"] = "1/2"


def b2(J):   # blind spot by design: one increment changed by 0.1 (resum still within 1/4 of 0 at l = 15)
    s = J["wind:15"]["pieces"][7]
    m = s["increment"].strip("[]").split("+/-")
    s["increment"] = f"[{float(m[0]) + 0.1} +/- {m[1].strip()}]"


SPEC = [
    ("M1", "P2-2 threshold raised (gap from the P2-2 side), l = 11", m1, ["C2.l11.a_identical"]),
    ("M2", "a piece removed from the l = 4 path", m2, ["C2.l4.path_traces_H"]),
    ("M3", "one increment off by pi, l = 9", m3, ["C2.l9.resum"]),
    ("M4", "l = 6 certificate in l = 7's slot", m4, ["C2.l7.wind_meta", "C2.l7.a_identical"]),
    ("M5", "lg certificate l = 12 missing", m5, ["C2.l12.present", "C2.all_l_0_to_19_present"]),
    ("M6", "lg criterion upper >= 2, l = 3", m6, ["C2.l3.lg_certified"]),
    ("M7", "lg strip depth 0.5 kappa_-, l = 5", m7, ["C2.l5.lg_strip_depth"]),
    ("M8", "winding certificate uncertified, l = 0", m8, ["C2.l0.wind_meta"]),
    ("M9", "an S tiling flag false", m9, ["C3.cover_pass", "C3.Im_range"]),
    ("M10", "L_PB_plus missing", m10, ["C3.regime_pass.L_PB_plus"]),
    ("M11", "L Im range [-1/2, 0]", m11, ["C3.design_L", "C3.Im_range"]),
    ("M12", "a T box not reaching -d_delta", m12, ["C3.design_T", "C3.Im_range"]),
    ("M13", "an oversized disc, l = 2", m13, ["C2.l2.discs_clear_poles"]),
    ("M14", "a 2^-40 gap in the l = 10 path", m14, ["C2.l10.path_traces_H"]),
    ("M15", "l = 13 b raised, path not moved", m15, ["C2.l13.path_traces_H"]),
    ("M16", "Krawczyk enclosure at -0.95 i kappa_-", m16, ["C6.beta_upper"]),
    ("M17", "params kappa_+ inconsistent (0.6 kappa_-)", m17, ["C1.overlap.kappa_plus"]),
    ("M18", "stored winding 0 at l = 0", m18, ["C2.l0.winding_value", "C2.l0.resum"]),
    ("M19", "cover zeta_SU = 1/10", m19, ["C3.S_T_overlap"]),
    ("M20", "cover h0 = 2/39", m20, ["C3.h0_is_2_41", "C3.no_l_gap"]),
    ("M21", "a P2-0 identity step failed", m21, ["C5.P2_0_identities"]),
    ("M22", "Z4 failed", m22, ["C3.symmetry_Z4"]),
    ("B1", "a_8 raised by 1 in lg, the box and the path coordinates equal to a_8", b1, []),
    ("B2", "BY DESIGN NOT CAUGHT: one increment +0.1 at l = 15", b2, []),
    ("B3", "BY DESIGN NOT CAUGHT: lg matching point changed at l = 6", b3, []),
]


def main():
    t0 = time.time()
    inp = assembly.load_inputs()
    pt = Point("PB", assembly.PREC)
    cov = assembly.replay_cover()
    base = assembly.check(inp, cov, pt)
    rows = []
    for mid, what, fn, expect in SPEC:
        mi = copy.deepcopy(inp)
        fn(mi["json"])
        try:
            res = assembly.check(mi, cov, pt)
            failed, crashed = res["failed"], None
        except Exception as exc:  # a crash is a failure to certify, but record it as such
            failed, crashed = [], f"{type(exc).__name__}: {exc}"
        if mid.startswith("B"):
            verdict = ("passes (not a bookkeeping error; only the ingredient's own replay can see it)"
                       if not failed and not crashed else f"caught incidentally by {failed or crashed}")
            ok = True
        else:
            ok = crashed is None and bool(failed) and all(x in failed for x in expect)
            verdict = "fails as expected" if ok else ("CRASH" if crashed else "NOT CAUGHT as expected")
        rows.append({"id": mid, "what": what, "expected_failures": expect, "failed": failed, "crashed": crashed,
                     "verdict": verdict, "ok": ok})
        print(f"{mid:4s} {what:72s} {verdict}  ({len(failed)} failed{'; ' + crashed if crashed else ''})")
    out = {"tool": "mutations.py (S5 review)", "baseline_pass": base["pass"], "baseline_checks": len(base["checks"]),
           "rows": rows, "all_ok": all(r["ok"] for r in rows), "timing_s": round(time.time() - t0, 2)}
    OUT.write_text(json.dumps(out, indent=1) + "\n")
    print(f"baseline pass = {base['pass']} ({len(base['checks'])} checks); all mutations as expected: {out['all_ok']}")
    return 0 if out["all_ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
