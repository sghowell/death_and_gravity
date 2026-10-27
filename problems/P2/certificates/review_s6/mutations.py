#!/usr/bin/env python
"""S6 review, check C11: extra negative controls (mutations) for the author's assembly check p2.beta_assembly.

NB6a-NB6i are the author's controls.  M1-M36 probe items they do not isolate: the traced path, the
discs, pole clearance, the re-sum, the stored winding, a certificate in the wrong slot, the band box
top, part (ii') per l (h0, W7, a gap at Re w = a'_l, a missing row), Lemma E (Wt_l, B5, positivity),
part (i) (T, A, S, L, cover), beta (P2-K6's ball below the band, or moved inside it), the symmetry,
the KAT/control/replay records and P2-4's pins.  B1-B3 document what the check does not catch by
design (it trusts each ingredient's own replay; A2 asks for +dirty tags to be listed, not failed).
Each mutation acts on an in-memory deep copy of the loaded inputs; nothing on disk is changed.

Run (from this directory):  python mutations.py      -> mutations_out.json
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

from p2 import beta_assembly as BA  # noqa: E402  (the code under test)
from p2 import beta_part1  # noqa: E402
from p2.params import Point  # noqa: E402


def shift_coord(bd, old, new):
    for s in bd["pieces"]:
        for k in ("p0", "p1", "c"):
            s[k] = [new if v == old else v for v in s[k]]


def inc_add(s, dx):
    m = s["increment"].strip("[]").split("+/-")
    s["increment"] = f"[{float(m[0]) + dx!r} +/- {m[1].strip()}]"


def piece_near(bd, x, y):
    best = None
    for s in bd["pieces"]:
        cx, cy = (float(Fr(v)) for v in s["c"])
        d = (cx - x) ** 2 + (cy - y) ** 2
        if best is None or d < best[0]:
            best = (d, s)
    return best[1]


def km_kp(inp):
    q = inp["json"]["params"]["points"]["PB"]["quantities"]
    return float(Fr(q["kappa_minus"]["approx"])), float(Fr(q["kappa_plus"]["approx"]))


# ---- mutations: (id, description, function, expected failing checks)
def m1(inp):
    p = inp["json"]["band:0"]["pieces"]
    del p[len(p) // 2]


def m2(inp):
    inc_add(inp["json"]["band:12"]["pieces"][7], 2 * 3.141592653589793)


def m3(inp):
    inp["json"]["band:8"] = copy.deepcopy(inp["json"]["band:9"])


def m4(inp):
    """band l = 2: the top edge moved down to 0.56 kappa_- (a gap below S_delta), path kept consistent"""
    km, _ = km_kp(inp)
    bd = inp["json"]["band:2"]
    old = bd["region"]["b0"]
    new = Fr(int(0.56 * km * 2 ** 30), 2 ** 30)
    bd["region"]["b0"] = str(new)
    shift_coord(bd, str(-Fr(old)), str(-new))


def m5(inp):
    """band l = 0: the disc nearest omega0 inflated until its box reaches -i kappa_+ (endpoints still inside)"""
    km, kp = km_kp(inp)
    bd = inp["json"]["band:0"]
    s = piece_near(bd, 0.0, -float(Fr(bd["region"]["b_prime"])))
    cx, cy = (Fr(v) for v in s["c"])
    d = ((float(cx)) ** 2 + (float(cy) + kp) ** 2) ** 0.5
    s["rho"] = str(Fr(d * 0.8).limit_denominator(2 ** 40))


def m6(inp):
    inp["json"]["band:5"]["certified"] = False


def m7(inp):
    inp["json"]["band:0"]["winding"] = -1


def m8(inp):
    p = inp["json"]["band:3"]["pieces"]
    p[3], p[4] = p[4], p[3]


def m9(inp):
    s = inp["json"]["band:11"]["pieces"][10]
    s["rho"] = str(Fr(s["rho"]) / 4)


def m10(inp):
    for r in inp["json"]["A_fixed"]["rows"]:
        if r["l"] == 3:
            r["box"]["h0"] = "2/41"


def m11(inp):
    for r in inp["json"]["A_fixed"]["rows"]:
        if r["l"] == 0:
            r["W7"]["pass"] = False


def m12(inp):
    """l = 6: a'_6 lowered below zeta_A Omega_P/h_6 in both the A certificate and the band box (a gap)"""
    rows = {r["l"]: r for r in inp["json"]["A_fixed"]["rows"]}
    r = rows[6]
    old = r["a_prime"]
    new = str(Fr(old) - Fr(1, 2 ** 12))
    r["a_prime"] = new
    bd = inp["json"]["band:6"]
    bd["region"]["a_prime"] = new
    shift_coord(bd, old, new)


def m13(inp):
    inp["json"]["A_fixed"]["rows"] = [r for r in inp["json"]["A_fixed"]["rows"] if r["l"] != 11]


def m14(inp):
    for r in inp["json"]["epole"]["Wt_l_le_19"]:
        if r["l"] == 0:
            r["Wt"] = "[+/- 0.05] + [+/- 1e-30]j"


def m15(inp):
    inp["json"]["epole"]["B5"]["rows"][3]["pass"] = False


def m16(inp):
    inp["json"]["epole"]["V_l_positive_for_l_ge_1"]["r_plus_minus_Q2_over_M"] = "[-0.01 +/- 1e-12]"


def m17(inp):
    del inp["json"]["T"]["boxes"][4]


def m18(inp):
    inp["json"]["T"]["boxes"][2]["nu_low"] = "[-0.02127099741 +/- 2.58e-12]"     # -kappa_-/2 only


def m19(inp):
    inp["json"]["A"]["monotonicity_W7"]["pass"] = False


def m20(inp):
    inp["json"]["L"]["design"]["Im_range_over_km"] = ["-9/10", "-11/20"]


def m21(inp):
    inp["json"]["L"]["min_Re_y_plus_lower"] = -0.01


def m22(inp):
    inp["json"]["S_plus"]["missing"] = [["box", "3/50", "2/25"]]


def m23(inp):
    inp["json"]["S_c"]["depth_over_kappa_minus"] = "11/20"


def m24(inp):
    k = next(iter(inp["json"]["cover"]["inequalities"]))
    inp["json"]["cover"]["inequalities"][k] = False


def m25(inp):
    inp["json"]["krawczyk"]["enclosure_over_kappa_minus"]["im"] = "[-0.930000000000000000000000000000 +/- 3.16e-31]"


def m26(inp):
    inp["json"]["krawczyk"]["enclosure_over_kappa_minus"]["im"] = "[-0.850000000000000000000000000000 +/- 3.16e-31]"


def m27(inp):
    inp["json"]["zeq"]["Z4"] = False


def m28(inp):
    inp["json"]["kat3"]["boxes"][1]["as_expected"] = False


def m29(inp):
    inp["json"]["replay3"]["pass"] = False


def m30(inp):
    inp["json"]["p2_4"]["pass"] = False


def m31(inp):
    k = next(iter(inp["p2_4_inputs_unchanged"]))
    inp["p2_4_inputs_unchanged"][k] = False


def m32(inp):
    inp["json"]["controls1"]["NB4_T_S4_design_at_e"]["fails_as_required"] = False


def m33(inp):
    p = inp["json"]["band:19"]["pieces"]
    p[-1]["p1"] = [p[-1]["p1"][0], str(Fr(p[-1]["p1"][1]) + Fr(1, 2 ** 40))]


def m34(inp):
    inp["json"]["epole"]["Wt_l_le_19"] = [r for r in inp["json"]["epole"]["Wt_l_le_19"] if r["l"] != 19]


def m35(inp):
    del inp["json"]["S_c"]


def m36(inp):
    inp["json"]["kat3"]["NB2"]["fails_as_required"] = False


# ---- not caught by design
def b1(inp):
    inc_add(inp["json"]["band:4"]["pieces"][5], 0.2)                  # within the 1/4 tolerance


def b2(inp):
    inp["json"]["T"]["commit"] = inp["json"]["T"]["commit"] + "+dirty"   # A2: listed, not failed


def b3(inp):
    """a disc far from E_+ u E_c inflated by 5/4 (the stored increment is not re-verified by the assembly)"""
    bd = inp["json"]["band:7"]
    a = float(Fr(bd["region"]["a_prime"]))
    s = piece_near(bd, a, -0.7 * km_kp(inp)[0])
    s["rho"] = str(Fr(s["rho"]) * Fr(5, 4))


SPEC = [
    ("M1", "a piece removed from the l = 0 path", m1, ["C5.l0.path_traces_H"]),
    ("M2", "an l = 12 increment off by 2 pi", m2, ["C5.l12.resum"]),
    ("M3", "the l = 9 band certificate in l = 8's slot", m3, ["C5.l8.meta"]),
    ("M4", "l = 2 band top moved to -0.56 kappa_- (gap below S_delta)", m4, ["C5.l2.b0"]),
    ("M5", "l = 0 disc near omega0 inflated to reach the pole", m5, ["C5.l0.discs_clear_E"]),
    ("M6", "l = 5 band certificate flagged uncertified", m6, ["C5.l5.meta"]),
    ("M7", "l = 0 stored winding -1 (ignoring omega_*)", m7, ["C5.l0.winding_value", "C5.l0.resum"]),
    ("M8", "two l = 3 pieces swapped", m8, ["C5.l3.path_traces_H"]),
    ("M9", "an l = 11 disc shrunk below its segment", m9, ["C5.l11.path_in_discs"]),
    ("M10", "part (ii') l = 3 certified at h0 = 2/41, not h_3", m10, ["C4.l3.pass"]),
    ("M11", "part (ii') l = 0 W7 fails", m11, ["C4.l0.pass"]),
    ("M12", "a'_6 below zeta_A Omega_P/h_6 (gap at Re w = a'_6)", m12, ["C4.l6.pass"]),
    ("M13", "part (ii') row l = 11 missing", m13, ["C4.l11.pass", "C5.l11.a_identical"]),
    ("M14", "Wt_0 not separated from 0", m14, ["C3.l0.Wt_nonzero", "C5.l0.pole_forced_consistent_C3"]),
    ("M15", "a B5 row fails", m15, ["C3.B5_residue_formula"]),
    ("M16", "V_l positivity input negative", m16, ["C3.V_l_positive_l_ge_1"]),
    ("M17", "a T box removed (zeta gap)", m17, ["C6.T_design"]),
    ("M18", "a T box certified only to -kappa_-/2", m18, ["C6.T_design"]),
    ("M19", "part (i) A: W7 fails", m19, ["C6.A_design"]),
    ("M20", "L rows only to -0.9 kappa_-", m20, ["C6.L_band_rows"]),
    ("M21", "L margin min Re y_+ negative", m21, ["C6.L_W6_prime_margin"]),
    ("M22", "an S+ job missing", m22, ["C6.S_design"]),
    ("M23", "S_c run at depth 11/20", m23, ["C6.S_c.pass"]),
    ("M24", "a cover inequality false", m24, ["C6.cover_pass", "C6.cover_replay_identical"]),
    ("M25", "omega_* moved below the band (-0.93 kappa_-)", m25, ["C8.omega_star_in_B0"]),
    ("M26", "omega_* moved inside the band (-0.85 kappa_-)", m26, ["C8.kat6_consistent"]),
    ("M27", "symmetry Z4 false", m27, ["C6.symmetry_Z4"]),
    ("M28", "KAT B3 not as expected", m28, ["C9.KAT_B2_B3_B4_B6"]),
    ("M29", "band 5% replay failed", m29, ["C9.band_replay"]),
    ("M30", "P2-4's assembly not passing", m30, ["C2.P2_4_pass"]),
    ("M31", "a P2-4 input changed", m31, ["C2.P2_4_inputs_unchanged"]),
    ("M32", "NB4 not failing", m32, ["C9.NB3_NB4_NB5"]),
    ("M33", "l = 19 path ends 2^-40 short of (0, -b0)", m33, ["C5.l19.path_traces_H"]),
    ("M34", "Lemma E row l = 19 missing", m34, ["C3.l19.Wt_nonzero"]),
    ("M35", "the S_c certificate missing", m35, ["C6.S_c.pass"]),
    ("M36", "NB2 certified (a winding through omega_*)", m36, ["C9.NB2"]),
]
BYDESIGN = [
    ("B1", "an l = 4 increment changed by 0.2 (< 1/4 of pi)", b1),
    ("B2", "a +dirty commit tag on T_PB.json", b2),
    ("B3", "an l = 7 disc far from E inflated by 5/4", b3),
]


def main():
    t0 = time.time()
    inp = BA.load_inputs()
    pt = Point("PB", BA.PREC)
    cov = beta_part1.cover()
    base = BA.check(inp, cov, pt)
    out = {"baseline_pass": base["pass"], "baseline_failed": base["failed"], "baseline_checks": len(base["checks"]),
           "mutations": [], "by_design": []}
    for mid, what, fn, expect in SPEC:
        mi = copy.deepcopy(inp)
        fn(mi)
        res = BA.check(mi, cov, pt)
        caught = (not res["pass"]) and all(x in res["failed"] for x in expect)
        out["mutations"].append({"id": mid, "what": what, "expected": expect, "failed": res["failed"][:12],
                                 "n_failed": len(res["failed"]), "caught_on_named_items": caught})
        print(mid, caught, res["failed"][:6], flush=True)
    for bid, what, fn in BYDESIGN:
        mi = copy.deepcopy(inp)
        fn(mi)
        res = BA.check(mi, cov, pt)
        out["by_design"].append({"id": bid, "what": what, "pass": res["pass"], "failed": res["failed"][:12],
                                 "dirty_tags_listed": res["summary"].get("dirty_tags")})
        print(bid, "passes" if res["pass"] else "fails", res["failed"][:6], res["summary"].get("dirty_tags"), flush=True)
    out["all_caught"] = all(m["caught_on_named_items"] for m in out["mutations"])
    out["wall_s"] = round(time.time() - t0, 1)
    (HERE / "mutations_out.json").write_text(json.dumps(out, indent=1) + "\n")
    print("all caught:", out["all_caught"], out["wall_s"])


if __name__ == "__main__":
    main()
