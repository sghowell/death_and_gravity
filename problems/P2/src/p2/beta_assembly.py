"""P2-6 at P_B (S6): a machine check of the assembly bookkeeping for the two-sided beta.

The written proof is notes/s6-proof.md; the plan is notes/s6-preregistration.md (section 3, Amendments
A2, A4, A5).  This module reads only committed certificates.  It pins each by sha256 and checks the
S2-S5 ones against their blobs on main 300c22a7 (A2).  With exact rationals and 256-bit Arb balls it
checks that the certificates cover S_e = {Im w >= -e kappa_-}, e = 23/25, for every degree l with no
gap, and that beta = -Im omega_*/kappa_- is enclosed by P2-K6's Krawczyk ball:

  C1  the pins; S_e meets E_+ u E_c only in omega0 = -i kappa_+; the band bounds b' >= e kappa_-, b0 <= kappa_-/2;
  C2  S_delta: P2-4's assembly certificate passes, and every input it pinned is unchanged;
  C3  Lemma E: the first Frobenius step, Wt_l(omega0) != 0 for l <= 19, the bound V_l > 0 (l >= 1), KAT B5;
  C4  part (ii'): A at fixed h_l passes for every l <= 19 (W7, Jost decay, depth e), a'_l >= zeta_A Omega_P/h_l;
  C5  part (iii): for every l <= 19 the band certificate's a'_l is C4's; the stored discs trace the right half
      of dB_l exactly and clear E_+ u E_c; the re-summed winding is the stored one, 0 (l = 0) or -1 (l >= 1);
  C6  part (i) at depth e (l >= 20): T, A, S (both sides), L and the cover pass at depth e and h0 = 2/41, the
      cover replays identically, the cover's regimes and the band rows;
  C7  the l-ranges: l0 = ceil(1/h0 - 1/2) = 20 follows l <= 19;
  C8  beta: P2-K6's ball (reused unchanged, A4), omega_* in B_0 and off E_+ u E_c, below S_delta, the width;
  C9  the KATs B1-B7 and the controls NB1-NB5 as recorded, and the 5% band replay;
  C10 the ingredient levels and their minimum; the +dirty commit tags of the new certificates (A2).

The windings count zeros minus poles: W_l is meromorphic on B_l with at most one pole, at omega0, of
order <= 1 (Lemma E, step 2).  l = 0: winding 0 and omega_* in B_0 give Z = P = 1; l >= 1: winding -1
gives P = 1 and Z = 0.  Negative controls NB6a-NB6i mutate a copy of the inputs and must each fail on
their named items.

    python -m p2.beta_assembly            # run, write certificates/beta/assembly_beta_PB.json
    python -m p2.beta_assembly verify     # recompute and compare with the committed file
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import subprocess
import sys
import time
from fractions import Fraction
from pathlib import Path

from flint import arb

from .arbutil import A, git_commit, lohi_out, precision
from .assembly import Rec, _discs_ok, _F, _path_ok, _poles_ok
from .beta_common import CERT, EDGE, H0, L_SMALL, MAIN_PIN, ROOT, band_bounds
from .params import Point

PREC = 256
OUT = CERT / "assembly_beta_PB.json"
NAME = "PB"
LEVEL_ORDER = ["HEURISTIC", "CONJECTURED", "VERIFIED_N", "WRITTEN_PROOF", "CERTIFIED", "FORMALIZED"]
# The ingredients of P2-6 and their levels (notes/s6-preregistration.md section 6; the ledger on main
# 300c22a7 for P2-4, P2-D1' and P2-K6).
INGREDIENTS = {
    "P2-4 (S_delta, every l)": "WRITTEN_PROOF",
    "P2-D1' (set clause, any depth)": "WRITTEN_PROOF",
    "part (i) at depth e, l >= 20 (W1, W2, W5, W6')": "WRITTEN_PROOF",
    "part (ii'), A at fixed h, l <= 19 (W1, W2, W4)": "WRITTEN_PROOF",
    "Lemma E, l >= 1 (positivity; l >= 20 also by W6')": "WRITTEN_PROOF",
    "Lemma E, l = 0 and l <= 19 (Arb)": "CERTIFIED",
    "part (iii) band counts, l <= 19": "CERTIFIED",
    "P2-K6 Krawczyk ball of omega_*": "CERTIFIED",
}
BETA_WIDTH_MAX = Fraction(1, 10 ** 8)                # FORMULATION P2-6: beta_hi - beta_lo <= 1e-8
P24_BETA_UPPER = Fraction(90938098, 10 ** 8)         # P2-4/P2-K6 in the ledger: beta <= 0.90938098
D_DELTA = Fraction(11, 20)                           # S_delta's depth over kappa_- (P2-4)


def _rel(p: Path) -> str:
    return str(p.relative_to(ROOT))


def s5_paths() -> dict:
    """The S2-S5 certificates S6 reads; each must equal its blob on main 300c22a7 (A2)."""
    c = ROOT / "certificates"
    return {"p2_4": c / "assembly" / "assembly_PB.json",            # S_delta (P2-4)
            "params": c / "params.json",                             # the point
            "krawczyk": c / "part3" / "krawczyk_PB_l0.json",         # omega_* (P2-K6, A4)
            "kat6": c / "kat6_winding_krawczyk.json",                # omega_* (P2-K6, S1 code)
            "zeq": c / "part3" / "zeq_identity.json",                # Z4: W(-conj w) = conj W(w)
            "bulk": c / "part1" / "bulk_PB.json"}                    # W_bound, read by S and L


def new_paths() -> dict:
    """The S6 certificates."""
    b = CERT
    p = {"epole": b / "epole_PB.json", "A_fixed": b / "part2" / "A_fixed_PB.json",
         "T": b / "part1" / "T_PB.json", "A": b / "part1" / "A_PB.json",
         "S_plus": b / "part1" / "S_PB_plus.json", "S_c": b / "part1" / "S_PB_c.json",
         "L": b / "part1" / "L_PB.json", "cover": b / "part1" / "cover_PB.json",
         "controls1": b / "part1" / "controls_PB.json", "B7": b / "part1" / "B7_PB.json",
         "kat3": b / "part3" / "kat_controls_PB.json", "replay3": b / "part3" / "replay_PB.json"}
    for l in L_SMALL:
        p[f"band:{l}"] = b / "part3" / f"band_PB_l{l}.json"
    return p


def _blob_sha(path: Path):
    try:
        out = subprocess.run(["git", "show", f"{MAIN_PIN}:problems/P2/{_rel(path)}"], capture_output=True,
                             cwd=ROOT, check=True)
        return hashlib.sha256(out.stdout).hexdigest()
    except (OSError, subprocess.CalledProcessError):
        return None


def load_inputs() -> dict:
    """{'json', 'sha256', 'main_blob_match', 'p2_4_inputs_unchanged'}; missing files are simply absent."""
    js, sha = {}, {}
    for role, p in {**s5_paths(), **new_paths()}.items():
        if p.exists():
            raw = p.read_bytes()
            js[role] = json.loads(raw)
            sha[_rel(p)] = hashlib.sha256(raw).hexdigest()
    main_ok = {_rel(p): (_rel(p) in sha and sha[_rel(p)] == _blob_sha(p)) for p in s5_paths().values()}
    pinned = (js.get("p2_4") or {}).get("inputs_sha256", {})
    p24 = {}
    for k, v in pinned.items():
        f = ROOT / k
        p24[k] = f.exists() and hashlib.sha256(f.read_bytes()).hexdigest() == v
    return {"json": js, "sha256": sha, "main_blob_match": main_ok, "p2_4_inputs_unchanged": p24}


def _m(x: arb) -> str:
    return lohi_out(x, 10)


def _re_ball(s: str) -> arb:
    """The real part of an acb string '[re +/- r] + [im +/- r]j'."""
    return arb(s.split("] + [")[0] + "]")


def _dirty_tags(obj, path="") -> list:
    """Every commit tag ending in '+dirty' under the keys commit, commits, computed_at, assembled_at."""
    out = []
    if isinstance(obj, dict):
        for k, v in obj.items():
            if k in ("commit", "commits", "computed_at", "assembled_at"):
                vals = v if isinstance(v, list) else [v]
                out += [f"{path}/{k}: {x}" for x in vals if isinstance(x, str) and x.endswith("+dirty")]
            elif isinstance(v, (dict, list)):
                out += _dirty_tags(v, f"{path}/{k}")
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            out += _dirty_tags(v, f"{path}[{i}]")
    return out


# ---------------------------------------------------------------------------------------------
# the check

def check(inp: dict, cover_replay: dict | None, pt: Point) -> dict:
    """Run C1-C10.  Returns {'checks', 'summary', 'pass', 'failed'}."""
    J = inp["json"]
    R = Rec()
    S: dict = {}
    with precision(PREC):
        km, kp, kc = pt.km, pt.kp, pt.kc
        e = A(EDGE) * km
        bpq, b0q = band_bounds(pt)
        bp, b0 = Fraction(int(bpq.p), int(bpq.q)), Fraction(int(b0q.p), int(b0q.q))

        # ---------------- C1: pins, the E-points, the band bounds
        mb = inp["main_blob_match"]
        R.add("C1.pins_main_blobs", bool(mb) and all(mb.values()),
              note=f"{sum(mb.values())}/{len(mb)} S2-S5 certificates equal their blobs on main {MAIN_PIN}")
        params = (J.get("params") or {}).get("points", {}).get(NAME, {})
        try:
            q = params["quantities"]
            ok = all(getattr(pt, a).overlaps(arb(q[k]["ball"])) for a, k in
                     (("km", "kappa_minus"), ("kp", "kappa_plus"), ("kc", "kappa_c"), ("OmegaP", "Omega_P")))
            R.add("C1.point_overlaps_params", ok and params.get("qfrac") == "999/1000")
        except (KeyError, TypeError):
            R.add("C1.point_overlaps_params", False)
        R.add("C1.only_omega0_in_S_e", bool(e < 2 * kp) and bool(e < kc),
              margin=_m((2 * kp - e) / km) + " kappa_- (2 kappa_+ - e kappa_-), " + _m((kc - e) / km) + " (kappa_c - e kappa_-)")
        R.add("C1.omega0_in_S_e", bool(kp < e), margin=_m((e - kp) / km) + " kappa_-")
        R.add("C1.b_prime_ge_e", bool(A(bp) >= e) and bool(A(bp) < 2 * kp) and bool(A(bp) < kc),
              margin=_m(A(bp) - e))
        R.add("C1.b0_le_half", bool(A(b0) <= km / 2) and b0 > 0 and bool(A(b0) < kp), margin=_m(km / 2 - A(b0)))
        R.add("C1.band_overlaps_S_delta", bool(A(b0) < A(D_DELTA) * km))
        S["edge_over_kappa_minus"] = str(EDGE)
        S["band"] = {"b_prime": str(bp), "b0": str(b0), "b_prime_over_kappa_minus": _m(A(bp) / km),
                     "b0_over_kappa_minus": _m(A(b0) / km), "omega0_over_kappa_minus": "-i " + _m(kp / km)}

        # ---------------- C2: S_delta (P2-4)
        p24 = J.get("p2_4") or {}
        R.add("C2.P2_4_pass", p24.get("pass") is True and not p24.get("result", {}).get("failed")
              and all(c.get("fails_as_expected") for c in p24.get("negative_controls", [])),
              note=f"{len(p24.get('result', {}).get('checks', {}))} checks")
        pin = inp["p2_4_inputs_unchanged"]
        R.add("C2.P2_4_inputs_unchanged", bool(pin) and all(pin.values()), note=f"{sum(pin.values())}/{len(pin)}")
        R.add("C2.P2_4_depth", (p24.get("point") or {}).get("delta") == "1/20" and D_DELTA == Fraction(1, 2) + Fraction(1, 20))

        # ---------------- C3: Lemma E
        ep = J.get("epole")
        wt_ok = {}
        if R.add("C3.present", ep is not None):
            R.add("C3.sympy_first_step", ep["sympy_first_step"]["pass"] is True
                  and ep["sympy_first_step"]["D1_at_omega0"] == "0"
                  and ep["sympy_first_step"]["equals_D1_a1_plus_C0"] is True)
            v = ep["V_l_positive_for_l_ge_1"]
            R.add("C3.V_l_positive_l_ge_1", v["pass"] is True and bool(arb(v["r_plus_minus_Q2_over_M"]) > 0)
                  and bool(arb(v["two_minus_2Lambda_rc2_over_3"]) > 0) and bool(arb(v["kappa_c_minus_kappa_plus"]) > 0))
            rows = {r["l"]: r for r in ep.get("Wt_l_le_19", [])}
            for l in L_SMALL:
                r = rows.get(l)
                ok = r is not None and r.get("Wt_nonzero") is True and not _re_ball(r["Wt"]).contains(0)
                wt_ok[l] = ok
                R.add(f"C3.l{l}.Wt_nonzero", ok, note=(r["Wt"].split(" + ")[0] if r else None))
            R.add("C3.B5_residue_formula", ep["B5"]["pass"] is True and len(ep["B5"]["rows"]) == 8
                  and all(x["pass"] for x in ep["B5"]["rows"]))
            R.add("C3.normalisation_stated", "a_0 = 1" in ep.get("normalisation", "")
                  and "t = r - r_+" in ep.get("normalisation", ""))
            S["lemma_E"] = {"Wt_l": {l: rows[l]["Wt"].split(" + ")[0] for l in L_SMALL if l in rows},
                            "samples_l_ge_20": [r["l"] for r in ep.get("Wt_samples_l_ge_20", []) if r.get("Wt_nonzero")]}

        # ---------------- C4: part (ii'), A at fixed h_l
        af = J.get("A_fixed")
        arows = {}
        if R.add("C4.present", af is not None):
            arows = {r["l"]: r for r in af.get("rows", [])}
            R.add("C4.depth_and_design", af.get("depth_over_kappa_minus") == str(EDGE)
                  and (af.get("design") or {}).get("phi0") == "-1/5")
            for l in L_SMALL:
                r = arows.get(l)
                h = Fraction(2, 2 * l + 1)
                ok = r is not None and r.get("pass") is True and r["W7"]["pass"] is True
                if ok:
                    bx = r["box"]
                    ok = (bx["pass"] is True and bx["jost_decay"] is True and bx["n_fails"] == 0
                          and bx["h0"] == str(h) and bool(arb(bx["mu_cert"]) > 0)
                          and [_F(z) for z in bx["zeta_over_OmegaP"]] == [_F(r["zeta_A_over_OmegaP"])] * 2
                          and arb(bx["nu_low"]).contains(-e)
                          and bool(A(_F(r["a_prime"])) >= A(_F(r["zeta_A_over_OmegaP"])) * pt.OmegaP / A(h)))
                R.add(f"C4.l{l}.pass", ok)

        # ---------------- C5: part (iii), the band counts
        per_l = {}
        for l in L_SMALL:
            bd = J.get(f"band:{l}")
            if not R.add(f"C5.l{l}.present", bd is not None):
                continue
            reg = bd["region"]
            a, bpl, b0l = _F(reg["a_prime"]), _F(reg["b_prime"]), _F(reg["b0"])
            R.add(f"C5.l{l}.meta", bd.get("point") == NAME and bd.get("l") == l and bd.get("certified") is True)
            R.add(f"C5.l{l}.a_identical", l in arows and a == _F(arows[l]["a_prime"]),
                  note=f"a'_l = {a}")
            R.add(f"C5.l{l}.b_prime_ge_e", bool(A(bpl) >= e) and bpl == bp, margin=_m(A(bpl) - e))
            R.add(f"C5.l{l}.b0", b0l == b0)
            R.add(f"C5.l{l}.omega0_interior", bool(A(bpl) > kp) and bool(kp > A(b0l)) and a > 0)
            pcs = bd.get("pieces", [])
            ok, msg = _path_ok(pcs, a, bpl, -b0l)
            R.add(f"C5.l{l}.path_traces_H", ok, note=msg)
            ok, msg = _discs_ok(pcs)
            R.add(f"C5.l{l}.path_in_discs", ok, note=msg)
            ok, msg = _poles_ok(pcs, kp, kc, bpl)
            R.add(f"C5.l{l}.discs_clear_E", ok, note=msg)
            tot = sum((arb(s["increment"]) for s in pcs), arb(0))
            wn = tot / arb.pi()
            k = bd.get("winding")
            R.add(f"C5.l{l}.resum", isinstance(k, int) and bool(abs(wn - k) < arb("0.25")), margin=wn.str(8))
            expected = 0 if l == 0 else -1
            R.add(f"C5.l{l}.winding_value", k == expected, note=f"stored {k}, expected {expected}")
            # P = 1 is forced (l = 0 with omega_* inside: Z >= 1, Z = P <= 1; l >= 1: P = Z + 1 >= 1),
            # so the residue at omega0 is nonzero: Lemma E's Wt_l != 0 again (C3)
            R.add(f"C5.l{l}.pole_forced_consistent_C3", k == expected and wt_ok.get(l, False))
            per_l[l] = {"a_prime_over_kappa_minus": _m(A(a) / km), "winding": k, "zeros": 1 if l == 0 else 0,
                        "poles": 1, "pieces": len(pcs), "resum_over_pi": wn.str(6),
                        "min_abs_W_lower": bd.get("min_abs_W_on_H_lower")}
        R.add("C5.all_l_present", all(f"band:{l}" in J for l in L_SMALL))
        S["band_counts"] = {"per_l": per_l, "pieces_total": sum(r["pieces"] for r in per_l.values())}

        # ---------------- C6: part (i) at depth e, l >= 20
        cv = J.get("cover") or {}
        names = ("T", "A", "S_plus", "S_c", "L")
        for n in names:
            c = J.get(n) or {}
            R.add(f"C6.{n}.pass", c.get("pass") is True and c.get("point") == NAME
                  and c.get("depth_over_kappa_minus") == str(EDGE) and c.get("h0") == cv.get("h0"))
        R.add("C6.cover_pass", cv.get("pass") is True and all(cv.get("inequalities", {}).values())
              and all(cv.get("certificates_pass", {}).values()) and all(cv.get("tiling", {}).values())
              and len(cv.get("tiling", {})) == 6 and len(cv.get("inequalities", {})) == 11)
        if cover_replay is not None:
            same = all(cover_replay.get(k) == cv.get(k) for k in
                       ("constants", "inequalities", "certificates_pass", "tiling", "cited", "pass", "h0"))
            R.add("C6.cover_replay_identical", same)
        try:
            T, Ac, L = J["T"], J["A"], J["L"]
            C = {k: _F(v) for k, v in cv["constants"].items()}
            R.add("C6.T_design", [_F(x) for x in T["design"]["zeta_over_OmegaP"]] == [C["zeta_TL"], C["zeta_TU"]]
                  and all(arb(bx["nu_low"]).contains(-e) and bx["pass"] for bx in T["boxes"]) and len(T["boxes"]) == 10)
            R.add("C6.A_design", _F(Ac["boxes"][0]["zeta_over_OmegaP"][0]) == C["zeta_A"]
                  and arb(Ac["boxes"][0]["nu_low"]).contains(-e) and Ac["monotonicity_W7"]["pass"] is True)
            R.add("C6.S_design", all(_F(x["constants"]["Omega_S_over_km"]) == C["Omega_S_over_km"]
                                     and _F(x["constants"]["zeta_SU_over_OmegaP"]) == C["zeta_SU"]
                                     and x.get("tail_region_in_strip") is True and not x.get("missing")
                                     for x in (J["S_plus"], J["S_c"])))
            d = L["design"]
            R.add("C6.L_band_rows", [_F(v) for v in d["Im_range_over_km"]] == [Fraction(-19, 20), Fraction(-11, 20)]
                  and [_F(v) for v in d["Re_range_over_km"]] == [Fraction(0), C["Omega_L_over_km"]]
                  and [_F(v) for v in d["h2_range"]] == [Fraction(0), _F(cv["h0"]) ** 2]
                  and L.get("tiles_ok") is True and not L.get("missing") and L.get("grid_boxes") == 1280
                  and "W6'" in d.get("series", ""),
                  note="Im w in [-19/20, -11/20] kappa_-: reaches below -e kappa_- and into S_delta")
            R.add("C6.L_rows_reach", Fraction(-19, 20) <= -EDGE and Fraction(-11, 20) >= -D_DELTA)
            R.add("C6.L_W6_prime_margin", bool(arb(L["min_Re_y_plus_lower"]) > 0) and bool(arb(L["max_Re_y_c_upper"]) < 0),
                  note=f"min Re y_+ >= {L['min_Re_y_plus_lower']:.4f}, max Re y_c <= {L['max_Re_y_c_upper']:.4f}")
            top_S = A(C["zeta_SU"]) * pt.OmegaP / A(_F(cv["h0"]))
            R.add("C6.S_reaches_above_L", bool(top_S > A(C["Omega_L_over_km"]) * km),
                  margin=_m(top_S / km - A(C["Omega_L_over_km"])) + " kappa_-")
            R.add("C6.regime_overlaps", C["Omega_S_over_km"] < C["Omega_L_over_km"] and C["zeta_TL"] < C["zeta_SU"]
                  and C["zeta_A"] < C["zeta_TU"])
            R.add("C6.Re_plus_Im_positive", bool(A(C["Omega_S_over_km"]) * km > e)
                  and bool(A(C["zeta_TL"]) * pt.OmegaP / A(_F(cv["h0"])) > e))
            S["part_i"] = {"h0": cv["h0"], "constants": cv["constants"],
                           "Re_cover": "[0, 8 km] (L, band rows) u [6 km, (9/10) OmP/h] (S) u [(4/5) OmP/h, "
                                       "(13/10) OmP/h] (T) u [(6/5) OmP/h, oo) (A) = [0, oo) for every h <= h0",
                           "T_mu_over_a_star": [bx["mu_over_a_star"] for bx in T["boxes"]],
                           "A_mu_over_a_star": Ac["boxes"][0]["mu_over_a_star"],
                           "L_margins": [L["min_Re_y_plus_lower"], L["max_Re_y_c_upper"]]}
        except (KeyError, IndexError, TypeError):
            R.add("C6.readable", False)
        zq = J.get("zeq") or {}
        R.add("C6.symmetry_Z4", zq.get("Z4") is True and zq.get("pass") is True)

        # ---------------- C7: the l-ranges
        try:
            h0 = _F(cv["h0"])
            l0 = math.ceil(1 / h0 - Fraction(1, 2))
            R.add("C7.no_l_gap", l0 == max(L_SMALL) + 1 and h0 == H0, note=f"l0 = ceil(1/h0 - 1/2) = {l0}; l <= 19 counted")
            S["l0"] = l0
        except (KeyError, TypeError, ZeroDivisionError):
            R.add("C7.no_l_gap", False)

        # ---------------- C8: beta
        kz = J.get("krawczyk") or {}
        try:
            im = arb(kz["enclosure_over_kappa_minus"]["im"])
            re_ = arb(kz["enclosure_over_kappa_minus"]["re"])
            R.add("C8.krawczyk_certified", kz.get("certified_unique_zero") is True and kz.get("point") == NAME
                  and kz.get("l") == 0)
            k6 = [c for c in (J.get("kat6") or {}).get("cases", []) if c.get("point") == NAME and c.get("l") == 0]
            R.add("C8.kat6_consistent", len(k6) == 1 and k6[0]["krawczyk"]["certified_unique_zero"] is True
                  and arb(k6[0]["krawczyk"]["enclosure_over_kappa_minus"]["im"]).overlaps(im))
            beta = -im
            width = 2 * arb(beta.rad())
            R.add("C8.beta_width", bool(width < A(BETA_WIDTH_MAX)), margin=width.str(3))
            b0band = J.get("band:0")
            if b0band is not None:
                a0 = A(_F(b0band["region"]["a_prime"]))
                ins = (bool(beta * km < A(_F(b0band["region"]["b_prime"]))) and bool(beta * km > A(_F(b0band["region"]["b0"])))
                       and bool(abs(re_) * km < a0))
                R.add("C8.omega_star_in_B0", ins, margin=_m(A(_F(b0band["region"]["b_prime"])) / km - beta) + " kappa_- below")
            else:
                R.add("C8.omega_star_in_B0", False, note="no l = 0 band certificate")
            R.add("C8.omega_star_off_E", bool(beta * km > kp) and bool(beta * km < 2 * kp) and bool(beta * km < kc))
            R.add("C8.omega_star_below_S_delta", bool(beta > A(D_DELTA)))
            R.add("C8.consistent_with_P2_4", bool(beta <= A(P24_BETA_UPPER)), margin=_m(A(P24_BETA_UPPER) - beta))
            S["beta"] = {"enclosure": beta.str(30, radius=True), "width": width.str(3),
                         "source": "certificates/part3/krawczyk_PB_l0.json (P2-K6, S3 re-derivation; unchanged, A4)",
                         "omega_star_over_kappa_minus": f"{re_.str(5)} + i {im.str(30)}"}
        except (KeyError, TypeError):
            R.add("C8.krawczyk_readable", False)

        # ---------------- C9: KATs and controls
        k3 = J.get("kat3") or {}
        boxes = k3.get("boxes", [])
        R.add("C9.KAT_B2_B3_B4_B6", len(boxes) == 6 and all(b.get("as_expected") for b in boxes),
              note=", ".join(f"{b.get('id')}(l={b.get('l', '?')}): {b.get('winding')}" for b in boxes))
        R.add("C9.NB1", (k3.get("NB1") or {}).get("as_expected") is True)
        R.add("C9.NB2", (k3.get("NB2") or {}).get("fails_as_required") is True)
        lh = k3.get("left_halves", [])
        R.add("C9.left_halves_symmetry", len(lh) == 2 and all(x.get("overlap") for x in lh))
        R.add("C9.B1_krawczyk_replay", (k3.get("B1") or {}).get("pass") is True)
        c1 = J.get("controls1") or {}
        R.add("C9.NB3_NB4_NB5", c1.get("pass") is True and all(
            (c1.get(k) or {}).get("fails_as_required") is True
            for k in ("NB3_L_pole_box_S4_series", "NB4_T_S4_design_at_e", "NB5_A_h0_2_zeta_5")))
        rp = J.get("replay3") or {}
        R.add("C9.band_replay", rp.get("pass") is True and rp.get("recomputed", 0) >= 3 * len(L_SMALL))
        b7 = J.get("B7") or {}
        R.add("C9.B7_float", b7.get("pass") is True, note="float evidence; not used by the proof")

        # ---------------- C10: levels and +dirty tags
        lmin = min(INGREDIENTS.values(), key=LEVEL_ORDER.index)
        R.add("C10.level_minimum", lmin == "WRITTEN_PROOF", note=lmin)
        dirty = []
        for role in new_paths():
            if role in J:
                dirty += [f"{role}{t}" for t in _dirty_tags(J[role])]
        S["dirty_tags"] = dirty
        S["levels"] = {"ingredients": INGREDIENTS, "minimum": lmin}

    failed = R.failed()
    return {"checks": R.items, "summary": S, "pass": not failed, "failed": failed}


# ---------------------------------------------------------------------------------------------
# negative controls (NB6 of the pre-registration, with A5, and more); each must fail on its named items

def _replace_coord(bd: dict, old: str, new: str) -> None:
    for s in bd["pieces"]:
        for key in ("p0", "p1", "c"):
            s[key] = [new if v == old else v for v in s[key]]


def _nb6a(inp):
    del inp["json"]["band:7"]
    return "band_PB_l7.json removed"


def _nb6b(inp):
    bd = inp["json"]["band:7"]
    old = bd["region"]["a_prime"]
    new = str(Fraction(old) - Fraction(1, 2 ** 24))
    bd["region"]["a_prime"] = new
    _replace_coord(bd, old, new)
    return f"l = 7: a' {old} -> {new}"


def _nb6c(inp):
    inp["json"]["band:3"]["winding"] = 0
    return "l = 3: stored winding -1 -> 0"


def _nb6d(inp, k):
    new = Fraction(k, 2 ** 30)
    for l in L_SMALL:
        bd = inp["json"][f"band:{l}"]
        old = bd["region"]["b_prime"]
        bd["region"]["b_prime"] = str(new)
        _replace_coord(bd, str(-Fraction(old)), str(-new))
    return f"every l: b' -> {new} (the dyadic just below e kappa_-; A5)"


def _nb6e(inp):
    for r in inp["json"]["epole"]["Wt_l_le_19"]:
        if r["l"] == 5:
            r["Wt"] = "[+/- 1.00e-10] + [+/- 1.00e-30]j"
    return "Lemma E: Wt_5 replaced by a ball containing 0"


def _nb6f(inp):
    for r in inp["json"]["A_fixed"]["rows"]:
        if r["l"] == 9:
            r["box"]["nu_low"] = "[-0.0300000 +/- 1e-8]"
    return "part (ii'), l = 9: the A box certified only to depth 0.705 kappa_-"


def _nb6g(inp):
    inp["json"]["cover"]["h0"] = "2/43"
    return "part (i): h0 2/41 -> 2/43 (l0 = 21)"


def _nb6h(inp):
    k = next(iter(inp["main_blob_match"]))
    inp["main_blob_match"][k] = False
    return f"{k} differs from its blob on main {MAIN_PIN}"


def _nb6i(inp):
    kz = inp["json"]["krawczyk"]["enclosure_over_kappa_minus"]
    kz["im"] = "[-0.9093809 +/- 1e-7]"
    return "P2-K6's ball widened to 2e-7"


def controls_spec(k_below: int) -> list:
    return [
        ("NB6a", "one band certificate removed", _nb6a, ["C5.l7.present", "C5.all_l_present"]),
        ("NB6b", "a'_7 lowered by 2^-24", _nb6b, ["C5.l7.a_identical"]),
        ("NB6c", "a winding set to 0 at l = 3", _nb6c, ["C5.l3.winding_value", "C5.l3.resum"]),
        ("NB6d", "b' just below e kappa_- (A5)", lambda i: _nb6d(i, k_below), ["C5.l0.b_prime_ge_e", "C5.l19.b_prime_ge_e"]),
        ("NB6e", "Lemma E: Wt_5 not certified nonzero", _nb6e, ["C3.l5.Wt_nonzero", "C5.l5.pole_forced_consistent_C3"]),
        ("NB6f", "part (ii') short of the edge at l = 9", _nb6f, ["C4.l9.pass"]),
        ("NB6g", "an l gap (h0 = 2/43)", _nb6g, ["C7.no_l_gap"]),
        ("NB6h", "an S5 certificate off its main blob", _nb6h, ["C1.pins_main_blobs"]),
        ("NB6i", "beta's enclosure too wide", _nb6i, ["C8.beta_width"]),
    ]


def run_controls(inp: dict, cover_replay: dict | None, pt: Point) -> list:
    with precision(PREC):
        e = A(EDGE) * pt.km
        k = math.floor(float((e * 2 ** 30).lower()))
        while not bool(A(Fraction(k, 2 ** 30)) < e):
            k -= 1
    out = []
    for cid, what, mutate, expect in controls_spec(k):
        mi = copy.deepcopy(inp)
        how = mutate(mi)
        res = check(mi, cover_replay, pt)
        failed = res["failed"]
        out.append({"id": cid, "what": what, "mutation": how, "expected_failures": expect, "failed": failed,
                    "fails_as_expected": (not res["pass"]) and all(x in failed for x in expect)})
    return out


# ---------------------------------------------------------------------------------------------
# driver

def assemble(write: bool = True) -> dict:
    from . import beta_part1
    t0, c0 = time.time(), time.process_time()
    inp = load_inputs()
    pt = Point(NAME, PREC)
    cov = beta_part1.cover()
    res = check(inp, cov, pt)
    ctr = run_controls(inp, cov, pt)
    ok = res["pass"] and all(c["fails_as_expected"] for c in ctr)
    out = {
        "claim": "P2-6 at P_B (assembly bookkeeping): the committed certificates cover S_e = {Im w >= -e kappa_-}, "
                 "e = 23/25, for every l with no gap; in S_e the only zeros of the W_l are 0 and omega_* (both l = 0, "
                 "simple), and the only point of E_+ u E_c is omega0 = -i kappa_+, where every W_l has a simple pole "
                 "(Lemma E); so beta = -Im omega_*/kappa_-, enclosed by P2-K6's Krawczyk ball",
        "point": {"name": NAME, "qfrac": "999/1000", "delta": "1/20", "Lambda": "3/50", "M": "1", "edge": str(EDGE)},
        "decomposition": [
            "S_delta = {Im w >= -(11/20) kappa_-}, every l: P2-4 (certificates/assembly/assembly_PB.json)",
            "Re w < 0: W_l(-conj w) = conj W_l(w) (zeq_identity Z4)",
            "band B = {-e kappa_- <= Im w <= -kappa_-/2}, Re w >= 0, l >= 20: part (i) at depth e (L, S, T, A; cover_PB.json)",
            "B, l <= 19: B_l = [-a'_l, a'_l] x [-b', -b0] (band counts) and Re w >= a'_l (A at fixed h_l)",
            "omega0 = -i kappa_+, every l: Lemma E (epole_PB.json and notes/s6-proof.md section 2)",
        ],
        "result": res,
        "negative_controls": ctr,
        "inputs_sha256": inp["sha256"],
        "main_blob_match": inp["main_blob_match"],
        "p2_4_inputs_unchanged": inp["p2_4_inputs_unchanged"],
        "pass": ok,
        "prec_bits": PREC,
        "commit": git_commit(),
        "timing_s": {"wall": round(time.time() - t0, 2), "cpu": round(time.process_time() - c0, 2)},
    }
    if write:
        OUT.parent.mkdir(parents=True, exist_ok=True)
        OUT.write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n")
    return out


def _comparable(o: dict) -> dict:
    return {k: o[k] for k in ("result", "negative_controls", "inputs_sha256", "main_blob_match",
                              "p2_4_inputs_unchanged", "pass")}


def verify() -> list[str]:
    """Recompute and compare with the committed file; returns a list of differences."""
    if not OUT.exists():
        return ["no committed assembly_beta_PB.json"]
    stored = json.loads(OUT.read_text())
    new = json.loads(json.dumps(assemble(write=False), ensure_ascii=False))
    a, b = _comparable(stored), _comparable(new)
    return [k for k in a if a[k] != b[k]]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m p2.beta_assembly", description=__doc__.split("\n")[0])
    ap.add_argument("cmd", nargs="?", default="run", choices=["run", "verify"])
    args = ap.parse_args(argv)
    if args.cmd == "verify":
        diff = verify()
        print("verify:", "identical" if not diff else f"differs in {diff}")
        return 0 if not diff else 1
    out = assemble(write=True)
    res = out["result"]
    print(f"assembly beta {NAME}: pass = {out['pass']}  (checks: {len(res['checks'])}, failed: {res['failed']})")
    for c in out["negative_controls"]:
        print(f"  {c['id']:5s} {c['what']:42s} fails as expected: {c['fails_as_expected']}  ({len(c['failed'])} failed)")
    print(f"written {_rel(OUT)}  ({out['timing_s']['cpu']} s CPU)")
    return 0 if out["pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
