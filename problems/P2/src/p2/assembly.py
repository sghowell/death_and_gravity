"""P2-4 at P_B (S5): a machine check of the assembly bookkeeping.

The written proof is notes/s5-assembly.md, and the plan is notes/s5-plan.md.  This module reads only
committed certificates and checks, with exact rationals and Arb balls, that they cover the strip

    S_delta = {Im w >= -d_delta},   d_delta = (1/2 + delta) kappa_-   (= (11/20) kappa_- at P_B),

for every degree l with no gap, and that the numerical hypotheses of [HV17] and [CF17] hold.

Decomposition, with U = {Im w > 0} and H_+- = {+-Re w >= 0, -d_delta <= Im w <= 0}:
  * U, every l:            P2-0 (W_l(w) = 0 with Im w >= 0 iff (l, w) = (0, 0)).
  * H_-, every l:          W_l(-conj w) = conj W_l(w) (FORMULATION 1.2; zeq_identity Z4) maps it onto H_+.
  * H_+, l >= 20:          P2-1, certificates/part1/cover_PB.json (h = 1/(l+1/2) <= h0 = 2/41).
  * H_+ u H_-, l <= 19:    R'_l = [-a_l, a_l] x [-b, e] (P2-3B counts) and {|Re w| >= a_l} x [-d_delta, 0] (P2-2).
  * E_+ u E_c lies below S_delta (P2-par), so W_l is analytic on S_delta and D1 applies there; the same
    inequality puts S_delta inside P2-D1''s multiplicity window Im w > -min(kappa_+, kappa_c).

Checks C1-C9 are listed in notes/s5-plan.md section 3, and the negative controls N1-N8 there too.  Each
control mutates a copy of the inputs and must make the check fail on its named item.

    python -m p2.assembly            # run, write certificates/assembly/assembly_PB.json
    python -m p2.assembly verify     # recompute and compare with the committed file
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import sys
import time
from fractions import Fraction
from pathlib import Path

from flint import arb

from .arbutil import A, git_commit, lohi_out, precision
from .params import Point

PREC = 256
ROOT = Path(__file__).resolve().parents[2]          # problems/P2
CERT = ROOT / "certificates"
OUT = CERT / "assembly" / "assembly_PB.json"
NAME = "PB"
L_SMALL = list(range(20))                            # l <= 19: P2-2 + P2-3B
REGIMES = ("bulk_PB", "T_PB", "A_PB", "S_PB_plus", "S_PB_c", "L_PB_plus", "L_PB_c")

LEVEL_ORDER = ["HEURISTIC", "CONJECTURED", "VERIFIED_N", "WRITTEN_PROOF", "CERTIFIED", "FORMALIZED"]
# Ledger levels on main 3dd0b2bc (CLAIMS.md; tests/test_assembly.py cross-checks this table), and the use
# each ingredient has in P2-4.
INGREDIENTS = {
    "P2-par": ("CERTIFIED", "E_+ u E_c below S_delta; S_delta inside D1''s window; alpha < gamma0; kappa_- < 2 kappa_+"),
    "P2-0": ("WRITTEN_PROOF", "Im w > 0 for every l; W_l(0) != 0 for l >= 1; D4 (W_0(0) = 0, simple)"),
    "P2-D1'": ("WRITTEN_PROOF", "zeros of W_l = resonances of the forward problem on Omega_23, multiplicity (2l+1) m_l"),
    "P2-1": ("WRITTEN_PROOF", "l >= 20 on the closed lower half-strip with Re w >= 0"),
    "P2-2": ("CERTIFIED", "l <= 19 on {|Re w| >= a_l} x [-d_delta, 0]"),
    "P2-3B": ("CERTIFIED", "l <= 19: counts on R'_l (1 at l = 0, 0 otherwise)"),
}
BETA_UPPER = Fraction(90938098, 10 ** 8)             # P2-K6, as stated in the ledger


# ---------------------------------------------------------------------------------------------
# inputs

def _rel(p: Path) -> str:
    return str(p.relative_to(ROOT))


def input_paths() -> dict:
    """Every committed file the check reads, keyed by role."""
    paths = {
        "params": CERT / "params.json",
        "cover": CERT / "part1" / "cover_PB.json",
        "replay_part1": CERT / "part1" / "replay_PB.json",
        "p2_0": CERT / "p2_0_identity.json",
        "zeq": CERT / "part3" / "zeq_identity.json",
        "krawczyk": CERT / "part3" / "krawczyk_PB_l0.json",
        "kat6": CERT / "kat6_winding_krawczyk.json",
    }
    for n in REGIMES:
        paths[f"regime:{n}"] = CERT / "part1" / f"{n}.json"
    for l in L_SMALL:
        paths[f"lg:{l}"] = CERT / "lg" / f"lg_{NAME}_l{l}.json"
        paths[f"wind:{l}"] = CERT / "part3" / f"winding_{NAME}_l{l}.json"
    return paths


def load_inputs() -> dict:
    """{'json': {role: parsed}, 'sha256': {relative path: digest}}; missing files are simply absent."""
    js, sha = {}, {}
    for role, p in input_paths().items():
        if p.exists():
            raw = p.read_bytes()
            js[role] = json.loads(raw)
            sha[_rel(p)] = hashlib.sha256(raw).hexdigest()
    return {"json": js, "sha256": sha}


def replay_cover() -> dict:
    """P2-1's cover check recomputed from the regime certificates on disk (part1.cover())."""
    from . import part1
    return part1.cover()


# ---------------------------------------------------------------------------------------------
# helpers

class Rec:
    """Named checks: name -> {'ok', optional 'margin', optional 'note'}."""

    def __init__(self):
        self.items: dict[str, dict] = {}

    def add(self, name: str, ok, margin: str | None = None, note: str | None = None) -> bool:
        ok = bool(ok)
        rec = {"ok": ok}
        if margin is not None:
            rec["margin"] = margin
        if note is not None:
            rec["note"] = note
        self.items[name] = rec
        return ok

    def failed(self) -> list[str]:
        return [k for k, v in self.items.items() if not v["ok"]]


def _F(s) -> Fraction:
    return Fraction(str(s))


def _ball(entry) -> arb:
    """An Arb ball from a params.json quantity ({'ball': '[m +/- r]'}) or a ball string."""
    return arb(entry["ball"] if isinstance(entry, dict) else entry)


def _m(x: arb) -> str:
    return lohi_out(x, 10)


def _pt(p) -> tuple[Fraction, Fraction]:
    return (_F(p[0]), _F(p[1]))


def _path_ok(pieces: list, a: Fraction, b: Fraction, e: Fraction) -> tuple[bool, str]:
    """The pieces trace the right half H of the boundary of [-a, a] x [-b, e] exactly, counterclockwise:
    (0, -b) -> (a, -b) -> (a, e) -> (0, e), contiguous, axis-parallel, monotone along each edge."""
    if not pieces:
        return False, "no pieces"
    if _pt(pieces[0]["p0"]) != (Fraction(0), -b):
        return False, "does not start at (0, -b)"
    if _pt(pieces[-1]["p1"]) != (Fraction(0), e):
        return False, "does not end at (0, e)"
    last_edge = 0
    for i, s in enumerate(pieces):
        (x0, y0), (x1, y1) = _pt(s["p0"]), _pt(s["p1"])
        ed = s["edge"]
        if ed < last_edge or ed not in (0, 1, 2):
            return False, f"piece {i}: edge order"
        last_edge = ed
        if ed == 0 and not (y0 == y1 == -b and 0 <= x0 < x1 <= a):
            return False, f"piece {i}: not on the bottom edge"
        if ed == 1 and not (x0 == x1 == a and -b <= y0 < y1 <= e):
            return False, f"piece {i}: not on the right edge"
        if ed == 2 and not (y0 == y1 == e and a >= x0 > x1 >= 0):
            return False, f"piece {i}: not on the top edge"
        if i + 1 < len(pieces) and _pt(pieces[i + 1]["p0"]) != (x1, y1):
            return False, f"gap after piece {i}"
    edges = [s["edge"] for s in pieces]
    if not (0 in edges and 1 in edges and 2 in edges):
        return False, "an edge is missing"
    return True, f"{len(pieces)} pieces"


def _discs_ok(pieces: list) -> tuple[bool, str]:
    """Each piece lies within rho/2 of its centre (exact)."""
    for i, s in enumerate(pieces):
        cx, cy = _pt(s["c"])
        r2 = (_F(s["rho"]) / 2) ** 2
        for p in (s["p0"], s["p1"]):
            x, y = _pt(p)
            if (x - cx) ** 2 + (y - cy) ** 2 > r2:
                return False, f"piece {i}: endpoint farther than rho/2 from the centre"
    return True, "all endpoints within rho/2"


def _poles_ok(pieces: list, kp: arb, kc: arb, b: Fraction) -> tuple[bool, str]:
    """The box [c +- rho]^2 of every piece (corner distance sqrt(2) rho) avoids E_+ u E_c = {-i n kappa}.
    n = 1, 2 are checked directly; n >= 3 by |c + i n kappa| >= n kappa + Im c >= 3 min(kappa) + Im c."""
    s2 = arb(2).sqrt()
    kmin = kp if bool(kp < kc) else kc
    for i, s in enumerate(pieces):
        cx, cy = (A(v) for v in _pt(s["c"]))
        rr = s2 * A(_F(s["rho"]))
        if not bool(cy >= -A(b)):
            return False, f"piece {i}: centre below the box"
        for kap in (kp, kc):
            for n in (1, 2):
                d2 = cx * cx + (cy + n * kap) ** 2
                if not bool(d2 > rr * rr):
                    return False, f"piece {i}: within sqrt(2) rho of -i {n} kappa"
        if not bool(3 * kmin + cy > rr):
            return False, f"piece {i}: the n >= 3 bound fails"
    return True, "all disc boxes clear of E_+ u E_c"


# ---------------------------------------------------------------------------------------------
# the check

def check(inp: dict, cover_replay: dict | None, pt: Point | None = None) -> dict:
    """Run C1-C9 on the inputs.  Returns {'checks', 'summary', 'pass', 'failed'}."""
    J = inp["json"]
    R = Rec()
    S: dict = {}
    with precision(PREC):
        # ---------------- C1: the point
        params = J.get("params")
        if params is None or "PB" not in params.get("points", {}):
            R.add("C1.params_present", False)
            return {"checks": R.items, "summary": S, "pass": False, "failed": R.failed()}
        P = params["points"]["PB"]
        Qn = P["quantities"]
        qfrac, delta = _F(P["qfrac"]), _F(P["delta"])
        conv = params.get("conventions", "")
        R.add("C1.qfrac", qfrac == Fraction(999, 1000), note=str(qfrac))
        R.add("C1.delta", delta == Fraction(1, 20), note=str(delta))
        R.add("C1.Lambda_M", "M = 1" in conv and "Lambda = 3/50" in conv, note=conv)
        km, kp, kc = _ball(Qn["kappa_minus"]), _ball(Qn["kappa_plus"]), _ball(Qn["kappa_c"])
        g0, OmP = _ball(Qn["gamma0"]), _ball(Qn["Omega_P"])
        rn, rm, rp, rc = (_ball(Qn[k]) for k in ("r_n", "r_minus", "r_plus", "r_c"))
        Q = _ball(Qn["Q"])
        if pt is not None:
            for key, mine, theirs in (("kappa_minus", pt.km, km), ("kappa_plus", pt.kp, kp), ("kappa_c", pt.kc, kc),
                                      ("gamma0", pt.gamma0, g0), ("Omega_P", pt.OmegaP, OmP), ("r_n", pt.rn, rn),
                                      ("r_minus", pt.rm, rm), ("r_plus", pt.rp, rp), ("r_c", pt.rc, rc), ("Q", pt.Q, Q)):
                R.add(f"C1.overlap.{key}", mine.overlaps(theirs))
        dfrac = Fraction(1, 2) + delta                  # d_delta / kappa_-
        d = A(dfrac) * km                               # d_delta
        S["d_delta_over_kappa_minus"] = str(dfrac)
        S["d_delta"] = _m(d)

        # ---------------- C2: l <= 19 (P2-2 + P2-3B)
        small = {}
        present_all = True
        for l in L_SMALL:
            lg, wd = J.get(f"lg:{l}"), J.get(f"wind:{l}")
            ok_p = R.add(f"C2.l{l}.present", lg is not None and wd is not None)
            present_all = present_all and ok_p
            if not ok_p:
                continue
            rec = {}
            R.add(f"C2.l{l}.lg_meta", lg.get("point") == "PB" and lg.get("l") == l and _F(lg.get("delta")) == delta
                  and _F(lg.get("qfrac")) == qfrac)
            crit = lg.get("criterion", {})
            R.add(f"C2.l{l}.lg_certified", crit.get("certified") is True and bool(arb(crit["value_upper"]) < 2),
                  margin=_m(2 - arb(crit["value_upper"])))
            R.add(f"C2.l{l}.lg_strip_depth", arb(lg["strip_depth_beta"]).overlaps(d))
            reg = wd.get("region", {})
            R.add(f"C2.l{l}.wind_meta", wd.get("point") == "PB" and wd.get("l") == l and wd.get("certified") is True
                  and reg.get("a_source") == f"certificates/lg/lg_PB_l{l}.json")
            a, b, e = _F(reg["a"]), _F(reg["b"]), _F(reg["e"])
            a_lg = _F(lg["W_l_exact"])
            R.add(f"C2.l{l}.a_identical", a == a_lg, note=f"a_l = {a} (P2-2: {a_lg})")
            R.add(f"C2.l{l}.b_ge_d", bool(A(b) >= d), margin=_m(A(b) - d))
            R.add(f"C2.l{l}.e_pos", e > 0, note=str(e))
            R.add(f"C2.l{l}.b_lt_kappa_plus", bool(A(b) < kp) and bool(A(b) < kc), margin=_m(kp - A(b)))
            pcs = wd.get("pieces", [])
            ok, msg = _path_ok(pcs, a, b, e)
            R.add(f"C2.l{l}.path_traces_H", ok, note=msg)
            ok, msg = _discs_ok(pcs)
            R.add(f"C2.l{l}.path_in_discs", ok, note=msg)
            ok, msg = _poles_ok(pcs, kp, kc, b)
            R.add(f"C2.l{l}.discs_clear_poles", ok, note=msg)
            tot = sum((arb(s["increment"]) for s in pcs), arb(0))
            wn = tot / arb.pi()
            k = wd.get("winding")
            R.add(f"C2.l{l}.resum", isinstance(k, int) and bool(abs(wn - k) < arb("0.25")), margin=wn.str(8))
            R.add(f"C2.l{l}.winding_value", k == (1 if l == 0 else 0), note=f"stored {k}")
            # containment of the lower strip: [-d_delta, 0] in [-b, e]; Re: [-a, a] u {|x| >= a_l} = R
            R.add(f"C2.l{l}.lower_strip_in_box_and_P2_2", bool(-A(b) <= -d) and e >= 0 and a == a_lg)
            rec.update({"a_l": str(a), "a_l_over_kappa_minus": _m(A(a) / km), "winding": k, "pieces": len(pcs),
                        "resum_over_pi": wn.str(6)})
            small[l] = rec
        R.add("C2.all_l_0_to_19_present", present_all and all(f"lg:{l}" in J and f"wind:{l}" in J for l in L_SMALL))
        bs = {J[f"wind:{l}"]["region"]["b"] for l in L_SMALL if f"wind:{l}" in J}
        es = {J[f"wind:{l}"]["region"]["e"] for l in L_SMALL if f"wind:{l}" in J}
        S["l_le_19"] = {"b": sorted(bs), "e": sorted(es), "pieces_total": sum(r["pieces"] for r in small.values()),
                        "per_l": small}

        # ---------------- C3: l >= 20 (P2-1)
        cov = J.get("cover")
        ok_c = R.add("C3.cover_present", cov is not None)
        if ok_c:
            R.add("C3.cover_pass", cov.get("pass") is True and all(cov.get("certificates_pass", {}).values())
                  and all(cov.get("inequalities", {}).values()) and all(cov.get("tiling", {}).values()))
            if cover_replay is not None:
                same = all(cover_replay.get(k) == cov.get(k)
                           for k in ("constants", "inequalities", "certificates_pass", "tiling", "cited", "pass"))
                R.add("C3.cover_replay_identical", same)
            C = {k: _F(v) for k, v in cov["constants"].items()}
            h0 = C["h0"]
            l0 = math.ceil(1 / h0 - Fraction(1, 2))
            R.add("C3.h0_is_2_41", h0 == Fraction(2, 41), note=str(h0))
            R.add("C3.no_l_gap", l0 == max(L_SMALL) + 1, note=f"l0 = ceil(1/h0 - 1/2) = {l0}; l <= {max(L_SMALL)} counted")
            S["l0"] = l0
            # regime certificates: pass, point, h0 and design constants agree with the cover
            for n in REGIMES:
                c = J.get(f"regime:{n}")
                R.add(f"C3.regime_pass.{n}", c is not None and c.get("pass") is True)
            for n in ("T_PB", "A_PB", "S_PB_plus", "S_PB_c", "L_PB_plus", "L_PB_c"):
                c = J.get(f"regime:{n}") or {}
                R.add(f"C3.regime_h0.{n}", c.get("point") == "PB" and c.get("h0") == str(h0))
            T, Ac = J.get("regime:T_PB") or {}, J.get("regime:A_PB") or {}
            Sp, Sc = J.get("regime:S_PB_plus") or {}, J.get("regime:S_PB_c") or {}
            Lp, Lc = J.get("regime:L_PB_plus") or {}, J.get("regime:L_PB_c") or {}
            try:
                R.add("C3.design_T", [_F(x) for x in T["design"]["zeta_over_OmegaP"]] == [C["zeta_TL"], C["zeta_TU"]]
                      and all(arb(bx["nu_low"]).contains(-d) for bx in T["boxes"]))
                R.add("C3.design_A", _F(Ac["boxes"][0]["zeta_over_OmegaP"][0]) == C["zeta_A"]
                      and arb(Ac["boxes"][0]["nu_low"]).contains(-d) and Ac["monotonicity_W7"]["pass"] is True)
                R.add("C3.design_S", all(_F(x["constants"]["Omega_S_over_km"]) == C["Omega_S_over_km"]
                                         and _F(x["constants"]["zeta_SU_over_OmegaP"]) == C["zeta_SU"]
                                         and x.get("tail_region_in_strip") is True for x in (Sp, Sc)))
                R.add("C3.design_L", all(_F(x["design"]["Omega_L_over_km"]) == C["Omega_L_over_km"]
                                         and [_F(v) for v in x["design"]["Im_range_over_km"]] == [-dfrac, Fraction(0)]
                                         and [_F(v) for v in x["design"]["h2_range"]] == [Fraction(0), h0 * h0]
                                         for x in (Lp, Lc)))
            except (KeyError, IndexError, TypeError):
                R.add("C3.design_readable", False)
            rp1 = J.get("replay_part1") or {}
            R.add("C3.part1_full_replay", rp1.get("pass") is True and rp1.get("stride") == 1
                  and rp1.get("mode") == "full" and rp1.get("problems") == [])
            # Re w coverage for every h in (0, h0]: L [0, 8 km], S [6 km, zSU OmP/h], T [zTL, zTU] OmP/h,
            # A [zA OmP/h, oo).  OmP/h decreases in h, so the worst case of the one non-exact inequality is h = h0.
            R.add("C3.L_S_overlap", C["Omega_S_over_km"] < C["Omega_L_over_km"])
            R.add("C3.S_T_overlap", C["zeta_TL"] < C["zeta_SU"])
            R.add("C3.T_A_overlap", C["zeta_A"] < C["zeta_TU"])
            top_S = A(C["zeta_SU"]) * OmP / A(h0)
            R.add("C3.S_reaches_above_L", bool(top_S > A(C["Omega_L_over_km"]) * km),
                  margin=_m(top_S / km - A(C["Omega_L_over_km"])) + " kappa_-")
            R.add("C3.S_tails_decay", bool(A(C["Omega_S_over_km"]) * km > d))
            s_tiles = [v for k, v in cov.get("tiling", {}).items() if k.startswith("S ")]
            im_ok = (all([_F(v) for v in x.get("design", {}).get("Im_range_over_km", [])] == [-dfrac, Fraction(0)]
                         for x in (Lp, Lc))
                     and all(arb(bx["nu_low"]).contains(-d) for bx in T.get("boxes", []) + Ac.get("boxes", []))
                     and len(s_tiles) == 2 and all(s_tiles))
            R.add("C3.Im_range", im_ok, note=f"[-d_delta, 0], d_delta/kappa_- = {dfrac}: L's design range, T's and A's "
                                             "nu_low balls contain -d_delta, S's eps-boxes tile [0, eps_max] (cover)")
            zq = J.get("zeq") or {}
            R.add("C3.symmetry_Z4", zq.get("Z4") is True and zq.get("pass") is True)
            S["l_ge_20"] = {"h0": str(h0), "l0": l0, "Re_cover": "[0, 8 km] u [6 km, (9/10) OmP/h] u "
                            "[(4/5) OmP/h, 2 OmP/h] u [(3/2) OmP/h, oo) = [0, oo) for every h <= h0",
                            "S_top_at_h0_over_kappa_minus": _m(top_S / km)}

        # ---------------- C4: poles and D1''s window
        R.add("C4.E_plus_below_strip", bool(d < kp), margin=_m(kp - d))
        R.add("C4.E_c_below_strip", bool(d < kc), margin=_m(kc - d))
        ineq = P.get("inequalities", {})
        R.add("C4.params_flags", ineq.get("strip_vs_kp", {}).get("certified") is True
              and ineq.get("strip_vs_kc", {}).get("certified") is True)

        # ---------------- C5: the zero mode
        w0 = J.get("wind:0")
        if w0 is not None:
            a0, b0, e0 = (_F(w0["region"][k]) for k in ("a", "b", "e"))
            R.add("C5.zero_interior_R0", -a0 < 0 < a0 and -b0 < 0 < e0)
        else:
            R.add("C5.zero_interior_R0", False, note="no l = 0 certificate")
        dW = (rp * rp + rc * rc) / (rp * rc)            # W_0'(0) = i (r_+^2 + r_c^2)/(r_+ r_c)
        R.add("C5.D4_derivative_nonzero", bool(dW > 0), margin="i " + _m(dW))
        p20 = J.get("p2_0") or {}
        R.add("C5.P2_0_identities", bool(p20.get("steps")) and all(s.get("ok") is True for s in p20["steps"]))
        S["zero_mode"] = {"W0_prime_0": "i " + _m(dW), "count_R0": (w0 or {}).get("winding")}

        # ---------------- C6: beta
        kz = J.get("krawczyk") or {}
        try:
            im = arb(kz["enclosure_over_kappa_minus"]["im"])   # omega_*/kappa_-, imaginary part
            R.add("C6.krawczyk_certified", kz.get("certified_unique_zero") is True and kz.get("point") == "PB"
                  and kz.get("l") == 0)
            R.add("C6.beta_upper", bool(-im <= A(BETA_UPPER)), margin=_m(A(BETA_UPPER) + im))
            R.add("C6.omega_star_outside_strip", bool(-im > A(dfrac)), margin=_m(-im - A(dfrac)))
            k6 = [c for c in (J.get("kat6") or {}).get("cases", []) if c.get("point") == "PB" and c.get("l") == 0]
            R.add("C6.kat6_consistent", len(k6) == 1 and k6[0]["krawczyk"]["certified_unique_zero"] is True
                  and arb(k6[0]["krawczyk"]["enclosure_over_kappa_minus"]["im"]).overlaps(im))
            S["beta"] = {"lower": str(dfrac), "upper": f"{BETA_UPPER.numerator / BETA_UPPER.denominator:.8f}",
                         "omega_star_over_kappa_minus_im": im.str(15)}
        except (KeyError, TypeError):
            R.add("C6.krawczyk_readable", False)
        R.add("C6.beta_lower_is_half_plus_delta", dfrac == Fraction(11, 20))

        # ---------------- C7: [HV17] numerical hypotheses
        R.add("C7.Lambda_M_Q_positive", "Lambda = 3/50" in conv and "M = 1" in conv and bool(Q > 0))
        roots_ok = bool(rn < 0) and bool(rm > 0) and bool(rm < rp) and bool(rp < rc)
        R.add("C7.three_simple_positive_roots", roots_ok
              and P.get("defining_polynomials", {}).get("all_roots_located") is True
              and bool(km > 0) and bool(kp > 0) and bool(kc > 0),
              note="Delta is a quartic with four distinct real roots r_n < 0 < r_- < r_+ < r_c, all simple")
        R.add("C7.alpha_positive", bool(d > 0))
        R.add("C7.alpha_lt_gamma0", bool(d < g0), margin=_m(g0 - d))
        expo = Fraction(1, 2) + dfrac                       # 1/2 + alpha/kappa_1, kappa_1 = kappa_-
        R.add("C7.exponent_gt_1", expo > 1, note=f"1/2 + alpha/kappa_- = {expo}")
        S["hv17"] = {"alpha": "(11/20) kappa_-", "alpha_over_gamma0": _m(d / g0), "exponent": str(expo)}

        # ---------------- C8: [CF17] numerical hypotheses, with p = alpha
        R.add("C8.subextremal", roots_ok and bool(Q > 0))
        R.add("C8.p_lt_kappa_plus", bool(d < kp), margin=_m(kp - d),
              note="so min{p, kappa_+} = p")
        R.add("C8.two_p_gt_km", 2 * dfrac > 1, note=f"2p/kappa_- = {2 * dfrac}: holds for every delta > 0 (tautology)")
        R.add("C8.km_lt_2kp", bool(km < 2 * kp), margin=_m(2 * kp - km))
        R.add("C8.condition", (2 * dfrac > 1) and bool(km < 2 * kp),
              note="2 min{p, kappa_+} > kappa_-  <=>  2p > kappa_- and 2 kappa_+ > kappa_-")
        S["cf17"] = {"p": "alpha = (11/20) kappa_-", "two_min_over_kappa_minus": str(2 * dfrac),
                     "two_kappa_plus_over_kappa_minus": _m(2 * kp / km)}

        # ---------------- C9: levels
        lv = {k: v[0] for k, v in INGREDIENTS.items()}
        lmin = min(lv.values(), key=LEVEL_ORDER.index)
        S["levels"] = {"ingredients": lv, "minimum": lmin}
        R.add("C9.level_minimum", lmin == "WRITTEN_PROOF", note=lmin)

    failed = R.failed()
    return {"checks": R.items, "summary": S, "pass": not failed, "failed": failed}


# ---------------------------------------------------------------------------------------------
# negative controls (each must fail on its named item)

def _replace_coord(wd: dict, old: str, new: str) -> None:
    for s in wd["pieces"]:
        for key in ("p0", "p1", "c"):
            s[key] = [new if v == old else v for v in s[key]]


def _n1(inp, d_lower_k):
    w = inp["json"]["wind:0"]
    old = w["region"]["b"]
    new = str(Fraction(d_lower_k, 2 ** 28))
    w["region"]["b"] = new
    _replace_coord(w, str(-Fraction(old)), str(-Fraction(new)))
    return f"l = 0: b {old} -> {new} (just below d_delta)"


def _n2(inp):
    del inp["json"]["wind:7"]
    return "winding_PB_l7.json removed"


def _n3(inp):
    w = inp["json"]["wind:5"]
    old = w["region"]["a"]
    new = str(Fraction(old) - Fraction(1, 2 ** 24))
    w["region"]["a"] = new
    _replace_coord(w, old, new)
    return f"l = 5: a {old} -> {new}"


def _n4(inp):
    inp["json"]["cover"]["constants"]["h0"] = "2/43"
    return "cover h0 2/41 -> 2/43 (l0 = 21)"


def _n5(inp):
    w = inp["json"]["wind:0"]
    old = w["region"]["e"]
    w["region"]["e"] = "0"
    _replace_coord(w, old, "0")
    return f"l = 0: e {old} -> 0"


def _n6(inp):
    inp["json"]["wind:3"]["winding"] = 1
    return "l = 3: stored winding 0 -> 1"


def _n7(inp):
    inp["json"]["regime:T_PB"]["pass"] = False
    inp["json"]["cover"]["certificates_pass"]["T_PB"] = False
    return "T_PB pass flag -> false (regime certificate and cover)"


def _n8a(inp):
    q = inp["json"]["params"]["points"]["PB"]["quantities"]
    with precision(PREC):
        q["gamma0"] = {"ball": (arb(q["kappa_minus"]["ball"]) / 2).str(40, radius=True)}
    return "gamma0 -> kappa_-/2"


def _n8b(inp):
    q = inp["json"]["params"]["points"]["PB"]["quantities"]
    with precision(PREC):
        q["kappa_plus"] = {"ball": (arb(q["kappa_minus"]["ball"]) * 45 / 100).str(40, radius=True)}
    return "kappa_+ -> 0.45 kappa_-"


def controls_spec(d_lower_k: int) -> list:
    return [
        ("N1", "b just below d_delta", lambda i: _n1(i, d_lower_k), ["C2.l0.b_ge_d"]),
        ("N2", "one l's certificate missing", _n2, ["C2.l7.present", "C2.all_l_0_to_19_present"]),
        ("N3", "a gap in Re w between R'_5 and P2-2", _n3, ["C2.l5.a_identical"]),
        ("N4", "an l gap (h0 = 2/43)", _n4, ["C3.no_l_gap"]),
        ("N5", "e = 0", _n5, ["C2.l0.e_pos", "C5.zero_interior_R0"]),
        ("N6", "a wrong count", _n6, ["C2.l3.winding_value", "C2.l3.resum"]),
        ("N7", "a failed regime", _n7, ["C3.regime_pass.T_PB", "C3.cover_pass"]),
        ("N8a", "HV17 alpha < gamma0 is exercised", _n8a, ["C7.alpha_lt_gamma0"]),
        ("N8b", "CF17 kappa_- < 2 kappa_+ is exercised", _n8b, ["C8.km_lt_2kp", "C8.condition"]),
    ]


def run_controls(inp: dict, cover_replay: dict | None, pt: Point | None) -> list:
    with precision(PREC):
        km = arb(inp["json"]["params"]["points"]["PB"]["quantities"]["kappa_minus"]["ball"])
        d = A(Fraction(1, 2) + _F(inp["json"]["params"]["points"]["PB"]["delta"])) * km
        k = math.floor(float((d * 2 ** 28).lower()))
        while not bool(A(Fraction(k, 2 ** 28)) < d):
            k -= 1
    out = []
    for cid, what, mutate, expect in controls_spec(k):
        mi = copy.deepcopy(inp)
        how = mutate(mi)
        res = check(mi, cover_replay, pt)
        failed = res["failed"]
        out.append({"id": cid, "what": what, "mutation": how, "expected_failures": expect,
                    "failed": failed, "fails_as_expected": (not res["pass"]) and all(x in failed for x in expect)})
    return out


# ---------------------------------------------------------------------------------------------
# driver

def assemble(write: bool = True) -> dict:
    t0, c0 = time.time(), time.process_time()
    inp = load_inputs()
    pt = Point(NAME, PREC)
    cov = replay_cover()
    res = check(inp, cov, pt)
    ctr = run_controls(inp, cov, pt)
    ok = res["pass"] and all(c["fails_as_expected"] for c in ctr)
    out = {
        "claim": "P2-4 at P_B (assembly bookkeeping): the committed certificates cover S_delta = {Im w >= -(1/2+delta)"
                 " kappa_-} for every l with no gap; E_+ u E_c lies below S_delta; the zero at 0 is the simple l = 0"
                 " zero; and the numerical hypotheses of [HV17] and [CF17] hold with alpha = p = (1/2+delta) kappa_-",
        "point": {"name": NAME, "qfrac": "999/1000", "delta": "1/20", "Lambda": "3/50", "M": "1"},
        "decomposition": [
            "U = {Im w > 0}, every l: P2-0",
            "H_- = {Re w <= 0, -d_delta <= Im w <= 0}, every l: symmetry W_l(-conj w) = conj W_l(w) onto H_+",
            "H_+ = {Re w >= 0, -d_delta <= Im w <= 0}, l >= 20: P2-1 (cover_PB.json, h <= h0 = 2/41)",
            "H_+ u H_-, l <= 19: R'_l = [-a_l, a_l] x [-b, e] (P2-3B) and {|Re w| >= a_l} x [-d_delta, 0] (P2-2)",
            "E_+ u E_c below S_delta (P2-par); S_delta inside P2-D1''s window Im w > -min(kappa_+, kappa_c)",
        ],
        "result": res,
        "negative_controls": ctr,
        "inputs_sha256": inp["sha256"],
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
    return {"result": o["result"], "negative_controls": o["negative_controls"],
            "inputs_sha256": o["inputs_sha256"], "pass": o["pass"]}


def verify() -> list[str]:
    """Recompute and compare with the committed file; returns a list of differences."""
    if not OUT.exists():
        return ["no committed assembly_PB.json"]
    stored = json.loads(OUT.read_text())
    new = json.loads(json.dumps(assemble(write=False), ensure_ascii=False))
    a, b = _comparable(stored), _comparable(new)
    return [k for k in a if a[k] != b[k]]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m p2.assembly", description=__doc__.split("\n")[0])
    ap.add_argument("cmd", nargs="?", default="run", choices=["run", "verify"])
    args = ap.parse_args(argv)
    if args.cmd == "verify":
        diff = verify()
        print("verify:", "identical" if not diff else f"differs in {diff}")
        return 0 if not diff else 1
    out = assemble(write=True)
    res = out["result"]
    print(f"assembly {NAME}: pass = {out['pass']}  (checks: {len(res['checks'])}, failed: {res['failed']})")
    for c in out["negative_controls"]:
        print(f"  {c['id']:4s} {c['what']:45s} fails as expected: {c['fails_as_expected']}  ({len(c['failed'])} failed)")
    print(f"written {_rel(OUT)}  ({out['timing_s']['cpu']} s CPU)")
    return 0 if out["pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
