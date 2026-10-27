"""S6: the two-sided beta at P_B (P2-6; notes/s6-preregistration.md, notes/s6-proof.md).

Fast (~1 min): Lemma E (the stored certificate; the sympy first step and Wt_0(omega0) recomputed live; the
closed form Wt_l = -(1/r_+) int (L + r F') Z2 dr/r of notes/s6-proof.md section 2.4 against the stored
balls, in floats); the part (i) cover replayed; one T box, the A box, part (ii') at l = 3 and three L grid
boxes (one contains omega0) recomputed; the assembly check passes, its negative controls fail on their
named items, and the stored certificate matches a recomputation; two stored l = 0 band pieces recomputed;
the ledger levels used for P2-6's level match CLAIMS.md.
P2_SLOW=1: one S job per side, the KAT boxes B2 and B3 (l = 0), the B1 Krawczyk replay, and a 5% replay of
the band pieces for l <= 2.
"""
import json
import math
import os
from fractions import Fraction

import pytest
from flint import arb

from p2 import beta_assembly, beta_common, beta_epole, beta_part1, beta_part3
from p2.arbutil import precision
from p2.params import Point

SLOW = pytest.mark.skipif(os.environ.get("P2_SLOW") is None, reason="P2_SLOW=1")
CERT = beta_common.CERT


def _load(rel):
    return json.loads((CERT / rel).read_text())


# ---------------------------------------------------------------------------------------------
# Lemma E

def test_lemma_E_certificate():
    ep = _load("epole_PB.json")
    assert ep["pass"] and ep["all_l_le_19_nonzero"] and ep["B5"]["pass"]
    assert ep["sympy_first_step"]["pass"] and ep["sympy_first_step"]["D1_at_omega0"] == "0"
    assert ep["V_l_positive_for_l_ge_1"]["pass"]
    with precision(128):
        for r in ep["Wt_l_le_19"] + ep["Wt_samples_l_ge_20"]:
            re_ = beta_assembly._re_ball(r["Wt"])
            assert not re_.contains(0), r["l"]
            assert (re_ > 0) if r["l"] == 0 else (re_ < 0), r["l"]       # the closed form's signs
    # B5 is O(eps): the relative difference scales by 100 between eps = 1e-5 and 1e-7
    rows = ep["B5"]["rows"]
    for a, b in zip(rows[0::2], rows[1::2]):
        ratio = float(arb(a["rel_diff_zsol"]).mid()) / float(arb(b["rel_diff_zsol"]).mid())
        assert 99 < ratio < 101, (a["l"], ratio)


def test_lemma_E_live():
    assert beta_epole.sympy_first_step()["pass"]
    stored = {r["l"]: r for r in _load("epole_PB.json")["Wt_l_le_19"]}
    pt = Point("PB", 256)
    with precision(128):
        w0, res, _ = beta_epole.wt_ball(pt, 0)
        assert not w0.real.contains(0)
        assert w0.real.overlaps(beta_assembly._re_ball(stored[0]["Wt"]))


def _closed_form(l, t0=0.01, dl=1e-7):
    """-(1/r_+) int_{r_+}^{r_c} (L + r F'(r)) Z2(r) dr/r in floats (scipy; Z2 = t + ... at omega0)."""
    integrate = pytest.importorskip("scipy.integrate")
    import mpmath as mp
    pt = Point("PB", 256)
    rp, rc, kp, kc = (float(x.mid()) for x in (pt.rp, pt.rc, pt.kp, pt.kc))
    lam, q2 = float(pt.Lam.mid()), float(pt.Q2.mid())
    D = lambda r: -(lam / 3) * r ** 4 + r * r - 2 * r + q2
    Dp = lambda r: -(4 * lam / 3) * r ** 3 + 2 * r - 2
    L = l * (l + 1)
    with mp.workdps(30):
        A_ = mp.taylor(lambda t: D(rp + t), 0, 4)
        B_ = mp.taylor(lambda t: Dp(rp + t) - 2 * kp * (rp + t) ** 2, 0, 4)
        C_ = [-(2 * kp * rp + L), -2 * kp]
        z = [mp.mpf(0), mp.mpf(1)]
        for n in range(1, 60):
            s = sum(A_[k] * (n + 2 - k) * (n + 1 - k) * z[n + 2 - k] for k in range(2, 5) if 0 <= n + 2 - k < len(z))
            s += sum(B_[k] * (n + 1 - k) * z[n + 1 - k] for k in range(1, 5) if 0 <= n + 1 - k < len(z))
            s += sum(C_[k] * z[n - k] for k in range(2) if 0 <= n - k < len(z))
            z.append(-s / ((n + 1) * n * A_[1]))
        g = lambda r: (L + Dp(r) / r - 2 * D(r) / r ** 2) / r
        Z0 = float(sum(z[n] * mp.mpf(t0) ** n for n in range(len(z))))
        Z0d = float(sum(n * z[n] * mp.mpf(t0) ** (n - 1) for n in range(1, len(z))))
        I0 = float(mp.quad(lambda t: g(rp + t) * sum(z[n] * t ** n for n in range(len(z))), [0, t0]))

    def f(r, y):
        zz, zd, _ = y
        return [zd, -((Dp(r) - 2 * kp * r * r) * zd - (2 * kp * r + L) * zz) / D(r), g(r) * zz]

    end = {}
    for d in (dl, 2 * dl):
        sol = integrate.solve_ivp(f, [rp + t0, rc - d], [Z0, Z0d, I0], method="DOP853", rtol=1e-13, atol=1e-30)
        end[d] = sol.y[:, -1]
    p = kp / kc                                           # Z2 ~ a (r_c - r)^(-p) + b near r_c
    a = (end[dl][0] - end[2 * dl][0]) / (dl ** -p - (2 * dl) ** -p)
    b = end[dl][0] - a * dl ** -p
    tail = (L + Dp(rc) / rc) / rc * (a * dl ** (1 - p) / (1 - p) + b * dl)
    return -(end[dl][2] + tail) / rp


def test_lemma_E_closed_form():
    stored = {r["l"]: r for r in _load("epole_PB.json")["Wt_l_le_19"]}
    for l in (0, 1, 2, 5):
        ref = float(beta_assembly._re_ball(stored[l]["Wt"]).mid())
        val = _closed_form(l)
        assert abs(val - ref) <= 1e-8 * abs(ref), (l, val, ref)


# ---------------------------------------------------------------------------------------------
# part (i) at depth e (l >= 20) and part (ii') (l <= 19)

def test_part1_cover_replay():
    stored = _load("part1/cover_PB.json")
    new = beta_part1.cover()
    for k in ("constants", "inequalities", "certificates_pass", "tiling", "cited", "pass", "h0"):
        assert new[k] == stored[k], k
    assert stored["pass"] and math.ceil(1 / Fraction(stored["h0"]) - Fraction(1, 2)) == 20


def test_part1_T_A_boxes_replay():
    from p2 import part1_T as P
    pt = beta_common.deep_point()
    T = _load("part1/T_PB.json")
    A6 = _load("part1/A_PB.json")
    with precision(beta_common.PREC):
        _, rays, tails = P.setup(pt)
        bx = T["boxes"][0]
        zb = tuple(Fraction(z) for z in bx["zeta_over_OmegaP"])
        r = P.certify_box(P.TChecker(pt, rays, tails, beta_common.T_DESIGN), pt, zb, beta_common.T_DESIGN,
                          h0=beta_common.H0)
        assert r["pass"] and r["mu_over_a_star"] == bx["mu_over_a_star"] and r["cells"] == bx["cells"]
        ba = A6["boxes"][0]
        za = tuple(Fraction(z) for z in ba["zeta_over_OmegaP"])
        r = P.certify_box(P.TChecker(pt, rays, tails, beta_common.A_DESIGN), pt, za, beta_common.A_DESIGN,
                          h0=beta_common.H0)
        assert r["pass"] and r["mu_over_a_star"] == ba["mu_over_a_star"]


def test_part2_l3_replay():
    from p2 import part1_T as P
    row = [r for r in _load("part2/A_fixed_PB.json")["rows"] if r["l"] == 3][0]
    pt = beta_common.deep_point()
    with precision(beta_common.PREC):
        _, rays, tails = P.setup(pt)
        z = Fraction(row["zeta_A_over_OmegaP"])
        r = P.certify_box(P.TChecker(pt, rays, tails, P.DESIGN_A), pt, (z, z), P.DESIGN_A, h0=Fraction(2, 7))
        assert r["pass"] and r["mu_over_a_star"] == row["box"]["mu_over_a_star"] and r["jost_decay"]


def test_part1_L_boxes_replay():
    from p2 import lowfreq as LF
    from p2.part1 import load_bulk
    L = _load("part1/L_PB.json")
    assert L["refined_grid_boxes"] == 0                   # one leaf per grid box
    W_bound, _ = load_bulk()
    pt = beta_common.deep_point()
    leaves = {tuple(x["box"]): x for x in L["leaves_detail"]}
    with precision(LF.PREC):
        kp = float((pt.kp / pt.km).mid())
        boxes = beta_part1.l_boxes()
        pole = [b for b in boxes if b[0] == 0 and float(b[2]) <= -kp <= float(b[3])]   # contains omega0
        assert len(pole) == 1
        ch = beta_common.make_lchecker(pt, W_bound)
        for b in (pole[0], boxes[1], boxes[-1]):
            key = tuple(str(x) for x in b)
            new = json.loads(json.dumps(LF.certify_box_refined(ch, b, beta_common.L_MIN_BOX)))
            assert len(new) == 1 and new[0]["pass"], key
            assert new[0] == leaves[key], key


# ---------------------------------------------------------------------------------------------
# the assembly check

@pytest.fixture(scope="module")
def run():
    return json.loads(json.dumps(beta_assembly.assemble(write=False), ensure_ascii=False))


def test_assembly_passes(run):
    r = run["result"]
    assert r["pass"] and r["failed"] == [], r["failed"]
    s = r["summary"]
    assert s["l0"] == 20 and s["levels"]["minimum"] == "WRITTEN_PROOF"
    per_l = s["band_counts"]["per_l"]
    assert sorted(int(k) for k in per_l) == list(range(20))
    assert all(per_l[str(l)]["winding"] == (0 if l == 0 else -1) for l in range(20))
    assert s["beta"]["enclosure"].startswith("[0.909380978994805455524933762622 +/- ")
    assert s["dirty_tags"] == []
    for name in ("C1.pins_main_blobs", "C2.P2_4_inputs_unchanged", "C6.cover_replay_identical", "C7.no_l_gap",
                 "C8.omega_star_in_B0", "C8.beta_width", "C9.B1_krawczyk_replay"):
        assert r["checks"][name]["ok"], name


def test_assembly_negative_controls(run):
    ctr = run["negative_controls"]
    assert [c["id"] for c in ctr] == ["NB6a", "NB6b", "NB6c", "NB6d", "NB6e", "NB6f", "NB6g", "NB6h", "NB6i"]
    for c in ctr:
        assert c["fails_as_expected"], c["id"]
        assert set(c["expected_failures"]) <= set(c["failed"]), c["id"]


def test_assembly_stored(run):
    assert beta_assembly.OUT.exists(), "run python -m p2.beta_assembly"
    stored = json.loads(beta_assembly.OUT.read_text())
    assert stored["pass"] is True
    assert stored["inputs_sha256"] == run["inputs_sha256"], "an input certificate changed since the assembly ran"
    assert stored["result"] == run["result"]
    assert stored["negative_controls"] == run["negative_controls"]


def test_band_pieces_l0_recomputed():
    from p2 import part3
    bd = _load("part3/band_PB_l0.json")
    pcs = bd["pieces"]
    for r in (pcs[0], pcs[len(pcs) // 2]):
        out = part3._replay_task(("PB", 0, r))
        assert out["overlap"], r["key"]


def test_levels_match_claims():
    claims = beta_common.ROOT.parents[1] / "CLAIMS.md"
    if not claims.exists():
        pytest.skip("CLAIMS.md not found")
    rows = {}
    for line in claims.read_text().splitlines():
        if line.startswith("| P2-"):
            cols = [c.strip() for c in line.split("|")]
            rows[cols[1]] = cols[4].split()[0].rstrip(".,")
    lv = beta_assembly.INGREDIENTS
    assert rows["P2-4"] == lv["P2-4 (S_delta, every l)"]
    assert rows["P2-D1′"] == lv["P2-D1' (set clause, any depth)"]
    assert rows["P2-K6"] == lv["P2-K6 Krawczyk ball of omega_*"]


# ---------------------------------------------------------------------------------------------
# slow

@SLOW
def test_part1_S_extended_grid_box():
    # the S+ box [0.06, 0.08] Omega_P needed the extended design grid (pre-registration O4)
    from p2 import throat as T
    from p2.part1 import load_bulk
    beta_common.use_extended_s_grid()
    W_bound, _ = load_bulk()
    st = _load("part1/S_PB_plus.json")
    stored = [x for x in st["leaves"] if x["zeta_over_OmegaP"] == ["3/50", "2/25"]]
    assert len(stored) == 1
    pt = beta_common.deep_point()
    with precision(T.PREC):
        new = T.certify_zbox_refined(pt, T.FGeo(pt), 1, Fraction(3, 50), Fraction(2, 25), W_bound)
    new = json.loads(json.dumps(new))
    assert len(new) == 1 and new[0]["pass"]
    # the float design margin (used only to pick the design) may differ in the last digits; the
    # certified record must be identical
    fm_new, fm_old = new[0].pop("float_margin"), dict(stored[0]).pop("float_margin")
    assert abs(fm_new - fm_old) <= 1e-9 * abs(fm_old)
    assert new[0] == {k: v for k, v in stored[0].items() if k != "float_margin"}


@SLOW
def test_kat_boxes_B2_B3():
    from p2 import part3
    for bx, expect in (((-0.03, 0.03, -0.93, -0.89), 1), ((-0.03, 0.03, -0.86, -0.79), -1)):
        r = part3.rectangle_winding("PB", 0, bx, workers=1)
        assert r.get("certified") and r.get("winding") == expect, (bx, r.get("winding"))


@SLOW
def test_B1_krawczyk_replay():
    from p2 import part3
    stored = json.loads((beta_common.ROOT / "certificates" / "part3" / "krawczyk_PB_l0.json").read_text())
    kr = part3.krawczyk_z("PB", 0, ("0", "-0.9093809789948054555249"))
    assert kr["certified_unique_zero"]
    with precision(256):
        assert arb(kr["enclosure_over_kappa_minus"]["im"]).overlaps(arb(stored["enclosure_over_kappa_minus"]["im"]))


@SLOW
def test_band_replay_l_le_2():
    from p2 import part3
    for l in (0, 1, 2):
        pcs = _load(f"part3/band_PB_l{l}.json")["pieces"]
        k = max(3, math.ceil(0.05 * len(pcs)))
        for j in sorted({(i * 7919 + l * 104729) % len(pcs) for i in range(k)}):
            assert part3._replay_task(("PB", l, pcs[j]))["overlap"], (l, j)
