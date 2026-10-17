"""S4-3 tests: the commutator strip-exclusion certificate (tc_commutator).

Fast by default (~2-3 min: exact identities, helper algebra, one sonic
closing-piece direct run, one centre piece).  The full-piece analyze and the
Lemma D rerun are gated behind P4_SLOW=1.

Measured 2026-10-16 (M4 Pro): the fast set takes ~3 min without
test_centre_first_piece; that test alone did not finish within 15 min.
test_commutation_value_consistency and test_kdiv_roundtrip fail on their
relative-width thresholds although the compared balls overlap (see
notes/s4-3-recovery.md section 3).
"""
import json
import os

import pytest
from flint import arb, acb

import p4.validated.tc_commutator as tc
import p4.validated.tc_energy as te
from p4.validated.arbseries import precision
from p4.validated.tc_tm import sonic_components, AtomKP, LC, TM, SPLIT_REL

SLOW = pytest.mark.skipif(os.environ.get("P4_SLOW") is None,
                          reason="P4_SLOW=1 for the long runs")
RESULTS = os.path.join(os.path.dirname(__file__), "..", "results", "theorem_c")


@pytest.fixture(scope="module")
def atoms():
    return te.Atoms(full=True)


@pytest.fixture(scope="module")
def sonic_pr(atoms):
    """Raised pairs + commutation on the sonic closing domain [0, 1e-6]."""
    ex = te.sonic_expansion_cached(320)
    with precision(320):
        comps = sonic_components(ex, 1e-6)
        ak = AtomKP(atoms, comps, sd_shift=1)
        pr = tc.sonic_raise_inputs(tc.build_pairs(ak))
        cm = tc.commute(pr, pr["dT"][-1][0].trim().inv())
    return pr, cm


# ---------------------------------------------------------- exact identities
def test_e1_u6_top_identity(atoms):
    assert tc.check_u6_top_identity(atoms)


def test_e3_sonic_zeros(atoms):
    assert tc.check_sonic_zeros(atoms) == (True, True, True)


def test_e5_sonic_vanishing(atoms):
    # every kappa-coefficient of U6 vanishes at u0; Aer_3(u0) = Ber_3(u0) = 0
    assert tc.check_sonic_vanishing(atoms) == (True, True, True)


@SLOW
def test_e2_layer_identity(atoms):
    # check_layer_identity returns (pmin, ok): pmin is the lowest raw centre
    # layer of Z (before A -> 1) and ok says that layers pmin, pmin + 1 cancel
    # with A -> 1 and no odd layers exist; with A - 1 = O(t^2) this gives
    # val(Z o u) >= pmin + 2.  run_lemma_d raises Z to offset 16, so the
    # requirement is pmin + 2 >= 16 (the recorded value is pmin = 14).  The
    # previous assertion `val >= 16` compared pmin itself with 16 and failed.
    val, ok = tc.check_layer_identity(atoms)
    assert ok and val + 2 >= 16


# ------------------------------------------------ certified centre background
def test_centre_expansion_encloses_certified_background():
    """_centre_expansion (Lemma D and the centre region) is built at (a, mu)
    balls containing the certified (a*, mu*) = a3 midpoint + K_minus_m (the
    balls of modecount.certified_centre), not at the bare midpoint, which
    misses mu* by 2.4e-9 (notes/s4-3-recovery.md section 3)."""
    from p4.validated import modecount
    ce = tc._centre_expansion()
    with precision(256):
        a = arb(modecount.A_STAR) + arb("7.031e-12") + arb(0, 8.87e-16)
        mu = arb(modecount.MU_STAR) + arb("2.3757e-9") + arb(0, 4.61e-14)
        assert ce.what.contains(mu * (2 * a).exp())
        assert ce.nhat.contains((-a).exp())
        a_m, mu_m = arb(modecount.A_STAR), arb(modecount.MU_STAR)
        assert not ce.what.overlaps(mu_m * (2 * a_m).exp())


def test_lemma_d_requires_layer_valuation_16(monkeypatch):
    """run_lemma_d shifts Z to offset 16 on the strength of E2; it must refuse
    when the E2 layer minimum only gives val(Z o u) = pmin + 2 < 16."""
    monkeypatch.setattr(tc, "check_layer_identity", lambda at: (13, True))
    with pytest.raises(RuntimeError, match="pmin"):
        tc.run_lemma_d(None)


def test_direct_only_needs_a_direct_range(sonic_pr):
    """tq=False with tstar_fac <= 1 would certify nothing: refused."""
    pr, cm = sonic_pr
    with precision(320):
        with pytest.raises(ValueError):
            tc.analyze_piece(pr, cm, 0.0, 1e-6, 14.0, [(0.0, 0.03)], 2,
                             (0, 0, 0), dp_raise=1, tq=False, tstar_fac=1.0)


# ------------------------------------------------------------ helper algebra
def test_bnd_arithmetic():
    b = tc.Bnd.of(2.0, -1) + tc.Bnd.of(3.0, 0)
    up = float(b.eval_up(0.5, 1.0, 1))
    # 2 ell^-1 + 3 evaluated with one ell of scaling on [1/2, 1]
    assert abs(up - (2.0 + 3.0)) < 1e-6


def test_lc_split_lead_exact():
    x = LC(TM([arb(2), arb(5), arb(7)], 1.0), -3)
    lead, rest = tc.lc_split_lead(x, 1)
    assert lead.off == -3 and rest.off == -2
    assert float(lead.tm.c[0]) == 2.0 and float(rest.tm.c[0]) == 5.0


def test_kp_helpers():
    mk = lambda v, off: LC(TM([arb(v)], 1.0), off)
    a = [mk(1, 0), mk(2, 0)]           # 1 + 2k
    b = [mk(3, 0), mk(4, 0)]           # 3 + 4k
    prod = tc.kp_mul(a, b)             # 3 + 10k + 8k^2
    assert [float(x.tm.c[0]) for x in prod] == [3.0, 10.0, 8.0]
    dk = tc.kp_dk(prod)                # 10 + 16k
    assert [float(x.tm.c[0]) for x in dk] == [10.0, 16.0]
    sub = tc.kp_sub(prod, prod)
    assert all(float(abs(x.tm.c[0]).upper()) == 0.0 for x in sub)


# ------------------------------------------- commutation consistency (values)
def test_commutation_value_consistency(sonic_pr):
    """A1 = A0 - (log B0)' and B1 = B0 + A0' - A0 (log B0)' hold as rational
    value identities: check nA1/dA1 and nB1/dB1 against jet chains of the
    input pairs at a sample (kappa, xi)."""
    pr, cm = sonic_pr
    with precision(320):
        xi_lo, xi_hi = 8e-7, 8.2e-7
        kap = acb(1.0, 37.0)

        def kp_val(kp, sl):
            acc = acb(0)
            for x in reversed([e[sl] for e in kp]):
                x = x.trim()
                acc = acc * kap + x.tm.eval(xi_lo, xi_hi) * arb(xi_lo) ** x.off
            return acc

        nT = [kp_val(pr["nT"], s) for s in (0, 1)]
        dT = [kp_val(pr["dT"], s) for s in (0, 1)]
        nB = [kp_val(pr["nB"], s) for s in (0, 1)]
        dB = [kp_val(pr["dB"], s) for s in (0, 1)]
        A0 = nT[0] / dT[0]
        A0p = (nT[1] * dT[0] - nT[0] * dT[1]) / (dT[0] * dT[0])
        B0 = nB[0] / dB[0]
        B0p = (nB[1] * dB[0] - nB[0] * dB[1]) / (dB[0] * dB[0])
        lB0 = B0p / B0
        A1 = A0 - lB0
        B1 = B0 + A0p - A0 * lB0
        A1c = kp_val(cm["nA1"], 0) / kp_val(cm["dA1"], 0)
        B1c = kp_val(cm["nB1"], 0) / kp_val(cm["dB1"], 0)
        assert abs(A1 - A1c).upper() / abs(A1).lower() < 1e-4
        assert abs(B1 - B1c).upper() / abs(B1).lower() < 1e-4


def test_kdiv_roundtrip(sonic_pr):
    """nB = sum q_j k^j dB + rem as values (division consistency)."""
    pr, _ = sonic_pr
    with precision(320):
        SPLIT_REL[0] = 1e-4
        dB = tc.slot0(pr["dB"])
        dBi = dB[-1].inv()
        dBtj = pr["dB"][-1][:3]
        dBi1 = -1 * (dBtj[1].trim() * (dBi * dBi))
        dBij = [dBi, dBi1, 2 * (dBtj[1].trim() * dBtj[1].trim()
                                * (dBi * dBi * dBi))
                - dBtj[2].trim() * (dBi * dBi)]
        uqj, remUj = tc.kdiv_jlc(pr["nB"], pr["dB"], dBij, 2, depth=2)
        xi_lo, xi_hi = 8e-7, 8.2e-7
        kap = acb(0.5, 29.0)

        def lval(x):
            x = x.trim()
            return x.tm.eval(xi_lo, xi_hi) * arb(xi_lo) ** x.off

        def pval(kp, sl=0):
            acc = acb(0)
            for e in reversed(kp):
                acc = acc * kap + lval(e[sl] if isinstance(e, list) else e)
            return acc

        lhs = pval(pr["nB"])
        q = sum(lval(uqj[j][0]) * kap ** j for j in range(3))
        rhs = q * pval(pr["dB"]) + pval([x[0] for x in remUj])
        assert abs(lhs - rhs).upper() / abs(lhs).lower() < 1e-8


# --------------------------------------------------------------- region runs
def test_sonic_closing_direct_phase(sonic_pr):
    """The closing piece [0, 1e-6]: direct phase (tau in [b1, tstar]) on two
    sample eps cells with the ell = 0-safe offset-aware chains."""
    pr, cm = sonic_pr
    with precision(320):
        dB = tc.slot0(pr["dB"])
        dBi = dB[-1].inv()
        dBtj = pr["dB"][-1][:3]
        dBi1 = -1 * (dBtj[1].trim() * (dBi * dBi))
        dBij = [dBi, dBi1, 2 * (dBtj[1].trim() * dBtj[1].trim()
                                * (dBi * dBi * dBi))
                - dBtj[2].trim() * (dBi * dBi)]
        uqj, remUj = tc.kdiv_jlc(pr["nB"], pr["dB"], dBij, 2, depth=2)
        ok, worst, fails = tc.direct_S_boxes(pr, cm, uqj, remUj, 0.0, 1e-6,
                                             14.0, 224.0,
                                             [(0.0, 0.03), (3.4, 4.6)], 2)
    assert ok and not fails
    assert worst < -0.05


def test_centre_first_piece(atoms):
    ce = tc._centre_expansion()
    rep = tc.centre_piece(atoms, ce, 1e-8, 1.6e-8, 14.0,
                          [(0.0, 0.03), (3.4, 4.6)], 320)
    assert rep["ok"], rep


@SLOW
def test_sonic_laurent_piece_full(atoms):
    r = tc.run_sonic(atoms, 14.0, piece_slice=(1, 2))
    assert r["ok"], r["detail"]
    assert r["detail"][0]["open"]  # tau > tstar tail marked open there


# ---------------------------------------------------------------- artifacts
def test_certificate_json_consistency():
    path = os.path.join(RESULTS, "tc_commutator.json")
    if not os.path.exists(path):
        pytest.skip("tc_commutator.json not yet merged")
    d = json.load(open(path))
    assert d["b1"] == 14.0 and d["m"] == 1
    assert d["checks"]["E1_U6top"] and all(d["checks"]["E3_sonic_zeros"])
    assert all(d["checks"]["E5_sonic_vanishing"])
    assert d["pieces"] == len(d["detail"])
    bad = [p for p in d["detail"] if not p["ok"]]
    assert not bad, bad[:3]


@SLOW
def test_lemma_d_rerun(atoms):
    res = tc.run_lemma_d(atoms)
    assert res["ok"] and res["mins"]["q2"] > 1.9


def test_lemma_d_recorded():
    path = os.path.join(RESULTS, "tc_energy.json")
    d = json.load(open(path))
    assert d["lemma_d"]["ok"] and d["all_ok"]
    assert d["lemma_d"]["mins"]["q2"] > 1.9
