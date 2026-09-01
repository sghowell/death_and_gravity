"""S4-2 tests: the direct-energy certificate (tc_energy) and the P4-3 tube sign
certificate (tc_monotone), on coarse subsamples plus known-answer guards.

Fast by default (~1-2 min: exact algebra + sonic-region pieces).  With
P4_TUBE_CACHE pointing to the saved certified tube, the tube subsamples run too;
the exact field identity QN_red(u0) = 0 (~2 min) is gated behind P4_SLOW.
"""
import math
import os

import pytest
from flint import arb

import p4.validated.tc_energy as te
import p4.validated.tc_monotone as tm
from p4.validated.arbseries import precision

C0 = 6.0
S0 = 0.090105970507920


@pytest.fixture(scope="module")
def atoms():
    return te.Atoms(full=True)


@pytest.fixture(scope="module")
def sonic_ex():
    return te.sonic_expansion_cached()


def test_exact_reductions_build(atoms):
    # Atoms() verifies at build time: SD = Dp SD_red, D2 = Dp^3 D2_red, Bc = Dp^2
    # Bc_red, Ac = Dp^2 Ac_red, Ae/Be = Dp * red, QNm = Dp^10 QN_red, TN1 = Dp^5
    # TN1_red, and the exact reduced identity U1 = Dp^6 (Bcr UN3r - Acr EBc D2r).
    for n in ("QNr", "TN1r", "SDr", "D2r", "Bcr", "Acr", "Aer", "Ber"):
        assert n in atoms.cmp


def test_closing_box_and_sigma_known_answer(atoms, sonic_ex):
    ok, m, extra = te.closing_box(atoms, sonic_ex, C0)
    assert ok, m
    # x^2 Q -> (sigma(kappa)^2 - 1)/4: all three kappa-coefficients contained
    assert extra["known_answer_sigma"] == [True, True, True]
    # E_D = D(x)/x -> D_1 = 7.5606 (s2-theorem-b.md section 3.8)
    assert abs(extra["ED"][0] - 7.5605725) < 1e-5
    # x T(x; kappa) -> 1 + sigma(kappa) = (1 - s0)(1 - kappa) - s0 + s0 = ...
    t0, t1 = extra["xT_affine"]
    assert abs(t0 - (1 - S0)) < 1e-6 and abs(t1 + (1 - S0)) < 1e-6
    assert extra["xT_rem_bound"] < 1e-3


def test_sonic_ladder_subsample(atoms, sonic_ex):
    # a few ladder boxes at the working widths (w ~ |x|/1500)
    with precision(256):
        for (m0, w) in ((-0.045, 2e-5), (-0.01, 5e-6), (-0.001, 5e-7), (-1e-5, 5e-9)):
            ev = te.BoxEval(atoms, sonic_ex.eval(arb(m0)), sonic_ex.eval(arb(m0, w)), w)
            QNj, QDj, v = te.qn_qd(atoms, ev)
            ok, mg, _ = te.conditions_jet(QNj, QDj, sonic_ex.eval(arb(m0, w)),
                                          v["Bcr"], C0, ev, fac=v["fac"], p3j=v["P3"])
            assert ok, (m0, {k: float(x.lower()) for k, x in mg.items()})
            # q2 -> (1 - s0)^2/4 / x^2 near the sonic point
            q2 = float(mg["q2"].lower())
            assert q2 * m0 * m0 < 0.25 and q2 > 0


def test_monotone_sonic_and_centre():
    rs = tm.run_sonic()
    assert rs["ok"], rs
    assert rs["mins"]["dvrel"] > 0.4          # v_rel'(0) = 0.49 (float diagnostics)
    rc = tm.run_centre()
    assert rc["ok"], rc
    assert rc["mins"]["neg_drho"] > 50        # -(rho^'/t^2) -> 2|w_2| + w0 (A-1)/t^2 ...


TUBE = os.environ.get("P4_TUBE_CACHE")
needs_tube = pytest.mark.skipif(not (TUBE and os.path.exists(TUBE or "")),
                                reason="needs P4_TUBE_CACHE (saved certified tube)")


@needs_tube
def test_tube_energy_subsample(atoms):
    from p4.validated import shootsys as ss
    from p4.validated.lintube import Tube
    tube = Tube.load(TUBE)
    with precision(256):
        for idx in (0, 40, 120, 244):
            sd = tube.steps[idx]
            if float(sd.x) - sd.h < -4.5:
                continue
            h = arb(sd.h)
            nsub = max(1, min(512, int(sd.h * 1500 / (2 * abs(float(sd.x))))))
            ok = False
            while not ok and nsub <= 4096:    # the first sub-box, doubling as run_tube
                w = sd.h / (2 * nsub)
                sm = (-h) / (2 * nsub)
                Zm = [zi + arb(0, sd.eps_z) for zi in ss.horner_vec(sd.co, sm)]
                Zb = [zi + arb(0, sd.eps_z)
                      for zi in ss.horner_vec(sd.co, sm + arb(0, w))]
                ev = te.BoxEvalScaled(atoms, Zm, Zb, w)
                QNj, QDj, v = te.qn_qd(atoms, ev)
                u = [ev.Ab, Zb[0] / Zb[3], Zb[1] * Zb[3] ** 2, Zb[2] * Zb[3]]
                ok, mg, _ = te.conditions_jet(QNj, QDj, u, v["Bcr"], C0, ev,
                                              fac=v["fac"], p3j=v["P3"])
                nsub *= 2
            assert ok, (idx, float(sd.x), {k: float(x.lower()) for k, x in mg.items()})


@needs_tube
def test_monotone_tube_subsample():
    from p4.validated.lintube import Tube
    tube = Tube.load(TUBE)
    tube.steps = tube.steps[::20]
    r = tm.run_tube(tube)
    assert r["ok"], r


@pytest.mark.skipif(os.environ.get("P4_SLOW") is None, reason="P4_SLOW=1 for the field identity")
def test_field_identity_QN_red_vanishes_at_u0(atoms):
    assert te.field_identity(atoms)
