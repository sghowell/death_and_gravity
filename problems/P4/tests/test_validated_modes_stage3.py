"""Theorem B, Stage 3: the 4D sonic series and its kappa-Taylor model (linsonic4), the Taylor-model
propagation (lintaylor), the winding number of E around a rectangle (modecount), the interior
analyticity certificates (analyticity), Krawczyk enclosures and the gauge mode (krawczyk_kappa).

Fast tests (< 1 min) use short tubes or no tube.  Slow tests need the full A3 tube to x = -3:
set P4_TUBE_CACHE to a tube saved by ``Tube.save`` (built otherwise, ~3 min) -- they recompute a
small rectangle around kappa_1 (winding 1) and the Krawczyk enclosure of kappa_1; the full-R
result (winding 2) is verified from the JSON written by ``modecount.winding_number`` when
P4_STAGE3_WINDING points to it (notes/s2-theorem-b.md section 3 records the run)."""
import json
import os

import pytest
from flint import acb, arb
from p4.validated import analyticity, krawczyk_kappa, linmatch, linprop, linscaled, linsonic, linsonic4, linstep
from p4.validated import lintaylor, lintube, modecount
from p4.validated.arbseries import precision
from p4.validated.linsys import abs_up, to_arb

V0_EC, W_V0 = modecount.V0_EC, modecount.W_V0
KAPPA1, KGAUGE = modecount.KAPPA1, modecount.KGAUGE
SLOW = pytest.mark.skipif(os.environ.get("P4_SLOW") is None and os.environ.get("P4_TUBE_CACHE") is None,
                          reason="needs the full tube (P4_TUBE_CACHE) or P4_SLOW=1")


@pytest.fixture(autouse=True)
def _prec():
    with precision(256):
        yield


@pytest.fixture(scope="module")
def bg():
    with precision(256):
        return linmatch.box_background(V0_EC, W_V0, K=41)


def rad(z):
    return float(abs_up(z - z.mid()))


# ---------------------------------------------------------------------------------------------
# fast
# ---------------------------------------------------------------------------------------------
def test_4d_sonic_series_matches_stage1_and_is_regular_at_A0(bg):
    """The 4D series (A_p kept) at kappa_1, the gauge value, 0 and 15+14i: every coefficient ball
    overlaps the Stage-1 (3D) one and the values at x0 = -0.05 overlap; D is kappa-free with a single
    exponent sigma(kappa) ~ -(kappa + 0.099)/1.099; the certificate also holds at kappa = A_0 = 1.861
    and 1.5 where the 3D form degenerates; the constraint is satisfied as a ball identity."""
    for kc in (acb(KAPPA1), acb(KGAUGE), acb(0), acb(15, 14)):
        e4 = linsonic4.linear_sonic_expansion4(bg, kc, K=40)
        c4 = e4.certify()
        e3 = linsonic.linear_sonic_expansion(bg, kc, K=40, check=False)
        c3 = e3.certify()
        assert c4.ok and c3.ok and float(c4.nu) >= float(c3.nu)
        assert all(e4.balls[n][i].overlaps(([e3.Ap[n]] + list(e3.balls[n]))[i]) for n in range(41) for i in range(4))
        if float(c3.nu) > 0.05:
            assert all(a.overlaps(b) and rad(b) < 1e-11 for a, b in zip(e3.eval(-0.05), e4.eval(-0.05)))
        sg = c4.details["sigma"]
        assert (sg + (kc + arb("0.099")) / arb("1.099")).abs_upper() < 0.02 * (1 + float(abs(kc)))
        cr = linsonic4.constraint_residual(e4, 20)
        assert all(z.contains(acb(0)) for z in cr)
    for kc in (acb("1.8614267226"), acb("1.5")):
        e4 = linsonic4.linear_sonic_expansion4(bg, kc, K=40)
        assert e4.certify().ok and e4.cert.nu > 0.05


def test_4d_delta_model_remainders(bg):
    """kappa-Taylor model of the sonic data on |delta| <= w: relative remainder <= 1e-6 at w = 0.25
    (m = 8) on the far sides of R, <= 1e-5 at w = 2^-7 near kappa = 0 (pole at -0.099), and <= 1e-18
    at Krawczyk scale; the box tail certificate covers x0 = -0.05."""
    for kc, w, m, tol in ((acb(15, 14), 0.25, 8, 1e-6), (acb(0, 14), 0.25, 8, 1e-6), (acb(15), 0.25, 8, 1e-6),
                          (acb(0), 2.0**-7, 8, 1e-5), (acb(KAPPA1), 2e-7, 3, 1e-18), (acb(KGAUGE), 2e-5, 3, 1e-15)):
        ex = linsonic4.linear_sonic_expansion4(bg, kc, width=w, m=m, K=40)
        cert = ex.certify()
        assert cert.ok and cert.nu > 0.05
        co, rem = ex.delta_model(-0.05, w)
        mag = max(float(abs_up(c)) for c in co[0])
        assert max(float(r) for r in rem) < tol * mag
        assert max(rad(c) for c in co[0]) < 1e-11 * mag


def test_gauge_mode_identities(bg):
    """kappa-bar = 2 - A_0 + 2 W_0/3 as a ball (radius < 1e-15); the pure-gauge generator solves the
    4D linearised system and the constraint order by order (balls contain 0) and coincides with the
    4D sonic series over a box containing kappa-bar (all balls overlap)."""
    g = krawczyk_kappa.gauge_checks(bg)
    kb = g["kappa_bar"]
    assert kb.rad() < 1e-15 and abs(kb.mid() - arb(KGAUGE)) < 1e-10
    assert g["residual_ok"] and g["constraint_ok"] and g["series_overlap"] and g["cert_ok"]


def test_interior_box_certificates(bg):
    """(S) and (C) on sample boxes of the cover of R (w = 0.5 far from the pole, 0.25 next to it)."""
    ce, S4 = modecount.certified_centre(), linscaled.full_system()
    for (cx, cy, w) in ((0.25, 0.25, 0.25), (0.5, 13.5, 0.5), (14.5, 13.5, 0.5), (2.0, 0.0, 0.5), (7.5, -7.5, 0.5)):
        ok, d = analyticity.sonic_box(bg, cx, cy, w)
        assert ok and d["nu"] > 0.05 and d["re_sigma_max"] < 1 and d["Z"] < 1
        ok, d = analyticity.centre_box(S4, ce, cx, cy, w)
        assert ok and d["Z"] < 1 and d["g"] < 51
    assert not analyticity.sonic_box(bg, 0.5, 0.5, 0.5)[0]           # too close to the pole at -0.099: bisected


def test_taylor_model_propagation_contains_point_runs(bg):
    """Fundamental-matrix Taylor model on a short certified tube (-0.05 -> -0.08), |delta| <= 0.25,
    m = 8 around 15+14i: composed with the sonic model it contains the Stage-2 point propagation of
    the 4D sonic data at delta = 0, w, iw and a generic delta; relative width < 1e-6."""
    tube = lintube.Tube.build(V0_EC, W_V0, x_c=-0.08)
    S4 = linscaled.full_system()
    kc, w, m = acb(15, 14), 0.25, 8
    ex = linsonic4.linear_sonic_expansion4(bg, kc, width=w, m=m, K=40)
    assert ex.certify().ok
    co, rem = ex.delta_model(-0.05, w)
    T0 = to_arb(-0.05).exp()
    sc = [1 / (T0 * T0), T0, 1 / (T0 * T0), 1 / T0]
    coefs = [[c[i] * sc[i] for i in range(4)] for c in co]
    rem0 = [abs_up(rem[i] * sc[i]) for i in range(4)]
    Phi, Rcols, log = lintaylor.propagate_tm(tube, S4, kc, w, m)
    cf, rm = lintaylor.tm_apply(Phi, Rcols, coefs, rem0, w)
    assert all(l["epsY"] < 1e-20 and l["bound"] < 1e-11 for l in log[:-1])
    for dk in (acb(0), acb(w), acb(0, w), acb(-0.7 * w, 0.3 * w)):
        e2 = linsonic4.linear_sonic_expansion4(bg, kc + dk, K=40)
        e2.certify()
        y0 = [v * s for v, s in zip(e2.eval(-0.05), sc)]
        st = linstep.LohnerSet(y0)
        linprop.propagate(tube, S4, kc + dk, [st])
        val = lintaylor.tm_eval(cf, rm, dk)
        assert all(a.overlaps(b) for a, b in zip(val, st.hull()))
        assert max(rad(z) for z in val) < 1e-6 * max(float(abs_up(z)) for z in val)   # sonic remainder 9e-8 dominates


def test_segment_increment_is_rigorous_and_subdivision_independent():
    """Scalar Taylor-model arithmetic and the argument increment on a synthetic model E(delta) =
    (delta - 0.3 i)^2 e^{delta} truncated: increments for 4 and 32 pieces overlap, and a model whose
    enclosure contains 0 is rejected."""
    r = 0.25
    e = modecount.TM([acb(1) / arb(k).fac() for k in range(9)], arb(1e-9), r)   # e^delta + O(1e-9)
    lin = modecount.TM([acb(0, "-0.3"), acb(1)] + [acb(0)] * 7, arb(0), r)
    E = lin * lin * e
    a, b = modecount.segment_increment(E, r, acb(1), 4), modecount.segment_increment(E, r, acb(1), 32)
    assert a is not None and b is not None and a[0].overlaps(b[0]) and a[0].rad() < 1e-6
    bad = modecount.TM([acb(0), acb(1)] + [acb(0)] * 7, arb(0), r)                # delta: vanishes on the segment
    assert modecount.segment_increment(bad, r, acb(1), 8) is None


# ---------------------------------------------------------------------------------------------
# slow: the full tube (P4_TUBE_CACHE) -- winding number of a small rectangle, Krawczyk at kappa_1
# ---------------------------------------------------------------------------------------------
@pytest.fixture(scope="module")
def tube_path(tmp_path_factory):
    cache = os.environ.get("P4_TUBE_CACHE")
    if cache and os.path.exists(cache):
        return cache
    path = str(tmp_path_factory.mktemp("tube") / "tube_full.json")
    with precision(256):
        lintube.Tube.build(V0_EC, W_V0, x_c=-3.0).save(path)
    return path


@SLOW
def test_winding_number_small_rectangle_is_one(tube_path):
    """[2.5, 3.1] x [-0.3, 0.3] around kappa_1: winding number 1, certified (~1 min with 4 workers)."""
    res = modecount.winding_number(tube_path, rect=(2.5, 3.1, -0.3, 0.3), workers=4, chunk=0.3, verbose=False)
    assert res["winding"] == 1 and res["certified"] and res["winding_ball"][1] < 1e-2
    for corner, direction, length in modecount.sides_of(res["rect"]):
        segs = sorted((s["s0"], s["s1"]) for s in res["segments"]
                      if abs(complex(*s["kc"]) - (complex(*corner) + complex(*direction) * (s["s0"] + s["w"]))) < 1e-12)
        assert segs[0][0] == 0.0 and abs(segs[-1][1] - length) < 1e-12
        assert all(abs(segs[i][1] - segs[i + 1][0]) < 1e-12 for i in range(len(segs) - 1))


@SLOW
def test_krawczyk_kappa1_and_gauge(tube_path):
    """Krawczyk on kappa_1 +/- 1e-7 and on the gauge value +/- 1e-5: contraction, unique simple zeros,
    gamma = 1/kappa_1, the gauge box contains the certified kappa-bar (~100 s)."""
    ctx = modecount.Context(lintube.Tube.load(tube_path))
    r = krawczyk_kappa.krawczyk(ctx, arb(KAPPA1), 2e-7, m=3)
    assert r["ok"] and abs(r["dE"]) > 0.02 and r["zero"].contains(acb("2.8105525488")) and rad(r["zero"]) < 1e-8
    gamma = 1 / r["zero"]
    assert gamma.real.contains(arb("0.355801922")) and gamma.imag.contains(arb(0))
    g = krawczyk_kappa.krawczyk(ctx, arb(KGAUGE), 2e-5, m=3)                  # the gauge zero (~20 s)
    kb = krawczyk_kappa.gauge_value(ctx.bg)
    assert g["ok"] and abs(g["dE"]) > 0.06 and g["zero"].contains(acb(kb)) and rad(g["zero"]) < 1e-8
    assert not r["zero"].overlaps(g["zero"])


def test_full_rectangle_result_if_available():
    """Verifies the JSON of the full-R run (P4_STAGE3_WINDING): winding 2, certified, contour tiled."""
    path = os.environ.get("P4_STAGE3_WINDING")
    if not path or not os.path.exists(path):
        pytest.skip("no full-rectangle result file")
    res = json.load(open(path))
    assert res["winding"] == 2 and res["certified"] and res["min_abs_E"] > 0
    with precision(256):
        total = sum((modecount.dec(s["inc"]) for s in res["segments"]), arb(0))
        assert (total / (2 * arb.pi()) - 2).abs_upper() < 0.5


# ---------------------------------------------------------------------------------------------
# constraint-surface invariance (exact) and its exponent conditions; stored full-R results
# ---------------------------------------------------------------------------------------------
RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "theorem_b")


def test_linearised_constraint_propagated_exactly_and_exponents(bg):
    """Exact fmpq_mpoly identity 16 S D c' = Lambda_lin c for the linearised momentum constraint along
    the 4D linearised flow (remainder 0, Lambda_lin kappa-free, 19 terms, D = Delta~/(4SW)); the sonic
    exponent gamma contains 0 (no positive integer), D_1 != 0, the centre exponent is -3, and e_W is not
    in Sigma(x_d): p~, r_1, r_2 lie on the constraint surface for every kappa (notes section 3.8)."""
    from p4.validated import linsys
    lam, D = linsys.linear_constraint_propagation()
    assert lam.degrees()[1] == 0 and lam.degrees()[2] == 3 and D.degrees()[2] == 2
    ce = modecount.certified_centre()
    r = analyticity.constraint_exponents(bg, ce, -3.0)
    assert r["ok"] and r["sonic_ok"] and r["centre_ok"] and r["eW_ok"]
    assert r["gamma_sonic"].contains(arb(0)) and r["gamma_sonic"].rad() < 1e-12 and r["D1"] > 7
    assert r["rho_centre"].contains(arb(-3)) and r["rho_centre"].rad() < 1e-12
    assert r["D0"].contains(arb(0)) and abs(r["lSigma_eW"]) > 0.5


def test_context_refuses_mismatched_tube():
    tube = lintube.Tube.build(V0_EC, W_V0, x_c=-0.06)
    tube.info["V0"] = "0.1124394013880984"
    with pytest.raises(ValueError):
        modecount.Context(tube)


def _sides_tiled(res):
    for corner, direction, length in modecount.sides_of(res["rect"]):
        segs = sorted((s["s0"], s["s1"]) for s in res["segments"]
                      if abs(complex(*s["kc"]) - (complex(*corner) + complex(*direction) * (s["s0"] + s["w"]))) < 1e-12)
        assert segs[0][0] == 0.0 and abs(segs[-1][1] - length) < 1e-12
        assert all(abs(segs[i][1] - segs[i + 1][0]) < 1e-12 for i in range(len(segs) - 1))


def test_full_rectangle_result_if_available():
    """Re-checks the stored JSON of the full-R run (results/theorem_b/winding_R.json, or P4_STAGE3_WINDING):
    winding 2 from the exactly encoded tile increments, contour tiled without gaps, constants as certified."""
    path = os.environ.get("P4_STAGE3_WINDING") or os.path.join(RESULTS, "winding_R.json")
    if not os.path.exists(path):
        pytest.skip("no full-rectangle result file")
    res = json.load(open(path))
    assert res["winding"] == 2 and res["certified"] and res["min_abs_E"] > 0 and res["n_tiles"] == len(res["segments"])
    assert res["constants"]["V0"] == V0_EC and res["tube_info"]["V0"] == V0_EC and float(res["tube_info"]["w"]) == W_V0
    assert all(s["minE"] > 0 and s["w"] <= res["wmax"] for s in res["segments"])
    _sides_tiled(res)
    with precision(256):
        total = sum((modecount.dec(s["inc"]) for s in res["segments"]), arb(0))
        assert (total / (2 * arb.pi()) - 2).abs_upper() < 0.5


def test_stored_analyticity_cover():
    """Re-checks the stored cover (results/theorem_b/analyticity_R.json, or P4_STAGE3_ANALYTICITY): the closed
    tiles partition R, each certified on a strictly larger box, nu > |x0|, Re sigma < 1, Z < 1, g < K+1, and the
    constraint-exponent certificate holds."""
    path = os.environ.get("P4_STAGE3_ANALYTICITY") or os.path.join(RESULTS, "analyticity_R.json")
    if not os.path.exists(path):
        pytest.skip("no analyticity result file")
    res = json.load(open(path))
    a, b, c, d = res["rect"]
    bx = res["boxes"]
    assert all(b_["w_cert"] > b_["w"] and a <= b_["cx"] - b_["w"] and b_["cx"] + b_["w"] <= b
               and c <= b_["cy"] - b_["w"] and b_["cy"] + b_["w"] <= d for b_ in bx)
    assert sum((2 * b_["w"]) ** 2 for b_ in bx) == (b - a) * (d - c)
    for i, p in enumerate(bx):
        for q in bx[i + 1:]:
            assert abs(p["cx"] - q["cx"]) >= p["w"] + q["w"] or abs(p["cy"] - q["cy"]) >= p["w"] + q["w"]
    assert all(b_["sonic"]["nu"] > 0.05 and b_["sonic"]["re_sigma_max"] < 1 and b_["sonic"]["Z"] < 1
               and b_["centre"]["Z"] < 1 and b_["centre"]["g"] < 51 for b_ in bx)
    assert res["tube_regular"] and res["constraint"]["ok"] and res["tube_info"]["V0"] == V0_EC
