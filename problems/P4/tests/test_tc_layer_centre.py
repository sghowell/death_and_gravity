"""S4-6 tests: the centre variant of the layer engine (tc_layer_centre).

* graded values (GV): products, folding and the inverse against exact
  arithmetic;
* pivot-form identities on synthetic Laurent polynomials: the graded small-
  and large-z forms equal P/(y^m z^p) and P/(y^m z^d) at points, in both
  regimes;
* the normalization identity: regimes L (Sigma = S/s^2) and S (Sigma = S)
  give the same Sigma_c = S/(1 + s^2), and it agrees with the direct chain
  (tc_layer.point_S_direct on point pairs) at interior points of c29 and c00;
* the negative grading exponents cancel exactly on boxes, and box enclosures
  contain point values.

The two piece fixtures build the centre pairs of c29 and c00 (about 50 s
each).
"""
import random

import pytest
from flint import acb, arb

import p4.validated.tc_commutator as tc
import p4.validated.tc_energy as te
from p4.validated import tc_layer as tl
from p4.validated import tc_layer_centre as tcc
from p4.validated.arbseries import precision
from p4.validated.ec_centre import certified_centre_expansion
from p4.validated.tc_layer_centre import GV, Grade
from p4.validated.tc_tm import LC, TM

INF = float("inf")


# ------------------------------------------------------------- graded values
def _close(a, b, tol):
    return abs(a - b).upper() <= tol * (1 + abs(b).upper())


def test_gv_product_and_value():
    with precision(256):
        g = arb("0.01")
        x = GV({0: acb(1), 1: acb(2, 1), -1: acb(3)}, g)
        y = GV({0: acb(-2), 2: acb(0, 5)}, g)
        p = (x * y).value(allow_negative=True)
        xv = acb(1) + 2 * g * acb(1) + g * acb(0, 1) + 3 / g
        yv = acb(-2) + g * g * acb(0, 5)
        assert _close(p, xv * yv, 1e-60)


def test_gv_folding_contains_high_powers():
    with precision(256):
        gh = arb(0).union(arb("0.05"))
        x = GV({1: acb(1)}, gh)
        p = x ** (tcc.KMAX + 4)
        assert set(p.c) == {tcc.KMAX + 1}
        for g in ("0", "0.01", "0.05"):
            gv = arb(g)
            assert (p.c[tcc.KMAX + 1] * gv ** (tcc.KMAX + 1)).contains(gv ** (tcc.KMAX + 4))


def test_gv_inverse():
    with precision(256):
        for g in ("0.003", "0.02", "0.07"):
            gh = arb(g)
            q = GV({0: acb(2, 1), 1: acb(-3), 2: acb(0, 4), 3: acb(1)}, gh)
            inv = q.inv().value()
            qv = q.value()
            assert (inv * qv - 1).contains(0) and abs(inv * qv - 1).upper() < 1e-50
        with pytest.raises(ArithmeticError):
            GV({0: acb(0, 1e-3).union(acb(0)), 1: acb(1)}, arb("0.01")).inv()


def test_gv_negative_exponent_needs_exact_zero():
    with precision(128):
        g = arb("0.01")
        assert GV({-1: acb(0), 0: acb(2)}, g).value().contains(2)
        with pytest.raises(ArithmeticError):
            GV({-1: acb(arb(0, 1e-30)), 0: acb(2)}, g).value()


# -------------------------------------------- pivot forms (synthetic Laurent)
def _kp(coefs, offs, R=1.0):
    return [LC(TM([arb(c) for c in cs], R), o) for cs, o in zip(coefs, offs)]


def _P(kp, y, k):
    acc = acb(0)
    for x in reversed(kp):
        acc = acc * k + x.tm.eval(y, y) * arb(y) ** x.off
    return acc


@pytest.mark.parametrize("offs", [(10, 8, 8, 8, 10, 10), (9, 7, 7, 7, 7, 7), (4, 2, 3, 5, 6, 7)])
@pytest.mark.parametrize("regime", ["L", "S"])
def test_pivot_forms_graded(offs, regime):
    coefs = [("1.5", "-0.25", "3"), ("-2", "0.5"), ("0.75", "1", "-1"), ("1.25", "0", "2"),
             ("-0.5", "0.125"), ("2", "-1", "0.5")]
    kp = _kp(coefs, offs)
    F = tl.PForm(kp, "syn")
    with precision(256):
        for y, eps, tau in ((0.004, 0.0, 60.0), (0.01, 2.5, 30.0), (3e-5, 6.0, 4e4)):
            s = y * tau
            if regime == "S" and s > tcc.S_SPLIT:
                continue
            if regime == "L" and s < 1e-3:
                continue
            box = tcc.point_box(y, eps, tau)
            g = tl.box_geometry(_Dummy(), box)
            gd = Grade(regime, g["L"], g["E"], g["s"], g["tau"])
            L = g["L"]
            k = acb(arb(eps), arb(tau))
            z = arb(y) * k
            P = _P(kp, y, k)
            yb = arb(y)
            # the box's s-range is the point s to 1e-15 (relative)
            small = tcc._form_small(F, L, gd, {}).value()
            ref = P / (yb ** F.m * z ** F.p)
            assert small.overlaps(ref) and _close(small, ref, 1e-12)
            if regime == "L":
                large = tcc._form_large(F, L, gd, {}).value()
                ref = P / (yb ** F.m * z ** F.d)
                assert large.overlaps(ref) and _close(large, ref, 1e-12)


class _Dummy:
    """Stands in for a Laurent LayerData in tl.box_geometry (xi = ell = y)."""
    laurent = True
    xik = None

    @staticmethod
    def xi_of_ell(x):
        return x

    @staticmethod
    def ell_of_xi(x):
        return x

    @staticmethod
    def xi_hull(L):
        return tl.iv(L[0], L[1])


# --------------------------------------------------------- the centre pieces
@pytest.fixture(scope="module")
def atoms_ce():
    return te.Atoms(full=True), certified_centre_expansion(256, K=30)


@pytest.fixture(scope="module")
def c29(atoms_ce):
    at, ce = atoms_ce
    lo, hi = tc.y_subranges()[-1]
    return lo, hi, tcc.centre_layer_data(at, ce, hi, label="c29")


@pytest.fixture(scope="module")
def c00(atoms_ce):
    at, ce = atoms_ce
    lo, hi = tc.y_subranges()[0]
    return lo, hi, tcc.centre_layer_data(at, ce, hi, label="c00")


def _sigma(ld, y, eps, tau, regime):
    box = tcc.point_box(y, eps, tau)
    tot, gd = tcc.sigma_graded(ld, box, regime=regime)
    return tot.value().real * gd.norm


def test_regimes_agree_and_normalization(c29):
    lo, hi, ld = c29
    with precision(320):
        for y in (lo, (lo * hi) ** 0.5, hi):
            for eps in (0.0, 1.0, 6.0):
                for s in (0.12, 0.3, 0.45):
                    tau = s / y
                    a = _sigma(ld, y, eps, tau, "L")
                    b = _sigma(ld, y, eps, tau, "S")
                    assert a.overlaps(b)
                    assert abs(a.mid() - b.mid()) <= max(a.rad(), b.rad()) + 1e-12


def test_agrees_with_direct_chain_c29(c29, atoms_ce):
    at, ce = atoms_ce
    lo, hi, ld = c29
    with precision(320):
        for y in (0.0075, 0.0105):
            pr = tcc.point_pairs_centre(at, ce, y)
            for eps in (0.0, 0.5, 4.0):
                for tau in (14.0, 60.0, 224.0, 1e4, 1e9):
                    d, S = tcc.point_sigma_direct(pr, y, eps, tau)
                    v = tcc.sigma_box_centre(ld, tcc.point_box(y, eps, tau))
                    assert v.overlaps(d)
                    assert abs(v.mid() - d.mid()) < 1e-6 * abs(d.mid()) + 1e-12


def test_agrees_with_direct_chain_c00(c00, atoms_ce):
    """Interior points of the closing piece where the direct chain's ball is
    tight (its cancellation grows like y^-10: tau large enough)."""
    at, ce = atoms_ce
    lo, hi, ld = c00
    with precision(320):
        y = 1.2e-8
        pr = tcc.point_pairs_centre(at, ce, y)
        n = 0
        for eps in (0.0, 1.0):
            for tau in (1e8, 1e9, 1e11):
                d, S = tcc.point_sigma_direct(pr, y, eps, tau)
                if not float(d.rad()) < 1e-3 * abs(float(d.mid())):
                    continue
                v = tcc.sigma_box_centre(ld, tcc.point_box(y, eps, tau))
                assert v.overlaps(d)
                assert abs(v.mid() - d.mid()) < 1e-3 * abs(d.mid())
                n += 1
        assert n >= 4


def test_negative_exponents_cancel_exactly_on_boxes(c29):
    lo, hi, ld = c29
    with precision(320):
        for box in ((lo, hi, 0.0, 0.03, 10.0, 20.0, 224.0, INF),
                    (0.009, 0.0092, 3.4, 4.6, 0.15, 0.16, 14.0, INF),
                    (lo, hi, 1.0, 1.6, 1e3, INF, 224.0, INF)):
            tot, gd = tcc.sigma_graded(ld, box)
            for e, v in tot.c.items():
                if e < 0:
                    assert v.is_zero(), (box, e)


def test_box_enclosure_contains_points(c29):
    lo, hi, ld = c29
    rnd = random.Random(7)
    with precision(320):
        for box in ((0.0091, 0.0093, 0.0, 0.03, 0.4, 0.42, 14.0, INF),
                    (0.0091, 0.0093, 3.4, 3.7, 0.15, 0.155, 14.0, INF),
                    (lo, hi, 1.0, 1.6, 10.0, 20.0, 224.0, INF),
                    (lo, hi, 0.0, 0.03, 1e3, INF, 224.0, INF)):
            v = tcc.sigma_box_centre(ld, box)
            g = tl.box_geometry(ld, box)
            (ya, yb), (sa, sb) = g["L"], g["s"]
            for _ in range(6):
                y = rnd.uniform(ya, yb)
                e = rnd.uniform(box[2], box[3])
                s = rnd.uniform(sa, min(sb, 1e6))
                tau = s / y
                if not g["tau"][0] <= tau <= g["tau"][1]:
                    continue
                pv = tcc.sigma_box_centre(ld, tcc.point_box(y, e, tau))
                assert v.contains(pv.mid())


def test_sigma_negative_with_margin_on_grid(c29):
    lo, hi, ld = c29
    with precision(320):
        worst = max(float(_sigma(ld, y, e, s / y, None).mid())
                    for y in (lo, hi) for e in (0.0, 3.0, 6.0)
                    for s in (14 * lo, 0.5, 2.0, 50.0, 1e5) if s / y >= 14.0)
        assert worst < -6.0
