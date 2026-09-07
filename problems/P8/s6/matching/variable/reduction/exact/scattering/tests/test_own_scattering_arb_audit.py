"""Independent univariate recurrence, exact domain cover and source audit."""
from fractions import Fraction
from itertools import product

import pytest
import sympy as sp
from flint import ctx, fmpq
from p8_own_scattering import action
from p8_own_scattering import arb_audit as audit


def contains(value, expected):
    expected = Fraction(expected)
    return value.contains(fmpq(expected.numerator, expected.denominator))


def test_independent_arb_arithmetic_at_exact_coefficients():
    with ctx.workprec(256):
        a = audit.Jet((1, 2, 3, 4))
        b = audit.Jet((2, -1, 0, 1))
        assert all(contains(value, target) for value, target in zip((a*b).a, (2, 3, 4, 6), strict=True))
        assert all(contains(value, target) for value, target in zip((a/a).a, (1, 0, 0, 0), strict=True))
        assert all(contains(value, target) for value, target in zip(a.derivative().a, (2, 6, 12), strict=True))
        cubic = a**3
        assert all(contains(value, target) for value, target in zip(cubic.root3().a, (1, 2, 3, 4), strict=True))
        reciprocal = audit.Jet((2, 1, 0, 0)).inverse()
        assert all(contains(value, target) for value, target in zip(reciprocal.a, (Fraction(1, 2), Fraction(-1, 4), Fraction(1, 8), Fraction(-1, 16)), strict=True))


def test_unavailable_taylor_derivatives_are_not_fabricated():
    with pytest.raises(ValueError):
        audit.Jet((1,)).derivative()
    with pytest.raises(ValueError):
        audit.Jet((1, 2)).truncate(2)
    with pytest.raises(ValueError):
        _ = audit.Jet((1, 2))+audit.Jet((1, 2, 3))
    with pytest.raises(ZeroDivisionError):
        audit.Jet((0, 1)).inverse()
    with pytest.raises(ValueError):
        audit.Jet((-1, 1)).root3()
    with pytest.raises(TypeError):
        _ = audit.Jet((1, 2))**True


def test_eight_boxes_cover_every_axis_without_gaps_or_domain_growth():
    boxes = audit.cover_boxes()
    assert len(boxes) == 8
    expected = []
    for name, (lo, hi) in audit.domain().items():
        midpoint = (lo+hi)/2
        expected.append((name, ((lo, midpoint), (midpoint, hi))))
    assert boxes == [dict(zip((item[0] for item in expected), intervals, strict=True))
                     for intervals in product(*(item[1] for item in expected))]


def test_independent_whole_cover_proves_both_sharp_derivative_bounds():
    before = ctx.prec
    report = audit.report()
    assert ctx.prec == before
    assert report["whole_domain_not_sampled_grid"]
    assert report["independent_univariate_implicit_recurrence"]
    assert all(value is True for record in report["boxes"] for value in record["checks"].values())


@pytest.mark.parametrize("precision", (128, 256, 384))
def test_independent_joint_center_implicit_jets(precision):
    box = {"v": (0, 0), "delta": (0, 0), "zeta": (12, 12)}
    data = audit.enclosure(precision, box)
    assert contains(data["zeta"][1], 204)
    assert contains(data["b"][0], 2)
    assert contains(data["b"][1], -8)
    assert contains(data["b"][2], -28)
    assert contains(data["N"][0], 2)
    assert contains(data["N"][1], 8)
    assert contains(data["A"][0], 16)
    assert contains(data["A_v"], -48)
    assert contains(data["ell"], -16)
    assert contains(data["pump"], -16)


@pytest.mark.parametrize("bad", (True, 0.1, "1/100", sp.Rational(1, 100), None))
def test_arb_domain_rejects_inexact_or_unsupported_inputs(bad):
    with pytest.raises(TypeError):
        audit.rational(bad)
    with pytest.raises(TypeError):
        audit.ball(bad, 1)
    with pytest.raises(TypeError):
        audit.enclosure(bad)


def test_arb_domain_cannot_be_enlarged_by_precision_or_tuple_shortcuts():
    for value in (0, 64, -1):
        with pytest.raises(ValueError):
            audit.enclosure(value)
    with pytest.raises(ValueError):
        audit.enclosure(box={"v": (0, 0)})
    box = audit.domain()
    box["v"] = (0, Fraction(1, 1000))
    with pytest.raises(ValueError):
        audit.enclosure(box=box)
    box["v"] = [0, Fraction(1, 10000)]
    with pytest.raises(TypeError):
        audit.enclosure(box=box)


def test_full_canonical_action_source_and_asinh_identities():
    assert len(action.identities()) == 9
    assert set(action.identities().values()) == {0}


def test_canonical_and_asinh_source_weights_are_both_required():
    c, a = action.canonical(), action.asinh_map()
    assert c["canonical_source_weight"] == c["spring"]/sp.sqrt(c["k"])
    assert a["source_weight_after_asinh"] == (action.U**2+action.DELTA/8)**sp.Rational(3, 4)
    assert c["canonical_source_weight"] != c["spring"]/c["k"]
    assert not action.calibration()["canonical_source_bound_or_zero_data_response_proved_here"]


def test_physical_and_dimensionless_endpoint_columns_remain_distinct():
    a = action.asinh_map()
    assert a["endpoint_from_Q_QT"][1, 1] == action.TAU*a["endpoint_from_Q_Qu"][1, 1]
    assert sp.simplify(a["endpoint_from_Q_Qu"].det()) == a["k"]/action.RHO
