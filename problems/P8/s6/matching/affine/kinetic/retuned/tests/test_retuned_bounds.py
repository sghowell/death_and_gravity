"""Warm-cache strict inputs and excluded charts are not silently accepted."""
import pytest
import sympy as sp
from p8_affine_retuned import bounds


@pytest.mark.parametrize("bad", (True, False, 1.0, sp.Float(1), "1", sp.I, sp.oo, sp.nan, sp.Symbol("unproved")))
def test_inexact_or_unproved_inputs_rejected(bad):
    bounds.require_domain(sp.Rational(1, 2), sp.Rational(1, 2000), 1, 1)
    with pytest.raises((TypeError, ValueError)):
        bounds.require_domain(sp.Rational(1, 2), bad, 1, 1)


@pytest.mark.parametrize("p,coupling,time,q", ((-sp.Rational(1, 2), sp.Rational(1, 2000), 1, 1),
    (sp.Rational(1, 3), sp.Rational(1, 2000), 1, 1),
    (sp.Rational(1, 2), -sp.Rational(1, 2000), 1, 1),
    (sp.Rational(1, 2), sp.Rational(1, 1000), 1, 1),
    (sp.Rational(1, 2), sp.Rational(1, 2000), 0, 1),
    (sp.Rational(1, 2), sp.Rational(1, 2000), 1, 0),
    (sp.Rational(1, 2), sp.Rational(1, 2000), 1, -1)))
def test_uncertified_domains_rejected(p, coupling, time, q):
    with pytest.raises(ValueError):
        bounds.require_domain(p, coupling, time, q)


def test_zero_curl_auxiliary_limit_and_first_order_center_are_separate():
    assert bounds.require_domain(sp.Rational(1, 2), 0, 0, 0)["branch"] == "exact_auxiliary_control"
    data = bounds.require_domain(sp.Rational(1, 2), sp.Rational(1, 2000), 0, 10, punctured=False)
    assert data["branch"] == "positive_principal_rolling_candidate"
    assert data["nonlinear_or_UV_health_claim"] is False


@pytest.mark.parametrize("values", ((0, 1, 1), (1, 0, 1), (1, 1, 0), (1, 1, -1), (1, 1.0, 1)))
def test_invalid_units_rejected(values):
    with pytest.raises((TypeError, ValueError)):
        bounds.units(*values)
