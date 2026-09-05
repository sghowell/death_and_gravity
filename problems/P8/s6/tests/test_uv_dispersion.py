import pytest
import sympy as sp
from p8_uv.dispersion import necessary_bound


@pytest.mark.parametrize("error,gravity", [(0, None), (None, 0), (None, None)])
def test_unknown_bounds_never_default_to_zero(error, gravity):
    assert necessary_bound(100, error=error, gravity_lower_radius=gravity)["status"] == "UNTESTED_MISSING_BOUND"


@pytest.mark.parametrize("b2,error,gravity,expected", [
    (-1, 0, 2, "NOT_EXCLUDED_BY_THIS_NECESSARY_BOUND"),
    (-3, 0, 2, "EXCLUDED_BY_SUPPLIED_CONDITIONAL_BOUND"),
    (0, 0, 0, "BOUNDARY_OR_UNRESOLVED"),
    (-1, 1, 1, "BOUNDARY_OR_UNRESOLVED"),
    (1, 2, 0, "BOUNDARY_OR_UNRESOLVED"),
    (sp.Rational(1, 3), sp.Rational(1, 6), 0, "NOT_EXCLUDED_BY_THIS_NECESSARY_BOUND")])
def test_supplied_interval_bounds(b2, error, gravity, expected):
    assert necessary_bound(b2, error=error, gravity_lower_radius=gravity)["status"] == expected


@pytest.mark.parametrize("b2,error,gravity", [(1, -1, 0), (1, 0, -1), (0.1, 0, 0),
                                            (1, sp.oo, 0), (sp.Symbol("c"), 0, 0)])
def test_inexact_or_invalid_bounds_rejected(b2, error, gravity):
    with pytest.raises(ValueError):
        necessary_bound(b2, error=error, gravity_lower_radius=gravity)
