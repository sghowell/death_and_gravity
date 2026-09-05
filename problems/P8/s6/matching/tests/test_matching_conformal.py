import pytest
import sympy as sp
from p8_matching import conformal as c


def test_conformal_identities_and_threshold():
    result = c.derive()
    assert set(result["exact_residuals"].values()) == {"0"}
    assert result["common_C2_error_necessary_lower_bound"] == "4/7"


@pytest.mark.parametrize("error", [0, sp.Rational(1, 10), sp.Rational(1, 2), sp.Rational(57, 100)])
def test_subthreshold_common_errors_are_excluded(error):
    assert c.error_bound(error, error, error)["status"] == "EXCLUDED_BY_CONFORMAL_NEC"


def test_threshold_is_not_a_sufficient_construction():
    result = c.error_bound(*([sp.Rational(4, 7)]*3))
    assert result["lower_margin"] == "0"
    assert result["status"] == "NEC_LOWER_BOUND_INCONCLUSIVE"


def test_missing_time_derivative_terms_would_change_answer():
    q0 = sp.Symbol("q0", positive=True)
    q1, q2 = sp.symbols("q1 q2", real=True)
    numerator = c.identities()["necessary_bounce_numerator"]
    jets = {q0: sp.Rational(2, 5), q1: sp.Rational(3, 5), q2: -sp.Rational(3, 5)}
    assert numerator.subs(jets) == -sp.Rational(7, 20)
    assert (4*q0+q2).subs(jets) == 1


def test_small_amplitude_alone_does_not_exclude():
    assert c.error_bound(sp.Rational(1, 100), 2, sp.Rational(1, 100))["status"] == "NEC_LOWER_BOUND_INCONCLUSIVE"


def test_ten_percent_exact_margin_and_geometry_tolerance():
    errors = [sp.Rational(1, 10)]*3
    assert c.error_bound(*errors)["lower_margin"] == "209/60"
    assert c.error_bound(*errors, geometry_error=sp.Rational(1, 10))["lower_margin"] == "991/300"


@pytest.mark.parametrize("errors,eta", [((1, 0, 0), 0), ((-1, 0, 0), 0), ((0, -1, 0), 0),
                                       ((0, 0, -1), 0), ((0, 0, 0), 2), ((0.1, 0, 0), 0)])
def test_invalid_bounds_rejected(errors, eta):
    with pytest.raises(ValueError):
        c.error_bound(*errors, geometry_error=eta)
