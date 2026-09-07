from fractions import Fraction

import pytest
import sympy as sp
from p8_variable_vacuum import bounds, independent, model, obstruction


def test_literal_full_metric_and_clock_variations():
    assert len(model.checks()) == 8
    assert set(model.checks().values()) == {0}
    assert model.controls()["same_metric_values_with_nonzero_clock_slope"] == 2


def test_actual_profile_residual_and_clock_chain_identities():
    assert len(obstruction.checks()) == 11
    assert set(obstruction.checks().values()) == {0}
    controls = obstruction.controls()
    assert controls["actual_center_f_and_clock_only"]["g_flat_residual"] == -2
    assert controls["actual_center_g_and_clock_only"]["f_flat_residual"] == -sp.Rational(469, 432)
    assert set(controls["separately_changed_b0_plus2_at_c4_center"].values()) == {0}


def test_continuous_bound_coefficients_and_strict_gaps():
    assert set(bounds.checks().values()) == {0}
    d = bounds.derive()
    assert d["endpoint_margin_Bernstein"] == [0, sp.Rational(177, 50), sp.Rational(3, 25)]
    assert d["constants"]["g_U_upper_c1"] == -sp.Rational(123, 25)
    assert d["constants"]["g_U_upper_positive_branch"] == -sp.Rational(459, 800)
    rebuilt = independent.coefficients()
    assert tuple(map(Fraction, d["endpoint_margin_Bernstein"])) == rebuilt["margin_Bernstein"]


@pytest.mark.parametrize("fixture", independent.fixtures())
def test_independent_fraction_values_and_clock_derivatives(fixture):
    exact = independent.profile(*fixture)
    primary = obstruction.evaluate(*fixture)
    assert set(exact) == set(primary)
    assert exact == {key: Fraction(value) for key, value in primary.items()}
    c = fixture[1]
    assert exact["b1"]*exact["b4"] < 0
    if c == 1:
        assert 0 < exact["r_star_cubed"] <= exact["y"]**3
        assert exact["metric_eliminant"] < 0
    else:
        assert exact["y"]**3 <= exact["r_star_cubed"] < Fraction(9, 4)**3
        assert exact["metric_eliminant"] > 0


@pytest.mark.parametrize("fixture", independent.coframe_fixtures())
def test_independent_literal_fraction_coframe_variations(fixture):
    ag, ng, ratio, betas, slopes = fixture
    exact = independent.coframe_fixture(*fixture)
    d = model.derive()
    derivatives = {sp.diff(beta, d["phi"]): sp.Rational(slope) for beta, slope in zip(d["betas"], slopes)}
    at = {d["a_g"]: sp.Rational(ag), d["N_g"]: sp.Rational(ng), d["r"]: sp.Rational(ratio)}
    at.update({beta: sp.Rational(value) for beta, value in zip(d["betas"], betas)})
    assert exact == {key: Fraction(d[key].subs(derivatives).subs(at)) for key in exact}


@pytest.mark.parametrize("args,error", [
    ((0.0, 4), TypeError), ((True, 4), TypeError), ((0, 4.0), TypeError),
    ((sp.Rational(101, 1000), 4), ValueError), ((-sp.Rational(101, 1000), 4), ValueError),
    ((0, 2), ValueError), ((0, sp.Rational(3, 2)), ValueError), ((0, 5), ValueError),
    ((0, 4, 0), ValueError), ((0, 4, -1), ValueError), ((0, 4, 2.0), TypeError),
    ((0, 4, True), TypeError), ((sp.oo, 4), TypeError),
])
def test_domain_validation_before_any_cached_profile(args, error):
    obstruction.evaluate(0, 4, 2)
    with pytest.raises(error):
        obstruction.evaluate(*args)


def test_action_parameter_and_vacuum_ratio_are_independent():
    d = obstruction.evaluate(0, 4, 2)
    forced_ratio = obstruction.evaluate(0, 4, 4)
    assert d["f_flat_residual"] == 0
    assert forced_ratio["f_flat_residual"] != 0
    assert forced_ratio["g_flat_residual"] != 0


def test_zero_vacuum_matter_is_not_the_rolling_profile_stress():
    # nbar labels coefficient functions inherited from reconstruction; it is
    # not an imposed nonzero scalar velocity in the candidate vacuum.
    d = obstruction.evaluate(0, 4, 2)
    assert d["nbar"] == 8 and d["kbar"] == sp.Rational(799, 100)
    assert d["g_flat_residual"] == -d["nbar"]/4
