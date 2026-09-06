from fractions import Fraction

import pytest
import sympy as sp
from p8_composite import background, bounce, branches, reconstruction


@pytest.mark.parametrize("check", [background.variation_checks, background.bianchi_checks,
                                   branches.checks, bounce.checks, reconstruction.checks,
                                   reconstruction.scale_checks])
def test_exact_background_and_reconstruction_identities(check):
    assert all(sp.simplify(value) == 0 for value in check().values())


def test_explicit_physical_metric_is_composite_not_g():
    data = background.equations()
    assert data["N_eff"] == background.ALPHA*background.Ng+background.BETA*background.Nf
    assert sp.factor(data["a_eff"]-background.a) != 0
    assert sp.factor(data["Q"]-data["P"]) != 0


def test_dynamic_and_intersection_coefficients_are_positive():
    assert all(value.is_positive is True for value in branches.positive_factors().values())


def test_free_bounce_has_positive_lapses_scalar_and_nonzero_branch_error():
    data = bounce.derive()
    y = data["y"]
    assert data["D"].is_positive is True
    assert data["Ng"].is_positive is True
    assert data["Nf"].is_positive is True
    assert data["varphi_prime"].is_positive is True
    assert sp.simplify((data["Y"]-data["X"]/y).subs(y, 1)) != 0
    assert all(sp.sympify(value) != 0 for value in bounce.negative_controls().values())


def test_double_pressure_root_is_not_a_zero_time_velocity():
    data = bounce.derive()
    y = data["y"]
    assert sp.diff(data["p"], y).subs(y, 1) == 0
    assert sp.diff(data["p"], y, 2).subs(y, 1) == -sp.Rational(1, 4)
    assert data["yprime"].subs(y, 1) != 0
    assert sp.simplify((sp.diff(data["H_eff"], y)*data["yprime"]).subs(y, 1)) == sp.Rational(7, 36)


def test_free_history_is_not_CD_even_after_positive_time_rescaling():
    u, k = sp.symbols("u k", real=True, positive=True)
    cd = (1+(u/k)**2)**2
    cd_invariant = sp.diff(cd, u, 4).subs(u, 0)-sp.Rational(3, 2)*sp.diff(cd, u, 2).subs(u, 0)**2
    assert sp.simplify(cd_invariant) == 0
    jets = bounce.time_jets()
    assert jets["a_fourth"]-sp.Rational(3, 2)*jets["a_second"]**2 == -sp.Rational(9, 32)
    assert jets["H_third"]+sp.Rational(3, 2)*jets["H_prime"]**2 == -sp.Rational(9, 32)


def test_finite_box_continuation_has_strict_rational_margins():
    box = reconstruction.positivity_box()
    assert all(value > 0 for value in box["strict_positive_margins"].values())
    assert box["ratio_speed_upper"] < 9
    assert box["density_speed_upper"] < sp.Rational(1, 4)
    assert box["null_lower"] > sp.Rational(1, 2)


def test_prescribed_potential_is_not_silently_the_free_action():
    data = reconstruction.derive()
    assert sp.simplify(data["potential_value"]) != 0
    assert sp.diff(data["potential_value"], data["rho"]) == sp.Rational(1, 2)


def test_unproved_larger_window_fails_the_sufficient_box_gate():
    with pytest.raises(ValueError, match="margin"):
        reconstruction.positivity_box(Fraction(1, 8))
    assert reconstruction.positivity_box(Fraction(1, 128))["delta_u"] == Fraction(1, 128)
    # Failure is of this sufficient enclosure, not a claim of physical instability.


@pytest.mark.parametrize("delta", [True, 0.01, float("inf"), float("nan")])
def test_inexact_window_radius_is_not_a_certified_input(delta):
    with pytest.raises(TypeError, match="exact"):
        reconstruction.positivity_box(delta)


@pytest.mark.parametrize("delta", [0, -1])
def test_nonpositive_window_radius_is_rejected(delta):
    with pytest.raises(ValueError, match="positive"):
        reconstruction.positivity_box(delta)
