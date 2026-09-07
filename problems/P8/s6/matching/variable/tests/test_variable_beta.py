from fractions import Fraction

import pytest
import sympy as sp
from p8_variable_beta import background, canonical, independent, inner, model, tensors


@pytest.mark.parametrize("checks", (model.checks, model.source_checks, background.checks,
                                    tensors.checks, tensors.spatial_curvature_check, tensors.vector_checks,
                                    canonical.checks, canonical.action_checks, canonical.center_checks,
                                    inner.checks))
def test_exact_symbolic_source_background_and_tensor_identities(checks):
    assert all(sp.simplify(value) == 0 for value in checks().values())


def test_independent_fraction_coefficients_and_true_second_jets():
    report = independent.checks()
    assert len(report["coefficientwise_identities"]) == 6
    assert len(report["actual_background_Fraction_fixtures"]) == 7
    assert len(report["canonical_second_jet_Fraction_fixtures"]) == 5
    x = independent.Jet(Fraction(2), Fraction(1))
    assert (x/x).value == 1
    assert (x/x).first == (x/x).second_coefficient == 0
    assert (1/x).second_coefficient == Fraction(1, 8)


def test_strict_uniform_clock_box_and_opposite_vector_controls():
    assert all(value.is_positive for value in background.domain_checks().values())
    assert background.evaluate(0, 1)["kbar"] == sp.Rational(5599, 100)
    assert background.evaluate(0, 4)["kbar"] == sp.Rational(799, 100)
    assert tensors.controls()["c1_center_shift_coefficient"] == -sp.Rational(128, 3)
    assert tensors.controls()["c4_center_shift_coefficient"] == sp.Rational(8, 3)


@pytest.mark.parametrize("bad", (True, 0.1, sp.Float("0.1"), sp.oo, 0, 2, 5))
def test_invalid_lapse_or_inexact_domain_is_not_silently_admitted(bad):
    with pytest.raises((ValueError, TypeError)):
        background.evaluate(0, bad)


def test_no_fixed_physical_window_or_gap_is_inferred_from_inner_omission_control():
    correction = inner.canonical_adiabatic_checks()
    assert correction["correction_at_center_limit"] == canonical.derive()["kbar"]**2-sp.Rational(18, 5)
    control = inner.controls()
    assert control["omitted_derivative_residual"] == sp.Rational(1, 5)
    assert control["particular_to_algebraic_ratio"] == sp.Rational(5, 6)
    assert "not a specified physical-metric source" in control["source_scope"]
