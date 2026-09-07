from fractions import Fraction

import pytest
import sympy as sp
from p8_a26_vacuum import coefficients, independent


def test_complete_coefficient_and_numerator_identities():
    assert coefficients.checks() and not any(coefficients.checks().values())


@pytest.mark.parametrize("a,g", [(Fraction(1, 7), Fraction(1, 5)), (3, Fraction(2, 3)),
                                 (Fraction(2, 5), Fraction(1, 5)), (-2, 0)])
def test_fraction_laurent_expansion_matches_symbolic_residue(a, g):
    actual = coefficients.residues(a, g)
    for printed in (False, True):
        assert independent.residue_fixture(a, g, printed=printed) == tuple(actual[k] for k in ("A3", "A4", "A5"))


def test_pole_cancellation_is_not_a_claim_about_the_next_order():
    actual = coefficients.residues(sp.Rational(2, 5), sp.Rational(1, 5))
    assert not actual["generic_obstruction"]
    assert tuple(actual[key] for key in ("A3", "A4", "A5")) == (0, 0, 0)


def test_tensor_sign_not_confused_with_nonzero_denominator():
    assert coefficients.residues(3, sp.Rational(2, 3))["F2_at_zero"] == -sp.Rational(1, 6)
    assert not coefficients.residues(3, sp.Rational(2, 3))["healthy_tensor_sign"]


def test_formal_stationary_equations_are_separate_from_derivative_chart():
    result = coefficients.stationary_conditions(12, 3, 2, 20, 5)
    assert result == {"zero_stress": 0, "phi": 0, "chi": 6}


def test_source_clock_and_vacuum_are_different_points():
    x = sp.Symbol("X")
    data = coefficients.literal(x, sp.Rational(1, 3), sp.Rational(1, 5))
    assert data["F2"].subs(x, -1) == sp.Rational(1, 2)
    assert data["A1"].subs(x, -1) == 0
    assert data["F2"].subs(x, 0) == sp.Rational(3, 10)


def test_independent_series_rejects_degenerate_leading_denominator():
    with pytest.raises(ValueError):
        independent.residue_fixture(1, Fraction(1, 2), printed=True)


def test_independent_local_numerator_does_not_divide_by_hubble_at_bounce():
    actual = independent.local_background()
    assert actual["Zprime"] == -Fraction(1, 20)
    assert actual["r_over_g1"] == -Fraction(3, 52)


def test_symbolic_parameters_use_literal_not_the_exact_numeric_verdict_api():
    with pytest.raises(ValueError):
        coefficients.residues(sp.Symbol("a", real=True), sp.Rational(1, 5))
