from fractions import Fraction

import pytest
import sympy as sp
from p8_variable_microlocal import (
    bounds,
    domain,
    energy,
    independent,
    rational,
    symbol,
    weighted,
)


def test_exact_quadratic_weighted_completion():
    assert set(weighted.checks().values()) == {0}


def test_generic_not_sampled_principal_symbol():
    d = symbol.checks()
    assert not d["determinant_leading_kinetic_identity"]
    assert d["det_C_degree"] == 0
    assert len(d["groups"]) == 5


def test_continuous_principal_matrix_bounds():
    d = domain.checks()
    assert [row["count"] for row in d.values()] == [87, 190, 112, 222]
    assert sum(row["count"] for row in d.values()) == 611
    assert all(row["minimum"] > 0 for row in d.values())


def test_constructive_uniform_cauchy_threshold():
    d = bounds.derive()
    assert d["rounded_coefficient_bounds"] == (16, 14, 1)
    assert d["relative_determinant_error_upper"] == Fraction(6944, 20625)
    assert bounds.threshold(Fraction(1, 100)) == 20_000_000_000
    assert bounds.threshold(2) == 500_000


def test_actual_normalization_and_uniform_energy_identity():
    assert energy.checks()["energy_skew_identity_count"] == 36


@pytest.mark.parametrize("value", [True, False, 0.01, 0, -1, 3])
def test_threshold_input_domain(value):
    with pytest.raises((TypeError, ValueError)):
        bounds.threshold(value)


def test_native_exact_arithmetic_against_sympy():
    u, c, K = sp.symbols("u c K", real=True)
    expressions = [
        (u+1)**2/(u+1), (u-c)/(K+u**2),
        (K**4*(u+c)+K**2)/(K**2+1),
        (u**2-c**2)/(u+c)+(K-1)/(K**2-1),
    ]
    for expression in expressions:
        native = rational.from_sympy(expression, (u, c, K))
        assert sp.cancel(native.sympy()-expression) == 0
        assert sp.cancel(native.derivative().sympy()-sp.diff(expression, u)) == 0
        num, den = sp.fraction(sp.cancel(expression))
        assert native.degree() == sp.degree(num, K)-sp.degree(den, K)
        assert sp.cancel(native.leading().sympy()-sp.LC(sp.Poly(num, K))/sp.LC(sp.Poly(den, K))) == 0


def test_native_matrix_inverse_against_sympy():
    u, c, K = sp.symbols("u c K", real=True)
    matrix = sp.Matrix([[u+1, c, 0], [K, u-c, 1], [1, 0, K+1]])
    native = [[rational.from_sympy(x, (u, c, K)) for x in row] for row in matrix.tolist()]
    inv = rational.inverse(native)
    recovered = sp.Matrix([[x.sympy() for x in row] for row in inv])
    assert (recovered*matrix-sp.eye(3)).applyfunc(sp.cancel) == sp.zeros(3)
    assert sp.cancel(rational.determinant(native).sympy()-matrix.det()) == 0


def test_independent_full_fraction_cauchy_bridge():
    d = symbol.derive()
    for f in independent.fixtures():
        at = (f["u"], f["c"], f["K"])
        for key in ("C", "inverse_C", "Mq", "Mv", "symplectic"):
            assert [[x.evaluate(at) for x in row] for row in d[key]] == f[key]
        assert d["det_C"].evaluate(at) == f["det_C"]


def test_independent_supplied_polynomial_bernstein_conversion():
    d = domain.derive()
    p = domain.leading.derive()
    for record in d["records"].values():
        terms = sp.Poly(record["numerator"], p["u"], p["c"]).as_dict()
        degrees, values = independent.bernstein_from_time_lapse_terms(
            {power: Fraction(coefficient) for power, coefficient in terms.items()})
        assert degrees == record["degrees"]
        assert values == record["coefficients"]
