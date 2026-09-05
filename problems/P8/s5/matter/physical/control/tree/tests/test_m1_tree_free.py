import sympy as sp
from p8_m1_tree import free


def test_two_species_initial_covariance_and_exact_canonical_flow():
    residuals = free.checks()
    assert len(residuals) == 56
    assert set(residuals.values()) == {sp.Integer(0)}


def test_connection_and_scalar_channel_omissions_are_detected():
    for value in free.negative_controls().values():
        entries = tuple(value) if isinstance(value, (list, sp.MatrixBase)) else (value,)
        assert any(sp.expand(entry) != 0 for entry in entries)


def test_independent_exponential_enclosure_has_slack():
    data = free.exponential_enclosure()
    assert 1 < sp.Rational(data["lower"]) < sp.Rational(data["upper"]) < sp.Rational(11, 7)
