"""Original-action, independent-premise and physical-units interfaces."""
import pytest
import sympy as sp
from p8_affine_kinetic import bridges, scalar


def test_all_independent_interfaces():
    for value in bridges.checks().values():
        assert all(entry == 0 for entry in value) if isinstance(value, sp.MatrixBase) else value == 0


def test_all_continuous_domain_premises():
    assert all(value is True for value in bridges.proof_checks().values())


def test_nonunit_values():
    data = bridges.calibration()
    assert data["normalized_zeta"] == sp.Rational(5, 12)
    assert data["isolated_proca_mass_squared"] == sp.Rational(8, 5)
    assert data["normalized_isolated_mass_squared"] == sp.Rational(32, 5)
    assert data["physical_q_threshold_at_center"] == sp.Rational(3597, 12800)
    assert data["threshold_is_squared_momentum_not_frequency_cutoff"] is True


@pytest.mark.parametrize("args", ((0, 2, 5), (3, 0, 5), (3, 2, 0), (True, 2, 5),
                                  (3, 2.0, 5), (3, 2, "5")))
def test_unit_controls(args):
    with pytest.raises((TypeError, ValueError)):
        scalar.units(*args)
