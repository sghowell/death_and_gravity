"""Frozen-action and continuous-domain connections, without report caching."""
from fractions import Fraction

import sympy as sp
from p8_own_scattering import action, bridges, connection, potential


def test_eighteen_exact_action_profile_endpoint_interfaces():
    assert len(bridges.identities()) == 18
    assert set(bridges.identities().values()) == {0}


def test_conditional_remainder_is_actually_discharged():
    checks = bridges.checks()
    assert len(checks) == 9
    assert all(value is True for value in checks.values())
    assert potential.calibration()["derived_potential_bound"] < connection.REMAINDER_BOUND
    assert potential.POTENTIAL_CAP == connection.REMAINDER_BOUND == 44


def test_precise_smaller_transfer_domain_is_not_the_full_background_box():
    assert connection.HALF_WIDTH == potential.RADIUS == Fraction(1, 100)
    assert 0 < connection.DELTA_MAX == Fraction(1, 10**6) < potential.DELTA_MAX


def test_source_measure_and_physical_normalization_are_explicit():
    canonical, clock = action.canonical(), action.asinh_map()
    assert action.calibration()["asinh_source"] == (
        clock["r"]**sp.Rational(3, 2)*canonical["spring"]*canonical["q"]/sp.sqrt(canonical["k"]))
    assert sp.simplify(clock["endpoint_from_Q_QT"].det()-action.TAU*clock["k"]/connection.RHO) == 0


def test_reference_error_margin_does_not_shrink_with_delta_in_this_proof():
    data = connection.calibration()
    assert data["transfer_error_upper"] == Fraction(4, 399) > 0
    assert data["effective_B_lower"] == Fraction(239, 15960) > Fraction(1, 80)


def test_nonzero_transfer_requires_nonzero_data_or_forcing():
    # A linear homogeneous transfer always maps the zero state to zero.
    aa, bb = sp.symbols("A B", complex=True)
    transfer = sp.Matrix([[sp.conjugate(aa), -sp.conjugate(bb)], [-bb, aa]])
    assert transfer*sp.zeros(2, 1) == sp.zeros(2, 1)
    assert "not a zero-data" in connection.calibration()["data_scope"]
