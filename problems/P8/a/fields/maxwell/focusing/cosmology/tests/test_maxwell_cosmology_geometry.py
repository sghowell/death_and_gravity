"""Full power interval, short clock and robust C3 history controls."""

import pytest
import sympy as sp
from p8a_maxwell_cosmology import geometry


def test_every_exact_clock_and_jet_identity():
    assert all(sp.simplify(value) == 0 for value in geometry.identities().values())


def test_radiation_and_matter_endpoints_match_known_geometric_jets():
    assert geometry.history_point(sp.Rational(1, 2), 0) == [-2, -8, -64, -768]
    assert geometry.history_point(sp.Rational(2, 3), 0) == [-2, -6, -36, -324]
    assert geometry.history_point(sp.Rational(1, 2), -sp.Rational(1, 100))[0] == -sp.Rational(25, 13)


def test_actual_tube_has_strict_margins_not_just_boundary_comparators():
    data = geometry.neighborhood()
    assert data["strict_cap_margins"] == [sp.Rational(9, 100), sp.Rational(1, 2), 3, 16]
    assert data["strict_contraction_margin"] == sp.Rational(17, 1300)
    assert data["log_scale_factor_error_upper"] == sp.Rational(1, 10000)


def test_anchored_and_unanchored_observer_clocks_are_not_conflated():
    assert geometry.neighborhood()["observed_H0_times_tau_range"] == [sp.Rational(199, 100), sp.Rational(201, 100)]
    assert geometry.neighborhood(anchored=True)["observed_H0_times_tau_range"] == [2, 2]
    with pytest.raises(TypeError):
        geometry.neighborhood(anchored=1)


def test_wrong_history_side_and_out_of_domain_power_are_rejected():
    for p, x in ((sp.Rational(1, 2), sp.Rational(1, 100)), (1, 0), (0, 0)):
        with pytest.raises(ValueError):
            geometry.history_point(p, x)


def test_loose_tube_cannot_keep_false_contraction_or_jet_claims():
    for errors in ((sp.Rational(1, 10), sp.Rational(1, 2), 3, 16),
                   (sp.Rational(1, 100), 2, 3, 16)):
        with pytest.raises(ValueError):
            geometry.neighborhood(errors)


@pytest.mark.parametrize("bad", [True, 0.01, sp.oo, sp.nan, sp.Symbol("p")])
def test_numeric_history_input_is_exact_finite_and_resolved(bad):
    with pytest.raises((TypeError, ValueError)):
        geometry.history_point(bad, 0)
