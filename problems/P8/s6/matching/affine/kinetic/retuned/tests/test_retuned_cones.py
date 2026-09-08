"""Exact physical front-cone obstruction, distinct from massive phase speed."""
import pytest
import sympy as sp
from p8_affine_retuned import cones
from p8_affine_retuned import dynamics as d


def test_every_exact_full_characteristic_identity():
    for value in cones.checks().values():
        assert all(item == 0 for item in value) if isinstance(value, sp.MatrixBase) else value == 0


@pytest.mark.parametrize("time", (sp.Rational(1, 100), sp.Rational(1, 2), -sp.Rational(1, 2), 3))
def test_positive_physical_characteristic_above_matter_cone(time):
    data = cones.characteristic()
    K, G = [data[key].subs(d.u, time) for key in ("K", "G")]
    kappa = data["kappa"].subs(d.u, time)
    assert kappa > 0
    assert (G-K).extract((0, 2), (0, 2)).det() < 0
    assert all(K[:i, :i].det() > 0 and G[:i, :i].det() > 0 for i in (1, 2, 3))
    assert data["speed_plus"].subs(d.u, time) > 1
    assert 0 < data["speed_minus"].subs(d.u, time) < 1


def test_front_limit_is_not_a_massive_finite_momentum_phase_speed():
    data = cones.characteristic()
    assert data["new_curl_absent_from_principal_speeds"] is True
    assert not data["K"].has(d.ZETA, d.q)
    assert not data["G"].has(d.ZETA, d.q)
    assert data["zero_curl_is_a_separate_rank_limit"] is True
    assert data["violation_below_an_EFT_cutoff_claim"] is False
    assert data["kappa"].subs(d.u, 0) == 0
