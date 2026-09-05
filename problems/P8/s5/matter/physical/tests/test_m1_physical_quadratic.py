import pytest
import sympy as sp
from p8_m1_physical import (
    background,
    lagrangian,
    momentum,
    quadratic,
    regressions,
    vertices,
)
from p8_m1_physical.vertices import Leg


@pytest.mark.parametrize("checks", [background.compact_checks, momentum.symplectic_projector_checks,
                                    regressions.boundary_checks, quadratic.identities,
                                    vertices.stationary_bridge_checks, lagrangian.stationary_checks])
def test_generic_identities(checks):
    assert all(sp.cancel(value) == 0 for value in checks().values())


def test_symbolic_time_full_quadratic_bridge():
    assert len(regressions.quadratic_checks()) == 35
    assert all(value == 0 for value in regressions.quadratic_checks().values())


@pytest.mark.parametrize("point", [0, sp.Rational(1, 3), -1, 1])
def test_literal_quadratic_at_bounce_interior_and_tail_limits(point):
    assert all(value == 0 for value in regressions.quadratic_checks(point).values())


def test_full_gamma_kinetic_is_not_principal_matrix():
    K = quadratic.response(0, "gamma", 8)["kinetic"]
    assert K == sp.Matrix([[24, sp.Rational(1, 5)], [sp.Rational(1, 5), sp.Rational(401, 800)]])
    wave, opposite = (2, 2, 0), (-2, -2, 0)
    for i, left in enumerate(("s_dot", "m_dot")):
        for j, right in enumerate(("s_dot", "m_dot")):
            assert lagrangian.kernel((Leg(wave, left), Leg(opposite, right)), 0, "gamma")["kernel"] == 2*K[i, j]


def test_gamma_velocity_pole_is_not_phase_singularity():
    with pytest.raises(ValueError, match="Singular scalar velocity chart"):
        quadratic.response(0, "gamma", 6)
    legs = (Leg((1, 1, 2), "p"), Leg((-1, -1, -2), "p"))
    assert not vertices.hamiltonian_kernel(legs, 0, "gamma")["kernel"].has(sp.zoo, sp.nan)
    with pytest.raises(ValueError, match="not positive"):
        quadratic.response(0, "gamma", 5)
    with pytest.raises(ValueError, match="Singular scalar velocity chart"):
        quadratic.response(0, "unitary", 8)


def test_boundary_and_constraint_omissions_fail():
    data = regressions.negative_controls()
    assert sp.Rational(data["linear_only_York_quartic_error"]) == -sp.Rational(3, 100)
    assert sp.Rational(data["missing_matter_clock_tadpole"]) != 0
    assert sp.Rational(data["missing_mixed_boundary_quadratic_error"]) != 0
