"""Independent joint constraints and the original physical matter cone."""
import pytest
import sympy as sp
from p8_affine_aligned import dynamics as d


def test_all_quadratic_constraint_and_cone_identities():
    for value in d.checks().values():
        assert all(item == 0 for item in value) if isinstance(value, sp.MatrixBase) else value == 0


@pytest.mark.parametrize("time", (sp.Rational(1, 2), -sp.Rational(1, 2), 3))
@pytest.mark.parametrize("momentum", (sp.Rational(1, 10), 17, 10000))
def test_original_temporal_variable_joint_solve(time, momentum):
    old = d.old
    bg = old.background()
    values = {symbol: bg[name].subs(d.u, time) for symbol, name in
              ((old.theta, "theta"), (old.ell, "ell"), (old.lam, "lam"), (old.J, "J"), (old.w, "w"))}
    values.update({d.q: momentum, d.ZETA: sp.Rational(1, 2000),
                   d.previous.DD: d.previous.source()["d"].subs(d.u, time)})
    L = d.scalar()["literal_before_temporal_shift"].subs(values)
    auxiliaries = (old.n, old.shift, old.temporal)
    solution = sp.solve([sp.diff(L, item) for item in auxiliaries], auxiliaries)
    assert len(solution) == 3
    K = (sp.hessian(L.subs(solution), (old.vd, old.sd, old.sigmad))/2).applyfunc(sp.factor)
    assert K == d.scalar()["expected_kinetic"].subs(values)
    assert all(K[:i, :i].det() > 0 for i in (1, 2, 3))
    assert K[0, 2] == K[1, 2] == 0


@pytest.mark.parametrize("time", (sp.Rational(1, 100), sp.Rational(1, 2), -3, 10))
def test_full_scalar_front_cones_are_exactly_the_matter_cone(time):
    data = d.cones()
    K, G = [matrix.subs(d.u, time) for matrix in (data["K"], data["G"])]
    assert (K-G).applyfunc(sp.factor) == sp.zeros(3)
    assert all(K[:i, :i].det() > 0 for i in (1, 2, 3))
    assert data["finite_q_massive_phase_speed_not_used"] is True


def test_first_order_crossing_is_the_old_system_plus_a_regular_vector():
    data = d.first_order()
    assert data["no_Theta_inverse"] is True
    assert not data["H"].has(d.previous.DD)
    vector = (1+1/(d.ZETA*d.q))*d.old.pi**2/2+d.q*d.old.sigma**2/2
    scalar = data["H"]-vector
    assert not scalar.has(d.ZETA, d.old.pi, d.old.sigma)
    assert sp.diff(data["H"], d.old.pi, 2).is_positive
