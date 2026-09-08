"""Joint scalar constraints, principal-energy square and crossing controls."""
import pytest
import sympy as sp
from p8_affine_retuned import bounds
from p8_affine_retuned import dynamics as d


def test_all_exact_source_constraint_and_canonical_identities():
    for value in d.checks().values():
        assert all(item == 0 for item in value) if isinstance(value, sp.MatrixBase) else value == 0


@pytest.mark.parametrize("time", (sp.Rational(1, 2), -sp.Rational(1, 2), 3))
@pytest.mark.parametrize("coupling", (sp.Rational(1, 2000), sp.Rational(1, 12000)))
@pytest.mark.parametrize("momentum", (sp.Rational(1, 10), 17, 10000))
def test_joint_lapse_shift_vector_temporal_solve(time, coupling, momentum):
    old = d.old
    bg = old.background()
    values = {symbol: bg[name].subs(d.u, time) for symbol, name in
              ((old.theta, "theta"), (old.ell, "ell"), (old.lam, "lam"), (old.w, "w"), (old.J, "J"))}
    values.update({old.q: momentum, d.ZETA: coupling, d.DD: d.source()["d"].subs(d.u, time)})
    L = d.action()["before"].subs(values)
    auxiliaries = (old.n, old.shift, old.temporal)
    solution = sp.solve([sp.diff(L, item) for item in auxiliaries], auxiliaries)
    assert len(solution) == 3
    kinetic = (sp.hessian(L.subs(solution), (old.vd, old.sd, old.sigmad))/2).applyfunc(sp.factor)
    assert kinetic == d.action()["kinetic"].subs(values)
    assert all(kinetic[:k, :k].det() > 0 for k in (1, 2, 3))
    energy = d.principal()["G"].subs({d.u: time, d.ZETA: coupling, d.q: momentum})
    assert all(energy[:k, :k].det() > 0 for k in (1, 2, 3))


def test_exact_rank_one_gradient_comparison_not_a_time_scan():
    data = d.principal()
    assert bounds.positive_even(data["bound_denominator"])
    assert bounds.positive_even(data["bound_witness"])
    assert sp.Poly(data["bound_witness"], d.u).degree() == 50
    assert min(sp.Poly(data["bound_witness"], d.u).coeffs()) == 10790325
    assert d.source()["f"].subs(d.u, 0) == 5
    assert sp.diff(d.source()["f"], d.u).subs(d.u, 0) == 0
    assert data["point_transformation_residual"] == data["principal_action_residual"] == 0


def test_independent_principal_Legendre_square_with_all_velocity_position_terms():
    data = d.principal()
    values = {d.u: sp.Rational(1, 2), d.q: 17, d.ZETA: sp.Rational(1, 2000)}
    K, G, linear = [data[key].subs(values) for key in ("K", "G", "linear_velocity_coefficient")]
    momentum = sp.Matrix(sp.symbols("P_v P_s P_rho", real=True))
    velocity = K.inv()*(momentum-linear)/2
    L = data["L"].subs(values)
    H = (momentum.T*velocity)[0]-L.subs(dict(zip(data["velocity"], velocity, strict=True)), simultaneous=True)
    expected = ((momentum-linear).T*K.inv()*(momentum-linear))[0]/4+17*(data["x"].T*G*data["x"])[0]
    assert sp.factor(H-expected) == 0


@pytest.mark.parametrize("momentum", (sp.Rational(61, 10), 7, 100))
def test_regular_center_includes_the_vector_and_no_Theta_inverse(momentum):
    data = d.center()
    K = data["kinetic"].subs({d.q: momentum, d.ZETA: sp.Rational(1, 12000)})
    assert all(K[:k, :k].det() > 0 for k in (1, 2, 3))
    assert not sp.denom(sp.together(d.first_order()["H"])).has(d.old.theta)
    assert not sp.denom(sp.together(d.first_order()["H"])).has(d.DD)
    assert d.source()["d"].subs(d.u, 0) == 0
    # The new first coefficient jet is not claimed identical to the old one.
    assert sp.diff(d.source()["d"], d.u).subs(d.u, 0) == 15


def test_transverse_vector_has_positive_actual_FLRW_Hamiltonian():
    a, momentum, field, velocity = sp.symbols("a P T T_dot", positive=True)
    L = a*(d.ZETA*velocity**2-(d.ZETA*d.q+1)*field**2)/2
    solved = momentum/(a*d.ZETA)
    H = momentum*solved-L.subs(velocity, solved)
    expected = momentum**2/(2*a*d.ZETA)+a*(d.ZETA*d.q+1)*field**2/2
    assert sp.factor(H-expected) == 0
    assert sp.diff(H, momentum, 2).is_positive
    assert sp.diff(H, field, 2).is_positive


def test_nonunit_scale_is_an_isolated_mass_not_a_coupled_gap():
    data = bounds.units(3, 2, sp.Rational(1, 1000))
    assert data["normalized_zeta"] == sp.Rational(1, 12000)
    assert data["isolated_Proca_mass_squared_physical"] == 3000
    assert data["isolated_Proca_mass_squared_normalized"] == 12000
    assert data["isolated_mass_is_not_a_coupled_gap_or_cutoff"] is True
