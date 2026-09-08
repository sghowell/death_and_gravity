"""Full-rank-one branch, all-constraint, and higher-derivative controls."""

from fractions import Fraction

import pytest
import sympy as sp
from p8_affine_traces import dynamics as d


def _zero(value):
    return all(item == 0 for item in value) if isinstance(value, sp.MatrixBase) else value == 0


def test_every_exact_identity_and_continuous_control():
    assert all(_zero(value) for value in d.checks().values())
    assert all(value is True for value in d.proof_checks().values())


def test_independent_two_trace_normalization_and_projective_kernel():
    result = d.geometry()
    eta = sp.diag(1, -1, -1, -1)
    expected = sp.Matrix.vstack(sp.Matrix.hstack(3*eta/8, -21*eta/8),
                               sp.Matrix.hstack(-21*eta/8, 3*eta/8))
    assert result["D_two"] == expected
    assert result["projective_residual"] == sp.zeros(8, 4)
    assert result["map"].rank() == 8


def test_full_rolling_U_is_not_merely_minus_V():
    rolling = d.rolling()
    assert rolling["parts"]["U"]["alpha"] == -rolling["parts"]["V"]["alpha"]
    assert rolling["parts"]["U"]["beta"] == -rolling["parts"]["V"]["beta"]
    assert sp.factor(rolling["parts"]["U"]["n"]+rolling["parts"]["V"]["n"]) != 0
    assert rolling["parts"]["U"]["n"] == -27*rolling["H"]/(8*rolling["h"])


@pytest.mark.parametrize("first,second,coupling,time", [
    (1, 0, sp.Rational(1, 3), sp.Rational(1, 100)),
    (0, 1, 2, -sp.Rational(1, 100)),
    (1, -1, sp.Rational(2, 5), sp.Rational(1, 2)),
    (1, 1, 1, -sp.Rational(1, 2)),
    (1, 2, sp.Rational(7, 3), sp.Rational(1, 10)),
])
def test_literal_joint_n_shift_W0_solve_has_one_negative_scalar_pivot(first, second, coupling, time):
    bg = d.old.background()
    gamma = sp.Rational(3, 8)*(first**2-14*first*second+second**2)
    h = bg["h"].subs(d.u, time)
    hubble = bg["H"].subs(d.u, time)
    J = bg["J"].subs(d.u, time)
    e = (first-second)/h
    momentum = (2*gamma*J/e**2 if gamma > 0 else -1/(gamma*coupling))+1
    values = {d.gamma: gamma, d.zeta: coupling, d.q: momentum, d.J: J, d.e: e,
              d.d: 3*(7*first+3*second)*hubble/(8*h),
              d.theta: bg["theta"].subs(d.u, time), d.ell: bg["ell"].subs(d.u, time),
              d.w: bg["w"].subs(d.u, time), d.lam: bg["lam"].subs(d.u, time)}
    L = d.regular()["before"].subs(values)
    aux = (d.n, d.shift, d.temporal)
    solution = sp.solve([sp.diff(L, field) for field in aux], aux)
    assert len(solution) == 3
    direct = sp.hessian(L.subs(solution), (d.vd, d.sd, d.sigmad))/2
    assert (direct-d.regular()["kinetic"].subs(values)).applyfunc(sp.factor) == sp.zeros(3)
    assert sp.factor(direct.det()) < 0
    branch = d.require_domain(first, second, coupling, time, momentum)
    assert branch["branch"] == ("negative_clock_pivot" if gamma > 0 else "negative_longitudinal_pivot")


def test_complete_square_not_a_frozen_metric_scalar_inference():
    K = d.regular()["kinetic"]
    assert K[1, 1] == sp.Rational(1, 2)
    assert K[0, 1] == d.w/(2*d.theta)
    assert sp.factor(K[0, 2]-d.regular()["C"]*d.d/d.theta) == 0
    assert d.regular()["shift_residual"] == 0


@pytest.mark.parametrize("ratio", (7+4*sp.sqrt(3), 7-4*sp.sqrt(3)))
@pytest.mark.parametrize("coupling", (-2, 3))
def test_both_exact_null_rays_and_both_coupling_signs(ratio, coupling):
    result = d.require_domain(ratio, 1, coupling, sp.Rational(1, 10), 2)
    assert result["gamma"] == 0
    assert result["branch"] == "null_schur_ostrogradsky"
    assert sp.simplify(ratio-1) != 0
    determinant = d.ostrogradsky()["highest_determinant"].subs(
        {d.a: 1, d.e: ratio-1, d.theta: 3, d.q: 2, d.zeta: coupling})
    assert sp.simplify(determinant).is_zero is False


def test_null_schur_forces_T_not_a_regular_Proca_limit():
    result = d.null_schur()
    assert result["M"].det() != 0
    assert result["projected_response"] == sp.zeros(1)
    assert result["trace_constraint"] == sp.zeros(1, 1)
    assert result["Euler_residual"] == sp.zeros(2, 1)
    # A formal gamma->0 in the REGULAR kinetic expression misses the constraint.
    assert d.regular()["C"].subs(d.gamma, 0) == d.zeta*d.q/2


def test_lapse_and_shift_are_both_accounted_for_in_the_higher_derivative_chart():
    nd = sp.Symbol("n_dot", real=True)
    raw = d.old.action()["base"]+d.zeta*d.q*(-d.e*nd+d.f*d.n)**2/2
    assert sp.factor(sp.diff(raw, d.shift)/(2*d.q)-d.theta*d.n+d.vd+d.ell*d.matter/2) == 0
    assert sp.diff(raw, d.n, d.shift) == 2*d.q*d.theta
    assert sp.diff(raw, nd, d.shift) == 0
    # Thus the lapse Euler equation reconstructs b even though it contains n_ddot.
    assert d.ostrogradsky()["lapse_reconstruction_coefficient"] == 2*d.q*d.theta


def test_independent_generic_highest_derivative_legendre_transform():
    volume, zz, qq, ae = sp.symbols("volume zz qq ae", nonzero=True, real=True)
    be, r, low0, low1, acceleration, velocity, p1, ps, p0, Q = sp.symbols(
        "be r low0 low1 acceleration velocity p1 ps p0 Q", real=True)
    L = volume*(low0+low1*velocity+velocity**2/2+zz*qq*(ae*acceleration+be*velocity+r)**2/2)
    highest = sp.hessian(L, (acceleration, velocity))
    assert sp.factor(highest.det()-volume**2*zz*qq*ae**2) == 0
    solution = sp.solve((sp.diff(L, acceleration)-p1, sp.diff(L, velocity)-ps),
                        (acceleration, velocity))
    H = sp.factor(p0*Q+p1*solution[acceleration]+ps*solution[velocity]-L.subs(solution))
    target = (p0*Q+(ps-(be/ae)*p1-volume*low1)**2/(2*volume)
              +p1**2/(2*volume*zz*qq*ae**2)-(r/ae)*p1-volume*low0)
    assert sp.factor(H-target) == 0
    assert sp.diff(H, p0) == Q


def test_actual_free_M1_makes_joint_highest_hessian_nondegenerate():
    result = d.ostrogradsky()
    assert result["highest_determinant"] == d.a**6*d.e**2*d.q*d.zeta/d.theta**2
    assert result["missing_free_matter_rank_control"] == 0
    assert result["P1_residual"] == result["Ps_residual"] == result["Legendre_residual"] == 0


@pytest.mark.parametrize("coupling", (-1, 1))
def test_affine_momentum_unbounded_for_either_sign(coupling):
    H = d.ostrogradsky()["Hamiltonian"].subs(d.zeta, coupling)
    offset = sp.Symbol("offset", real=True)
    assert sp.factor(H.subs(d.P0, d.P0+offset)-H-offset*d.vd) == 0
    assert sp.diff(H, d.P0, 2) == 0
    assert sp.diff(H, d.P0).subs(d.vd, 1) == 1


def test_nonzero_null_direction_cannot_remove_acceleration_by_equal_traces():
    gamma = d.coefficients()["gamma"]
    assert sp.factor(gamma.subs(d.B, d.A)) == -9*d.A**2/2
    assert gamma.subs({d.A: 1, d.B: 1}) != 0
    assert d.ostrogradsky()["highest_determinant"].subs(d.e, 0) == 0


@pytest.mark.parametrize("arguments", ((0, 0, 3, 0, 1), (2, 5, 0, 0, 1)))
def test_unchanged_auxiliary_controls(arguments):
    assert d.require_domain(*arguments)["branch"] == "unchanged_auxiliary"


@pytest.mark.parametrize("first,second", ((1, 0), (1, 1), (1, -1)))
def test_negative_coupling_nonnull_uses_only_the_unconstrained_transverse_test(first, second):
    result = d.require_domain(first, second, -1, 0, 1)
    assert result["branch"] == "negative_transverse_kinetic"
    assert result["gamma"] != 0


def test_rescaling_same_rank_one_action_preserves_threshold():
    original = d.require_domain(1, 0, 2, sp.Rational(1, 100), 8)
    scaled = d.require_domain(3, 0, sp.Rational(2, 9), sp.Rational(1, 100), 8)
    assert original["branch"] == scaled["branch"]
    assert original["threshold"] == scaled["threshold"]
    assert scaled["gamma"] == 9*original["gamma"]


def test_nonunit_physical_scales_and_null_mass_nonclaim():
    values = d.units(3, 2, 5)
    assert values["normalized_zeta"] == sp.Rational(5, 12)
    assert values["isolated_mass_squared_physical"] == sp.Rational(8, 5)
    assert values["isolated_mass_squared_normalized"] == sp.Rational(32, 5)
    null = d.units(3, 2, -5, 7+4*sp.sqrt(3), 1)
    assert null["gamma"] == 0
    assert null["isolated_mass_squared_physical"] is None


@pytest.mark.parametrize("value", [True, False, 1.0, sp.Float(1), "1", sp.true, sp.false, sp.I,
                                    sp.oo, sp.nan, sp.Symbol("unproved"), [1]])
def test_strict_exact_inputs_after_warmup(value):
    d.require_domain(1, 0, 1, sp.Rational(1, 100), 8)
    with pytest.raises((TypeError, ValueError)):
        d.require_domain(value, 0, 1, sp.Rational(1, 100), 8)


@pytest.mark.parametrize("arguments", [(1, 0, 1, 0, 8), (1, 0, 1, 1, 0),
                                       (1, 0, 1, 1, -1), (1, 1, 1, 1, sp.Rational(2, 9))])
def test_wrong_chart_zero_mode_and_temporal_pole_rejected(arguments):
    with pytest.raises(ValueError):
        d.require_domain(*arguments)


def test_strict_clock_threshold_and_exact_fraction_support():
    bg = d.old.background()
    time = sp.Rational(1, 100)
    threshold = sp.factor(3*bg["J"].subs(d.u, time)*bg["h"].subs(d.u, time)**2/4)
    with pytest.raises(ValueError):
        d.require_domain(1, 0, 1, time, threshold)
    assert d.require_domain(Fraction(1), 0, Fraction(1), time, 8)["gamma"] == sp.Rational(3, 8)


@pytest.mark.parametrize("arguments", [(0, 1, 1), (1, 0, 1), (-1, 1, 1), (1, 1.0, 1)])
def test_units_exclude_invalid_gravitational_scales(arguments):
    with pytest.raises((TypeError, ValueError)):
        d.units(*arguments)
