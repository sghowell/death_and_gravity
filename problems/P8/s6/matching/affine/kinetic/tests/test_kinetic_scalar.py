"""Independent unreduced-action, constraint, and failed-inference controls."""

from fractions import Fraction

import pytest
import sympy as sp
from p8_affine_kinetic import scalar as s


def _zero(value):
    return all(item == 0 for item in value) if isinstance(value, sp.MatrixBase) else value == 0


def test_every_exact_identity():
    assert all(_zero(value) for value in s.checks().values())


def test_every_continuous_margin_and_control():
    assert all(bool(value) for value in s.proof_checks().values())


def test_literal_rolling_vector_has_no_spatial_metric_velocity():
    result = s.vector_reconstruction()
    assert result["trace_perturbation_residual"] == 0
    assert result["lapse_coefficient"] == -3*result["H"]/(8*result["h"])
    assert result["d"] == 21*result["H"]/(8*result["h"])


def test_true_gradient_shift_includes_coefficient_derivative():
    alpha, alphad, ndot, ni, nidot, w0i, widot = sp.symbols(
        "alpha alphad ndot ni nidot w0i widot", real=True)
    # V_i=W_i-alpha*n_i and V_0=W_0-alpha*n_dot-alpha_dot*n.
    direct = widot-alphad*ni-alpha*nidot-(w0i-alpha*nidot-alphad*ni)
    assert sp.expand(direct-(widot-w0i)) == 0
    naive = widot-alphad*ni-alpha*nidot-(w0i-alpha*nidot)
    assert sp.expand(naive-(widot-w0i)) == -alphad*ni
    assert sp.diff(naive, ndot) == 0


def test_raw_lapse_velocity_is_not_loss_of_degeneracy():
    nvel, sigmavel, v0, q, coupling = sp.symbols("nvel sigmavel v0 q coupling", positive=True)
    raw = sp.Rational(4, 3)*(v0+3*nvel/2)**2+coupling*q*(sigmavel-v0)**2/2
    assert sp.diff(raw, nvel, 2) == 6
    v0_solution = sp.solve(sp.diff(raw, v0), v0)[0]
    reduced = sp.factor(raw.subs(v0, v0_solution))
    hessian = sp.hessian(reduced, (nvel, sigmavel))
    assert sp.factor(hessian.det()) == 0
    assert hessian.rank() == 1
    assert hessian*sp.Matrix([2, -3]) == sp.zeros(2, 1)


@pytest.mark.parametrize("time,coupling", [
    (sp.Rational(1, 100), sp.Rational(1, 5)),
    (-sp.Rational(1, 100), sp.Rational(3, 2)),
    (sp.Rational(1, 2), sp.Integer(7)),
    (-sp.Rational(3, 2), sp.Rational(1, 13)),
])
def test_literal_three_auxiliary_elimination_matches_configuration_kinetic(time, coupling):
    """Solve n, shift and W0 together, not a copied kinetic matrix."""
    bg = s.background()
    point = {s.theta: bg["theta"].subs(s.u, time), s.J: bg["J"].subs(s.u, time),
             s.lam: bg["lam"].subs(s.u, time), s.ell: bg["ell"].subs(s.u, time),
             s.w: bg["w"].subs(s.u, time), s.dd: bg["d"].subs(s.u, time),
             s.e: bg["e"].subs(s.u, time), s.zeta: coupling}
    point[s.q] = sp.factor(3*point[s.J]/(4*point[s.e]**2)+1)
    lagrangian = s.action()["full_before_temporal"].subs(point)
    auxiliary = (s.n, s.shift, s.temporal)
    solution = sp.solve([sp.diff(lagrangian, item) for item in auxiliary], auxiliary)
    assert len(solution) == 3
    reduced = sp.factor(lagrangian.subs(solution))
    direct = sp.hessian(reduced, (s.vd, s.sd, s.sigmad))/2
    expected = s.unitary()["kinetic"].subs(point)
    assert (direct-expected).applyfunc(sp.factor) == sp.zeros(3)
    # Positive matter and vector pivots leave a strictly negative clock Schur.
    assert direct[1, 1] == sp.Rational(1, 2)
    assert direct[2, 2] > 0
    assert sp.factor(direct.det()) < 0


def test_square_completion_retains_the_actual_free_matter_mixing():
    K = s.unitary()["kinetic"]
    assert K[0, 1] == s.w/(2*s.theta)
    assert K[1, 1] == sp.Rational(1, 2)
    assert K[1, 2] == 0
    assert s.unitary()["pivots"][-1] == (s.J-sp.Rational(4, 3)*s.q*s.e**2)/s.theta**2
    assert sp.factor((s.background()["a"]**3*s.background()["ell"])-sp.Rational(1, 10)) == 0


def test_isolated_positive_proca_control_does_not_prove_coupled_health():
    extra = s.action()["extra"].subs(s.n, 0)
    assert sp.factor(sp.diff(extra, s.sigmad, 2)/2-s.action()["C"]) == 0
    assert s.action()["C"].is_positive
    assert s.center()["kinetic"][:2, :2].subs(s.q, 8).det() < 0


@pytest.mark.parametrize("qvalue", [8, 16, 100, sp.Rational(801, 100)])
def test_center_high_momentum_corrobation_and_baseline_control(qvalue):
    result = s.center()
    K = result["kinetic"].subs({s.q: qvalue, s.zeta: sp.Rational(1, 7)})
    assert K[0, 0] < 0 < K[1, 1]
    assert K[:2, :2].det() < 0 < K[2, 2]
    baseline = sp.Matrix([[6*s.q/(s.q-6), s.q/(20*(s.q-6))],
                          [s.q/(20*(s.q-6)), (200*s.q-1199)/(400*(s.q-6))]])
    baseline = baseline.subs(s.q, qvalue)
    assert baseline[0, 0] > 0 and baseline.det() > 0


def test_center_chart_divisors_are_distinct_and_not_a_mass_cutoff():
    ce = s.center()
    assert ce["lapse_pole"] != ce["velocity_chart_pole"] != ce["baseline_velocity_chart_pole"]
    assert ce["lapse_pole"] == sp.Rational(3597, 3200)
    assert ce["velocity_chart_pole"] == sp.Rational(18, 19)
    assert s.calibration()["EFT_low_frequency_ghost_or_cutoff_claim"] is False


def test_exact_center_coefficient_jets_and_cauchy_relation():
    jets = s.center_jets()
    assert (jets["H_1"], jets["theta_1"], jets["d_1"]) == (4, 3, sp.Rational(21, 2))
    assert jets["e_1"] == jets["alpha_1"] == 0
    result = s.center_pencil()
    assert result["cauchy_residual"] == sp.zeros(3, 6)
    assert result["omitted_coefficient_jets"] != sp.zeros(3, 6)
    omega = sp.zeros(6)
    omega[:3, 3:] = sp.eye(3)
    omega[3:, :3] = -sp.eye(3)
    assert result["phase"].T*omega+omega*result["phase"] == sp.zeros(6)


def test_zero_zeta_is_a_rank_change_not_a_continuous_healthy_proca_limit():
    raw = s.action()["extra_before_temporal"].subs(s.zeta, 0)
    solution = sp.solve([sp.diff(raw, field) for field in (s.temporal, s.sigma)],
                        (s.temporal, s.sigma))
    assert solution == {s.temporal: -s.dd*s.n, s.sigma: -s.e*s.n}
    assert sp.factor(raw.subs(solution)) == 0
    assert s.action()["C"].subs(s.zeta, 0) == 0


@pytest.mark.parametrize("time,momentum,coupling", [
    (sp.Rational(1, 100), 8, Fraction(1, 2)),
    (-sp.Rational(1, 100), 8, sp.sqrt(2)),
])
def test_exact_unhealthy_domain(time, momentum, coupling):
    result = s.require_domain(time, momentum, coupling, unhealthy=True)
    assert result["q"] > result["threshold"] > 0


@pytest.mark.parametrize("value", [True, False, 1.0, sp.Float(1), "1", sp.true, sp.false, sp.I,
                                    sp.oo, sp.nan, sp.Symbol("unproved"), [1]])
def test_guarded_inputs_reject_before_any_cached_exact_result(value):
    s.require_domain(sp.Rational(1, 100), 8, 1)
    with pytest.raises((TypeError, ValueError)):
        s.require_domain(sp.Rational(1, 100), 8, value)


@pytest.mark.parametrize("time,momentum,coupling", [(0, 8, 1), (1, 0, 1), (1, -1, 1),
                                                     (1, 8, 0), (1, 8, -1)])
def test_domain_excludes_crossing_zero_mode_and_wrong_maxwell_sign(time, momentum, coupling):
    with pytest.raises(ValueError):
        s.require_domain(time, momentum, coupling)


def test_strict_threshold_equality_and_wrong_flag_are_rejected():
    time = sp.Rational(1, 10)
    threshold = s.require_domain(time, 1, 1)["threshold"]
    with pytest.raises(ValueError):
        s.require_domain(time, threshold, 1, unhealthy=True)
    with pytest.raises(TypeError):
        s.require_domain(time, 8, 1, punctured=1)


def test_nonunit_normalization_is_not_physical_zeta():
    result = s.units(3, 2, 5)
    assert result == {"normalized_zeta": sp.Rational(5, 12),
                      "isolated_proca_mass_squared": sp.Rational(8, 5),
                      "normalized_isolated_mass_squared": sp.Rational(32, 5)}
    assert result["normalized_zeta"] != 5


@pytest.mark.parametrize("arguments", [(0, 1, 1), (1, 0, 1), (1, 1, -1), (1, 1, 1.0)])
def test_physical_units_have_exact_positive_domain(arguments):
    with pytest.raises((TypeError, ValueError)):
        s.units(*arguments)
