"""Independent full-source, constraint and scope regressions."""

import copy

import pytest
import sympy as s
from p8_vacuum_affine_canonical_boundary_corrected_hybrid import audit as previous
from p8_vacuum_affine_global_lapse_chart_obstruction import (
    audit,
    auxiliary,
    branch,
    fields,
    fold,
)

RESIDUALS = audit.residuals()
GATES = audit.gates()


@pytest.mark.parametrize("name", tuple(RESIDUALS))
def test_every_exact_full_source_residual(name):
    value = RESIDUALS[name]
    assert (
        value == s.zeros(*value.shape)
        if isinstance(value, s.MatrixBase)
        else value == 0
    )


@pytest.mark.parametrize("name", tuple(GATES))
def test_every_stated_proof_gate(name):
    assert GATES[name] is True


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_every_unsupported_input_is_rejected(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("order", range(4))
def test_clock_jets_from_independent_explicit_primitive_Hamiltonian(order):
    N, m, sh, c = fold.N, fold.M, fold.SHEAR, fold.CURVATURE
    rho, p = fold.RHO, fold.PRESSURE
    # This is only a representative of exact finite CLOCK jets. The nonzero
    # primitive is included; the complete off-clock source is not replaced.
    literal = (
        24
        - s.Rational(72, 25) * N ** (-s.Rational(3, 2))
        - s.Rational(2847, 200) * s.sqrt(N)
        + s.Rational(28, 25) * N ** s.Rational(5, 2)
        + (m * m - c) / (2 * s.sqrt(N))
        + 2 * sh * N ** s.Rational(3, 2)
        + (rho + p) * s.sqrt(N) / 2
        + (p - rho) * N ** s.Rational(5, 2) / 2
    )
    independently_differentiated = s.diff(literal, N, order + 1).subs(N, 1)
    assert s.factor(independently_differentiated - fold.constraint_jet(order)) == 0


def test_primitive_cannot_be_silently_dropped():
    N = fold.N
    Iu = 24 - 6 * N ** (-s.Rational(3, 2)) - 18 * s.sqrt(N)
    assert Iu.subs(N, 1) == 0
    assert s.diff(Iu, N).subs(N, 1) == 0
    assert s.diff(Iu, N, 2).subs(N, 1) == -18
    assert s.diff(Iu, N, 3).subs(N, 1) == 72
    assert s.factor(s.diff(Iu, N) - fold.original.tree(fold.original.PRIMITIVE_N)) == 0


@pytest.mark.parametrize(
    "rho_value,pressure_value",
    (
        (0, 0),
        (s.Rational(1, 10**500), -s.Rational(2, 10**501)),
        (-s.Rational(3, 10**500), s.Rational(1, 10**500)),
    ),
)
def test_literal_profile_fold_solves_both_equations(rho_value, pressure_value):
    # Independent small algebraic diagnostics, NOT the actual profile choice.
    N, m, sh, c = fold.N, fold.M, fold.SHEAR, fold.CURVATURE
    rho, p = s.sympify(rho_value), s.sympify(pressure_value)
    R = N**-2
    F = -(624 * N**-4 + 753 * N**-2 + 224) / 200 - p - (rho + p) * (N**-2 - 1) / 2
    Iu = 24 - 6 * N ** (-s.Rational(3, 2)) - 18 * s.sqrt(N)
    H = (
        N
        * (
            -(R ** (-s.Rational(3, 4))) * F
            + 2 * sh * R ** (-s.Rational(1, 4))
            + (m * m - c) * R ** s.Rational(3, 4) / 2
        )
        + Iu
    )
    solved_s = s.Rational(81, 160) - 2 * p / 3 + 7 * rho / 12
    assert not H.has(s.Float)
    assert not solved_s.has(s.Float)
    rule = {
        N: 1,
        m * m: s.Rational(1, 100) + 6 * p - 4 * rho,
        sh: solved_s,
        c: -12 * solved_s,
    }
    for degree in (1, 2):
        assert s.factor(s.diff(H, N, degree).subs(rule, simultaneous=True)) == 0
    assert (
        s.factor(
            s.diff(H, N, 3).subs(rule, simultaneous=True)
            - (11391 + 1400 * p - 1600 * rho) / 400
        )
        == 0
    )


@pytest.mark.parametrize("axis", range(3))
def test_axis_input_is_in_the_actual_band_and_preserves_all_sources(axis):
    row = fields.fixture(axis)
    assert row["whole_original_frequency_binding"] == fold.fifth.LOW == 10**64
    assert len(set(row["whole_input_axes"])) == 3
    assert row["whole_full_source_momentum_generator"] == s.zeros(3, 1)
    assert row["whole_metric"] != s.eye(3)
    assert row["whole_full_raw_metric_momentum"] != s.zeros(3)
    assert row["whole_M1_momentum_input"].has(s.cos)
    assert row["whole_shape_input"] == s.zeros(3)


def test_general_off_diagonal_density_generator_contact():
    x = s.symbols("independent_x0:3", real=True)
    v = s.Function("arbitrary_conformal_field")(*x)
    kap = s.Symbol("positive_kappa", positive=True)
    a, b, c, d, e = (
        s.Function(name)(*x) for name in ("tt_a", "tt_b", "tt_c", "tt_d", "tt_e")
    )
    tt = s.Matrix([[a, c, d], [c, b, e], [d, e, -a - b]])
    metric, pi = s.exp(2 * v) * s.eye(3), -kap * s.exp(-2 * v) * tt
    for i in range(3):
        density_generator = -2 * sum(
            s.diff((metric * pi)[i, j], x[j]) for j in range(3)
        )
        density_generator += sum(
            pi[j, k] * s.diff(metric[j, k], x[i]) for j in range(3) for k in range(3)
        )
        assert (
            s.simplify(
                density_generator
                - 2 * kap * sum(s.diff(tt[i, j], x[j]) for j in range(3))
            )
            == 0
        )


@pytest.mark.parametrize("value", (s.Rational(3, 5), 1, (4 + s.sqrt(7)) / 3, 3, 10))
def test_independent_mixed_stationarity_at_and_away_from_sequential_crossing(value):
    K, T, p, G = s.symbols(
        "independent_K independent_T independent_p independent_G", real=True
    )
    aa, U, r = -(value ** s.Rational(1, 4)) / 3, value ** (-s.Rational(3, 4)), value - 1
    B, F, lower = s.Rational(3, 7), s.Rational(5, 11), -s.Rational(2, 9)
    L = aa * K * K + B * K + F + U * (T - r * K - lower) ** 2 / 2
    kstar = (p - B - r * G) / (2 * aa)
    tstar = r * kstar + lower - G / U
    for residual in (s.diff(L, K) - p, s.diff(L, T) + G):
        assert s.simplify(residual.subs({K: kstar, T: tstar}, simultaneous=True)) == 0
    assert s.simplify(s.hessian(L, (K, T)).det() + 2 / (3 * s.sqrt(value))) == 0


def test_actual_source_profiles_and_remainders_remain_nonzero_objects():
    packet = fold.data()
    assert len(packet["whole_actual_fixed_profile_bindings"]) == 2
    assert all(value != 0 for value in packet["whole_actual_source_remainders"])
    assert packet["whole_all_three_fixed_vacuum_constants"] != 0
    assert fold.whole()["H"].has(fold.original.Iu)
    assert auxiliary.data()["whole_sequential_crossing_R"] != 1


@pytest.mark.parametrize("jj", (-3, 0, 7))
def test_second_class_rank_change_cannot_be_cancelled_by_secondary_brackets(jj):
    D = s.diag(0, -1)
    J = s.Matrix([[0, jj], [-jj, 0]])
    block = s.zeros(2).row_join(-D).col_join(D.row_join(J))
    assert block.det() == 0 and block.rank() == 2
    for a in (-2, 3):
        regular = s.diag(a, -1)
        test = s.zeros(2).row_join(-regular).col_join(regular.row_join(J))
        assert test.det() == a * a and test.rank() == 4


def test_transverse_variation_differs_from_the_tangent_path():
    packet = branch.data()
    sh = fold.SHEAR
    assert s.diff(packet["whole_connected_path_C_N"], sh) == 6
    assert fold.clock_jet(s.diff(fold.whole()["C"], sh), 0) == 3
    # Changing BOTH curvature and shear along the path keeps C zero; holding
    # the metric fixed for the transverse derivative has nonzero C_s.
    assert (
        s.factor(
            fold.constraint_jet(0).subs(
                {fold.M**2: fold.MATTER_SQUARE, fold.CURVATURE: -12 * sh},
                simultaneous=True,
            )
        )
        == 0
    )


@pytest.mark.parametrize("sign", (-1, 1))
def test_local_fold_control_has_two_branches_and_unbounded_envelope(sign):
    delta = s.Symbol("positive_distance_to_fold", positive=True)
    D = s.Symbol("positive_second_jet", positive=True)
    n = sign * s.sqrt(6 * delta / D)
    assert s.simplify(-3 * delta + D * n * n / 2) == 0
    envelope = -9 / (D * n)
    assert s.simplify(envelope * envelope - 27 / (2 * D * delta)) == 0
    assert s.limit(abs(envelope), delta, 0, dir="+") == s.oo


@pytest.mark.parametrize("call,good", ((fields.fixture, 1), (fold.constraint_jet, 1)))
def test_warm_cache_does_not_accept_boolean_or_float_aliases(call, good):
    call(good)
    call(s.Integer(good))
    for bad in (True, 1.0, s.Float(1)):
        with pytest.raises(TypeError):
            call(bad)


def test_full_previous_frontiers_and_physical_qualifications_are_unchanged():
    assert audit.frontier() == previous.frontier()
    assert audit.matching()[:-1] == previous.matching()
    assert audit.qualifications() == previous.qualifications()
    assert len(audit.qualifications()) == 6
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 132
    assert audit.validate_scope(audit.frontier(), audit.matching())


@pytest.mark.parametrize("which", ("primitive", "matching"))
def test_changed_status_without_changed_count_is_rejected(which):
    original, matched = copy.deepcopy(audit.frontier()), copy.deepcopy(audit.matching())
    rows = original if which == "primitive" else matched
    rows[0]["status"] = "UNSUPPORTED_ORIGINAL_P8_COMPLETION"
    with pytest.raises(ValueError):
        audit.validate_scope(original, matched)


@pytest.mark.parametrize("observable", audit.OBSERVABLES)
def test_only_scoped_observables_are_admitted(observable):
    assert audit.require_observable(observable) == observable
    assert audit.require_model(audit.MODEL) == audit.MODEL


def test_affine_profile_enclosure_does_not_accept_quadratic_shortcut():
    with pytest.raises(ValueError):
        fold.affine_error(fold.RHO**2 + fold.PRESSURE)
