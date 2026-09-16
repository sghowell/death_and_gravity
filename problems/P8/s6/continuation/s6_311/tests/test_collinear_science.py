"""Independent off-shell, angular-endpoint and original-bound regression tests."""

import pytest
import sympy as s
from p8_vacuum_affine_complete_two_graviton_tree import checks as samples
from p8_vacuum_affine_complete_two_graviton_tree import trees as e
from p8_vacuum_affine_massive_scalar_graviton_bremsstrahlung import tree as null_only
from p8_vacuum_affine_minimal_gravity_radiation import ward as null_ward
from p8_vacuum_affine_two_real_collinear_current import (
    audit,
    bounds,
    collinear,
    current,
    source,
)

ROWS, GATES = audit.residuals(), audit.gates()


@pytest.mark.parametrize("name", tuple(ROWS))
def test_all_exact_residuals(name):
    value = ROWS[name]
    assert all(
        v == 0 for v in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    assert not s.sympify(value).has(s.Float)


@pytest.mark.parametrize("name", tuple(GATES))
def test_written_proof_gates(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_original_scope_rejections(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_null_root_restriction_recovers_independent_frozen_channel():
    assert (
        s.cancel(
            current.rank_one_channel().subs(current.ell, 0)
            - null_ward.rank_one_channel()
        )
        == 0
    )


def test_timelike_root_cannot_use_the_null_scalar_propagators():
    points, qs, _ = samples.nonopposite_state()
    Q = qs[0] + qs[1]
    correct = current.matter_current(
        points, Q, 128, s.Rational(2, 3), s.Rational(-1, 7)
    )
    # Negative control only: the frozen null-only function is used outside its domain.
    wrong = null_only.whole_tensor(
        points,
        Q,
        mass=1,
        heavy=128,
        cubic=s.Rational(2, 3),
        contact=s.Rational(-1, 7),
        kappa=1,
    )
    assert (Q.T * e.ETA * correct).applyfunc(s.factor) == s.zeros(1, 4)
    assert (Q.T * e.ETA * wrong).applyfunc(s.factor) != s.zeros(1, 4)


@pytest.mark.parametrize("pair", ((0, 0), (0, 1), (1, 0), (1, 1)))
@pytest.mark.parametrize(
    "fraction", (s.Rational(1, 1000), s.Rational(1, 2), s.Rational(999, 1000))
)
@pytest.mark.parametrize(
    "angle", (s.S.Zero, s.Rational(1, 10), s.S.One, s.Integer(10), s.Integer(10) ** 40)
)
def test_literal_coefficients_obey_all_angle_budgets(pair, fraction, angle):
    a, b, r, records, _ = collinear.coefficients()
    values = [
        s.factor(v.subs({a: fraction, b: 1 - fraction, r: angle}, simultaneous=True))
        for v in records[pair]
    ]
    norm = fraction * (1 - fraction) * sum(abs(v) for v in values)
    expected = {(0, 0): 530, (0, 1): 128, (1, 0): 88, (1, 1): 268}
    assert norm <= expected[pair]


@pytest.mark.parametrize("pair", ((0, 0), (0, 1), (1, 0), (1, 1)))
def test_antiparallel_angular_endpoint_has_no_coefficient_pole(pair):
    a, b, r, records, _ = collinear.coefficients()
    for value in records[pair]:
        limit = s.limit(value.subs({a: s.Rational(1, 3), b: s.Rational(2, 3)}), r, s.oo)
        assert limit.is_finite is True


def test_unconserved_hard_tensor_cannot_drop_the_angular_pole():
    r = s.Symbol("r", positive=True)
    c = (1 - r * r) / (1 + r * r)
    d = 2 * r / (1 + r * r)
    q1 = e.imm([1, 0, 0, 1])
    q2 = e.imm([1, d, 0, c])
    Q = q1 + q2
    A = e.imm(s.diag(0, 1, -1, 0))
    v = s.Matrix([c, 0, -d])
    y = s.Matrix([0, 1, 0])
    B = s.zeros(4)
    B[1:, 1:] = v * v.T - y * y.T
    B = e.imm(B)
    U = s.zeros(4)
    U[0, 0] = 1
    assert Q.T * e.ETA * U != s.zeros(1, 4)
    H = e.imm(e.ETA * U * e.ETA - e.ETA * s.trace(e.ETA * U) / 2)
    wrong = s.factor(-e.cubic_gravity((A, B, H), (q1, q2, -Q)) / e.old.dot(Q, Q))
    assert s.limit(r * r * wrong, r, 0, dir="+") != 0


@pytest.mark.parametrize(
    "coefficients",
    (
        (s.Rational(3, 5) / s.sqrt(2), s.Rational(4, 5) / s.sqrt(2)),
        (s.Rational(1, 2), s.I / 2),
        (1 / s.sqrt(2), 0),
    ),
)
def test_unit_TT_basis_coefficients_have_l1_norm_at_most_one(coefficients):
    a, b = coefficients
    assert s.simplify(2 * (abs(a) ** 2 + abs(b) ** 2) - 1) == 0
    assert s.simplify(abs(a) + abs(b)) <= 1


def test_spatial_frobenius_norm_is_preserved_by_the_rational_rotation():
    R = s.Matrix(
        [
            [s.Rational(3, 5), -s.Rational(4, 5), 0],
            [s.Rational(4, 5), s.Rational(3, 5), 0],
            [0, 0, 1],
        ]
    )
    x = s.symbols("x0:6", real=True)
    U = s.Matrix([[x[0], x[1], x[2]], [x[1], x[3], x[4]], [x[2], x[4], x[5]]])
    V = R * U * R.T
    assert R.T * R == s.eye(3)
    assert s.expand(sum(v * v for v in V) - sum(v * v for v in U)) == 0


@pytest.mark.parametrize(
    "invalid",
    (
        True,
        False,
        1.0,
        s.Float(1),
        "1",
        None,
        s.oo,
        s.I,
        s.nan,
        s.Symbol("missing"),
        0,
        -1,
        s.Rational(1, 4),
    ),
)
def test_energy_contract_rejects_inexact_or_outside_inputs(invalid):
    with pytest.raises((TypeError, ValueError)):
        bounds.require_energy(invalid)
    with pytest.raises((TypeError, ValueError)):
        bounds.hard_current_relative_upper(invalid)
    with pytest.raises((TypeError, ValueError)):
        bounds.selected_pair_upper(invalid, s.Rational(1, 64))


@pytest.mark.parametrize(
    "first,second",
    (
        (s.Rational(1, 8), s.Rational(1, 8)),
        (s.Rational(1, 10), s.Rational(1, 10)),
        (s.Rational(1, 8), s.Rational(1, 10**40)),
    ),
)
def test_total_energy_cut_is_not_two_individual_cuts(first, second):
    with pytest.raises(ValueError):
        bounds.selected_pair_upper(first, second)


@pytest.mark.parametrize(
    "first,second",
    (
        (s.Rational(1, 16), s.Rational(1, 16)),
        (s.Rational(1, 1000), s.Rational(1, 10)),
        (s.Rational(1, 10**800), s.Rational(1, 20)),
    ),
)
def test_selected_pair_bound_keeps_energy_symmetry_and_soft_poles(first, second):
    value = bounds.selected_pair_upper(first, second)
    assert value == bounds.selected_pair_upper(second, first)
    assert value == 9 * (1 / first + 1 / second) / s.sqrt(source.KAPPA)
    assert bounds.selected_pair_upper(first / 2, second / 2) == 2 * value


def test_source_hierarchy_is_original_not_a_new_matching_choice():
    assert source.HEAVY_MASS2 == s.Rational(10**200, 512) + 2
    assert source.KAPPA == 10**800 and source.CUBIC == s.Rational(1, 8192)
    assert source.CONTACT == -(source.CUBIC**2) * (
        3 / (source.HEAVY_MASS2 - 2) - 2 / (source.HEAVY_MASS2 - 2) ** 2
    )


def test_selected_class_bound_does_not_claim_ir_finite_detector_error():
    scope = audit.observable()["not_established"]
    assert "integrated two-real" in scope and "all-N" in scope
    assert len(audit.matching()) == 167 and len(audit.frontier()) == 9
