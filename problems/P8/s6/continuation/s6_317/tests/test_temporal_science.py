"""Exact full-source identities, nonlinear controls and the bounded planar domain."""

import pytest
import sympy as s
from p8_vacuum_affine_temporal_tree_reorganization import (
    audit,
    bounds,
    chart,
    source,
    trees,
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
def test_all_written_proof_gates(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_original_scope_rejections(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("original", (False, True))
def test_complete_original_and_diagnostic_amplitudes(original):
    result = trees.complete_calibration(original)
    assert result["graph_count"] == 5116 and result["nonzero_amplitude"]
    assert all(value == 0 for value in result["checks"].values())
    if original:
        assert result["parameters"] == source.original_parameters()


@pytest.mark.parametrize("invalid", (0, 1, None, "original", s.S.One))
def test_calibration_selector_validated_before_cache(invalid):
    with pytest.raises(TypeError):
        trees.complete_calibration(invalid)


@pytest.mark.parametrize(
    "entry",
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
        s.zoo,
        s.Symbol("unknown"),
    ),
)
def test_inexact_or_unspecified_formal_array_entries_rejected(entry):
    with pytest.raises((TypeError, ValueError)):
        chart.formal_array([[entry]])


def test_exact_real_rational_function_array_is_not_numeric_point_api():
    r = s.Symbol("r", real=True)
    assert chart.formal_array([1 / (1 + r * r), s.sqrt(2)]) == s.Matrix(
        [1 / (1 + r * r), s.sqrt(2)]
    )
    with pytest.raises(ValueError):
        chart.formal_array([s.sin(r)])


@pytest.mark.parametrize(
    "Q",
    (
        [0, 1, 0, 0],
        [-1, 0, 0, 0],
        [1, 0, 0],
        [[1, 0, 0, 0]],
        [s.Symbol("W", real=True), 0, 0, 0],
    ),
)
def test_temporal_chart_requires_positive_energy_four_vector(Q):
    with pytest.raises((TypeError, ValueError)):
        chart.require_energy(Q)


@pytest.mark.parametrize("Q", ([3, 1, 1, 1], [2, 0, 0, 2], [5, 0, 0, 0]))
def test_projection_temporal_idempotent_and_covariant_gauge_difference(Q):
    Q = s.Matrix(Q)
    H = s.Matrix([[2, 1, 3, 4], [1, 5, 6, 7], [3, 6, 8, 9], [4, 7, 9, 10]])
    T = chart.project(H, Q)
    assert T[0, :] == s.zeros(1, 4) and T == T.T
    assert chart.project(T, Q) == T
    Z = trees.pair_shift(H, Q)
    q = chart.ETA * Q
    assert (H - T - (q * (chart.ETA * Z).T + (chart.ETA * Z) * q.T) / 2).applyfunc(
        s.factor
    ) == s.zeros(4)


def test_nonlinear_compensation_is_essential():
    result = trees.cluster_calibration()
    assert result["whole_linear_only_pole"] != s.zeros(4)
    assert result["checks"][
        "nonlinear_temporal_projection_removes_known_pole"
    ] == s.zeros(4)
    assert result["whole_nonlinear_limit"][2, 2] == s.Rational(71, 2)


def test_general_literal_chart_checks_all_four_leaf_subsets():
    metric, momenta, shifts, mapped = chart.four_literal()
    assert len(shifts) == 11 and len(mapped) == 16
    assert all(A[0, :] == s.zeros(1, 4) for mask, A in mapped.items() if mask)
    assert mapped[15] != 2 * chart.project(metric[15] / 2, momenta[15])


def test_literal_chart_requires_complete_additive_data():
    q = s.Matrix([1, 0, 0, 1])
    metric = {0: chart.ETA, 1: s.zeros(4), 2: s.zeros(4), 3: s.zeros(4)}
    momenta = {0: s.zeros(4, 1), 1: q, 2: q, 3: 2 * q}
    _, mapped = chart.temporal_chart(metric, momenta, 2)
    assert mapped[3] == s.zeros(4)
    with pytest.raises(ValueError):
        chart.temporal_chart({0: chart.ETA}, momenta, 2)
    bad = dict(momenta)
    bad[3] = 3 * q
    with pytest.raises(ValueError, match="additive"):
        chart.temporal_chart(metric, bad, 2)
    bad = dict(metric)
    bad[0] = s.eye(4)
    with pytest.raises(ValueError, match="background"):
        chart.temporal_chart(bad, momenta, 2)


@pytest.mark.parametrize("invalid", (True, 1.0, "3", None, -1, s.Rational(3, 2)))
def test_finite_chart_multiplicity_rejections(invalid):
    with pytest.raises((TypeError, ValueError)):
        chart.temporal_chart({}, {}, invalid)


@pytest.mark.parametrize("a,b", ((1, 1), (1, 2), (3, 5), (s.Rational(1, 100), 2)))
def test_pair_energy_hierarchy_is_explicit(a, b):
    assert bounds.pair_majorant(a, b, 16) == 530 * (s.sympify(a) + b) ** 2 / (4 * a * b)
    assert bounds.pair_majorant(a, b, 16) == bounds.pair_majorant(b, a, 16)


@pytest.mark.parametrize(
    "invalid", (True, False, 1.0, s.Float(1), "1", None, 0, -1, s.oo, s.I)
)
def test_norm_bound_parameter_rejections(invalid):
    with pytest.raises((TypeError, ValueError)):
        bounds.pair_majorant(invalid, 1)
    with pytest.raises((TypeError, ValueError)):
        bounds.planar_majorant(invalid)


@pytest.mark.parametrize("index,target", tuple(enumerate(bounds.PLANAR_BUDGETS)))
def test_all_eight_exact_planar_budgets(index, target):
    label = format(index, "03b")
    result = bounds.planar_data()["whole_planar_coefficient_budgets"][label]
    assert result["budget"] == target < 120000
    assert sum(row.get("budget", 0) for row in result["entries"].values()) == target


def test_original_coupling_restoration_is_exact():
    assert bounds.planar_majorant(source.KAPPA) * source.KAPPA == s.Rational(
        25982722, 225
    )
    assert bounds.planar_majorant(16) == bounds.planar_majorant(1) / 16


def test_planar_support_checker_rejects_missing_angular_vanishing_and_new_poles():
    r, t = s.symbols("r t", real=True)
    D = 8 * r * r - 12 * r * t + 9 * t * t + 5 * r * r * t * t
    value, _ = bounds.coefficient_budget((r * r + t * t) / D, D, r, t)
    assert value == 1
    with pytest.raises(ValueError, match="quadratic"):
        bounds.coefficient_budget(1 / D, D, r, t)
    with pytest.raises(ValueError, match="denominator"):
        bounds.coefficient_budget(r * r / ((r - t) ** 2 * D), D, r, t)


@pytest.mark.parametrize(
    "bits", ((True, 0, 0), (1.0, 0, 0), (0, 0), (0, 0, 2), "000", None)
)
def test_planar_basis_selector_rejections(bits):
    with pytest.raises(ValueError):
        bounds.planar_configuration(s.Rational(1, 2), -s.Rational(1, 3), bits)


def test_common_energy_and_canonical_coupling_homogeneity_at_a_generic_point():
    hs = bounds.planar_configuration(s.Rational(1, 2), -s.Rational(1, 3), (0, 1, 0))
    legs = [("h", q, A) for q, A in hs]
    H, n = trees.TemporalSoft(legs).current(7, "h")
    scaled = [("h", s.Rational(3, 4) * q, A) for q, A in hs]
    G, m = trees.TemporalSoft(scaled, kappa=16).current(7, "h")
    assert n == m == 4 and (16 * G - H).applyfunc(s.factor) == s.zeros(4)


def test_restricted_bound_keeps_original_frontiers_open():
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.matching()) == 173
    assert len(audit.frontier()) == 9 and len(audit.qualifications()) == 6
    text = audit.observable()["not_established"]
    assert all(
        word in text
        for word in ("three-dimensional", "energy-hierarchy", "all-N", "Regge")
    )
