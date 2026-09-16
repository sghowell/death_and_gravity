"""Exact tree, independent action/topology and original physical scope checks."""

from itertools import permutations

import pytest
import sympy as s
from p8_vacuum_affine_complete_two_graviton_tree import (
    audit,
    checks,
    source,
    topology,
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
def test_written_proof_gates(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_original_scope_rejections(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("count,bell", ((1, 1), (2, 2), (3, 5), (4, 15), (5, 52)))
def test_literal_set_partitions_are_complete_and_disjoint(count, bell):
    mask = (1 << count) - 1
    rows = trees.partitions(mask)
    assert len(rows) == bell and len({frozenset(r) for r in rows}) == bell
    for row in rows:
        assert all(x > 0 for x in row) and sum(row) == mask
        assert all(a & b == 0 for i, a in enumerate(row) for b in row[i + 1 :])


@pytest.mark.parametrize(
    "number,expected", ((0, (1, 3, 3)), (1, (5, 21, 21)), (2, (38, 198, 198)))
)
def test_independent_three_sector_graph_counts(number, expected):
    value = topology.counts()[number]
    assert (
        value.coeff(topology.C),
        value.coeff(topology.g, 2),
        value.subs({topology.C: 0, topology.g: 0}),
    ) == expected


@pytest.mark.parametrize("permutation", tuple(permutations(range(3))))
def test_mixed_scalar_metric_jet_is_fully_symmetric(permutation):
    fields = (
        trees.imm(s.diag(1, 2, 0, -1)),
        trees.imm([[0, 1, 0, 0], [1, 0, 1, 0], [0, 1, 0, 1], [0, 0, 1, 2]]),
        trees.imm(s.diag(2, 0, 1, 1)),
    )
    p, q = trees.imm([2, 1, 0, 1]), trees.imm([3, -1, 2, 0])
    shuffled = tuple(fields[i] for i in permutation)
    assert trees.scalar_vertex(p, q, fields, 7) == trees.scalar_vertex(
        p, q, shuffled, 7
    )
    assert trees.determinant_coefficient(fields) == trees.determinant_coefficient(
        shuffled
    )


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
        -s.oo,
        s.nan,
        s.I,
        s.Symbol("missing"),
    ),
)
def test_inexact_or_unknown_physical_entries_rejected(invalid):
    with pytest.raises((TypeError, ValueError)):
        trees.exact_real(invalid)
    with pytest.raises((TypeError, ValueError)):
        trees.exact_array([1, invalid, 0, 0])


@pytest.mark.parametrize(
    "mutation",
    (
        "wrong_mass",
        "wrong_null",
        "negative_radiation",
        "wrong_conservation",
        "trace",
        "time_component",
        "nontransverse",
        "nonsymmetric",
        "three_real",
        "missing_scalar",
        "inexact_momentum",
        "string_momentum",
    ),
)
def test_original_physical_point_api_rejects_outside_contract(mutation):
    ps, qs, eps = checks.nonopposite_state()
    ps = [s.Matrix(p) for p in ps]
    qs = [s.Matrix(q) for q in qs]
    eps = [s.Matrix(a) for a in eps]
    if mutation == "wrong_mass":
        ps[2][0] += 1
    elif mutation == "wrong_null":
        qs[0][1] += 1
    elif mutation == "negative_radiation":
        qs[0] = -qs[0]
    elif mutation == "wrong_conservation":
        qs[0] = qs[0] / 2
    elif mutation == "trace":
        eps[0][2, 2] += 1
    elif mutation == "time_component":
        eps[0][0, 0] = 1
    elif mutation == "nontransverse":
        eps[0][1, 2] = eps[0][2, 1] = 1
    elif mutation == "nonsymmetric":
        eps[0][1, 2] = 1
    elif mutation == "three_real":
        qs.append(qs[0])
        eps.append(eps[0])
    elif mutation == "missing_scalar":
        ps.pop()
    elif mutation == "inexact_momentum":
        ps[0][0] = s.Float(ps[0][0])
    elif mutation == "string_momentum":
        ps[0] = ps[0].tolist()
        ps[0][0][0] = str(ps[0][0][0])
    with pytest.raises((TypeError, ValueError)):
        trees.original_amplitude(ps, list(zip(qs, eps)))


def test_exactly_collinear_two_null_point_is_not_silently_evaluated():
    E = s.Rational(5, 4)
    a = s.Rational(19, 10)
    Ep = (a + 1 / a) / 2
    w = (E - Ep**2 / E) / 2
    q = trees.imm([w, 0, 0, w])
    u = s.Matrix([s.Rational(3, 5), s.Rational(4, 5), 0])
    points, _, _ = source.recoil.momenta(E, [q, q], u)
    eps = trees.imm(s.diag(0, 1, -1, 0))
    with pytest.raises(ValueError, match="collinear"):
        trees.original_amplitude(points, [(q, eps), (q, eps)])


def test_original_four_point_full_Born_calibration():
    _points, _q, born = trees.old.sample(0)
    value, count = trees.original_amplitude(born, [])
    from p8_vacuum_affine_massive_scalar_graviton_bremsstrahlung import tree

    expected = (
        tree.born_continuation(born, source.HEAVY_MASS2, source.CUBIC, source.CONTACT)
        + trees.old.born(born) / source.KAPPA
    )
    assert count == 7 and s.factor(value - expected) == 0


def test_current_cache_is_owned_by_one_finite_evaluator():
    points, qs, eps = checks.nonopposite_state()
    legs = [("phi", p, 1) for p in points[1:]] + [
        ("h", qs[0], eps[0]),
        ("h", qs[1], eps[1]),
    ]
    first, second = trees.TreeEngine(legs), trees.TreeEngine(legs)
    assert first._memo is not second._memo
    first.momentum(1)
    assert "momentum" in first._memo and second._memo == {}


def test_wrong_required_vertices_are_rejected():
    _rows, gates = checks.six_points()
    assert gates["wrong_Einstein_quartic_breaks_Ward"]
    assert gates["wrong_scalar_third_metric_breaks_Ward"]


def test_full_marked_state_not_massive_only_soft_current():
    _rows, gates = checks.original_soft()
    assert gates["original_simultaneous_two_soft_calibration"]
    assert gates["original_state_correct_hierarchical_calibration"]
    assert gates["massive_only_hierarchical_current_negative_control"]
