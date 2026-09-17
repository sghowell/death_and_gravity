"""Whole finite-source identities, graph accounting and failed shortcuts."""

import pytest
import sympy as s
from p8_vacuum_affine_all_multiplicity_tree_source import topology
from p8_vacuum_affine_complete_pair_factorization import (
    audit,
    forests,
    obstruction,
    pairs,
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
def test_all_written_proof_gates(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_original_scope_rejections(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("n,target", tuple(enumerate((7, 47, 387, 3814, 44050))))
def test_closed_and_independent_no_pair_counts(n, target):
    value = forests.no_pair_count(n)
    assert value == forests.independent_no_pair_count(n)
    assert value.subs({forests.C: 1, forests.g: 1}) == target


@pytest.mark.parametrize("n", range(9))
def test_matching_enumeration_and_full_inventory(n):
    rows = forests.labeled_matchings(range(n))
    assert len(rows) == len(set(rows))
    for row in rows:
        flat = [label for pair in row for label in pair]
        assert len(flat) == len(set(flat))
    for k in range(n // 2 + 1):
        assert sum(len(row) == k for row in rows) == forests.matching_count(n, k)
    assert forests.matching_count(n, n + 1) == 0
    assert forests.no_pair_count(n).subs(
        {forests.C: 1, forests.g: 1}
    ) <= topology.topology_closed(n).subs({forests.C: 1, forests.g: 1})


@pytest.mark.parametrize(
    "invalid",
    (
        True,
        False,
        1.0,
        s.Float(1),
        "1",
        None,
        -1,
        s.Rational(1, 2),
        s.oo,
        s.I,
        s.nan,
        s.Symbol("unknown"),
        [],
        {},
        (1,),
    ),
)
def test_multiplicity_rejections(invalid):
    for call in (
        source.require_multiplicity,
        forests.no_pair_count,
        forests.independent_no_pair_count,
    ):
        with pytest.raises((TypeError, ValueError)):
            call(invalid)
    with pytest.raises((TypeError, ValueError)):
        forests.matching_count(4, invalid)


@pytest.mark.parametrize("original", (False, True))
def test_complete_original_and_diagnostic_pair_factorizations(original):
    result = pairs.three_real_calibration(original)
    assert result["whole_counts"] == (5116, 3814, 434, 434, 434)
    assert result["nonzero_pair_pieces"] and result["nonzero_whole_coefficient"]
    assert result["all_currents_symmetric"]
    if original:
        assert result["parameters"] == source.original_parameters()


@pytest.mark.parametrize("invalid", (0, 1, None, "original", s.S.One))
def test_explicit_calibration_selector(invalid):
    with pytest.raises(TypeError):
        pairs.three_real_calibration(invalid)


def test_four_real_intersections_cannot_be_dropped():
    assert 73444 - 6 * 5116 + 3 * 434 == 44050
    assert 73444 - 6 * 5116 != 44050


def test_symbolic_isolated_cluster_pole_and_other_finite_endpoint():
    result = obstruction.cluster_data()
    r = result["parameter"]
    value = result["rational_cluster"]
    assert s.limit(s.factor((r - 1) ** 2 * value), r, 1) == 2
    assert s.limit(value, r, 0) == -s.Rational(32, 3)
    assert result["nested_Ward_defect"] != s.zeros(1, 4)


def test_actual_original_conserved_source_does_not_bound_isolated_cluster():
    result = obstruction.actual_source_data()
    for sample in result["samples"].values():
        assert sample["hard_graph_count"] == 47
        assert sample["isolated_cluster_pole"] != 0
        assert not sample["isolated_cluster_pole"].has(s.Float)


def test_one_offshell_root_not_two_offshell_roots():
    result = obstruction.four_soft_data()
    assert result["checks"]["complete_four_soft_root_Ward"] == s.zeros(1, 4)
    assert result["whole_multi_offshell_Ward_defect"] != s.zeros(1, 4)


def test_conditional_bound_preserves_every_open_frontier():
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.matching()) == 172
    assert len(audit.frontier()) == 9 and len(audit.qualifications()) == 6
    scope = audit.observable()
    assert "completed-remainder" in scope["not_established"]
    assert "all-N" in scope["not_established"]
    assert "Regge" in scope["not_established"]
