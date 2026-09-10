"""Independent topology, exact polynomial and compact point controls."""

from fractions import Fraction
from itertools import combinations

import pytest
import sympy as sp
from p8_vacuum_two_loop_denom import (
    audit,
    calibration,
    forests,
    graphs,
    positive,
    subgraphs,
)


@pytest.mark.parametrize(
    "name,value", list(audit.residuals().items()), ids=list(audit.residuals())
)
def test_exact_graph_or_polynomial_identity(name, value):
    assert value == 0, name


@pytest.mark.parametrize(
    "name,value", list(audit.gates().items()), ids=list(audit.gates())
)
def test_continuous_domain_or_scope_gate(name, value):
    assert bool(value), name


@pytest.mark.parametrize(
    "name,call,args",
    calibration.bad_cases(),
    ids=[row[0] for row in calibration.bad_cases()],
)
def test_reject_unsupported_two_loop_input(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize(
    "name,value",
    [(k, v) for k, v in audit.controls().items() if k != "rejected_inputs"],
    ids=[k for k in audit.controls() if k != "rejected_inputs"],
)
def test_nontrivial_negative_control(name, value):
    assert bool(value), name


@pytest.mark.parametrize(
    "kind,choices", graphs.cases(), ids=[audit.label(*case) for case in graphs.cases()]
)
def test_every_refinement_independent_positive_witness_and_subgraphs(kind, choices):
    d = forests.data(kind, choices)
    p = positive.data(kind, choices)
    alpha = d["parameters"]
    assert sp.expand(d["matrix_tree_U"] - d["U"]) == 0
    witness = (
        p["explicit_tadpole_base"]
        + sum(
            w * alpha[k] * (alpha[i] - alpha[j]) ** 2
            for w, k, i, j in p["pair_square_terms"]
        )
        + p["nonnegative_monomial_remainder"]
    )
    assert (
        sp.expand(d["coarse_mass_sixteen_F"] - d["U"] * sum(alpha) / 4 - witness) == 0
    )
    assert all(w > 0 for w, k, i, j in p["pair_square_terms"])
    assert all(
        c >= 0 for _, c in sp.Poly(p["nonnegative_monomial_remainder"], *alpha).terms()
    )
    assert all(
        c >= 0
        for _, c in sp.Poly(p["nonnegative_base_monomial_remainder"], *alpha).terms()
    )
    assert all(
        c >= 0 for _, c in sp.Poly(d["derivative_majorant_remainder"], *alpha).terms()
    )
    g = d["graph"]
    assert len(g["edges"]) - len(g["vertices"]) + 1 == 2
    assert g["heavy_edge_count"] == sum(value != 0 for value in choices)
    uv = subgraphs.data(kind, choices)
    for row in uv["UV_subgraphs"]:
        assert (
            row["superficial_momentum_degree"]
            == 4 * row["loops"] - 2 * len(row["edges"])
            >= 0
        )
        assert row["light_external_legs"] >= 0
        assert row["heavy_external_legs"] >= 0
    for forest in uv["restricted_forest_index_sets"]:
        for i, j in combinations(forest, 2):
            a, b = uv["UV_subgraphs"][i], uv["UV_subgraphs"][j]
            assert (
                set(a["vertices"]).isdisjoint(b["vertices"])
                or set(a["edges"]) < set(b["edges"])
                or set(b["edges"]) < set(a["edges"])
            )


@pytest.mark.parametrize(
    "kind,choices", graphs.cases(), ids=[audit.label(*case) for case in graphs.cases()]
)
def test_every_refinement_complex_disc_diagnostic_mass_boundary(kind, choices):
    size = len(graphs.refinement(kind, choices)["edges"])
    d = calibration.point(
        kind,
        choices,
        (sp.Rational(1, size),) * size,
        Fraction(3, 5),
        Fraction(4, 5),
        16,
    )
    assert d["U"] > 0
    assert d["real_F"] >= d["proved_real_F_lower"] > 0
    assert 0 <= d["relative_crossing_derivative_upper"] <= 4
    assert d["actual_parent_mass"] is False


@pytest.mark.parametrize("kind", graphs.KINDS)
@pytest.mark.parametrize("choices", ((0, 0, 0), (3, 3, 3)))
def test_actual_mass_compact_points_not_just_diagnostic_mass(kind, choices):
    size = len(graphs.refinement(kind, choices)["edges"])
    d = calibration.point(kind, choices, (sp.Rational(1, size),) * size, -1, 0)
    assert bool(d["actual_parent_mass"])
    assert d["real_F"] >= d["proved_real_F_lower"] > 0
    assert 0 <= d["relative_crossing_derivative_upper"] <= 4


@pytest.mark.parametrize("kind", graphs.KINDS)
def test_approach_parameter_boundary_without_asserting_boundary_convergence(kind):
    d = calibration.point(
        kind,
        (0, 0, 0),
        (sp.Rational(1, 10**12),) * 3 + (1 - sp.Rational(3, 10**12),),
        0,
        1,
        16,
    )
    assert d["real_F"] >= d["proved_real_F_lower"] > 0


def test_exhaustive_graph_family_and_label_weights():
    rows = graphs.skeletons()
    assert len(rows) == 9
    assert all(sum(r["external_counts"]) == 4 for r in rows)
    assert sum(r["external_label_assignments"] for r in rows) == 72
    assert graphs.data()["family_Wick_weights"] == {
        "double_bubble": sp.Rational(3, 4),
        "wineglass": sp.Integer(3),
        "tadpole_insertion": sp.Rational(3, 2),
    }


def test_double_bubble_vertex_overlap_is_not_independent():
    d = subgraphs.data("double_bubble", (0, 0, 0))
    a, b = d["UV_subgraphs"][:2]
    assert set(a["edges"]).isdisjoint(b["edges"])
    assert not set(a["vertices"]).isdisjoint(b["vertices"])
    assert d["restricted_forest_index_sets"] == ((), (0,), (1,), (2,), (0, 2), (1, 2))


def test_full_UV_type_inventory():
    types = {
        row for d in subgraphs.summary() for row in d["UV_external_types_and_degrees"]
    }
    assert types == {
        (0, 1, 2, 1),
        (0, 2, 0, 1),
        (0, 2, 0, 2),
        (2, 0, 0, 1),
        (2, 0, 2, 1),
        (2, 1, 0, 1),
        (2, 1, 0, 2),
        (4, 0, 0, 1),
        (4, 0, 0, 2),
    }
    assert {d["restricted_forest_count"] for d in subgraphs.summary()} == {1, 2, 4, 6}


@pytest.mark.parametrize("module", (graphs, forests, positive, subgraphs))
def test_validation_not_bypassed_by_equal_cached_boolean_or_sympy_integer(module):
    call = module.refinement if module is graphs else module.data
    call("double_bubble", (1, 0, 0))
    for bad in (True, 1.0, sp.Integer(1)):
        with pytest.raises((TypeError, ValueError)):
            call("double_bubble", (bad, 0, 0))


def test_denominator_bounds_are_not_an_integrated_two_loop_result():
    assert (
        "not an integrated"
        in calibration.point("double_bubble", (0, 0, 0), (sp.Rational(1, 4),) * 4)[
            "scope"
        ]
    )
    assert "No subtraction integral" in subgraphs.data("wineglass", (0, 0, 0))["scope"]


def test_exact_counts_and_original_scope():
    assert len(audit.residuals()) == 1449
    assert len(audit.gates()) == 31
    assert len(audit.controls()) == 11
    assert audit.rejected_inputs() == 50
