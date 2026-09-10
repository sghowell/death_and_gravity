"""Exact selected graphs, legal forests, finite factors and actual bounds."""

from fractions import Fraction

import pytest
import sympy as sp
from p8_vacuum_two_loop_denom import subgraphs
from p8_vacuum_two_loop_double_bubble import (
    audit,
    bubble,
    calibration,
    forests,
    selection,
    triangle,
)


@pytest.mark.parametrize(
    "name,value", list(audit.residuals().items()), ids=list(audit.residuals())
)
def test_exact_selection_forest_factor_or_calibration_identity(name, value):
    assert value == 0, name


@pytest.mark.parametrize(
    "name,value", list(audit.gates().items()), ids=list(audit.gates())
)
def test_analytic_subtraction_or_original_scope_gate(name, value):
    assert bool(value), name


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[row[0] for row in audit.bad_cases()]
)
def test_reject_unsupported_graph_or_channel_input(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize(
    "name,value",
    [(k, v) for k, v in audit.controls().items() if k != "rejected_inputs"],
    ids=[k for k in audit.controls() if k != "rejected_inputs"],
)
def test_negative_control_or_scope(name, value):
    assert bool(value), name


@pytest.mark.parametrize("choices", selection.cases())
def test_every_selected_graph_has_only_the_claimed_local_logarithmic_cores(choices):
    selection.require_group(choices)
    d = subgraphs.data("double_bubble", choices)
    assert d["UV_subgraphs"]
    assert all(s["superficial_momentum_degree"] == 0 for s in d["UV_subgraphs"])
    assert all(
        (s["light_external_legs"], s["heavy_external_legs"]) in ((0, 2), (2, 1), (4, 0))
        for s in d["UV_subgraphs"]
    )
    assert len(d["restricted_forest_index_sets"]) in (2, 4, 6)


@pytest.mark.parametrize("choices", selection.cases())
def test_every_actual_forest_sum_and_no_illegal_overlap_pair(choices):
    d = forests.data()
    row = next(r for r in d["all_actual_forest_terms"] if r["choices"] == choices)
    original = subgraphs.data("double_bubble", choices)
    assert (
        tuple(r["forest"] for r in row["terms"])
        == original["restricted_forest_index_sets"]
    )
    assert sp.expand(sum(r["term"] for r in row["terms"]) - row["sum"]) == 0
    I, I0 = d["I"], d["I0"]
    if row["mode"] != "single":
        assert sp.expand(row["sum"] - (I - I0) ** 2) == 0
    if row["mode"] == "overlap":
        assert (0, 1) not in original["restricted_forest_index_sets"]
        # The two overlapping subtractions alone are not the finite square.
        assert sp.expand((I - I0) ** 2 - (I**2 - 2 * I0 * I)) == I0**2
    if row["mode"] == "disjoint":
        assert (0, 1) in original["restricted_forest_index_sets"]


@pytest.mark.parametrize("s", (-3, -1, 0, Fraction(1, 2), 1, 2, 3))
def test_exact_diagnostic_channel_api_and_nonzero_finite_factor_bounds(s):
    d = calibration.point(s)
    p = calibration.data()
    assert d["channel_invariant"] == sp.Rational(s)
    C = -p["actual_quartic"] + p["actual_cubic_squared"] / (
        p["actual_heavy_mass_squared"] - s
    )
    assert d["actual_external_pair_vertex"] == C
    assert abs(C) < d["external_pair_vertex_modulus_upper"]
    assert d["finite_bubble_modulus_upper"] == sp.Rational(1, 72)
    assert d["each_summed_heavy_triangle_modulus_upper"] > 0


def test_finite_bubble_zero_channel_and_logarithm_anchor():
    d = bubble.data()
    assert d["finite_light_bubble"].subs(d["s"], 0).doit() == 0
    assert d["complete_complex_logarithm_upper"] == 2
    assert d["positive_exp_two_partial_sum"] > 4
    assert all(v == 0 for v in d["checks"].values())


def test_triangle_primitive_including_infinite_endpoint():
    d = triangle.data()
    y, M, delta = d["y"], d["M"], d["delta"]
    primitive = d["anchored_radial_primitive"]
    assert sp.factor(sp.diff(primitive, y) - y / ((y + delta) ** 2 * (y + M))) == 0
    assert sp.limit(primitive, y, sp.oo) == 0
    assert (
        sp.simplify(-primitive.subs(y, 0) - d["exact_zero_to_infinity_radial_integral"])
        == 0
    )
    assert d["checks"]["conservative_mass_squared_ratio_margin_at_thirty_two"] == 0


def test_actual_partial_bound_and_all_remaining_normalization_terms():
    d = calibration.data()
    E = d["actual_double_bubble_group_b2_absolute_upper"]
    assert E == sp.Rational(19224, 20736) * d["actual_quartic"] ** 3
    assert 0 < E < sp.Rational(2, 10**614)
    assert E / d["actual_tree_b2"] < sp.Rational(3, 10**15)
    combined = d["three_integrated_groups_combined_upper"]
    assert combined == E + d["previous_two_integrated_groups_upper"]
    assert combined < sp.Rational(3, 10**607)
    assert combined / d["actual_tree_b2"] < sp.Rational(1, 10**7)
    assert d["represented_raw_refinement_count"] == 176
    assert d["remaining_raw_wineglass_refinement_count"] == 16
    assert "additional finite-potential counterterm insertions" in d["scope"]


def test_exact_input_validation_after_equal_cached_or_previous_inputs():
    calibration.point(1)
    selection.require_group((0, 0, 0))
    with pytest.raises(TypeError):
        calibration.point(True)
    with pytest.raises((TypeError, ValueError)):
        selection.require_group((False, 0, 0))


def test_selection_counts_and_finite_refinements_not_double_counted():
    d = selection.data()
    assert d["selected_refinement_count"] == 24
    assert d["forest_class_counts"] == {"single": 16, "disjoint": 4, "overlap": 4}
    assert d["already_finite_separable_count"] == 8
    assert set(selection.cases()).isdisjoint(
        d["already_finite_separable_refinements_excluded"]
    )
    assert len(d["remaining_wineglass_refinements"]) == 16


def test_exact_counts_and_original_problem_open():
    assert len(audit.residuals()) == 142
    assert len(audit.gates()) == 34
    assert len(audit.controls()) == 9
    assert audit.rejected_inputs() == 26
    assert audit.controls()["full_two_loop_and_original_P8_not_closed"]
