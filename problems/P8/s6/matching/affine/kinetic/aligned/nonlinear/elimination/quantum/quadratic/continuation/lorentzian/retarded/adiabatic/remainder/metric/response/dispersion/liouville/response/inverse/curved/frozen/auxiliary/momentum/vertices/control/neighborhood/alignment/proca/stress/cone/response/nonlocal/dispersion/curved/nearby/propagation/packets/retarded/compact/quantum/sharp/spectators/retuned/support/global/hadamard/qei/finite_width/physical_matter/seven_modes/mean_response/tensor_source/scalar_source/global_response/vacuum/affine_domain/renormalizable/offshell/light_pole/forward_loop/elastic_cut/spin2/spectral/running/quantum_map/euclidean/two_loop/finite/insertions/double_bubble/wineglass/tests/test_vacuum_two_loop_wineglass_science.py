"""Independent all-routing, subtraction, integral and actual-parameter tests."""

from fractions import Fraction

import pytest
import sympy as sp
from p8_vacuum_two_loop_denom import subgraphs
from p8_vacuum_two_loop_insertions import routing
from p8_vacuum_two_loop_wineglass import (
    audit,
    calibration,
    integrals,
    kernel,
    selection,
    subtraction,
)


@pytest.mark.parametrize(
    "name,value", list(audit.residuals().items()), ids=list(audit.residuals())
)
def test_exact_selection_subtraction_routing_integral_or_bound_identity(name, value):
    assert value == 0, name


@pytest.mark.parametrize(
    "name,value", list(audit.gates().items()), ids=list(audit.gates())
)
def test_common_analytic_subtraction_or_original_scope_gate(name, value):
    assert bool(value), name


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[row[0] for row in audit.bad_cases()]
)
def test_reject_unsupported_raw_graph_or_radial_input(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize(
    "name,value",
    [(k, v) for k, v in audit.controls().items() if k != "rejected_inputs"],
    ids=[k for k in audit.controls() if k != "rejected_inputs"],
)
def test_negative_control_or_completion_scope(name, value):
    assert bool(value), name


@pytest.mark.parametrize("choices", selection.cases())
def test_every_selected_graph_exact_local_core_inventory(choices):
    selection.require_group(choices)
    d = subgraphs.data("wineglass", choices)
    overall = choices[0] in (0, 1) and choices[1:] == (0, 0)
    assert len(d["UV_subgraphs"]) == 1 + int(overall)
    assert all(r["superficial_momentum_degree"] == 0 for r in d["UV_subgraphs"])
    assert len(d["restricted_forest_index_sets"]) == (4 if overall else 2)


@pytest.mark.parametrize("choices", selection.cases())
def test_each_forest_replayed_and_nested_overall_term_retained(choices):
    row = next(
        r
        for r in subtraction.data()["all_actual_forest_terms"]
        if r["choices"] == choices
    )
    actual = subgraphs.data("wineglass", choices)
    assert (
        tuple(r["forest"] for r in row["forest_terms"])
        == actual["restricted_forest_index_sets"]
    )
    if row["has_overall_core"]:
        assert (0, 1) in actual["restricted_forest_index_sets"]
        nested = next(r["term"] for r in row["forest_terms"] if r["forest"] == (0, 1))
        assert nested == sp.Symbol("same_regulated_reference") ** 2
        assert row["renormalized_form"] == (
            sp.Symbol("inner_subtracted_wineglass")
            - sp.Symbol("fixed_inner_subtracted_zero_reference")
        )


@pytest.mark.parametrize("row", kernel.data()["all_external_label_routings"])
def test_every_actual_bilinear_routing_and_convex_shift_identity(row):
    old = routing.data()
    ps, x = old["momenta"], old["x"]
    pair = row["outer_external_pair"]
    j = row["inner_external_label"]
    other = next(i for i in range(4) if i not in pair and i != j)
    r = row["inner_shift"]
    assert all(sp.expand(v) == 0 for v in r + (1 - x) * ps[j] - x * ps[other])
    assert sp.simplify(r.dot(r) + row["Delta"]) == 0
    P = ps[pair[0]] + ps[pair[1]]
    assert sp.simplify(P.dot(P) + row["channel_invariant"]) == 0


@pytest.mark.parametrize("y", (0, Fraction(1, 2), 1, 4, 10**200, 10**400))
def test_all_radius_diagnostic_enclosures_without_cutoff(y):
    d = calibration.point(y)
    yy = sp.Rational(y)
    assert d["radial_invariant"] == yy
    assert d["shifted_inner_bubble_modulus_upper"] == (sp.log(1 + yy) + 3) / 144
    assert d["decaying_inner_bubble_difference_modulus_upper"] == (
        8 * sp.sqrt(yy) + 4
    ) / ((yy + 4) * 144)
    assert d["centered_outer_light_inverse_real_lower"] == yy + sp.Rational(1, 4)
    assert (
        d["routed_heavy_inverse_real_lower"]
        == (yy + calibration.data()["actual_heavy_mass_squared"]) / 2
    )


def test_logarithm_difference_has_the_needed_decay_not_just_a_constant_bound():
    d = kernel.data()
    y = d["y"]
    D = d["decaying_logarithm_difference_modulus_upper"]
    assert sp.limit(D, y, sp.oo) == 0
    assert sp.limit(sp.sqrt(y) * D, y, sp.oo) == 8
    assert d["uniform_logarithm_difference_modulus_upper"] == 3


def test_both_independent_radial_primitives_and_endpoints():
    d = integrals.data()
    primitive = d["half_power_radial_primitive_after_sqrt_substitution"]
    t = next(iter(primitive.free_symbols))
    assert (
        sp.factor(sp.diff(primitive, t) - 2 * t**4 / (t**2 + sp.Rational(1, 4)) ** 3)
        == 0
    )
    assert primitive.subs(t, 0) == 0
    assert sp.limit(primitive, t, sp.oo) == 3 * sp.pi / 4
    p = d["integer_radial_primitive"]
    y = next(iter(p.free_symbols))
    assert sp.factor(sp.diff(p, y) - y / (y + sp.Rational(1, 4)) ** 3) == 0
    assert -p.subs(y, 0) == 2
    assert sp.limit(p, y, sp.oo) == 0


def test_all_channel_weights_are_exhaustively_rebuilt():
    d = selection.data()
    assert d["explicit_external_assignment_channel_weights"] == {"s": 1, "t": 1, "u": 1}
    assert d["total_family_Wick_weight"] == 3
    assert d["selected_refinement_count"] == 16
    assert d["overall_subtraction_refinement_count"] == 2
    assert d["finite_after_inner_subtraction_count"] == 14


def test_actual_four_disjoint_raw_groups_and_nonzero_error():
    d = calibration.data()
    E = d["actual_wineglass_group_b2_absolute_upper"]
    assert E == sp.Rational(14496240, 20736) * d["actual_quartic"] ** 3
    assert 0 < E < sp.Rational(1, 10**611)
    assert E / d["actual_tree_b2"] < sp.Rational(3, 10**12)
    combined = d["four_integrated_raw_graph_groups_combined_upper"]
    assert combined == E + d["previous_three_integrated_groups_upper"]
    assert combined < sp.Rational(3, 10**607)
    assert combined / d["actual_tree_b2"] < sp.Rational(1, 10**7)
    assert d["represented_raw_refinement_count"] == 192
    assert d["remaining_raw_refinement_count"] == 0


def test_exact_input_validation_after_equal_previous_calls():
    calibration.point(1)
    selection.require_group((0, 0, 0))
    with pytest.raises(TypeError):
        calibration.point(True)
    with pytest.raises((TypeError, ValueError)):
        selection.require_group((False, 0, 0))


def test_raw_graph_completion_does_not_erase_additional_counterterm_and_LSZ_tasks():
    assert (
        "fixed finite-potential counterterm insertions" in calibration.data()["scope"]
    )
    assert audit.gates()["fixed_finite_potential_counterterm_insertions_still_required"]
    assert audit.gates()["complete_two_loop_light_pole_residue_LSZ_still_required"]
    assert audit.controls()["full_two_loop_and_original_P8_not_closed"]


def test_exact_counts_and_all_external_routings():
    assert len(audit.residuals()) == 215
    assert len(audit.gates()) == 34
    assert len(audit.controls()) == 9
    assert audit.rejected_inputs() == 25
    assert kernel.data()["routing_count"] == 12
