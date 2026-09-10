"""Independent full-label cographs, canonical products, LSZ and actual bound tests."""

from itertools import combinations

import pytest
import sympy as sp
from p8_vacuum_two_loop_canonical import (
    assembly,
    audit,
    calibration,
    contractions,
    fields,
    ledger,
)
from p8_vacuum_two_loop_denom import graphs as parent


@pytest.mark.parametrize(
    "name,value", list(audit.residuals().items()), ids=list(audit.residuals())
)
def test_exact_complete_order_identity(name, value):
    assert value == 0, name


@pytest.mark.parametrize(
    "name,value", list(audit.gates().items()), ids=list(audit.gates())
)
def test_actual_bound_or_complete_order_scope(name, value):
    assert bool(value), name


@pytest.mark.parametrize(
    "name,call,args",
    calibration.bad_cases(),
    ids=[r[0] for r in calibration.bad_cases()],
)
def test_unsupported_loop_order_rejected(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize(
    "name,value",
    [(k, v) for k, v in audit.controls().items() if k != "rejected_inputs"],
    ids=[k for k in audit.controls() if k != "rejected_inputs"],
)
def test_negative_control_or_scope(name, value):
    assert bool(value), name


SINGLE = contractions.data()["all_single_core_cograph_weight_matches"]
PAIRS = contractions.data()["all_disjoint_tree_cograph_weight_matches"]


@pytest.mark.parametrize("row", SINGLE, ids=[str(r["index"]) for r in SINGLE])
def test_each_single_core_cograph_preserves_full_counterterm_valences(row):
    _check_vertices(row, 1)
    assert (
        row["contracted_raw_Wick_weight"]
        == row["independent_counterterm_Wick_weight"]
        > 0
    )
    assert row["contracted_signed_coefficient"] == row["independent_signed_coefficient"]
    assert (
        row["contracted_signed_coefficient"]
        == (-1) ** (row["coupling_powers_L_G"][0] + 1)
        * row["contracted_raw_Wick_weight"]
    )


@pytest.mark.parametrize("row", PAIRS, ids=[str(r["index"]) for r in PAIRS])
def test_each_disjoint_cograph_is_actual_heavy_tree_with_two_CT_vertices(row):
    _check_vertices(row, 2)
    assert all(t == "H" for _, _, t in row["typed_cograph_edges"])
    assert (
        row["contracted_raw_Wick_weight"]
        == row["independent_counterterm_Wick_weight"]
        > 0
    )
    assert row["contracted_signed_coefficient"] == row["independent_signed_coefficient"]
    assert (
        row["contracted_signed_coefficient"]
        == (-1) ** row["coupling_powers_L_G"][0] * row["contracted_raw_Wick_weight"]
    )


def _check_vertices(row, core_count):
    vertices = row["vertex_types_and_external_labels"]
    edges = row["typed_cograph_edges"]
    assert sorted(label for _, labels in vertices for label in labels) == [0, 1, 2, 3]
    assert sum(color.startswith("C") for color, _ in vertices) == core_count
    assert len(edges) - len(vertices) + 1 == 2 - core_count
    nL, nG = row["coupling_powers_L_G"]
    assert nL + sp.Rational(nG, 2) == 3
    expected = {"L": (4, 0), "G": (2, 1), "C4": (4, 0), "C3": (2, 1), "C2": (0, 2)}
    for i, (color, labels) in enumerate(vertices):
        phi, H = len(labels), 0
        for u, v, t in edges:
            count = int(u == i) + int(v == i)
            if t == "P":
                phi += count
            else:
                assert t == "H"
                H += count
        assert (phi, H) == expected[color]


def test_internal_relabelling_does_not_change_canonical_key():
    row = {"edges": ((0, 1), (0, 1)), "external_counts": (2, 2)}
    g = contractions.refine(row, ((0, 1), (2, 3)), (0, 2))
    colors, ext, edges = g
    rename = {v: 7 + 3 * v for v in colors}
    moved = (
        {rename[v]: c for v, c in colors.items()},
        {rename[v]: labels for v, labels in ext.items()},
        tuple((rename[v], rename[u], t) for u, v, t in reversed(edges)),
    )
    assert contractions.key(g, (2, 2)) == contractions.key(moved, (2, 2))


def test_external_labels_and_propagator_types_are_not_quotiented_out():
    tree = ({0: "G", 1: "G"}, {0: (0, 1), 1: (2, 3)}, ((0, 1, "H"),))
    crossed = (tree[0], {0: (0, 2), 1: (1, 3)}, tree[2])
    wrong_mass = (tree[0], tree[1], ((0, 1, "P"),))
    assert contractions.key(tree, (0, 2)) != contractions.key(crossed, (0, 2))
    assert contractions.key(tree, (0, 2)) != contractions.key(wrong_mass, (0, 2))


def test_identical_mass_insertions_do_not_supply_an_extra_ordering_factor():
    tree = ({0: "G", 1: "G"}, {0: (0, 1), 1: (2, 3)}, ((0, 1, "H"),))
    direct = contractions.insert_mass(tree, 0, 2)
    sequential = contractions.insert_mass(contractions.insert_mass(tree, 0), 0)
    assert contractions.key(direct, (0, 6)) == contractions.key(sequential, (0, 6))
    rows = [r for r in PAIRS if r["coupling_powers_L_G"] == (0, 6)]
    assert len(rows) == 3
    assert all(r["contracted_signed_coefficient"] == sp.Rational(1, 4) for r in rows)


def test_shared_vertex_proper_cores_are_not_disjoint_pairs():
    row = next(r for r in parent.skeletons() if r["kind"] == "double_bubble")
    g = contractions.refine(
        row, contractions.assignments(row["external_counts"])[0], (0, 0, 0)
    )
    cores = contractions.cores(g)
    assert len(cores) == 2
    assert not set(cores[0][1]).isdisjoint(cores[1][1])
    assert not [
        pair
        for pair in combinations(cores, 2)
        if set(pair[0][1]).isdisjoint(pair[1][1])
    ]


def test_all_actual_raw_refinements_have_one_certified_integrated_owner():
    rows = assembly.data()["exact_raw_refinement_ownership"]
    assert len(rows) == 192
    assert {(r["kind"], r["choices"]) for r in rows} == set(parent.cases())
    assert assembly.data()["raw_group_sizes"] == {
        "finite": 88,
        "insertion": 64,
        "double_bubble": 24,
        "wineglass": 16,
    }


def test_cubic_square_and_mixed_mass_terms_are_both_retained():
    d = ledger.data()
    first = d["first_order_canonical_vertex_counterterms"]
    second = d["second_order_canonical_vertex_counterterms"]
    G = d["G"]
    assert (
        sp.expand(
            d["derived_total_squared_coupling_second_order"]
            - 2 * G * second["G"]
            - first["G"] ** 2
        )
        == 0
    )
    assert d["diagnostic_omitted_cubic_product_b2"] != 0
    assert (
        sp.expand(
            d["full_second_order_tree_counterterm_amplitude"]
            - d["linear_second_order_tree_counterterms"]
            - d["actual_disjoint_first_order_tree_products"]
        )
        == 0
    )


def test_canonical_LSZ_cancels_bare_field_factor_at_both_orders():
    d = fields.data()
    tau = d["formal_loop_parameter"]
    assert (
        fields.truncate(
            d["bare_four_external_leg_LSZ_factor"]
            * d["bare_amputated_amplitude_through_two_loops"],
            tau,
        )
        == d["canonical_amputated_amplitude_through_two_loops"]
    )
    assert d["diagnostic_nonzero_extra_LSZ_error"] != 0
    assert d["actual_inherited_light_pole"]["light_pole_residue_through_two_loops"] == 1


@pytest.mark.parametrize("order", (0, 1, 2))
def test_actual_complete_canonical_interval_at_each_available_order(order):
    d = calibration.point(order)
    tree = calibration.data()["actual_tree_b2"]
    assert d["canonical_b2_lower"] > 0
    assert d["canonical_b2_lower"] <= tree <= d["canonical_b2_upper"]
    assert d["canonical_b2_lower"] + d["canonical_b2_upper"] == 2 * tree
    assert d["strictly_positive_canonical_b2_at_this_order"]
    assert "gravitationally decoupled" in d["scope"]


def test_complete_two_loop_bound_is_actual_exact_rational():
    d = calibration.data()
    E2 = d["actual_complete_two_loop_b2_absolute_upper"]
    assert not E2.free_symbols
    assert 0 < E2 < sp.Rational(3, 10**607)
    assert E2 / d["actual_tree_b2"] < sp.Rational(1, 10**7)
    assert d["actual_one_plus_two_loop_relative_error_upper"] < sp.Rational(1, 10**6)
    assert d["actual_two_loop_canonical_b2_lower"] > 0


def test_equal_inexact_loop_orders_rejected_after_exact_point():
    calibration.point(1)
    for bad in (True, 1.0, sp.Integer(1), sp.Float(1)):
        with pytest.raises((TypeError, ValueError)):
            calibration.point(bad)


def test_complete_scalar_order_two_is_not_source_or_original_P8_completion():
    assert audit.gates()["two_loop_derivative_coordinate_source_matching_not_computed"]
    assert audit.gates()["no_all_higher_loop_or_all_energy_error_bound"]
    assert audit.controls()["full_V_G_B_and_original_P8_not_closed"]


def test_exact_counts():
    assert len(audit.residuals()) == 921
    assert len(audit.gates()) == 37
    assert len(audit.controls()) == 9
    assert audit.rejected_inputs() == 15
