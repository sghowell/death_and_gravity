"""Independent inventory, denominator, integration, subtraction and pole tests."""

from fractions import Fraction
from itertools import combinations

import pytest
import sympy as sp
from p8_vacuum_two_loop_light_pole import (
    analytic,
    audit,
    calibration,
    denominators,
    graphs,
    sectors,
    subtraction,
)


@pytest.mark.parametrize(
    "name,value", list(audit.residuals().items()), ids=list(audit.residuals())
)
def test_exact_two_loop_identity(name, value):
    assert value == 0, name


@pytest.mark.parametrize(
    "name,value", list(audit.gates().items()), ids=list(audit.gates())
)
def test_actual_bound_or_written_scope_gate(name, value):
    assert bool(value), name


BAD = graphs.bad_cases() + calibration.bad_cases()


@pytest.mark.parametrize("name,call,args", BAD, ids=[r[0] for r in BAD])
def test_unsupported_input_rejected(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize(
    "name,value",
    [(k, v) for k, v in audit.controls().items() if k != "rejected_inputs"],
    ids=[k for k in audit.controls() if k != "rejected_inputs"],
)
def test_negative_control(name, value):
    assert bool(value), name


@pytest.mark.parametrize("kind,choice", graphs.cases())
def test_independent_denominator_at_exact_positive_parameters(kind, choice):
    d = denominators.point(kind, choice)
    a = d["parameters"]
    for values in (range(1, len(a) + 1), tuple(2 * i + 1 for i in range(len(a)))):
        assignment = dict(zip(a, values))
        U, P = [d[k].subs(assignment) for k in ("U", "P")]
        total = sum(values)
        assert U > 0 and P >= 0
        assert U * total / 4 - P >= 0
        assert (
            sp.expand(
                (U**2 * total / 4 - U * P)
                - d["resistance_square_polynomial"].subs(assignment)
            )
            == 0
        )


@pytest.mark.parametrize("kind,choice", graphs.cases())
def test_every_reported_forest_has_only_nested_or_disjoint_cores(kind, choice):
    g = graphs.refine(kind, choice)
    uv = graphs.uv(g)
    for forest in graphs.forests(g):
        for i, j in combinations(forest, 2):
            first, second = set(uv[i][0]), set(uv[j][0])
            va = {v for e in first for v in g["edges"][e]}
            vb = {v for e in second for v in g["edges"][e]}
            assert first <= second or second <= first or va.isdisjoint(vb)
    assert () in graphs.forests(g)


@pytest.mark.parametrize("kind,choice", graphs.cases())
def test_two_point_external_labels_and_full_graph_loop_number(kind, choice):
    g = graphs.refine(kind, choice)
    assert sorted(v for row in g["external"].values() for v in row) == [0, 1]
    assert len(g["edges"]) - len(g["vertices"]) + 1 == 2
    assert len(g["edges"]) == 3 + sum(c > 0 for c in choice)


def test_explicit_local_sunset_polynomials_and_overlapping_cores():
    d = denominators.point("sunset", (0, 0))
    a, b, c = d["parameters"]
    assert d["U"] == a * b + a * c + b * c
    assert d["P"] == a * b * c
    g = graphs.refine("sunset", (0, 0))
    assert len(graphs.uv(g)) == 4
    assert len(graphs.forests(g)) == 8
    assert (0, 1) not in graphs.forests(g)


def test_all_finite_sunset_orders_and_local_twice_differentiated_integral():
    d = sectors.data()
    rows = d["six_UV_finite_sunset_ordered_bounds"]
    assert len(rows) == 6
    assert sum(r["orders"] for r in rows) == 720
    assert all(r["unit_cube_integral_upper"] == 14 for r in rows)
    assert all(r["minimum_tail_exponent"] == 1 for r in rows)
    local = d["local_twice_differentiated_sunset_bound"]
    assert local["orders"] == 6 and local["unit_cube_integral_upper"] == 1
    assert local["ordered_integral_patterns"] == (((-2, -2, 2), (1, 2, 3)),)


@pytest.mark.parametrize("radius", (Fraction(1, 100), Fraction(1, 2), 1))
def test_actual_subdisc_unique_pole_and_unit_residue(radius):
    d = calibration.point(radius)
    bound = calibration.data()["actual_one_plus_two_loop_quadratic_coefficient_upper"]
    assert d["inverse_propagator_factored_error_upper"] == sp.Rational(radius) * bound
    assert d["positive_inverse_factor_modulus_lower"] > 0
    assert d["unique_mass_one_pole_through_two_loops"]
    assert d["light_pole_residue_through_two_loops"] == 1


def test_concrete_actual_six_group_bound():
    d = calibration.data()
    assert d["actual_cubic_squared"] == sp.Rational(1, 67108864)
    groups = d["actual_two_loop_group_quadratic_coefficient_uppers"]
    assert len(groups) == 6 and all(v > 0 for v in groups.values())
    B2 = d["actual_two_loop_quadratic_coefficient_upper"]
    assert B2 == sum(groups.values())
    assert sp.Rational(3, 10**19) < B2 < sp.Rational(4, 10**19)
    assert d["actual_one_plus_two_loop_quadratic_coefficient_upper"] < sp.Rational(
        1, 10**18
    )
    assert not B2.free_symbols


def test_actual_proper_counterterms_cannot_be_dropped():
    d = subtraction.data()
    bare = d["sunset_proper_UV_piece_before_subtraction"]
    ct = d["complete_momentum_dependent_parameter_counterterm_piece"]
    assert bare != 0 and ct != 0
    assert sp.expand(bare + ct) == 0
    assert d["complete_local_parameter_counterterm_piece"] != 0


def test_local_on_shell_operation_has_both_required_zeros():
    s = sp.Symbol("s")
    expr = 7 + 3 * s + s**2 + 2 * s**3
    projected = subtraction.on_shell(expr, s)
    assert projected.subs(s, 1) == 0
    assert sp.diff(projected, s).subs(s, 1) == 0
    assert sp.expand(projected - (s - 1) ** 2 * (2 * s + 5)) == 0


def test_anchored_log_difference_decays_and_has_real_positive_bound():
    d = analytic.data()
    y = d["y"]
    D = d["anchored_log_difference_majorant"]
    assert sp.limit(D, y, sp.oo) == 0
    for r in (0, Fraction(1, 4), 1, 4, 10**6):
        assert 0 < D.subs(y, sp.Rational(r)) < 7


def test_entire_routing_remains_on_exact_complex_mass_shell():
    d = analytic.data()
    s = d["s"]
    p = d["entire_external_routing"]
    for value in (1, 3, -1, 1 + 2 * sp.I):
        pv = p.subs(s, value)
        assert sp.expand((pv.T * pv)[0] + value) == 0
        assert sp.simplify((sp.conjugate(pv).T * pv)[0]) <= 3


def test_equal_inexact_graph_choices_remain_rejected_after_cache_population():
    graphs.refine("sunset", (0, 1))
    denominators.point("sunset", (0, 1))
    for bad in (True, 1.0, sp.Integer(1)):
        with pytest.raises((TypeError, ValueError)):
            graphs.refine("sunset", (0, bad))
        with pytest.raises((TypeError, ValueError)):
            denominators.point("sunset", (0, bad))


def test_equal_inexact_radius_remains_rejected_after_exact_point():
    calibration.point(1)
    for bad in (True, 1.0, sp.Float(1)):
        with pytest.raises((TypeError, ValueError)):
            calibration.point(bad)


def test_finite_order_pole_does_not_close_four_point_matching_or_original_P8():
    assert "fixed-order" in subtraction.data()["scope"]
    assert audit.gates()["complete_two_loop_four_point_ledger_not_inferred_here"]
    assert audit.gates()["two_loop_source_aware_derivative_map_not_computed"]
    assert audit.controls()["full_V_G_B_and_original_P8_not_closed"]


def test_exact_counts():
    assert len(audit.residuals()) == 743
    assert len(audit.gates()) == 32
    assert len(audit.controls()) == 9
    assert audit.rejected_inputs() == 33
