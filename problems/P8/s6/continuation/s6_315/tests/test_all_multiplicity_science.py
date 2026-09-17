"""Exact source, all-order majorant and finite physical-domain tests."""

from itertools import permutations

import pytest
import sympy as s
from p8_vacuum_affine_all_multiplicity_tree_source import (
    audit,
    jets,
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
def test_all_written_proof_gates(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_original_scope_rejections(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize(
    "n,expected",
    (
        (0, (1, 3, 3)),
        (1, (5, 21, 21)),
        (2, (38, 198, 198)),
        (3, (388, 2364, 2364)),
        (4, (4972, 34236, 34236)),
    ),
)
def test_closed_and_independent_labeled_topology(n, expected):
    got = topology.topology_closed(n)
    assert got == topology.count_root_cuts(n)
    assert (
        got.coeff(topology.C),
        got.coeff(topology.g, 2),
        got.subs({topology.C: 0, topology.g: 0}),
    ) == expected
    assert got.subs({topology.C: 1, topology.g: 1}) < topology.graph_majorant(n)


@pytest.mark.parametrize("n", range(13))
def test_exact_positive_vertex_majorants(n):
    b = jets.budgets(n)
    assert all(v.is_Integer is True and v >= 0 for v in b.values())
    assert b["determinant"] == s.factorial(n) * 2**n * (n + 1)
    assert b["density_inverse_operator"] == s.factorial(n) * 2**n * s.binomial(n + 2, 2)
    if n >= 3:
        assert b["Einstein_per_L_squared"] <= s.factorial(n) * 32**n


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
        -s.oo,
        s.I,
        s.nan,
        s.Symbol("unknown"),
        [],
        {},
        (1,),
    ),
)
def test_multiplicity_domain_rejected(invalid):
    for call in (
        source.require_multiplicity,
        topology.topology_closed,
        topology.count_root_cuts,
        topology.graph_majorant,
        jets.budgets,
    ):
        with pytest.raises((TypeError, ValueError)):
            call(invalid)


@pytest.mark.parametrize("n", range(9))
def test_nonnegative_multiplicity_includes_zero(n):
    assert source.require_multiplicity(n) == n
    assert source.require_multiplicity(s.Integer(n)) == n


@pytest.mark.parametrize("order", tuple(permutations(range(4))))
def test_new_metric_coefficient_has_full_permutation_symmetry(order):
    fs, _ = jets.calibration_fields()
    shuffled = tuple(fs[i] for i in order)
    assert jets.MetricJet(shuffled).det(15) == 620


@pytest.mark.parametrize("r", range(1, 8))
def test_no_finite_graviton_order_cutoff_in_matter_dispatch(r):
    engine = trees.TreeEngine([], cubic=1, contact=1)
    assert engine.permitted(("phi", "phi", *(("h",) * r)))
    assert engine.permitted(("H", "H", *(("h",) * r)))
    assert engine.permitted(("phi", "phi", "H", *(("h",) * r)))
    assert engine.permitted(("phi", "phi", "phi", "phi", *(("h",) * r)))


@pytest.mark.parametrize("r", range(3, 9))
def test_no_finite_Einstein_vertex_cutoff(r):
    assert trees.TreeEngine([]).permitted(("h",) * r)


def test_no_kinetic_vertex_double_counting():
    engine = trees.TreeEngine([], cubic=1, contact=1)
    assert not engine.permitted(("phi", "phi"))
    assert not engine.permitted(("H", "H"))
    assert not engine.permitted(("h", "h"))
    assert engine.permitted(("phi", "phi", "H"))
    assert engine.permitted(("phi",) * 4)


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
        "missing_scalar",
        "inexact_momentum",
        "string_momentum",
        "not_list",
    ),
)
def test_original_three_real_domain_rejected(mutation):
    ps, hs = trees.configuration()
    ps = [s.Matrix(p) for p in ps]
    qs = [s.Matrix(q) for q, _ in hs]
    eps = [s.Matrix(e) for _, e in hs]
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
    elif mutation == "missing_scalar":
        ps.pop()
    elif mutation == "inexact_momentum":
        ps[0][0] = s.Float(ps[0][0])
    elif mutation == "string_momentum":
        ps[0] = ps[0].tolist()
        ps[0][0][0] = str(ps[0][0][0])
    rays = list(zip(qs, eps)) if mutation != "not_list" else iter(zip(qs, eps))
    with pytest.raises((TypeError, ValueError)):
        trees.original_amplitude(ps, rays)


def test_exactly_collinear_point_has_no_assigned_finite_value():
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


def test_original_three_real_source_is_not_a_diagnostic_substitution():
    value = trees.data()["whole_three_real_calibration"]
    assert value["graph_count"] == 5116
    assert value["original_amplitude"].is_Rational
    assert value["original_amplitude"] != 0
    assert value["original_amplitude"] != value["diagnostic_amplitude"]


def test_required_new_vertices_cannot_be_dropped():
    gates = trees.data()["gates"]
    assert gates["omitted_EH5_breaks_Ward"]
    assert gates["omitted_scalar_fourth_metric_breaks_Ward"]
    assert gates["two_independent_three_real_gauge_directions_vanish"]


def test_current_and_metric_caches_are_owned_by_finite_instances():
    ps, hs = trees.configuration()
    legs = [("phi", p, 1) for p in ps[1:]] + [("h", q, e) for q, e in hs]
    first, second = trees.TreeEngine(legs), trees.TreeEngine(legs)
    assert first._memo is not second._memo
    first.momentum(1)
    assert "momentum" in first._memo and second._memo == {}
    fs, _ = jets.calibration_fields()
    j1, j2 = jets.MetricJet(fs), jets.MetricJet(fs)
    j1.det(1)
    assert j1._memo and j2._memo == {}


def test_graph_generating_function_radius_is_not_a_scattering_claim():
    scope = audit.observable()
    assert "all-N" in scope["not_established"]
    assert "Regge" in scope["not_established"]
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.matching()) == 171
    assert len(audit.frontier()) == 9
    assert len(audit.qualifications()) == 6
