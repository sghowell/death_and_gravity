"""Whole-graph topology, canonical budgets, energy limits and scope controls."""

import pytest
import sympy as s
from p8_vacuum_affine_complete_two_graviton_tree import trees as e
from p8_vacuum_affine_uniform_two_real_tree_bound import (
    audit,
    bounds,
    gaps,
    source,
    topology,
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
def test_scope_rejections(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("key", tuple(topology.EXPECTED_PROFILES))
def test_every_regular_topology_profile_has_exact_multiplicity(key):
    records = topology.data()["whole_regular_graph_profiles"]
    rows = [
        r
        for r in records
        if (
            r["sector"],
            r["scalar_propagators"],
            r["hard_gravitons"],
            r["heavy_metric_vertices"],
            r["Einstein_cubics"],
            r["Einstein_quartics"],
            r["light_metric_degree"],
        )
        == key
    ]
    assert len(rows) == 1 and rows[0]["count"] == topology.EXPECTED_PROFILES[key]


def test_no_regular_hard_cut_is_mistaken_for_the_pure_soft_pair():
    for _, paired, _, edges in topology.classified():
        if paired:
            assert sum(k == "h" and p == 0 and q == 3 for k, p, q in edges) == 1
        else:
            assert all(p.bit_count() == 2 for k, p, q in edges if k in ("h", "H"))


@pytest.mark.parametrize(
    "first,second",
    (
        (s.Rational(1, 16), s.Rational(1, 16)),
        (s.Rational(1, 10**400), s.Rational(1, 20)),
        (s.Rational(1, 10**40), s.Rational(1, 10**60)),
    ),
)
def test_bare_full_bound_retains_both_soft_poles(first, second):
    v = bounds.full_upper(first, second)
    assert v == bounds.full_upper(second, first)
    assert bounds.full_upper(first / 2, second / 2) == 4 * v
    assert bounds.full_upper(first / 2, second) == 2 * v
    assert v == 2 * bounds.regular_upper(first, second)


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
def test_energy_domain_rejects_inexact_or_outside_inputs(invalid):
    for call in (bounds.energies, bounds.regular_upper, bounds.full_upper):
        with pytest.raises((TypeError, ValueError)):
            call(invalid, s.Rational(1, 64))


@pytest.mark.parametrize(
    "first,second",
    (
        (s.Rational(1, 8), s.Rational(1, 8)),
        (s.Rational(1, 10), s.Rational(1, 10)),
        (s.Rational(1, 8), s.Rational(1, 10**40)),
    ),
)
def test_total_energy_cut_not_individual_cuts(first, second):
    with pytest.raises(ValueError):
        bounds.full_upper(first, second)


@pytest.mark.parametrize("r", (1, 2, 3))
def test_general_non_TT_scalar_metric_vertices_fit_written_budget(r):
    fields = (
        e.imm([[1, 1, -1, 1], [1, 0, 1, -1], [-1, 1, -1, 1], [1, -1, 1, 0]]),
        e.imm([[0, 1, 1, -1], [1, -1, 0, 1], [1, 0, 1, 1], [-1, 1, 1, -1]]),
        e.imm([[1, -1, 0, 1], [-1, 1, 1, 0], [0, 1, 1, -1], [1, 0, -1, 1]]),
    )[:r]
    p = e.imm([2, 1, 0, -1])
    q = e.imm([2, -1, 1, 0])
    assert abs(e.determinant_coefficient(fields)) <= bounds.DET[r]
    assert max(abs(x) for x in e.inverse_coefficient(fields)) <= bounds.INV[r]
    assert (
        max(abs(x) for x in e.density_inverse_coefficient(fields)) <= bounds.DENSITY[r]
    )
    assert abs(e.scalar_vertex(p, q, fields, s.S.One)) < 1024**r


@pytest.mark.parametrize("row", range(len(gaps.SAMPLES)))
def test_every_sample_has_all_four_cluster_assignments(row):
    label, *args = gaps.SAMPLES[row]
    points, qs, _, _ = gaps.configuration(*args)
    Q = qs[0] + qs[1]
    for mask in range(4):
        Qa = sum((qs[j] for j in range(2) if mask >> j & 1), e.VECTOR_ZERO)
        A = points[2] + Qa
        B = points[3] + Q - Qa
        assert A + B + points[0] + points[1] == s.zeros(4, 1)
        assert e.old.dot(A, A) >= 1 and e.old.dot(B, B) >= 1
    assert gaps.physical_samples()[1][label + "_all_cluster_and_scalar_gaps"]


def test_cluster_mass_cannot_be_silently_replaced_by_mass_one():
    points, qs, _, _ = gaps.configuration()
    A = points[2] + qs[0]
    assert e.old.dot(A, A) > 1
    E = -points[0][0]
    u = e.old.dot(A, A) - 1
    B = points[3] + qs[1]
    v = e.old.dot(B, B) - 1
    r2 = sum(A[j] ** 2 for j in range(1, 4))
    assert s.factor(E * E - 1 - r2 - (u + v) / 2 + (u - v) ** 2 / (16 * E * E)) == 0
    assert s.factor(E * E - 1 - r2) != 0


def test_small_prefactor_does_not_make_the_soft_envelope_integrable():
    mu, x = s.symbols("mu x", positive=True)
    inside_square = 4 * s.log(x / (2 * mu)) ** 2 / s.Integer(10) ** 790
    assert s.limit(inside_square, mu, 0, dir="+") == s.oo
    assert "infrared-finite" in audit.observable()["not_established"]


def test_original_source_and_all_frontiers_preserved():
    assert (
        source.KAPPA == 10**800 and source.HEAVY_MASS2 == s.Rational(10**200, 512) + 2
    )
    assert source.CUBIC == s.Rational(1, 8192)
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 168
    assert len(audit.qualifications()) == 6
    assert audit.matching()[:-1] == audit.previous.matching()
