"""Independent massive vertex, full tensor, complex bounds and scope tests."""

import pytest
import sympy as s
from p8_affine import verify as base
from p8_vacuum_affine_massive_matter_graviton_vertex import (
    audit,
    matching,
    source,
    vertex,
    ward,
)

ROWS = audit.residuals()
GATES = audit.gates()


@pytest.mark.parametrize("name", tuple(ROWS))
def test_every_exact_identity(name):
    value = ROWS[name]
    assert all(
        x == 0 for x in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    assert not value.atoms(s.Float)


@pytest.mark.parametrize("name", tuple(GATES))
def test_every_proof_gate(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_every_unsupported_input(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("mass,heavy", ((1, 3), (2, 7), (s.Rational(3, 2), 5)))
@pytest.mark.parametrize("order", (1, 2, 3, 4))
def test_independent_simplex_beta_moments(mass, heavy, order):
    mass, heavy = map(s.sympify, (mass, heavy))
    z, y = s.symbols("z y", real=True)
    actual = 0
    for active, spectator in ((mass, heavy), (heavy, mass)):
        A = (1 - z) * active + z * spectator - z * (1 - z) * mass
        actual += s.integrate(
            z * z * (y * (1 - z - y)) ** order, (y, 0, 1 - z)
        ) / A ** (order + 1)
    assert s.factor(actual - vertex.F1_coefficient(order, z, mass, heavy)) == 0


@pytest.mark.parametrize("mass,Q,L", ((1, 1, 2), (2, 2, 1), (s.Rational(3, 2), 1, 3)))
@pytest.mark.parametrize(
    "F1,F2", ((1, s.Rational(-1, 2)), (2, 3), (s.Rational(5, 4), s.Rational(-2, 3)))
)
def test_independent_literal_complete_metric_contraction(mass, Q, L, F1, F2):
    mass, Q, L, F1, F2 = map(s.sympify, (mass, Q, L, F1, F2))
    eta = s.diag(1, -1, -1, -1)
    E = s.sqrt(mass + Q * Q + L * L)
    p = s.Matrix([E, -Q, L, 0])
    r = s.Matrix([E, Q, L, 0])
    other_p = s.Matrix([E, Q, -L, 0])
    other_r = s.Matrix([E, -Q, -L, 0])
    q = r - p
    t = (q.T * eta * q)[0]
    ss = ((p + other_p).T * eta * (p + other_p))[0]
    A = 2 * ((p + r) / 2) * ((p + r) / 2).T * F1 + (q * q.T - eta * t) * F2
    B = (
        2 * ((other_p + other_r) / 2) * ((other_p + other_r) / 2).T * F1
        + (q * q.T - eta * t) * F2
    )
    contraction = s.trace(eta * A * eta * B) - s.trace(eta * A) * s.trace(eta * B) / 2
    assert (
        s.simplify(contraction - matching.metric_contraction(ss, t, F1, F2, mass)) == 0
    )


@pytest.mark.parametrize("mass,heavy", ((1, 3), (2, 7), (s.Rational(3, 2), 5)))
@pytest.mark.parametrize("transfer_ratio", (-2, 1, 2 * s.I))
@pytest.mark.parametrize(
    "z,v", ((s.Rational(1, 3), s.Rational(1, 2)), (s.Rational(3, 4), s.Rational(1, 3)))
)
def test_exact_complex_integrand_remainder_controls(mass, heavy, transfer_ratio, z, v):
    mass, heavy = map(s.sympify, (mass, heavy))
    t = mass * transfer_ratio
    w = (1 - z) ** 2 * (1 - v * v) / 4
    slope = 0
    for a, b in ((mass, heavy), (heavy, mass)):
        A = (1 - z) * a + z * b - z * (1 - z) * mass
        slope += (1 - z) * z * z * w / (A * A)
    value = s.cancel(vertex.F1_integrand(t, z, v, mass, heavy))
    r = s.Abs(t) / (4 * mass)
    assert s.Abs(value) <= s.Abs(t) * slope / (1 - r)
    assert s.Abs(value - t * slope) <= s.Abs(t) ** 2 * slope / (4 * mass * (1 - r))


@pytest.mark.parametrize(
    "p,r",
    (
        ((3, 1, 0, 1), (4, -1, 2, 0)),
        ((5, 2, 1, -1), (3, 0, -2, 1)),
        ((2, 1, 1, 0), (7, 2, -1, 3)),
    ),
)
def test_independent_polynomial_inverse_Ward(p, r):
    p, r = map(s.Matrix, (p, r))
    eta = s.diag(1, -1, -1, -1)
    q = r - p
    x = (p.T * eta * p)[0]
    y = (r.T * eta * r)[0]
    t = (q.T * eta * q)[0]
    K = lambda a: (
        a - 1 + s.Rational(1, 11) * (a - 1) ** 2 + s.Rational(1, 13) * (a - 1) ** 3
    )
    F = s.cancel((K(y) - K(x)) / (y - x))
    tensor = (p * r.T + r * p.T) * F - eta * (K(x) + K(y) - t * F) / 2
    assert tensor * eta * q == p * K(y) - r * K(x)


@pytest.mark.parametrize("coefficient", (-7, -1, 1, 5))
def test_transverse_Ricci_addition_changes_slope_without_Ward(coefficient):
    eta = s.diag(1, -1, -1, -1)
    P = s.Matrix([3, 0, 1, 0])
    q = s.Matrix([0, 2, 0, 0])
    t = (q.T * eta * q)[0]
    tensor = coefficient * ward.transverse_ricci(P, q)
    assert tensor * eta * q == s.zeros(4, 1)
    assert tensor == coefficient * t * P * P.T
    assert tensor[0, 2] / (2 * P[0] * P[2]) == coefficient * t / 2
    assert s.Rational(coefficient, 2) != 0


def test_actual_hierarchy_bound_is_not_a_full_P8_verdict():
    data = matching.data()
    bounds = data["actual_exact_arithmetic_upper_bounds"]
    assert 0 < bounds["Pi_second_upper"] < s.Rational(1, 10**405)
    assert 0 < bounds["whole_abs_t_le_2_b2_magnitude_upper"] < s.Rational(1, 10**1205)
    assert "not the whole" in data["t_channel_matching"]
    assert "not adopted" in data["primary_comparison"]
    assert "Ricci-derivative" in audit.observable()["not_established"]


def test_whole_fixed_source_and_light_curvature_frontier_retained():
    data = source.data()
    assert data["same_original_vacuum_parameters"]["g"] == s.Rational(1, 8192)
    assert "stays unmatched" in data["curved_boundary"]
    assert "not a new" in data["unchanged_OS_scope"]
    assert source.data()["checks"] is not source.previous.data()["checks"]


def test_heavy_metric_mixing_is_not_discarded_from_full_F2():
    data = vertex.data()
    assert "H-metric mixing" in data["complete_projected_graph_count"]
    assert "not discarded" in data["complete_projected_graph_count"]
    assert "complete physical F2" in data["matching_boundary"]
    assert "internal gravitons" in ward.data()["whole_background_scope"]


def test_parent_warmup_does_not_mutate_packet_contract():
    before = base.serialize(audit.packets())
    audit.previous.packets()
    assert base.serialize(audit.packets()) == before
    assert len(audit.gates()) == 36


def test_all_original_and_historical_frontiers_preserved():
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.matching()) == 142 and len(audit.frontier()) == 9
    assert audit.qualifications() == audit.previous.qualifications()
    assert len(audit.qualifications()) == 6
    assert any(row["status"].startswith("REJECTED_") for row in audit.matching())
    assert audit.controls()["rejected_inputs"] == 73
    assert audit.require_parameters(audit.parameters()) == audit.parameters()


def test_both_massive_insertions_and_whole_OS_are_required():
    z = vertex.Z
    F = vertex.MU * (1 - z) ** 2 + vertex.N * z
    light = s.integrate(vertex.raw_F1(vertex.MU, vertex.N, 0), (vertex.V, 0, 1))
    assert s.factor(light - z * (1 - z) / F) != 0
    assert vertex.F1_integrand(0) == 0
    assert "including its finite value" in vertex.data()["source_OS_normalization"]


def test_original_P8_not_closed_by_this_piece():
    statement = audit.observable()["not_established"]
    for text in ("finite crossed amplitude", "Regge", "V/G/B/P8", "massless", "S285"):
        assert text in statement
    with pytest.raises(ValueError):
        audit.require_observable("full_finite_amplitude")
