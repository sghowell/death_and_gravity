"""Independent dimensional tree, complete master and threshold tests."""

import pytest
import sympy as s
from p8_affine import verify as base
from p8_vacuum_affine_massive_dimensional_cut_completion import (
    audit,
    dimensional,
    polarizations,
    source,
    threshold,
)

ROWS = audit.residuals()
GATES = audit.gates()


@pytest.mark.parametrize("name", tuple(ROWS))
def test_each_exact_identity(name):
    value = ROWS[name]
    assert all(
        x == 0 for x in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    assert not value.atoms(s.Float)


@pytest.mark.parametrize("name", tuple(GATES))
def test_each_written_proof_gate(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_every_unsupported_input(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("dimension", (2, 3, 4, 5, 6))
def test_full_physical_projector_on_independent_vectors(dimension):
    n = dimension
    P = polarizations.tt_projector(n)
    L = s.Matrix([1, 2] + [0] * (n - 2))
    R = s.Matrix([3, -1] + [0] * (n - 2))
    QL = 2 * L * L.T - 7 * s.eye(n)
    QR = 2 * R * R.T - 11 * s.eye(n)
    assert P * P == P
    assert P.T == P
    assert s.trace(P) == n * (n + 1) // 2 - 1
    actual = s.trace(P * s.kronecker_product(QL, QL) * P * s.kronecker_product(QR, QR))
    expected = polarizations.projector_numerator(
        s.Integer(5),
        s.Integer(10),
        s.Integer(1),
        s.Integer(7),
        s.Integer(11),
        s.Integer(n),
    )
    assert s.factor(actual - expected) == 0


@pytest.mark.parametrize("dimension", (5, 6, 7, 8))
def test_literal_full_Einstein_and_scalar_graphs_independent_dimensions(dimension):
    for row in polarizations.literal_controls(dimension, dimension):
        actual = row["contact"] + row["scalar1"] + row["scalar2"] + row["EH_exchange"]
        assert s.factor(actual - row["expected"]) == 0
        assert s.factor(actual - row["actual"]) == 0


@pytest.mark.parametrize("mass", (s.Integer(1), s.Rational(3, 2), s.Integer(7)))
@pytest.mark.parametrize(
    "ratio,transfer_ratio", ((5, -s.Rational(1, 3)), (10, -2), (17, -5))
)
def test_entire_four_dimensional_master_limit_at_independent_kinematics(
    mass, ratio, transfer_ratio
):
    energy = ratio * mass
    transfer = transfer_ratio * mass
    row = dimensional.master_coefficients(energy, transfer, mass, s.Integer(0))
    prior = source.reconstruction.coefficients(energy, transfer, mass)
    for key in ("C00mu", "B00", "C0mumu", "Bmm"):
        assert s.factor(row[key] - prior[key]) == 0
    v = transfer**2 - 4 * mass * transfer + 2 * mass**2
    assert s.factor(row["ordered_box"] - v * v) == 0


@pytest.mark.parametrize(
    "epsilon", (s.Rational(-1, 10), s.Integer(0), s.Rational(1, 10), s.Rational(1, 3))
)
@pytest.mark.parametrize(
    "energy,transfer", ((5, -s.Rational(1, 3)), (10, -2), (17, -5))
)
def test_full_dimensional_both_species_coefficients(epsilon, energy, transfer):
    ss, tt, mu = map(s.sympify, (energy, transfer, 1))
    q = ss - 4 * mu
    z = 1 + 2 * tt / q
    n = 2 + 2 * epsilon
    d = n + 1
    vd = ss * ss - 4 * mu * ss + 2 * mu * mu + 2 * mu * mu * epsilon / (1 + epsilon)
    vt = tt * tt - 4 * mu * tt + 2 * mu * mu + 2 * mu * mu * epsilon / (1 + epsilon)
    P = 4 * vd / q
    b = -q * q / (4 * ss)
    a = (
        -2 * ss
        + 6 * mu
        - 2 * mu * mu / ss
        + q * q / (4 * ss)
        - 2 * mu * mu * epsilon / (ss * (1 + epsilon))
    )
    row = dimensional.master_coefficients(ss, tt, mu, epsilon)
    assert s.factor(row["ordered_box"] - vt * vt) == 0
    assert s.factor(row["C0mumu"] + 2 * vd * (a + b * z * z)) == 0
    reg = a * a + 2 * a * b / d + b * b * (1 + 2 * z * z) / (d * (d + 2))
    interference = 2 * P * b * (-z * z + (1 - z * z) / n)
    assert s.factor(2 * row["Bmm"] - reg - interference) == 0


@pytest.mark.parametrize(
    "energy,transfer", ((5, -s.Rational(1, 3)), (10, -2), (17, -5))
)
def test_finite_bubble_rational_is_actual_EP_coefficient(energy, transfer):
    ep = source.EP
    row = dimensional.master_coefficients(energy, transfer, s.Integer(1), ep)
    actual = -s.diff(row["B00"] + row["Bmm"], ep).subs(ep, 0)
    expected = dimensional.evanescent_rational(energy, transfer, s.Integer(1))
    assert s.factor(actual - expected) == 0
    assert actual != 0


@pytest.mark.parametrize("mass", (s.Integer(1), s.Rational(3, 2), s.Integer(7)))
def test_entire_threshold_and_tadpole_dictionary(mass):
    ep = source.EP
    scale = s.Integer(3)
    B = threshold.bubble_threshold(mass, ep, scale)
    A0 = threshold.tadpole_finite(mass, ep, scale)
    assert s.factor(A0 / mass + 1 - B) == 0
    C = threshold.triangle_threshold(mass)
    assert s.re(C) == s.log(2) / mass
    assert s.im(C) == -s.pi / (2 * mass)
    assert s.im(threshold.triangle_threshold_slope(mass)) == s.pi / (12 * mass**2)
    energy = 5 * mass
    transfer = -mass / 3
    q = energy - 4 * mass
    u = -q - transfer
    old = threshold.old_gram_choice(energy, transfer, mass, B)
    ev = threshold.evanescent_gram_choice(energy, transfer, mass)
    full = threshold.full_gram_choice(energy, transfer, mass, B)
    assert s.factor(full - old - ev) == 0
    expected = (
        -(12 * mass * A0 + 4 * mass * mass) * transfer * u / q**2
        + (5 * A0 + mass) * transfer * u / q
    )
    assert s.factor(full - expected) == 0
    assert s.factor(full - threshold.full_gram_choice(energy, u, mass, B)) == 0


@pytest.mark.parametrize(
    "energy,transfer", ((5, -s.Rational(1, 3)), (10, -2), (17, -5))
)
def test_compact_crossed_evans_preserves_all_physical_pole_terms(energy, transfer):
    ss, tt, mu = map(s.sympify, (energy, transfer, 1))
    u = 4 * mu - ss - tt
    raw = sum(
        dimensional.evanescent_rational(a, b, mu)
        - threshold.evanescent_gram_choice(a, b, mu)
        for a, b in ((ss, tt), (tt, ss), (u, ss))
    )
    assert s.factor(raw - dimensional.compact_crossed_evans(ss, tt, mu)) == 0
    assert s.factor(raw - dimensional.compact_crossed_evans(tt, ss, mu)) == 0


def test_parent_warmup_does_not_change_whole_packet_contract():
    before = base.serialize(audit.packets())
    audit.previous.packets()
    assert base.serialize(audit.packets()) == before
    assert source.data()["checks"] is not source.previous.data()["checks"]
    assert (len(ROWS), audit.scalar_entry_count(), len(GATES)) == (114, 114, 36)


def test_whole_prior_frontier_and_historical_qualifications():
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.matching()) == 140 and len(audit.frontier()) == 9
    assert audit.qualifications() == audit.previous.qualifications()
    assert len(audit.qualifications()) == 6
    assert any(row["status"].startswith("REJECTED_") for row in audit.matching())
    assert audit.controls()["rejected_inputs"] == 73


def test_full_physical_and_dimensional_boundaries_are_explicit():
    scope = audit.observable()
    assert "tadpole" in scope["not_established"]
    assert "physical finite pole" in scope["not_established"].lower()
    assert "V/G/B/P8 closure" in scope["not_established"]
    assert "noninteger" in scope["domain"]
    assert (
        "not claimed to be a positive" in polarizations.data()["dimensional_boundary"]
    )
    assert "threshold-supported" in dimensional.data()["completion_boundary"]
    assert "not a choice" in threshold.data()["prescription"]
    assert (
        source.data()["gates"]["leading_flat_reference_not_exact_quantum_vacuum"]
        is True
    )


def test_source_frozen_parameters_are_not_retuned():
    assert audit.require_parameters(audit.parameters()) == audit.parameters()
    assert source.data()["gates"]["fixed_vacuum_constants_not_retuned"] is True
    assert "all_three_vacuum_constants" in source.data()
