"""Independent common-master, literal self-energy and full-cut tests."""

import pytest
import sympy as s
from p8_affine import verify as base
from p8_vacuum_affine_massive_common_gravity_masters import (
    audit,
    legs,
    masters,
    reconstruction,
    source,
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


@pytest.mark.parametrize(
    "energy,angle", ((5, 0), (10, s.Rational(1, 3)), (17, s.Rational(-2, 5)))
)
def test_entire_source_and_coefficients_at_independent_kinematics(energy, angle):
    energy, angle = map(s.sympify, (energy, angle))
    mass = s.Integer(1)
    k = s.Integer(11)
    t = -(energy - 4 * mass) * (1 - angle) / 2
    u = 4 * mass - energy - t
    row = reconstruction.coefficients(energy, t, mass)
    p, _, _, a0, a2 = source.elastic.coefficients(
        energy=energy, mass=mass, kappa=k, cubic2=s.Integer(0), contact=s.Integer(0)
    )
    assert s.factor(row["P_times_kappa"] - k * p) == 0
    assert s.factor(row["H0"] - k * a0) == 0
    assert s.factor(row["H2"] - k * a2) == 0
    assert (
        s.factor(
            row["H"]
            - (-2 * energy + 6 * mass - 2 * mass * mass / energy + t * u / energy)
        )
        == 0
    )
    expected = p / (1 - angle * angle) + a0 + a2 * s.legendre(2, angle)
    assert s.factor(source.gravity_tree(energy, t, mass, k) - expected) == 0
    assert source.require_physical(energy, angle) == (energy, angle, mass)


@pytest.mark.parametrize(
    "energy,transfer", ((5, -s.Rational(1, 3)), (10, -2), (17, -5))
)
def test_graviton_cut_coefficient_basis_independently(energy, transfer):
    energy, transfer = map(s.sympify, (energy, transfer))
    mass = s.Integer(1)
    u = 4 * mass - energy - transfer
    row = reconstruction.coefficients(energy, transfer, mass)
    mt, mu_cut, ms = map(s.Rational, ("3/7", "11/13", "2/5"))
    vt = source.invariant.polynomials(transfer, mass)[0]
    vu = source.invariant.polynomials(u, mass)[0]
    actual = (
        4 * (vt * vt * mt + vu * vu * mu_cut) / energy
        - 2 * row["C00mu"] * ms
        + row["B00"]
    )
    expected = 4 * source.invariant.algebraic_cut(
        energy, transfer, mass, mt, mu_cut, ms
    )
    assert s.factor(actual - expected) == 0


@pytest.mark.parametrize(
    "energy,transfer", ((5, -s.Rational(1, 3)), (10, -2), (17, -5))
)
def test_full_evanescent_coefficients_not_discarded(energy, transfer):
    energy, transfer = map(s.sympify, (energy, transfer))
    ep = source.EP
    row = reconstruction.coefficients(energy, transfer, s.Integer(1))
    box, tri = reconstruction.evanescent_coefficients(energy, transfer, s.Integer(1))
    assert s.factor(s.diff(box, ep).subs(ep, 0) - 4 * row["V"]) == 0
    assert (
        s.factor(s.diff(tri, ep).subs(ep, 0) + 4 * row["H"] - 4 * row["V"] / energy)
        == 0
    )
    assert s.diff(box, ep).subs(ep, 0) != 0


@pytest.mark.parametrize("mass", (s.Integer(1), s.Rational(3, 2), s.Integer(7)))
def test_complete_scalar_box_zero_transfer_limit(mass):
    ep = source.EP
    actual = masters.box_finite(-4 * mass, s.Integer(0), mass, ep, s.Integer(3)).doit()
    expected = (-1 / ep + s.log(9 * s.pi / mass) - s.EulerGamma) / (-4 * mass * mass)
    assert s.simplify(actual - expected) == 0


@pytest.mark.parametrize("mass", (s.Integer(1), s.Rational(3, 2), s.Integer(7)))
def test_complete_massive_triangle_zero_energy_limit(mass):
    ep = source.EP
    actual = masters.massive_triangle_finite(
        s.Integer(0), mass, ep, s.Integer(3)
    ).doit()
    expected = (-1 / ep + s.log(36 * s.pi / mass) - s.EulerGamma) / (2 * mass)
    assert s.simplify(actual - expected) == 0


@pytest.mark.parametrize("mass", (s.Integer(1), s.Rational(3, 2), s.Integer(7)))
def test_massive_bubble_zero_energy_keeps_raw_constant(mass):
    actual = masters.massive_bubble_finite(
        s.Integer(0), mass, source.EP, s.Integer(3)
    ).doit()
    expected = -1 / source.EP - s.EulerGamma + s.log(36 * s.pi / mass)
    assert s.simplify(actual - expected) == 0


@pytest.mark.parametrize(
    "eps", (s.Rational(1, 10), s.Rational(-1, 10), s.Rational(1, 20))
)
def test_all_D_zero_self_mass_not_a_four_dimensional_numerator(eps):
    mass = s.Integer(3)
    common = s.Rational(7, 5)
    b0 = common / (1 - 2 * eps)
    a0 = mass * common / (1 - eps)
    assert s.factor(legs.normalized_self_energy(mass, b0, a0, mass, eps)) == 0
    dropped = 2 * mass * mass * b0 - 2 * mass * a0
    assert dropped != 0


@pytest.mark.parametrize("mass", (s.Integer(1), s.Rational(3, 2), s.Integer(7)))
def test_raw_residue_has_both_IR_and_UV_origins(mass):
    k = s.Integer(11)
    ep = source.EP
    expression = legs.raw_residue_derivative(mass, k, ep, s.Integer(3))
    assert (
        s.factor(s.limit(ep * expression, ep, 0) + 3 * mass / (16 * s.pi**2 * k)) == 0
    )
    assert expression.has(s.EulerGamma)
    assert not expression.has(s.Float)


def test_six_ordered_boxes_not_twelve_copies_for_two_cuts():
    rows = reconstruction.crossed_basis()["six_ordered_boxes"]
    pairs = {(row["massless_channel"], row["massive_channel"]) for row in rows}
    assert len(rows) == len(pairs) == 6
    assert all((b, a) in pairs for a, b in pairs)
    assert len(reconstruction.crossed_basis()["three_triangle_bubble_channels"]) == 3


def test_every_master_has_an_explicit_Feynman_limit():
    for function in (
        masters.box_finite,
        masters.massless_triangle,
        masters.massive_triangle_finite,
        masters.massless_bubble_finite,
        masters.massive_bubble_finite,
    ):
        value = function()
        assert isinstance(value, s.Limit)
        assert not value.has(s.Float)


def test_literal_offshell_tensor_and_exact_Ward_scope():
    data = legs.data()
    assert data["whole_offshell_stress_tensor"].shape == (4, 4)
    assert data["whole_general_D_contracted_numerator"].has(source.EPS)
    assert data["checks"]["all_four_offshell_Ward_components"] == s.zeros(4, 1)
    assert "not zero total" not in data["on_shell_mass_statement"]
    assert "does not assert zero total mass shift" in data["on_shell_mass_statement"]


def test_full_source_never_mutates_parent():
    before = dict(source.previous.data()["checks"])
    source.data()
    assert before == source.previous.data()["checks"]
    assert source.data()["checks"] is not source.previous.data()["checks"]
    assert len(before) == 43


def test_full_parent_first_order_preserves_whole_packet():
    before = base.serialize(audit.packets())
    audit.previous.packets()
    source.data.cache_clear()
    audit.packets.cache_clear()
    audit.residuals.cache_clear()
    audit.gates.cache_clear()
    assert base.serialize(audit.packets()) == before
    assert audit.residuals() == ROWS and audit.gates() == GATES


def test_all_historical_frontiers_and_rejected_packet_retained():
    assert audit.frontier() == audit.previous.frontier()
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.matching()) == 139 and len(audit.qualifications()) == 6
    assert any(row["status"].startswith("REJECTED_") for row in audit.matching())
    assert audit.validate_scope(audit.frontier(), audit.matching())
    assert audit.require_parameters(audit.parameters()) == audit.parameters()


def test_cut_matched_representation_does_not_fix_rational_anchor():
    assert "not automatically" in masters.data()["boundary"]
    assert "NOT fixed" in reconstruction.crossed_basis()["scope"]
    assert "not a renormalized physical pole" in legs.data()["finite_boundary"]
    assert "full_finite_amplitude" not in audit.OBSERVABLES


def test_exact_counts_and_controls():
    assert (len(ROWS), audit.scalar_entry_count(), len(GATES)) == (94, 97, 36)
    assert len(audit.controls()) == 8 and audit.rejected_inputs() == 73
    assert len(audit.frontier()) == 9
