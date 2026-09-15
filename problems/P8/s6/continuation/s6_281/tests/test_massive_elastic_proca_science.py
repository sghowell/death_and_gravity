"""Independent original-source, whole-angle, dimensional and Proca checks."""

import pytest
import sympy as s
from p8_vacuum_affine_massive_elastic_proca_infrared import (
    audit,
    elastic,
    infrared,
    proca,
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
def test_each_proof_gate(name):
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
def test_entire_tree_at_independent_massive_kinematics(energy, angle):
    mass = s.Integer(1)
    n = s.Integer(50)
    g2 = s.Rational(3, 17)
    c = -s.Rational(1, 16)
    k = s.Integer(11)
    t = -(energy - 4 * mass) * (1 - angle) / 2
    u = -(energy - 4 * mass) * (1 + angle) / 2
    inv = [s.Integer(energy), t, u]
    grav = (
        -sum(
            (2 * mass**2 - 2 * mass * v - inv[(i + 1) % 3] * inv[(i + 2) % 3]) / v
            for i, v in enumerate(inv)
        )
        / k
    )
    matter = c + g2 * sum(1 / (n - v) for v in inv)
    assert (
        s.factor(elastic.whole_tree(angle, energy, mass, k, n, g2, c) - grav - matter)
        == 0
    )
    assert elastic.require_domain(energy, angle, 1, 50) == (energy, angle, 1, 50)


@pytest.mark.parametrize("value", (s.Integer(2), s.Integer(3), s.Rational(7, 4)))
def test_continuous_master_discriminant_zero(value):
    assert elastic.master_I(value, value, s.Integer(1)) == 1 / (value * value - 1)


@pytest.mark.parametrize("z", (s.Integer(0), s.Rational(1, 3), s.Rational(-2, 5)))
def test_complete_angular_master_parity(z):
    w, v = s.Integer(2), s.Integer(3)
    assert elastic.master_J(w, v, z) == elastic.master_J(w, v, -z)
    assert s.simplify(elastic.master_J(w, v, z) - elastic.master_J(v, w, z)) == 0


def test_whole_original_source_does_not_mutate_parent_check_contract():
    before = dict(source.previous.data()["checks"])
    source.data()
    assert before == source.previous.data()["checks"]
    assert source.data()["checks"] is not source.previous.data()["checks"]
    assert len(before) == 21


def test_order_warmup_cannot_change_current_source_checks():
    before = dict(source.data()["checks"])
    source.previous.previous.data()
    source.data.cache_clear()
    assert source.data()["checks"] == before


def test_original_mass_and_heavy_threshold_not_rounded():
    assert source.VECTOR_MASS2 == 10**6
    assert source.HEAVY_MASS2 == s.Integer(10) ** 200 / 512 + 2
    assert 4 * source.VECTOR_MASS2 < 10**196 < source.HEAVY_MASS2
    assert source.HEAVY_MASS2 != 10**200


@pytest.mark.parametrize("mass", (s.Integer(1), s.Rational(3, 2), s.Integer(7)))
def test_proca_exact_threshold_spin_weights(mass):
    mass2 = mass * mass
    assert proca.spin_weights(4 * mass2, mass2) == (12 * mass2**2, 16 * mass2**2)
    assert (
        proca.whole_cut(4 * mass2, s.Rational(1, 3), s.Integer(1), mass2, s.Integer(11))
        == 0
    )


@pytest.mark.parametrize("ratio", (5, 9, 17))
def test_proca_forward_spin_weights_positive_on_full_open_channel(ratio):
    mass2 = s.Rational(7, 3)
    energy = ratio * mass2
    w0, w2 = proca.spin_weights(energy, mass2)
    assert w0 > 0 and w2 > 0
    assert (
        proca.sewn_shape(energy, s.Integer(1), s.Integer(1), mass2, s.Integer(13)) > 0
    )


def test_nine_pairs_are_not_a_four_polarization_sum():
    tensors = proca.data()["all_nine_literal_physical_stress_tensors"]
    assert len(tensors) == 9
    assert all(B.shape == (4, 4) for B in tensors)
    assert proca.spin_weights()[1].coeff(proca.S, 2) == s.Rational(13, 30)


def test_evanescent_massive_trace_is_nonzero_and_not_the_soft_tree_piece():
    point = {source.S: 10, source.MU: 1, source.K: 13, source.Z: s.Rational(1, 3)}
    p1 = (8 * source.MU**2 / (source.K * (source.S - 4 * source.MU))).subs(point)
    assert p1 > 0
    assert infrared.evanescent_tree().subs(point) != 0
    assert infrared.whole_stripped_cut().has(s.EulerGamma)
    assert all(
        infrared.whole_stripped_cut().has(x)
        for x in (infrared.E, elastic.C, elastic.G2, source.MU)
    )


@pytest.mark.parametrize(
    "energy,angle", ((5, 0), (10, s.Rational(1, 3)), (17, s.Rational(-2, 5)))
)
def test_complete_finite_subtraction_uses_dimensional_tree(energy, angle):
    eps = infrared.EP
    p = elastic.coefficients()[0]
    A = elastic.whole_tree()
    A1 = infrared.evanescent_tree()
    pref = s.sqrt(1 - 4 * source.MU / source.S) / (32 * s.pi)
    soft = pref * p * (1 / eps + 2 * s.log(infrared.E / infrared.NU)) * (A + eps * A1)
    expression = s.expand(
        infrared.whole_bare_cut() - soft - infrared.whole_stripped_cut()
    )
    point = {
        source.S: energy,
        source.MU: 1,
        source.K: 13,
        source.Z: angle,
        elastic.N: 50,
        elastic.G2: s.Rational(3, 17),
        elastic.C: -s.Rational(1, 16),
    }
    finite = s.limit(s.expand_log(expression.subs(point), force=True), eps, 0, dir="+")
    assert s.simplify(finite) == 0


def test_physical_resolution_flow_contains_full_original_tree():
    derivative = s.diff(infrared.whole_stripped_cut(), infrared.E)
    assert (
        s.simplify(
            derivative
            + s.sqrt(1 - 4 * source.MU / source.S)
            * elastic.coefficients()[0]
            * elastic.whole_tree()
            / (16 * s.pi * infrared.E)
        )
        == 0
    )
    assert not infrared.whole_stripped_cut().has(infrared.NU, infrared.EP)


def test_all_four_species_are_distinct_and_not_double_counted():
    species = audit.whole_species_cut()
    physical = [
        k
        for k in species
        if k not in ("domain", "no_double_count", "not_an_optical_positivity_theorem")
    ]
    assert len(physical) == 4 and len(set(physical)) == 4
    assert isinstance(species["two_physical_TT_gravitons"], s.Expr)
    assert species["two_physical_TT_gravitons"].has(s.Integral)
    assert (
        species["Phi_Phi_after_specified_soft_stripping"]
        == infrared.whole_stripped_cut()
    )
    assert "second time" in species["no_double_count"]


def test_all_original_obligations_and_rejected_packet_status_remain():
    assert audit.frontier() == audit.previous.frontier()
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.matching()) == 137 and len(audit.qualifications()) == 6
    assert any(r["status"].startswith("REJECTED_") for r in audit.matching())
    assert audit.validate_scope(audit.frontier(), audit.matching())
    assert audit.require_parameters(audit.parameters()) == audit.parameters()


def test_exact_counts_and_controls():
    assert (len(ROWS), audit.scalar_entry_count(), len(GATES)) == (118, 310, 33)
    assert len(audit.controls()) == 8 and audit.rejected_inputs() == 75
    assert len(audit.frontier()) == 9
