"""Independent source, Born, rate and forward-boundary tests."""

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_one_newton_inclusive_assembly import (
    audit,
    forward,
    inclusive,
    matter,
    source,
)

ROWS, GATES = audit.residuals(), audit.gates()


@pytest.mark.parametrize("name", tuple(ROWS))
def test_every_exact_identity(name):
    value = ROWS[name]
    assert all(
        x == 0 for x in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    assert not s.sympify(value).has(s.Float)


@pytest.mark.parametrize("name", tuple(GATES))
def test_every_written_proof_gate(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_every_scope_guard(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize(
    "name", tuple(forward.data()["whole_explicit_positive_margins"])
)
def test_each_forward_margin(name):
    assert forward.data()["whole_explicit_positive_margins"][name] > 0


@pytest.mark.parametrize("E", (s.Rational(5, 4), s.Rational(3, 2), s.Integer(2)))
@pytest.mark.parametrize("z", (s.Rational(-3, 5), s.Integer(0), s.Rational(4, 5)))
def test_independent_full_stress_tensor_gravity_Born(E, z):
    eta = s.diag(1, -1, -1, -1)
    r = s.sqrt(E * E - 1)
    sn = s.sqrt(1 - z * z)
    k = [
        s.Matrix([-E, 0, 0, -r]),
        s.Matrix([-E, 0, 0, r]),
        s.Matrix([E, r * sn, 0, r * z]),
        s.Matrix([E, -r * sn, 0, -r * z]),
    ]

    def dot(a, b):
        return (a.T * eta * b)[0]

    def stress(a, b):
        return a * b.T + b * a.T - eta * (dot(a, b) + 1)

    actual = 0
    for j in (1, 2, 3):
        other = [i for i in range(4) if i not in (0, j)]
        left = stress(k[0], k[j])
        right = stress(k[other[0]], k[other[1]])
        numerator = (
            s.trace(eta * left * eta * right)
            - s.trace(eta * left) * s.trace(eta * right) / 2
        )
        actual -= numerator / dot(k[0] + k[j], k[0] + k[j])
    assert s.factor(actual - forward.gravity_born(4 * E * E, z, 1)) == 0
    assert actual > 0


@pytest.mark.parametrize("n", (s.Integer(128), s.Integer(257), s.Integer(4096)))
@pytest.mark.parametrize("energy", (s.Rational(25, 4), s.Integer(9), s.Integer(16)))
@pytest.mark.parametrize("z", (s.Integer(0), s.Rational(3, 5), s.Rational(99, 100)))
def test_full_compact_Born_bounds(n, energy, z):
    am = forward.matter_born(energy, z, n, 1)
    ag = forward.gravity_born(energy, z, 1)
    xi = n**3
    w = 1 - z * z
    assert 4 / n**3 < am < 470 / n**3
    assert xi / (33 * w) < ag / am < 17 * xi / w


@pytest.mark.parametrize("exponent", (190, 200, 204, 220))
def test_original_forward_crossover_without_inexact_hierarchy(exponent):
    with mp.workdps(300):
        n = mp.mpf(int(source.HEAVY_MASS2.p)) / int(source.HEAVY_MASS2.q)
        g = mp.mpf(1) / 8192
        kap = mp.mpf(10) ** 800
        energy = mp.mpf(9)
        q = energy - 4
        w = mp.mpf(10) ** (-exponent)
        z = mp.sqrt(1 - w)
        channels = (energy, -q * (1 - z) / 2, -q * (1 + z) / 2)
        am = g * g / (n - 2) ** 2 * sum((a - 2) ** 2 / (n - a) for a in channels)
        ag = (
            (4 * q + 16 + 8 / q) / w
            - 2 * energy
            + 6
            - 2 / energy
            + q * q * w / (4 * energy)
        ) / kap
        xi = n**3 / (kap * g * g)
        ratio = ag / am
        assert xi / (33 * w) < ratio < 17 * xi / w
        if exponent == 190:
            assert ratio < mp.mpf("1.7e-9")
        if exponent >= 204:
            assert ratio > 150


@pytest.mark.parametrize(
    "ratio",
    (
        s.Rational(1, 10**100),
        s.Rational(1, 100),
        s.Integer(1),
        s.Integer(100),
        s.Integer(10) ** 100,
    ),
)
def test_interference_normalization_has_global_maximum(ratio):
    assert 0 < matter.interference_weight(ratio) <= s.Rational(1, 2)
    assert s.factor(s.Rational(1, 2) - matter.interference_weight(ratio)) == (
        ratio - 1
    ) ** 2 / (2 * (ratio + 1) ** 2)


def test_missing_interference_not_silently_zeroed():
    value = inclusive.normalized_one_Newton_correction(2, 3, 5, 7, 11, 13)
    assert value == s.Rational(154, 25)
    omitted = inclusive.normalized_one_Newton_correction(2, 3, 0, 7, 11, 13)
    assert value - omitted == s.Rational(6, 5)


def test_unmatched_hard_data_remain_in_complete_rate():
    assert (
        s.factor(s.diff(inclusive.normalized_one_Newton_correction(), inclusive.HD))
        != 0
    )
    assert (
        s.factor(s.diff(inclusive.normalized_one_Newton_correction(), inclusive.LM))
        != 0
    )
    assert inclusive.reference_error() < s.Rational(2, 10**199)


def test_same_complete_matter_renormalization_not_a_new_gravity_value():
    packet = source.data()
    assert len(packet["whole_finite_kappa_matter_couplings"]) == 7
    assert "S239" in packet["whole_renormalization_scope"]
    assert "No new finite value" in packet["whole_renormalization_scope"]
    assert "not interchange" in packet["whole_graph_equality_proof"]


def test_rate_not_used_as_an_analytic_amplitude():
    assert "not an analytic" in inclusive.data()["whole_matching_scope"]
    assert "before t tends to zero" in forward.data()["whole_forward_boundary"]
