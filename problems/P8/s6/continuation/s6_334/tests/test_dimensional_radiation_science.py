"""Independent component, exact regulator and preserved-source checks."""

from itertools import product

import pytest
import sympy as s
from p8_vacuum_affine_dimensional_gravity_radiation import (
    audit,
    sew,
    source,
    tensor,
    vertices,
    ward,
)
from p8_vacuum_affine_minimal_gravity_radiation import vertices as old_vertices

ROWS, GATES = audit.residuals(), audit.gates()


@pytest.mark.parametrize("name", tuple(ROWS))
def test_every_exact_residual(name):
    value = ROWS[name]
    assert all(
        v == 0 for v in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    assert not s.sympify(value).has(s.Float)


@pytest.mark.parametrize("name", tuple(GATES))
def test_every_written_proof_gate(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[row[0] for row in audit.bad_cases()]
)
def test_rejected_inputs_and_scope(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize(
    "bad",
    (
        None,
        True,
        False,
        "1",
        1.0,
        s.Float(1),
        [],
        {},
        s.Matrix([1]),
        s.I,
        s.oo,
        s.nan,
        s.Symbol("unknown"),
    ),
)
def test_exact_real_domains(bad):
    for call, args in (
        (sew.regulator, (bad,)),
        (sew.ratio, (bad,)),
        (sew.regulator_rate_upper, (bad, 1, 0)),
        (sew.regulator_rate_upper, (1, bad, 0)),
        (sew.regulator_rate_upper, (1, 1, bad)),
        (sew.dimension_sew, (s.eye(2), bad)),
    ):
        with pytest.raises((TypeError, ValueError)):
            call(*args)


@pytest.mark.parametrize("bad", (-1, -s.Rational(1, 100), s.Rational(1, 7), 1))
def test_closed_regulator_interval(bad):
    for call in (sew.regulator, sew.ratio):
        with pytest.raises(ValueError):
            call(bad)


@pytest.mark.parametrize("bad", (-3, 0, 1, 2, 3, s.Rational(7, 2)))
def test_polarization_dimension_floor(bad):
    with pytest.raises(ValueError):
        sew.dimension_sew(s.eye(2), bad)


@pytest.mark.parametrize("bad", (None, True, 4.0, s.Integer(4), 3, 0, -1, "4"))
def test_literal_component_integer_dimension(bad):
    with pytest.raises(ValueError):
        vertices.component_engine(bad)


@pytest.mark.parametrize(
    "bad",
    (
        None,
        True,
        [[1, 0], [0, 1]],
        s.eye(3),
        s.Matrix([[1, 2], [3, 4]]),
        s.Matrix([[s.Float(1), 0], [0, 1]]),
        s.Matrix([[s.oo, 0], [0, 1]]),
        s.Matrix([[s.Symbol("unfixed"), 0], [0, 1]]),
    ),
)
def test_canonical_core_domain(bad):
    for call, args in (
        (sew.core, (bad,)),
        (sew.coefficients, (bad, s.zeros(2))),
        (sew.coefficients, (s.zeros(2), bad)),
        (sew.dimension_sew, (bad, 4)),
    ):
        with pytest.raises((TypeError, ValueError)):
            call(*args)


PAIRS = (
    (s.eye(2), s.zeros(2)),
    (s.Matrix([[1, 2], [2, -1]]), s.Matrix([[2, 1], [1, 0]])),
    (s.Matrix([[s.I, 2], [2, 1]]), s.Matrix([[1, s.I], [s.I, -2]])),
    (s.zeros(2), s.Matrix([[1, s.I], [s.I, -1]])),
)


@pytest.mark.parametrize(
    "index,epsilon",
    tuple(
        product(
            range(len(PAIRS)),
            (s.S.Zero, s.Rational(1, 10**40), s.Rational(1, 32), s.Rational(1, 8)),
        )
    ),
)
def test_exact_full_complex_sew_and_quantitative_bound(index, epsilon):
    A, B = PAIRS[index]
    t = sew.ratio(epsilon)
    actual = sew.continued_rate(A, B, epsilon)
    assert s.expand(actual - sew.dimension_sew(A + t * B, 4 + 2 * epsilon)) == 0
    assert actual.is_nonnegative is True
    a = sum(abs(v) for v in A)
    b = sum(abs(v) for v in B)
    upper = sew.regulator_rate_upper(a, b, epsilon)
    assert (
        s.simplify(upper - abs(actual - sew.dimension_sew(A, 4))).is_nonnegative is True
    )
    assert not actual.has(s.Float) and not upper.has(s.Float)


@pytest.mark.parametrize(
    "epsilon", (s.Rational(1, 10**40), s.Rational(1, 32), s.Rational(1, 8))
)
def test_missing_extra_trace_is_detected(epsilon):
    assert sew.dimension_sew(s.eye(2), 4) == 0
    assert sew.continued_rate(s.eye(2), s.zeros(2), epsilon) == 2 * epsilon / (
        1 + epsilon
    )
    assert sew.continued_rate(s.eye(2), s.zeros(2), epsilon) > 0


@pytest.mark.parametrize("dim", (4, 5, 6))
def test_complete_literal_polarization_inventory(dim):
    rows = tensor.data()["whole_literal_fixture_inventory"]
    assert [row["polarizations"] for row in rows if row["dimension"] == dim] == [
        dim * (dim - 3) // 2
    ] * 2


def test_every_dimension_change_uses_both_trace_occurrences():
    assert ward.data()["wrong_closed_trace_nonzero_witness"] != 0
    assert ward.data()["checks"]["general_D_channel_Ward"] == 0
    assert ward.data()["checks"]["general_TT_auxiliary_decomposition"] == 0


@pytest.mark.parametrize(
    "mutation",
    (
        "missing_leg",
        "zero_emitted",
        "wrong_mass",
        "wrong_conservation",
        "above_resolution",
        "negative_outgoing",
    ),
)
def test_physical_kinematics_guard(mutation):
    old, k, _ = old_vertices.sample(0)
    ps = [p.copy() for p in old]
    k = k.copy()
    if mutation == "missing_leg":
        ps.pop()
    elif mutation == "zero_emitted":
        k = s.zeros(4, 1)
    elif mutation == "wrong_mass":
        ps[0][0] -= 1
    elif mutation == "wrong_conservation":
        k *= s.Rational(1, 2)
    elif mutation == "above_resolution":
        k *= 100
    elif mutation == "negative_outgoing":
        ps[2] = -ps[2]
    with pytest.raises(ValueError):
        tensor.kinematics(ps, k)


@pytest.mark.parametrize(
    "bad",
    (
        None,
        s.zeros(4, 1),
        s.Matrix([1, 0, 0, 0]),
        s.Matrix([-1, 0, 0, -1]),
        s.Matrix([1, 0, 0]),
        s.Matrix([s.Float(1), 0, 0, 1]),
    ),
)
def test_exact_transverse_frame_guard(bad):
    with pytest.raises((TypeError, ValueError)):
        tensor.frame(bad)


def test_original_source_and_separate_physical_calibration():
    assert source.KAPPA == s.Integer(10) ** 800
    assert source.HEAVY_MASS2 == s.Integer(10) ** 200 / 512 + 2
    assert source.CUBIC == s.Rational(1, 8192)
    assert (
        tensor.data()["checks"]["original_physical_source_D5_all47_all5_polarizations"]
        == 0
    )
    assert tensor.data()["checks"]["original_regulator_zero_rate"] == 0


@pytest.mark.parametrize("norms", ((-1, 1), (1, -1), (-1, -1)))
def test_norm_majorants_must_be_nonnegative(norms):
    with pytest.raises(ValueError):
        sew.regulator_rate_upper(*norms, s.Rational(1, 8))


def test_historical_frontier_is_not_closed():
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 190
    assert len(audit.qualifications()) == 6
    assert audit.matching()[:-1] == audit.previous.matching()
    assert "integrated" in audit.observable()["not_established"].lower()
    assert "loop" in audit.observable()["not_established"].lower()
