"""Curvature matching, explicit coefficient guards and conditional rate bounds."""

import pytest
import sympy as s
from p8_vacuum_affine_minimal_gravity_radiation import vertices as old
from p8_vacuum_affine_radiative_curvature_matching import (
    audit,
    bounds,
    calibration,
    contact,
    matching,
    source,
)

ROWS, GATES = audit.residuals(), audit.gates()


@pytest.mark.parametrize("name", tuple(ROWS))
def test_every_exact_residual(name):
    value = ROWS[name]
    assert all(
        v == 0 for v in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    assert not s.sympify(value).has(s.Float)


@pytest.mark.parametrize("name", tuple(GATES))
def test_every_proof_gate(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[row[0] for row in audit.bad_cases()]
)
def test_all_rejected_scopes(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize(
    "bad",
    (
        None,
        True,
        False,
        "0",
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
def test_explicit_exact_coefficient_required(bad):
    with pytest.raises((TypeError, ValueError)):
        contact.exact_coefficient(bad)
    ps, k, _ = old.sample(0)
    for call, args in (
        (matching.linear_interference, (ps, k, bad)),
        (matching.contact_square, (ps, k, bad)),
        (bounds.low_interference_upper, (9, 0, s.Rational(1, 1000), bad)),
        (bounds.contact_square_upper, (9, 0, s.Rational(1, 8), bad)),
    ):
        with pytest.raises((TypeError, ValueError)):
            call(*args)


@pytest.mark.parametrize("chi", (-3, s.Rational(-1, 7), 0, s.Rational(1, 5), 2))
def test_supplied_comparison_values_not_defaults(chi):
    assert contact.exact_coefficient(chi) == chi
    ps, k, _ = old.sample(0)
    linear = matching.linear_interference(ps, k, chi)
    square = matching.contact_square(ps, k, chi)
    assert linear == chi * matching.linear_interference(ps, k, 1)
    assert square == chi**2 * matching.contact_square(ps, k, 1)
    assert not linear.has(s.Float) and not square.has(s.Float)
    assert square >= 0
    if chi == 0:
        assert linear == square == 0


@pytest.mark.parametrize("index", (0, 1, 2))
def test_whole_original_tree_interference_and_Ward(index):
    ps, k, _ = old.sample(index)
    T = contact.core(ps, k)
    assert T * contact.ETA * k == s.zeros(4, 1)
    assert matching.linear_interference(ps, k, 1) != 0
    assert matching.linear_interference(ps, k, -1) == -matching.linear_interference(
        ps, k, 1
    )
    assert matching.contact_square(ps, k, 1) > 0
    assert (
        matching.contact_square(ps, k, 1)
        <= bounds.tensor_upper(k[0]) ** 2 / source.KAPPA
    )


@pytest.mark.parametrize(
    "ss,z", ((s.Rational(25, 4), 0), (9, s.Rational(9, 10)), (16, s.Rational(-99, 100)))
)
def test_conditional_rate_bounds_and_soft_limits(ss, z):
    delta = bounds.born.transfer_gap(ss, z)
    x = delta / 384
    for chi in (-2, 0, 3):
        first = bounds.low_interference_upper(ss, z, x, chi)
        second = bounds.contact_square_upper(ss, z, x, chi)
        assert first >= 0 and second >= 0
        assert first == bounds.low_interference_upper(ss, z, x, -chi)
        assert second == bounds.contact_square_upper(ss, z, x, -chi)
        assert bounds.low_interference_upper(ss, z, 0, chi) == 0
        assert bounds.contact_square_upper(ss, z, 0, chi) == 0
        assert not first.has(s.Float) and not second.has(s.Float)
        if chi:
            assert bounds.low_interference_upper(ss, z, x / 2, chi) <= first / 8
            assert bounds.contact_square_upper(ss, z, x / 2, chi) == second / 64
    with pytest.raises(ValueError):
        bounds.low_interference_upper(ss, z, delta / 96, 1)


@pytest.mark.parametrize("bad", (-1, s.Rational(1, 7), 1, s.Float(0), True, "0"))
def test_actual_resolution_domain(bad):
    for call, args in (
        (bounds.tensor_upper, (bad,)),
        (bounds.contact_square_upper, (9, 0, bad, 1)),
    ):
        with pytest.raises((TypeError, ValueError)):
            call(*args)


@pytest.mark.parametrize("ss,z", ((6, 0), (17, 0), (9, 1), (9, -1)))
def test_original_nonforward_Born_domain(ss, z):
    with pytest.raises(ValueError):
        bounds.low_interference_upper(ss, z, 0, 1)
    with pytest.raises(ValueError):
        bounds.contact_square_upper(ss, z, 0, 1)


def test_no_coefficient_zero_default():
    ps, k, _ = old.sample(0)
    for call, args in (
        (matching.linear_interference, (ps, k)),
        (matching.contact_square, (ps, k)),
        (bounds.low_interference_upper, (9, 0, 0)),
        (bounds.contact_square_upper, (9, 0, 0)),
    ):
        with pytest.raises(TypeError):
            call(*args)


def test_original_state_counts_and_boundary():
    assert [
        row["original_tree_interference_sign"]
        for row in calibration.data()["whole_original_physical_calibrations"]
    ] == [-1, 1, 1]
    assert source.KAPPA == s.Integer(10) ** 800 and source.CUBIC == s.Rational(1, 8192)
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 192
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.qualifications()) == 6
    assert "physical" in audit.observable()["not_established"].lower()
    assert "curvature" in audit.observable()["not_established"].lower()
