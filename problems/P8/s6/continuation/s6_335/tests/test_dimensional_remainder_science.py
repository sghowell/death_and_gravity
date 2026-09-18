"""Exact API guards, two-window budgets and unchanged-source checks."""

from itertools import product

import pytest
import sympy as s
from p8_vacuum_affine_dimensional_real_remainder import (
    audit,
    born,
    bounds,
    calibration,
    soft,
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
def test_every_written_proof_and_calibration_gate(name):
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
def test_exact_domains(bad):
    for call, args in (
        (born.domain, (bad, 0)),
        (born.domain, (9, bad)),
        (born.original, (9, 0, bad)),
        (bounds.resolution, (bad,)),
        (bounds.gap, (bad,)),
        (bounds.folded_weight, (bad, 0)),
        (bounds.folded_weight, (0, bad)),
        (bounds.pair_density, (bad, 0, 0)),
        (bounds.pair_density, (s.Rational(5, 4), bad, 0)),
        (bounds.pair_density, (s.Rational(5, 4), 0, bad)),
    ):
        with pytest.raises((TypeError, ValueError)):
            call(*args)


@pytest.mark.parametrize("ss,z", ((6, 0), (17, 0), (9, -1), (9, 1), (9, 2), (9, -2)))
def test_nonforward_original_hard_domain(ss, z):
    for call, args in (
        (born.domain, (ss, z)),
        (born.auxiliary, (ss, z)),
        (born.original, (ss, z, 0)),
        (bounds.real_remainder_upper, (ss, z, s.Rational(1, 8))),
    ):
        with pytest.raises(ValueError):
            call(*args)


@pytest.mark.parametrize("bad", (-1, -s.Rational(1, 100), s.Rational(1, 7), 1))
def test_regulator_interval(bad):
    for call, args in (
        (born.original, (9, 0, bad)),
        (bounds.folded_weight, (s.Rational(1, 2), bad)),
        (bounds.pair_density, (s.Rational(5, 4), s.Rational(1, 128), bad)),
    ):
        with pytest.raises(ValueError):
            call(*args)


@pytest.mark.parametrize("bad", (0, -1, s.Rational(101, 100), 2))
def test_transfer_gap_must_be_positive_and_bounded(bad):
    for call in (
        bounds.gap,
        bounds.unit_D6_remainder_upper,
        bounds.canonical_gravity_remainder_upper,
        bounds.canonical_gravity_high_upper,
    ):
        with pytest.raises(ValueError):
            call(bad)


@pytest.mark.parametrize(
    "ss,z,e",
    tuple(
        product(
            (s.Rational(25, 4), s.Integer(9), s.Integer(16)),
            (
                -s.Rational(999, 1000),
                s.S.Zero,
                s.Rational(999, 1000),
                1 - s.Rational(1, 10**220),
            ),
            (s.S.Zero, s.Rational(1, 32), s.Rational(1, 8)),
        )
    ),
)
def test_same_D_Born_positivity_and_uniform_reference_ratio(ss, z, e):
    delta = born.transfer_gap(ss, z)
    aux = born.auxiliary(ss, z)
    g0 = born.stripped_gravity(ss, z, 0)
    ratio = born.ratio_to_four(ss, z, e)
    assert 0 < delta <= 1 and aux > 0 and g0 > 8 / delta
    assert 1 <= ratio <= s.Rational(19, 18)
    assert born.original(ss, z, e) > 0
    assert s.factor(born.stripped_gravity(ss, z, e) - g0 - 2 * e * aux / (1 + e)) == 0
    assert source.HEAVY_MASS2 == s.Integer(10) ** 200 / 512 + 2
    if e == 0:
        assert ratio == 1


@pytest.mark.parametrize(
    "ss,z",
    tuple(
        product(
            (s.Rational(25, 4), s.Integer(9), s.Integer(16)),
            (-s.Rational(99, 100), s.S.Zero, s.Rational(99, 100)),
        )
    ),
)
def test_low_high_budget_split_and_zero_limit(ss, z):
    delta = born.transfer_gap(ss, z)
    a = delta / 192
    last = s.S.Zero
    for x in (s.S.Zero, a / 2, a, s.Rational(1, 8)):
        value = bounds.real_remainder_upper(ss, z, x)
        low = bounds.low_real_remainder_upper(ss, z, min(x, a))
        high = (
            0
            if x <= a
            else s.Integer(10) ** 35 * (x * x - a * a) / (source.KAPPA * delta**2)
        )
        assert s.simplify(value - low - high) == 0
        assert not value.has(s.Float)
        assert value >= last
        coarse = (
            s.Integer(10) ** 13 + s.Integer(10) ** 35 * x * x / delta**2
        ) / source.KAPPA
        # Replacing pi by3 only increases this positive upper bound.
        assert s.factor(coarse - value.subs(s.pi, 3)) > 0
        if x == 0:
            assert value == 0
        if x <= a:
            assert value == bounds.lower_cutoff_error_upper(ss, z, x)
        last = value
    with pytest.raises(ValueError):
        bounds.low_real_remainder_upper(ss, z, a * 2)
    assert bounds.real_remainder_upper(ss, z, a) == bounds.low_real_remainder_upper(
        ss, z, a
    )


@pytest.mark.parametrize(
    "E,w,e",
    tuple(
        product(
            (s.Rational(5, 4), s.Integer(2)),
            (s.S.Zero, s.Rational(1, 128), s.Rational(1, 8)),
            (s.S.Zero, s.Rational(1, 32), s.Rational(1, 8)),
        )
    ),
)
def test_actual_continued_pair_phase(E, w, e):
    value = bounds.pair_density(E, w, e)
    assert value > 0 and value <= 1 and 1 - value <= 3 * w
    assert not value.has(s.Float)
    if w == 0:
        assert value == 1


@pytest.mark.parametrize(
    "v,e",
    tuple(
        product(
            (s.S.Zero, s.Rational(1, 100), s.Rational(1, 2), s.S.One),
            (s.S.Zero, s.Rational(1, 32), s.Rational(1, 8)),
        )
    ),
)
def test_actual_positive_folded_measure(v, e):
    value = bounds.folded_weight(v, e)
    assert 0 <= value <= s.Rational(5, 4)
    assert not value.has(s.Float)
    if e == 0:
        assert value == 1


@pytest.mark.parametrize("bad", (-1, s.Rational(1, 7), 1))
def test_resolution_domain(bad):
    with pytest.raises(ValueError):
        bounds.real_remainder_upper(9, 0, bad)


@pytest.mark.parametrize("bad", (-1, s.Rational(101, 100), 2))
def test_folded_angular_domain(bad):
    with pytest.raises(ValueError):
        bounds.folded_weight(bad, 0)


@pytest.mark.parametrize("bad", (0, 1, s.Rational(201, 100), 3))
def test_pair_energy_domain(bad):
    with pytest.raises(ValueError):
        bounds.pair_density(bad, s.Rational(1, 128), 0)


def test_canonical_trace_and_both_windows_are_present():
    values = soft.component_budgets()
    assert values["D6_unit_remainder"] == 13396352288
    assert values["canonical_remainder"] == 26792749792
    assert bounds.canonical_gravity_remainder_upper(1) == 30000000000
    assert bounds.canonical_gravity_high_upper(1) == s.Integer(9000000000) * 30000**2
    inventory = calibration.data()["whole_calibration_domain"]
    assert inventory["low_recoil_states"] == 4 and inventory["high_recoil_states"] == 6
    assert inventory["all_original_source_parameters_retained"] is True


def test_source_and_original_frontier_are_not_changed():
    assert source.KAPPA == s.Integer(10) ** 800 and source.CUBIC == s.Rational(1, 8192)
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 191
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.qualifications()) == 6
    assert "virtual" in audit.observable()["not_established"].lower()
    assert "forward" in audit.observable()["not_established"].lower()
