"""Exact Gramian engine and literal source-chart checks."""
from copy import deepcopy
from fractions import Fraction as Q

import pytest
from flint import ctx
from p8_preparation_cost import core, gramian, independent


def test_literal_generator_and_continuous_bounds():
    checks = core.checks()
    assert len(checks["residuals"]) == 44
    assert all(value == 0 for value in checks["residuals"].values())
    assert all(value > 0 for value in checks["strict_margins"].values())


def test_separate_fraction_energy_arithmetic():
    primary, separate = core.calibration(), independent.checks()
    assert primary.keys() == separate.keys()
    assert all(primary[key] == value for key, value in separate.items())


def test_literal_rational_to_arb_jet_bridges():
    rows = core.coefficient_bridges()
    assert len(rows) == 4
    assert sum(row["comparisons"] for row in rows) == 256


def test_validated_gramian_has_full_remainder_and_pivots():
    result = gramian.certificate()
    assert result["steps"] == 64
    assert result["order"] == 18
    assert result["step_matrix_norm_upper"] < 26
    assert 0 < result["local_entry_remainder_upper"] < Q(1, 10**27)
    assert all(0 < lo <= hi for lo, hi in result["shifted_LDL_pivot_intervals"])
    assert result["full_L2_du_scaling"] == Q(4, 10**6)


def test_independent_rational_schur_proof():
    result = independent.positive_gramian(gramian.certificate()["scaled_outer_gramian"])
    assert result["strict_positive_pivots"] == 4
    assert all(value > 0 for value in result["pivot_lower_bounds"])


def test_changed_step_and_order_enclosures_overlap():
    primary = gramian.certificate()
    other = gramian.enclosure(order=16, steps=128, precision=192)
    for key in ("normalized_physical_gramian", "scaled_outer_gramian", "endpoint_transform"):
        for i in range(4):
            for j in range(4):
                a, b = primary[key][i][j], other[key][i][j]
                assert max(a[0], b[0]) <= min(a[1], b[1])


def test_process_precision_is_restored():
    old = ctx.prec, ctx.cap
    try:
        ctx.prec, ctx.cap = 137, 7
        gramian.enclosure(order=16, steps=128, precision=192)
        assert (ctx.prec, ctx.cap) == (137, 7)
    finally:
        ctx.prec, ctx.cap = old


@pytest.mark.parametrize("kwargs", [
    {"order": True}, {"order": 18.0}, {"order": 11}, {"order": 25},
    {"steps": True}, {"steps": 64.0}, {"steps": 0}, {"steps": 257},
    {"precision": True}, {"precision": 256.0}, {"precision": 127}, {"precision": 513},
])
def test_invalid_numerical_controls_rejected_before_replay(kwargs):
    with pytest.raises((TypeError, ValueError)):
        gramian.enclosure(**kwargs)


def test_tampered_gramian_is_not_certified():
    matrix = deepcopy(gramian.certificate()["scaled_outer_gramian"])
    matrix[0][0] = (Q(0), Q(0))
    with pytest.raises(ValueError):
        independent.positive_gramian(matrix)


@pytest.mark.parametrize("args", [(True,), (0.0,), (Q(2), Q(1))])
def test_exact_interval_guards(args):
    with pytest.raises((TypeError, ValueError)):
        independent.Interval(*args)
