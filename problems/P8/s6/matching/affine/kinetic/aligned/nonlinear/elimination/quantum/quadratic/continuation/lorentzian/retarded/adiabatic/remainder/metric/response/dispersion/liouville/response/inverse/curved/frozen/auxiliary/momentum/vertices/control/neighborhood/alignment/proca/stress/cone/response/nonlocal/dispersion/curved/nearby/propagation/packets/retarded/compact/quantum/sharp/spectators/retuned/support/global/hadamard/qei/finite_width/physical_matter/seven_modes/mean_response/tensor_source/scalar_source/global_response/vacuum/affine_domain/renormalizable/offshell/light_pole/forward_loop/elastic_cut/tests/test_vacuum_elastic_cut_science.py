"""Independent physical cut, threshold and one-loop improved-margin tests."""

from fractions import Fraction

import pytest
import sympy as sp
from p8_vacuum_elastic_cut import (
    amplitude,
    angular,
    audit,
    calibration,
    cut,
    normalization,
)
from p8_vacuum_forward_loop import calibration as loop


@pytest.mark.parametrize(
    "name,value", list(audit.residuals().items()), ids=list(audit.residuals())
)
def test_exact_physical_amplitude_angular_cut_or_phase_identity(name, value):
    assert value == 0, name


@pytest.mark.parametrize(
    "name,value", list(audit.gates().items()), ids=list(audit.gates())
)
def test_continuous_domain_bound_and_written_proof_gate(name, value):
    assert bool(value), name


@pytest.mark.parametrize(
    "name,call,args",
    calibration.bad_cases(),
    ids=[r[0] for r in calibration.bad_cases()],
)
def test_reject_unsupported_invariant(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


@pytest.mark.parametrize(
    "name,value",
    [(k, v) for k, v in audit.controls().items() if k != "rejected_inputs"],
    ids=[k for k in audit.controls() if k != "rejected_inputs"],
)
def test_negative_control_and_scope(name, value):
    assert bool(value), name


@pytest.mark.parametrize("s", (4, 5, 6, Fraction(9, 2), sp.Rational(11, 2)))
def test_exact_invariant_enclosures(s):
    d = calibration.point(s)
    assert d["physical_invariant"] == sp.Rational(s)
    assert sp.factor(d["two_body_phase_beta"] ** 2 - (1 - 4 / sp.Rational(s))) == 0
    if s == 4:
        assert d["cut_density_lower"] == d["cut_density_upper"] == 0
        assert d["threshold_density_exactly_zero"]
    else:
        assert 0 < d["cut_density_lower"] < d["cut_density_upper"]
        assert d["cut_density_strictly_positive"]


def test_validation_not_bypassed_by_equal_cached_float():
    calibration.point(5)
    for v in (True, 5.0, sp.Float(5)):
        with pytest.raises((ValueError, TypeError)):
            calibration.point(v)


def test_actual_amplitude_extrema_are_strict_not_rounded():
    d = amplitude.data()
    l = d["actual_fixed_lambda"]
    assert 23 * l < d["actual_threshold_minimum"] < 24 * l
    assert 42 * l < d["actual_positive_subwindow_minimum"] < 44 * l
    assert 71 * l < d["actual_full_window_maximum"] < 73 * l


@pytest.mark.parametrize(
    "s,z", ((4, -1), (5, 0), (6, 1), (sp.Rational(9, 2), sp.Rational(-1, 2)))
)
def test_independent_physical_s_monotonicity(s, z):
    d = amplitude.data()
    A = d["physical_tree_amplitude"]
    value = sp.diff(A, d["s"]).subs(
        {d["s"]: s, d["cosine"]: z, d["M"]: 10, d["g"]: 1, d["lambda4"]: 1}
    )
    assert value > 0


def test_full_angular_threshold_limit_and_branch():
    d = angular.data()
    assert (
        sp.factor(
            sp.limit(d["exact_even_squared_vertex_average"], d["k"], 0, dir="+")
            - d["threshold_average"]
        )
        == 0
    )
    assert d["anchored_squared_vertex_primitive"].subs(d["cosine"], 0) == 0


def test_angular_interference_retained_in_closed_average():
    d = angular.data()
    Q = d["exact_even_squared_vertex_average"]
    assert sp.diff(Q, d["C"], d["g"]) == 4 * sp.atanh(d["k"] / d["B"]) / d["k"]


def test_phase_space_identical_and_two_im_factors_are_distinct():
    d = normalization.data()
    assert (
        sp.factor(
            d["two_body_phase_space_per_cosine"]
            / d["one_loop_imaginary_part_per_cosine_squared_vertex"]
        )
        == 4
    )


def test_computed_low_cut_strictly_positive_and_retained():
    d = cut.data()
    l = amplitude.data()["actual_fixed_lambda"]
    assert d["actual_strict_cut_lower"] == l * l / 20 > 0
    assert d["actual_cut_upper"] == 3 * l * l
    assert d["selected_low_energy_cut"].limits[0][1:] == (4, 6)


def test_actual_complete_loop_error_and_cut_both_subtracted():
    d = cut.data()
    p = loop.point()
    assert (
        d["actual_tree_plus_one_loop_minus_cut_lower"]
        == p["actual_tree_b2"]
        - p["total_one_loop_b2_error_upper"]
        - d["actual_cut_upper"]
    )
    assert (
        0
        < d["actual_tree_plus_one_loop_minus_cut_lower"]
        < p["tree_plus_one_loop_b2_lower"]
    )


def test_exact_counts_and_original_scope():
    assert len(audit.residuals()) == 33
    assert len(audit.gates()) == 32
    assert len(audit.controls()) == 10
    assert audit.rejected_inputs() == 16
