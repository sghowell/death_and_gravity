"""Independent tensor, reference-flow and actual-window enclosure tests."""

from fractions import Fraction

import pytest
import sympy as sp
from p8_polynomial_vacuum_running import audit, calibration, flow, tensors


@pytest.mark.parametrize(
    "name,value", list(audit.residuals().items()), ids=list(audit.residuals())
)
def test_exact_tensor_reference_or_flow_identity(name, value):
    assert value == 0, name


@pytest.mark.parametrize(
    "name,value", list(audit.gates().items()), ids=list(audit.gates())
)
def test_continuous_bound_or_explicit_order_gate(name, value):
    assert bool(value), name


@pytest.mark.parametrize(
    "name,call,args",
    calibration.bad_cases(),
    ids=[r[0] for r in calibration.bad_cases()],
)
def test_reject_unsupported_reference_time(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


@pytest.mark.parametrize(
    "name,value",
    [(k, v) for k, v in audit.controls().items() if k != "rejected_inputs"],
    ids=[k for k in audit.controls() if k != "rejected_inputs"],
)
def test_negative_control_and_unclosed_obligation(name, value):
    assert bool(value), name


@pytest.mark.parametrize("ell", (0, 1, 9, Fraction(1, 2), sp.Rational(9, 2)))
def test_exact_reference_window_enclosures(ell):
    d = calibration.point(ell)
    assert d["logarithmic_reference_time"] == sp.Rational(ell)
    assert sp.Rational(3, 4) < d["flow_root_cubed"] <= 1
    assert 0 < d["cubic_squared_lower"] <= d["cubic_squared_upper"]
    assert 0 < d["heavy_mass_squared_lower"] <= d["heavy_mass_squared_upper"]
    assert 0 < d["completed_square_margin_lower"] <= d["completed_square_margin_upper"]
    assert 0 <= d["margin_relative_increase_upper"] < sp.Rational(1, 10**5)
    if ell:
        assert d["completed_square_margin_lower"] < d["completed_square_margin_upper"]


def test_validation_not_bypassed_by_equal_cached_float():
    calibration.point(1)
    for value in (True, 1.0, sp.Float(1)):
        with pytest.raises((TypeError, ValueError)):
            calibration.point(value)


def test_tensor_component_counts_and_nonzero_species():
    d = tensors.data()
    B4, B3, B2 = (
        d["literal_quartic_beta_tensor"],
        d["literal_cubic_beta_tensor"],
        d["literal_mass_beta_tensor"],
    )
    assert len(B4) == 16 and len(B3) == 8 and len(B2) == 4
    assert sum(v != 0 for v in B4.values()) == 1
    assert sum(v != 0 for v in B3.values()) == 3
    assert sum(v != 0 for v in B2.values()) == 2


def test_cubic_squared_requires_factor_two():
    d = tensors.data()
    G = d["G"]
    assert (
        2 * G * d["literal_cubic_beta_tensor"][(0, 0, 1)]
        == d["one_loop_beta_in_log_reference_over_16pi2"]["cubic_squared"]
    )
    assert (
        G * d["literal_cubic_beta_tensor"][(0, 0, 1)]
        != d["one_loop_beta_in_log_reference_over_16pi2"]["cubic_squared"]
    )


def test_light_UV_mass_counterterm_retained_but_not_physical_running():
    d = tensors.data()
    assert (
        d["literal_mass_beta_tensor"][(0, 0)]
        == d["L"] * d["light_mass_squared"] + 2 * d["G"] ** 2
    )
    assert "kept one" in d["light_mass_convention"]


@pytest.mark.parametrize(
    "q", (sp.Rational(1, 2), sp.Rational(3, 4), sp.Rational(9, 10), 1)
)
def test_general_closed_flow_positive_margin_diagnostic(q):
    d = flow.data()
    sub = {
        d["q"]: q,
        d["L0"]: sp.Rational(1, 100),
        d["g0"]: sp.Rational(1, 1000),
        d["M0"]: 1,
    }
    value = d["running_completed_square_margin"].subs(sub)
    assert value >= sp.Rational(7, 1000) > 0


def test_exact_initial_parameter_enclosure():
    d = calibration.point(0)
    p = calibration.data()
    assert (
        d["exact_running_polynomial_quartic"] == p["actual_initial_polynomial_quartic"]
    )
    assert (
        d["heavy_mass_squared_lower"]
        == d["heavy_mass_squared_upper"]
        == p["actual_initial_heavy_mass_squared"]
    )
    assert (
        d["completed_square_margin_lower"]
        == d["completed_square_margin_upper"]
        == p["actual_initial_positive_margin"]
    )
    assert d["margin_relative_increase_upper"] == 0


def test_actual_uniform_increments_are_nonzero_and_bounded():
    d = calibration.data()
    assert 0 < d["whole_window_margin_relative_increase_upper"] < sp.Rational(1, 10**5)
    assert (
        0 < d["whole_window_heavy_mass_squared_increase_upper"] < sp.Rational(3, 10**7)
    )
    assert (
        0 < d["whole_window_quartic_relative_increase_upper"] < sp.Rational(1, 10**202)
    )
    assert (
        0
        < d["whole_window_cubic_squared_relative_increase_upper"]
        < sp.Rational(1, 10**202)
    )


def test_formal_flow_boundary_not_quantum_Landau_claim():
    assert "not a proven quantum Landau pole" in flow.data()["scope"]


def test_exact_counts_and_original_scope():
    assert len(audit.residuals()) == 61
    assert len(audit.gates()) == 32
    assert len(audit.controls()) == 10
    assert audit.rejected_inputs() == 16
