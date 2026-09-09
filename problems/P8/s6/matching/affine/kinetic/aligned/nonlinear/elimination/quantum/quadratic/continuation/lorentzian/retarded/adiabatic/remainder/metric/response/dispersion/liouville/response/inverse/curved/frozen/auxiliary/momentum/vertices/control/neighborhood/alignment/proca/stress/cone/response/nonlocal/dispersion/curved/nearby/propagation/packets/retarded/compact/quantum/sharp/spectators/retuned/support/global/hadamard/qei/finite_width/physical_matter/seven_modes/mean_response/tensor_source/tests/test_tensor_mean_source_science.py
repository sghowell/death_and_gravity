"""Actual tensor action/source, nonlinear chart, clock exchange and global controls."""

from fractions import Fraction

import pytest
import sympy as sp
from p8_tensor_mean_source import (
    audit,
    bounds,
    calibration,
    chart,
    operator,
    response,
    state,
    verify,
)

ROWS = audit.residuals()


@pytest.mark.parametrize("name", list(ROWS))
def test_native_exact_identity(name):
    value = ROWS[name]
    assert all(
        v == 0 for v in (list(value) if isinstance(value, sp.MatrixBase) else [value])
    ), name


@pytest.mark.parametrize("name", list(audit.gates()))
def test_exact_and_written_proof_gate(name):
    assert bool(audit.gates()[name]), name


@pytest.mark.parametrize(
    "name,call,args",
    calibration.bad_cases(),
    ids=[v[0] for v in calibration.bad_cases()],
)
def test_invalid_tensor_covariance_amplitude(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_source_and_scientific_manifest():
    assert len(ROWS) == 53
    assert (
        sum(
            v.rows * v.cols if isinstance(v, sp.MatrixBase) else 1
            for v in ROWS.values()
        )
        == 89
    )
    assert len(audit.gates()) == 37
    assert audit.rejected_inputs() == 16
    assert len(verify.source_files()) == 20


@pytest.mark.parametrize("value", (0, Fraction(1, 10**20), sp.Rational(1, 10**21)))
def test_exact_valid_calibration_and_zero_response(value):
    d = calibration.calibrated(value)
    assert d["eta"] == sp.Rational(value)
    assert all(v >= 0 for v in d.values())
    if value == 0:
        assert all(v == 0 for v in d.values())


@pytest.mark.parametrize("bad", (0.1, sp.Float("0.1"), sp.oo, sp.nan))
def test_inexact_or_nonfinite_report_value_rejected(bad):
    with pytest.raises((TypeError, ValueError)):
        verify.serialize({"bad": bad})


@pytest.mark.parametrize(
    "name", [k for k in audit.controls() if k != "rejected_inputs"]
)
def test_actual_negative_and_scope_control(name):
    assert bool(audit.controls()[name]), name


def test_source_is_neither_minimal_Proca_nor_positive_proxy():
    d = operator.data()
    r, s = d["positive_on_clock_tensor_energy"], d["on_clock_spatial_pressure"]
    assert sp.factor(d["actual_hat_lapse_source"] - (r - 3 * s / (2 * operator.h))) != 0
    assert sp.factor(d["fixed_physical_metric_lapse_source"] - r) != 0


def test_two_polarization_anchor_normalization_is_not_one_or_nine_modes():
    d = state.data()
    assert sp.factor(d["two_polarization_anchor_density"] - (operator.q + 4) / 2) == 0
    assert sp.factor(d["two_polarization_anchor_pressure"] - (12 - operator.q) / 6) == 0
    T, S = sp.symbols("T S", real=True)
    gamma = sp.Matrix([[T, S, 0], [S, -T, 0], [0, 0, 0]])
    assert sp.expand(sp.trace(gamma * gamma) - 2 * (T * T + S * S)) == 0


def test_center_chart_correction_is_nonzero_and_has_required_sign():
    d = chart.data()
    assert d["induced_linear_homogeneous_scalar_density"] != 0
    assert (
        sp.factor(
            d["additive_tensor_density_second_variation"]
            + d["induced_linear_homogeneous_scalar_density"]
            - d["actual_exponential_chart_density_response"]
        )
        == 0
    )


def test_selected_tensor_band_has_actual_clock_exchange():
    c = response.center()
    q = operator.q
    assert c["anchor_clock_source_first_derivative"].subs(q, 1) == 82
    assert c["anchor_clock_source_first_derivative"].subs(q, 4) == 64
    d = bounds.data()
    assert d[
        "absolute_volume_weighted_clock_source_integral_upper_per_eta"
    ] == sp.Rational(1944, 7)
    assert d["whole_line_volume_weighted_clock_source_integral"] == 0
    future = d["future_volume_weighted_clock_source_integral_per_eta"]
    mu = next(iter(future.free_symbols))
    assert future.subs(mu, sp.Rational(5, 2)) == sp.Rational(3, 2)


def test_new_negative_center_acceleration_does_not_import_Proca_positivity():
    d = response.center()
    q = operator.q
    assert d["anchor_lapse_response_second_derivative"].subs(q, 1) < 0
    assert d["anchor_proper_Hubble_derivative_response"].subs(q, 1) > -15
    assert d["anchor_proper_Hubble_derivative_response"].subs(q, 4) < -9
    assert (
        bounds.data()["center_lapse_second_derivative_absolute_coefficient_bound"] < 30
    )


def test_actual_forcing_keeps_constraint_lapse_and_induced_matter_field():
    d = response.data()
    assert d["actual_tensor_induced_lapse"] != 0
    assert d["induced_matter_field_derivative"] != 0
    assert response.ward()["retained_background_matter_connection_variation"] != 0


def test_global_and_clock_bounds_use_new_tensor_energy_not_Proca_mass():
    d = bounds.data()
    assert d["source_specific_proxy_density_anchor_upper"] == 4
    assert d["global_weighted_mean_phase_upper_per_eta"] == 6400
    assert d["global_lapse_upper_per_eta"] == 1500
    assert d["global_actual_frame_log_scale_upper_per_eta"] == 8000
    assert (
        d["rounded_absolute_volume_weighted_clock_source_integral_upper_per_eta"] == 280
    )


def test_report_science_serializes_with_native_exact_values():
    for data in (
        operator.data(),
        operator.clock(),
        chart.data(),
        chart.metric(),
        state.data(),
        state.variational_clock(),
        response.data(),
        response.center(),
        response.ward(),
        bounds.data(),
    ):
        verify.serialize(verify.payload(data))
    verify.serialize(audit.controls())
    verify.serialize(calibration.calibrated())


def test_geometric_representative_not_promoted_to_full_quantum_solution():
    assert bounds.data()["complete_first_order_frame_representative_not_a_SEE_solution"]
    assert audit.controls()[
        "no_finite_amplitude_SEE_or_coevolved_quantum_state_claimed"
    ]
    assert audit.controls()[
        "no_full_covariant_graviton_or_BRST_renormalized_stress_claimed"
    ]
