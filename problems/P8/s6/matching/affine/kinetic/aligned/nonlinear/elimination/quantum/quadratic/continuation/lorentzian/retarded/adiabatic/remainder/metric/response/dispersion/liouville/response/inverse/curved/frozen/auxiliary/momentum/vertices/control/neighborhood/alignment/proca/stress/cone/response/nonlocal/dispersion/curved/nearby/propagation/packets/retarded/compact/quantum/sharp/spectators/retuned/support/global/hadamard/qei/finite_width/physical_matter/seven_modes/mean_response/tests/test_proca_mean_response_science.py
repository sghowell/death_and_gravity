"""Actual stress, full mean source, global bounds and strict scope controls."""

from fractions import Fraction

import pytest
import sympy as sp
from p8_proca_mean_response import anchor, audit, band, mean, stress, tails, verify

ROWS = audit.residuals()


@pytest.mark.parametrize("name", list(ROWS))
def test_native_exact_identity(name):
    v = ROWS[name]
    assert all(x == 0 for x in (list(v) if isinstance(v, sp.MatrixBase) else [v])), name


@pytest.mark.parametrize("name", list(audit.gates()))
def test_written_and_exact_proof_gate(name):
    assert bool(audit.gates()[name]), name


@pytest.mark.parametrize(
    "name,call,args", band.bad_cases(), ids=[v[0] for v in band.bad_cases()]
)
def test_invalid_relative_covariance_amplitude(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_source_and_scientific_manifest():
    assert len(ROWS) == 65
    assert (
        sum(
            v.rows * v.cols if isinstance(v, sp.MatrixBase) else 1
            for v in ROWS.values()
        )
        == 451
    )
    assert len(audit.gates()) == 47
    assert audit.rejected_inputs() == 16
    assert len(verify.source_files()) == 20


@pytest.mark.parametrize("value", (0, Fraction(1, 10**20), sp.Rational(1, 10**21)))
def test_exact_valid_calibration_and_zero_response(value):
    d = band.calibrated(value)
    assert d["eta"] == sp.Rational(value)
    assert all(v >= 0 for v in d.values())
    if value == 0:
        assert all(v == 0 for v in d.values())


@pytest.mark.parametrize("bad", (0.1, sp.Float("0.1"), sp.oo, sp.nan))
def test_inexact_or_nonfinite_report_values_rejected(bad):
    with pytest.raises((TypeError, ValueError)):
        verify.serialize({"bad": bad})


def test_positive_density_pieces_have_different_pressure_weights():
    d = stress.data()
    assert len(d["four_positive_density_pieces"]) == 4
    assert d["physical_isotropic_pressure_Hessian"] != d["physical_density_Hessian"]


def test_lapse_is_spatial_only_and_pressure_cannot_be_omitted():
    c = audit.controls()
    assert c["wrong_four_dimensional_lapse_rescaling_extra_source_nonzero"]
    assert c["omitting_physical_pressure_changes_lapse_source"]


def test_center_mean_lapse_is_not_silently_reset_to_zero():
    d = mean.center_control()
    assert d["initial_induced_lapse"] != 0
    assert d[
        "initial_combined_density_lower_per_positive_vector_density"
    ] == sp.Rational(134, 135)


def test_full_matter_response_keeps_induced_field_and_connection():
    d = mean.data()
    w = mean.ward()
    assert d["matter_clock_field_response_derivative"] != 0
    assert w["retained_connection_variation_term"] != 0
    assert audit.controls()["omitting_connection_variation_changes_Ward_identity"]


def test_global_positive_pivot_is_coefficient_not_sampling_proof():
    d = tails.data()
    assert all(c >= 0 for c in d["positive_J_remainder_polynomial_coefficients"])
    assert d["positive_J_remainder_polynomial_coefficients"][-1] > 0


def test_global_weighted_growth_and_source_have_finite_exact_integrals():
    d = tails.data()
    assert d["weighted_generator_integrable_entry_sum_upper"] < 2
    assert d["weighted_forcing_L1_upper_per_initial_density"] < 170
    assert d["weighted_global_Gronwall_factor_upper"] == 9


def test_independent_center_acceleration_retains_massive_stress_jet():
    d = anchor.data()
    assert d["anchor_pressure_second_derivative_per_unit_band_covariance"] != 0
    assert (
        d["absolute_anchor_proper_Hubble_derivative_response_upper_per_eta"]
        < 5 * 10**11
    )
    assert d["all_anchor_acceleration_response_polynomial_coefficients_positive"]


def test_global_physical_lapse_and_exact_frame_calibration():
    d = tails.data()
    eta = sp.Rational(1, 10**20)
    assert d["global_lapse_upper_per_eta"] * eta < sp.Rational(1, 10**10)
    assert d["global_actual_frame_log_scale_upper_per_eta"] * eta < sp.Rational(
        1, 10**10
    )


def test_complete_metric_representative_is_not_promoted_to_quantum_solution():
    assert tails.data()["representative_is_not_an_exact_semiclassical_solution"]
    assert audit.controls()[
        "no_quantum_state_transported_on_representative_metric_claimed"
    ]
    assert audit.controls()["no_nonlinear_or_higher_loop_remainder_bound_inferred"]


def test_absolute_common_tadpole_not_replaced_by_zero():
    assert audit.controls()["no_absolute_vacuum_source_claimed_zero"]
    assert audit.controls()["reference_scalar_tensor_state_and_tadpole_not_reassigned"]


def test_actual_source_sign_has_independent_GR_constraint_control():
    H, r, s, l = sp.symbols("H r s ell", real=True)
    # Hcal_GR=-3p^2/4+ell^2/2+r; p=-2H.
    # On the lapse constraint, p_dot=3H^2+ell^2/2+s.
    pdot = (3 * H * H + l * l / 2 + s).subs(H * H, (l * l / 2 + r) / 3)
    assert sp.factor(-pdot / 2 + (l * l + r + s) / 2) == 0


def test_matter_field_response_is_odd_under_source_reflection():
    d = mean.data()
    derivative = d["matter_clock_field_response_derivative"]
    reflected = derivative.subs(mean.u, -mean.u).subs(mean.dp, -mean.dp)
    assert sp.factor(reflected - derivative) == 0
