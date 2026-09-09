"""Native whole-line phase, compensated physical response and joint bounds."""

from fractions import Fraction

import pytest
import sympy as sp
from p8_scalar_global_response import (
    audit,
    bounds,
    calibration,
    identities,
    joint,
    phase,
    source,
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
    ids=[row[0] for row in calibration.bad_cases()],
)
def test_invalid_exact_amplitude_integral_or_envelope_rejected(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_source_and_scientific_manifest():
    assert len(ROWS) == 47
    assert (
        sum(
            v.rows * v.cols if isinstance(v, sp.MatrixBase) else 1
            for v in ROWS.values()
        )
        == 203
    )
    assert len(audit.gates()) == 48
    assert audit.rejected_inputs() == 51
    assert len(verify.source_files()) == 20


@pytest.mark.parametrize("value", (0, Fraction(1, 10**80), sp.Rational(1, 10**81)))
def test_exact_scalar_amplitude_and_calibration(value):
    d = calibration.calibrated(value)
    assert d["eta"] == sp.Rational(value)
    assert all(v >= 0 for v in d.values())
    if value == 0:
        assert all(v == 0 for v in d.values())


@pytest.mark.parametrize(
    "power,expected",
    (
        (1, sp.Rational(11, 7)),
        (2, sp.Rational(11, 14)),
        (sp.Rational(3, 2), 1),
        (sp.Rational(5, 2), sp.Rational(2, 3)),
        (sp.Rational(7, 2), sp.Rational(8, 15)),
    ),
)
def test_exact_integer_and_half_integer_integrals(power, expected):
    assert phase.half_line_integral(power) == expected


@pytest.mark.parametrize("bad", (0.1, sp.Float("0.1"), sp.oo, sp.nan))
def test_inexact_or_nonfinite_serialization_and_fingerprint_rejected(bad):
    with pytest.raises((TypeError, ValueError)):
        verify.serialize({"bad": bad})
    with pytest.raises((TypeError, ValueError)):
        verify.rational_signature(bad)


@pytest.mark.parametrize(
    "name", [key for key in audit.controls() if key != "rejected_inputs"]
)
def test_nonzero_negative_and_scope_control(name):
    assert bool(audit.controls()[name]), name


def test_phase_weight_and_anchor_map_are_exact():
    d = phase.data()
    t = 1 + phase.u**2
    assert d["old_to_weighted_phase"] == sp.diag(
        1, 1, sp.Rational(2, 5) / t**3, sp.Rational(2, 5) / t**3
    )
    assert d["natural_anchor_initial_infinity_norm"] == 1
    assert d["integrable_entry_sum_upper"] < 59


def test_only_negative_damping_is_removed_and_then_retained_for_momenta():
    d = phase.data()
    u = phase.u
    t = 1 + u * u
    assert (
        d["negative_damping_removed_generator"]
        - d["weighted_generator"]
        - sp.diag(0, 0, 6 * u / t, 6 * u / t)
    ).applyfunc(sp.factor) == sp.zeros(4)
    assert all(1 + C < 40 for C in d["momentum_remainder_row_bounds_per_one_over_t"])


def test_elementary_coefficientwise_envelope_is_not_sampled():
    u, k = phase.u, phase.k
    t = 1 + u * u
    d = phase.envelope((2 + 3 * u * u + k * k) / t**2)
    assert d["constant"] == 6 and d["parity"] == 0 and d["power"] == 1
    assert d["integral_upper"] == sp.Rational(66, 7)
    assert d["coefficientwise_remainder_is_nonnegative"]


def test_extra_momentum_half_order_is_retained_in_integral():
    u = phase.u
    t = 1 + u * u
    assert phase.envelope(1 / t, 1)["integral_upper"] == 1
    assert phase.envelope(u / t, 1)["integral_upper"] == 1
    assert phase.envelope(1 / t, 2)["integral_upper"] == sp.Rational(11, 14)
    with pytest.raises(ValueError):
        phase.envelope(u / t)


def test_pointwise_envelopes_need_not_be_integrable():
    u = phase.u
    t = 1 + u * u
    d = phase.envelope(u / t, require_integrable=False)
    assert d["integral_upper"] is None
    assert phase.uniform_upper(d) == sp.Rational(1, 2)
    assert phase.uniform_upper(phase.envelope(1, require_integrable=False)) == 1
    with pytest.raises(TypeError):
        phase.envelope(1 / t, require_integrable=1)


def test_zero_envelope_keeps_same_explicit_proof_structure():
    for integrable in (True, False):
        d = phase.envelope(0, 2, require_integrable=integrable)
        assert d["constant"] == 0
        assert d["coefficientwise_remainder_is_nonnegative"]
        assert d["integral_upper"] == (0 if integrable else None)


def test_all_actual_compensated_kernels_and_off_diagonal_entries_retained():
    d = source.envelopes()
    assert len(d) == 4
    for rows in d.values():
        assert len(rows["entries"]) == 16
        assert (
            rows["half_entry_L1_upper"]
            == sum(row["integral_upper"] for row in rows["entries"]) / 2
        )


def test_exact_tail_cancellation_does_not_assert_physical_divergence():
    d = identities.asymptotic_controls()
    assert d["actual_limits"] == {
        "uncompensated_scale_diagonal_times_u": 396,
        "uncompensated_weighted_trace_diagonal_over_u_four": -79200,
        "compensated_scale_diagonal_times_u": 0,
        "compensated_weighted_trace_diagonal_over_u_four": 0,
    }
    assert d["no_physical_divergence_inferred_from_uncompensated_forcing"]


def test_original_physical_trace_restored_with_weaker_tail_bound():
    d = bounds.data()
    u = phase.u
    t = 1 + u * u
    Q = sp.Integer(3) ** 59
    assert (
        d["global_original_trace_bound_per_eta"]
        == d["weighted_mean_phase_upper_per_eta"] / t**3 + 36 * sp.Abs(u) * Q**2 / t
    )
    assert source.data()["compensating_trace_coefficient"] == -36 * u / t


def test_intrinsic_and_mean_observable_bound_is_absolute_not_fractional():
    d = bounds.data()
    for row in d["global_physical_observable_bounds"].values():
        assert len(row["scaled_kernel_entry_envelopes"]) == 16
        assert row["absolute_coefficient_per_t_minus_four"] < 10**64
    assert audit.controls()[
        "no_uniform_density_fraction_claim_from_absolute_decay_bound"
    ]


def test_joint_global_domain_and_bounds_retain_all_three_blocks():
    d = joint.data()
    assert d["joint_maximum_lapse_bound"] < sp.Rational(1, 10**10)
    assert d["joint_maximum_exact_frame_log_scale_bound"] < sp.Rational(1, 10**10)
    for key in d["global_response_upper"]:
        assert d["global_response_upper"][key].free_symbols == {
            joint.eta_s,
            joint.eta_t,
            joint.eta_p,
        }
    assert all(v == 0 for v in joint.calibrated(0, 0, 0).values())


def test_joint_amplitudes_are_validated_independently():
    bad = sp.Rational(1, 10**19)
    with pytest.raises(ValueError):
        joint.calibrated(0, bad, 0)
    with pytest.raises(ValueError):
        joint.calibrated(0, 0, bad)
    with pytest.raises(ValueError):
        joint.calibrated(sp.Rational(1, 10**79), 0, 0)


def test_joint_exact_frame_center_bounce_uses_total_lapse():
    d = joint.data()
    value = d["center_exact_log_scale_second_derivative_lower"]
    assert value.subs(d["exact_amplitude_box"]) > 3
    assert joint.eta_s in value.free_symbols and joint.eta_t in value.free_symbols


def test_zero_scalar_limit_recovers_actual_tensor_and_proca_coefficients():
    value = joint.data()["global_response_upper"]["lapse"]
    assert value.subs(joint.eta_s, 0) == 1500 * joint.eta_t + 550000000 * joint.eta_p


def test_complete_kernel_fingerprints_rebuilt_without_rounding():
    for matrix in source.data()["actual_weighted_response_kernels"].values():
        result = verify.matrix_signature(matrix)
        assert result["rows"] == 4 and result["columns"] == 4
        assert len(result["row_major_entries"]) == 16
    assert verify.rational_signature(phase.u + phase.k) != verify.rational_signature(
        phase.u - phase.k
    )


def test_scientific_payloads_are_exactly_serializable():
    for d in (
        phase.data(),
        source.data(),
        source.envelopes(),
        bounds.data(),
        identities.asymptotic_controls(),
        joint.data(),
        calibration.calibrated(),
        audit.controls(),
    ):
        verify.serialize(verify.payload(d))


def test_global_relative_result_does_not_claim_original_P8_closure():
    c = audit.controls()
    assert c["no_absolute_or_finite_amplitude_SEE_solution_claim"]
    assert c["no_corrected_cone_or_higher_loop_transfer_from_geometric_completeness"]
    assert c["no_original_P8_or_V_G_B_closure_inferred"]
