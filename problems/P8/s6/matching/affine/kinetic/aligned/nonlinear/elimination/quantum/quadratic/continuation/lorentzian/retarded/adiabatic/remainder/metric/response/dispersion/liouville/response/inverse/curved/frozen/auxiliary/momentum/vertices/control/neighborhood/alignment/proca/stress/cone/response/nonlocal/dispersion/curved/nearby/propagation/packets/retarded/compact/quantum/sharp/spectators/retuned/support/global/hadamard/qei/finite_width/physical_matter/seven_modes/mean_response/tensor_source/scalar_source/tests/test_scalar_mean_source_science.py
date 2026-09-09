"""Native scalar mean variations, physical observable and whole-strip checks."""

from fractions import Fraction

import pytest
import sympy as sp
from p8_scalar_mean_source import (
    audit,
    bounds,
    calibration,
    geometry,
    model,
    observable,
    response,
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
def test_invalid_amplitude_or_rational_bound_structure(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_source_and_scientific_manifest():
    assert len(ROWS) == 33
    assert (
        sum(
            v.rows * v.cols if isinstance(v, sp.MatrixBase) else 1
            for v in ROWS.values()
        )
        == 166
    )
    assert len(audit.gates()) == 43
    assert audit.rejected_inputs() == 32
    assert len(verify.source_files()) == 21


@pytest.mark.parametrize("value", (0, Fraction(1, 10**20), sp.Rational(1, 10**21)))
def test_exact_valid_calibration_and_zero_response(value):
    d = calibration.calibrated(value)
    assert d["eta"] == sp.Rational(value)
    assert all(v >= 0 for v in d.values())
    assert all(not key.endswith("_per_eta") for key in d)
    if value == 0:
        assert all(v == 0 for v in d.values())


@pytest.mark.parametrize("bad", (0.1, sp.Float("0.1"), sp.oo, sp.nan))
def test_inexact_or_nonfinite_serialization_and_kernel_fingerprint_rejected(bad):
    with pytest.raises((TypeError, ValueError)):
        verify.serialize({"bad": bad})
    with pytest.raises((TypeError, ValueError)):
        verify.rational_signature(bad)


@pytest.mark.parametrize(
    "name", [key for key in audit.controls() if key != "rejected_inputs"]
)
def test_nonzero_negative_and_scope_control(name):
    assert bool(audit.controls()[name]), name


def test_existing_generator_requires_both_shifts_and_time_boundary():
    B = geometry.canonical()["old_to_natural_phase"]
    assert B[2, 1] != 0 and B[3, 0] != 0 and B[2, 0] != 0
    assert model.bridge()["generating_time_derivative"] != 0
    assert model.bridge()["difference"] == 0


def test_all_four_actual_source_matrices_are_nonzero():
    for value in source.data()["source_Hessians"].values():
        assert value != sp.zeros(4)


def test_actual_clock_pivot_not_frozen_before_mean_variation():
    g = model.generic()
    assert g["actual_mean_lapse_source"] != 0
    assert any(
        g["symbols"][name][3] in g["actual_mean_lapse_source"].free_symbols
        for name in "abcd"
    )


def test_scalar_time_reversal_has_odd_matter_coordinate():
    assert response.parity()["natural_phase_reversal"] == sp.diag(1, -1, -1, 1)


def test_zero_anchor_means_do_not_force_zero_lapse_or_matter_velocity():
    c = response.center()
    assert c["anchor_induced_lapse"] != 0
    assert c["initial_source_traces"]["mean_matter_field_direct_source"] != 0
    assert c["initial_source_traces"]["mean_hat_scale_direct_source"] == 0


def test_intrinsic_physical_density_is_not_discarded_or_lapse_counted_twice():
    d = observable.data()
    assert d["intrinsic_density_Hessian"] != sp.zeros(4)
    assert d["mean_lapse_not_counted_again_inside_intrinsic_quadratic_observable"]
    assert (
        observable.center()["checks"][
            "all_center_physical_density_pairs_match_independent_nonlinear_observable"
        ]
        == 0
    )


def test_rational_strip_bound_is_exact_and_coefficient_based():
    u, k = model.u, model.k
    assert bounds.exact_rational_bound(
        (2 + 3 * u * u + k * k) / (1 + 5 * u * u)
    ) == 6 + sp.Rational(3, 10000)
    assert bounds.exact_rational_bound(-7) == 7
    assert bounds.exact_rational_bound(0) == 0


def test_whole_strip_sources_and_observables_keep_the_selected_band():
    d = bounds.response_bounds()
    assert d["exact_generator_infinity_norm_upper"] == 90
    assert d["added_covariance_each_entry_absolute_upper_per_eta"] == 9
    assert d["mean_scale_and_trace_upper_per_eta"] == 5
    assert d["mean_lapse_upper_per_eta"] == 120
    assert d["complete_physical_density_and_pressure_absolute_upper_per_eta"] == 40


def test_center_acceleration_is_not_given_an_unproved_uniform_sign():
    c = response.center()["anchor_proper_Hubble_derivative_response"]
    assert c.subs(model.k, 1) > 0
    assert c.subs(model.k, 2) < 0
    assert (
        bounds.response_bounds()["exact_center_coefficient_bounds"][
            "anchor_proper_Hubble_derivative_response"
        ]
        < 200
    )


def test_kernel_fingerprint_preserves_all_exact_coefficients():
    u, k = model.u, model.k
    assert verify.rational_signature(
        (2 * u + 2) / (2 * u + 2)
    ) == verify.rational_signature(1)
    assert verify.rational_signature(u + k) != verify.rational_signature(u - k)
    matrix = verify.matrix_signature(sp.Matrix([[0, u], [k, 1]]))
    assert matrix["rows"] == 2 and matrix["columns"] == 2
    assert len(matrix["row_major_entries"]) == 4


def test_all_actual_source_and_observable_kernels_have_complete_exact_fingerprints():
    for matrix in list(source.data()["source_Hessians"].values()) + [
        source.data()["natural_phase_generator"],
        observable.data()["intrinsic_density_Hessian"],
        observable.data()["intrinsic_pressure_Hessian"],
    ]:
        encoded = verify.matrix_signature(matrix)
        assert len(encoded["row_major_entries"]) == 16
        for entry in encoded["row_major_entries"]:
            for side in ("numerator", "denominator"):
                assert len(entry[side]["sorted_exact_rational_monomial_sha256"]) == 64


def test_scientific_payloads_are_exactly_serializable():
    for d in (
        model.generic(),
        model.coefficients(),
        geometry.data(),
        geometry.canonical(),
        source.density_center(),
        response.data(),
        response.center(),
        observable.center(),
        bounds.response_bounds(),
        calibration.calibrated(),
        audit.controls(),
    ):
        verify.serialize(verify.payload(d))


def test_local_representative_not_promoted_to_global_or_exact_quantum_bounce():
    c = audit.controls()
    assert c["no_uniform_infinite_tail_bound_inferred_from_local_interval"]
    assert c[
        "no_positive_lapse_for_every_unbounded_Gaussian_field_configuration_claimed"
    ]
    assert c["no_finite_amplitude_SEE_or_full_quantum_Ward_identity_claimed"]
