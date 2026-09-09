"""Native constrained dynamics, physical covariance and global-state controls."""

from fractions import Fraction

import pytest
import sympy as sp
from p8_proca_seven_modes import (
    audit,
    bridge,
    covariance,
    deformation,
    model,
    projectors,
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
def test_written_and_exact_proof_gate(name):
    assert bool(audit.gates()[name]), name


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[row[0] for row in audit.bad_cases()]
)
def test_invalid_preparation_strip_inputs(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_manifest_and_exact_counts():
    assert len(ROWS) == 69
    assert (
        sum(
            v.rows * v.cols if isinstance(v, sp.MatrixBase) else 1
            for v in ROWS.values()
        )
        == 1360
    )
    assert len(audit.gates()) == 22
    assert audit.rejected_inputs() == 16
    assert len(verify.source_files()) == 21


@pytest.mark.parametrize("radius", (3, 4, Fraction(7, 2), sp.Rational(13, 3), 100))
def test_positive_all_frequency_compact_strip_bounds(radius):
    d = deformation.strip(radius)
    assert d["Proca_energy_lower"] > 0
    assert d["tensor_modified_energy_lower"] > 0
    assert d["actual_and_auxiliary_scale_upper"] > 1
    assert d["tensor_low_frequency_radial_integral_majorant"] == sp.Rational(7, 6)


@pytest.mark.parametrize("bad", (0.1, sp.Float("0.1"), sp.oo, sp.nan))
def test_inexact_or_nonfinite_report_values_rejected(bad):
    with pytest.raises((TypeError, ValueError)):
        verify.serialize({"bad": bad})


def test_full_state_retains_seven_modes_and_fourteen_CCR_channels():
    d = model.data()
    assert d["actual_full_fourteen_phase_generator"].shape == (14, 14)
    assert d["actual_full_symplectic_form"].rank() == 14
    assert state.data()["full_physical_mode_count"] == 7
    assert state.data()["additional_physical_tensor_and_vector_modes"] == 5


def test_massive_polarization_is_three_positive_columns_not_four_ghosts():
    d = covariance.data()
    assert d["positive_Proca_Gram_factor_without_sqrt_hbar_over_kappa"].shape == (6, 3)
    assert d["physical_Proca_polarization_factor"].shape == (4, 3)
    assert (
        audit.controls()[
            "naive_four_independent_covariant_polarizations_have_negative_time_component"
        ]
        == -1
    )


def test_nonzero_Proca_symbol_does_not_imply_ellipticity():
    d = covariance.data()["nonelliptic_nonzero_null_Proca_projection_symbol"]
    assert d != sp.zeros(4)
    assert d.rank() == 1 and d.det() == 0 and d * d == sp.zeros(4)


def test_no_angular_polarization_frame_or_fourth_vector_mode():
    d = projectors.data()
    assert sp.factor(sp.trace(d["TT"])) == 2
    assert sp.factor(sp.trace(d["longitudinal"])) == 1
    assert sp.factor(sp.trace(d["transverse"])) == 2
    assert audit.controls()["no_projector_direction_or_zero_frequency_atom_is_inserted"]


def test_reference_preparation_does_not_replace_actual_background():
    d = deformation.data()
    assert d["reference_preparation_initial_time"] == -3
    assert d["actual_anchor_time"] == 0
    assert d["auxiliary_metric_equals_actual_for_u_ge_minus_one"]
    assert state.data()["same_original_global_background_on_both_legs_of_final_state"]
    assert audit.controls()["no_instantaneous_Minkowski_reset_at_actual_bounce"]


def test_old_scalar_marginal_and_subtraction_boundary():
    d = state.data()
    assert d["scalar_restriction_is_exactly_the_preexisting_chosen_scalar_state"]
    assert d["all_covariances_retain_common_hbar_over_kappa"]
    assert audit.controls()[
        "no_old_tadpole_or_stress_subtraction_reassigned_to_new_product_state"
    ]


def test_physical_stress_bridge_keeps_both_scalar_boundary_shifts():
    d = bridge.data()["full_density_to_S6_102_spatial_phase"]
    assert bridge.data()["physical_stress_bridge_time"] == 0
    assert d.shape == (14, 14)
    assert d[1, 1] != 0 and d[3, 0] != 0
    assert d[4:, 4:] == sp.eye(10)


def test_tensor_realization_keeps_actual_factor_two():
    S = state.data()["scalar_to_each_tensor_polarization_canonical_map"]
    assert S == sp.diag(sp.sqrt(2), 1 / sp.sqrt(2))
    assert S * model.J2 * S.T == model.J2


def test_complete_preparation_energy_bounds_not_high_frequency_only():
    d = deformation.data()
    assert d["transition_absolute_scale_derivative_upper"] == 232
    assert d["Proca_energy_exponent_from_initial_to_anchor_upper"] == 702
    assert d["tensor_modified_energy_exponent_from_initial_to_anchor_upper"] == 705
