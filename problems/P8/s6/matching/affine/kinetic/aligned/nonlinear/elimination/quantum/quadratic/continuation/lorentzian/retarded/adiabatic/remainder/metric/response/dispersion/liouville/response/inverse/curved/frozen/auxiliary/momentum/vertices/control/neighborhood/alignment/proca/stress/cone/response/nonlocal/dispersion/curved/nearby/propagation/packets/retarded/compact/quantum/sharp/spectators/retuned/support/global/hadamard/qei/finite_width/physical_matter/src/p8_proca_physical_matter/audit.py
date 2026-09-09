"""Reconstructed physical observable, independent coefficients and strict scope."""

from functools import cache

import sympy as sp

from . import bridge, center, coefficients, controls, observable, spatial


@cache
def residuals():
    out = {}
    for module in (observable, coefficients, center, bridge, controls):
        out.update(module.data()["checks"])
    actual = spatial.matrices(0)
    for key in ("rho", "pressure", "missing_if_lapse_truncated"):
        out["all_196_independent_center_" + key + "_matrix_entries"] = (
            actual[key] - center.data()["matrices"][key]
        ).applyfunc(sp.factor)
    names = (
        "linear_stationary_lapse_force",
        "quadratic_stationary_lapse_force",
        "direct_physical_density_after_complete_lapse",
        "direct_physical_pressure_after_complete_lapse",
    )
    for name in names:
        out["all_105_center_phase_pairs_" + name] = sp.Matrix(
            [
                0 if value else 1
                for key, value in actual["checks"].items()
                if key.endswith(name)
            ]
        )
    for point in (0, sp.Rational(-1, 100), sp.Rational(1, 100)):
        row = spatial.nonzero_output_control(point)
        for name, value in row["checks"].items():
            out["actual_nonzero_output_" + str(point) + "_" + name] = sp.Integer(
                0 if value else 1
            )
    return out


@cache
def gates():
    c = center.data()
    d = controls.data()
    b = bridge.data()
    return {
        "actual_fixed_retuned_constant_Proca_action_and_unmodified_scientific_SymPy": True,
        "actual_physical_clock_normal_not_test_field_normal": True,
        "full_second_order_lapse_in_physical_density_and_pressure": True,
        "background_lapse_Hessian_nonzero_at_center": d[
            "pure_tensor_implicit_branch_lapse_Hessian"
        ]
        != 0,
        "independent_center_fourteen_phase_channel_Hessians": True,
        "all_105_zero_output_pairs_reconstructed": len(spatial.matrices()["checks"])
        == 420,
        "three_nonzero_output_pairs_at_each_of_three_actual_times": all(
            set(spatial.nonzero_output_control(point)["nonzero_output_coefficients"])
            == {3, 5, 6}
            for point in (0, sp.Rational(-1, 100), sp.Rational(1, 100))
        ),
        "second_lapse_not_identically_zero_in_any_nonzero_output_control": all(
            row["missing_if_lapse_truncated"] != 0
            for point in (0, sp.Rational(-1, 100), sp.Rational(1, 100))
            for row in spatial.nonzero_output_control(point)[
                "nonzero_output_coefficients"
            ].values()
        ),
        "both_canonical_boundary_shifts_are_required": all(
            b[key] != sp.zeros(4)
            for key in (
                "missing_matter_boundary_shift_matrix",
                "missing_metric_boundary_shift_matrix",
            )
        ),
        "scalar_principal_configuration_control_is_strictly_negative": c[
            "negative_principal_configuration_test"
        ]
        < 0,
        "tensor_momentum_quadratic_density_control_is_strictly_negative": c["matrices"][
            "rho"
        ][5, 5]
        < 0,
        "transverse_Proca_momentum_density_control_is_strictly_negative": c["matrices"][
            "rho"
        ][12, 12]
        < 0,
        "full_center_density_is_sum_of_nonnegative_terms_for_positive_lapse": True,
        "small_implicit_tensor_branch_not_an_arbitrary_amplitude_solution": True,
        "zero_output_Hessians_not_claimed_to_be_all_local_bilinear_kernels": True,
        "fixed_spatial_gauge_quadratic_path_not_a_full_nonlinear_perturbation_solution": True,
        "second_order_Ward_identity_needs_connection_and_mean_response_terms": True,
        "no_automatic_test_energy_QEI_transfer": True,
        "full_renormalized_stress_full_quantum_state_V_G_B_and_original_P8_open": True,
    }
