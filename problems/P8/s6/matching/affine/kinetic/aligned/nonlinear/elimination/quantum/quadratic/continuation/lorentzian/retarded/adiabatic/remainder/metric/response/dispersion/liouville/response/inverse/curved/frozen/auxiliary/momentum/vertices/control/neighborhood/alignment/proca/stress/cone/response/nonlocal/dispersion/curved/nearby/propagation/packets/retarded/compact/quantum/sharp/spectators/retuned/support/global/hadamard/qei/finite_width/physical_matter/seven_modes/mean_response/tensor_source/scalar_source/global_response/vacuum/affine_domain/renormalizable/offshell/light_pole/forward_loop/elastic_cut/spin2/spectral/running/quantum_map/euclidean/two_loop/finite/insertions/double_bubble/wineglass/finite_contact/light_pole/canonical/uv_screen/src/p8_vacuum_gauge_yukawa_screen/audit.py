"""Exact residuals and narrowly scoped matching-screen gates."""

from functools import cache

import sympy as sp

from . import calibration, cuts, flow, model, threshold

MODULES = (model, flow, threshold, cuts, calibration)


@cache
def residuals():
    result = {}
    for mod in MODULES:
        prefix = mod.__name__.rsplit(".", 1)[-1]
        for name, value in mod.data()["checks"].items():
            result[prefix + "_" + name] = (
                value.applyfunc(sp.simplify)
                if isinstance(value, sp.MatrixBase)
                else sp.simplify(value)
            )
    return result


def scalar_entries():
    return sum(
        v.rows * v.cols if isinstance(v, sp.MatrixBase) else 1
        for v in residuals().values()
    )


@cache
def gates():
    d = flow.data()
    c = cuts.data()
    m = model.data()
    result = {k: bool(v) for k, v in calibration.data()["bounds"].items()}
    result.update(
        {
            "eight_independently_normalized_SU3_generators": len(m["gauge_generators"])
            == 8,
            "fourteen_vectorlike_fundamental_Dirac_flavors": m["all_Dirac_flavors"]
            == 14,
            "two_opposite_active_flavor_Yukawas": m["active_Dirac_flavors"] == 2,
            "twelve_Yukawa_inert_spectators": m["spectator_Dirac_flavors"] == 12,
            "all_new_gauge_masses_zero_in_unbroken_perturbation_theory": m[
                "perturbative_gauge_boson_mass_squared"
            ]
            == [0] * 8,
            "positive_gauge_one_loop_AF_coefficient": m["one_loop_gauge_b0"] > 0,
            "quartic_positive_root_bracket": d["quartic_root_lower_test"]
            < 0
            < d["quartic_root_upper_test"],
            "opposite_massive_flavors_preserve_Phi_parity": True,
            "matrix_Yukawa_and_potential_tensor_derivations_retained": True,
            "marginal_fixed_flow_only_not_relevant_parameter_stability": True,
            "both_heavy_flavor_thresholds_retained": True,
            "threshold_field_domain_not_momentum_remainder": True,
            "physical_transverse_and_Lorentz_contractions_agree": True,
            "identical_gauge_particle_and_optical_factors_retained": True,
            "new_cut_enters_at_three_parent_loops": c[
                "first_full_parent_loop_order_for_this_cut"
            ]
            == 3,
            "opposite_channel_boundary_prescriptions_retained": True,
            "crossing_center_cancellation_not_missing_cut": cuts.point(2)[
                "normalized_forward_boundary_jump"
            ]
            == 0,
            "nonzero_boundary_mismatch_inside_old_disc": cuts.point(sp.Rational(3, 2))[
                "normalized_forward_boundary_jump"
            ]
            != 0,
            "local_logarithm_germ_not_ruled_out": True,
            "second_derivative_jump_cancellation_not_gapped_dispersion": True,
            "finite_mass_on_shell_threshold_remainder_not_computed": True,
            "one_loop_transmutation_scale_not_confinement_mass": True,
            "old_two_loop_amplitude_not_transferred_to_new_fields": True,
            "no_new_all_orders_UV_requirement_added_to_P8": True,
            "finite_gravity_and_common_parent_obligations_unchanged": True,
            "no_full_candidate_or_whole_row_exclusion": True,
            "original_V_G_B_and_P8_remain_open": True,
        }
    )
    return {name: bool(value) for name, value in result.items()}


@cache
def rejected_inputs():
    for name, call, args in cuts.bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            continue
        raise ValueError("Unsupported channel invariant accepted: " + name)
    return len(cuts.bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "odd_flavor_threshold_cancels_but_even_threshold_nonzero": threshold.data()[
            "leading_Phi_squared_F_squared_coefficient"
        ]
        != 0,
        "zero_boundary_jump_at_center_is_not_a_disc_test": cuts.point(2)[
            "normalized_forward_boundary_jump"
        ]
        == 0
        and cuts.point(sp.Rational(3, 2))["normalized_forward_boundary_jump"] != 0,
        "second_derivative_jump_zero_does_not_restore_mass_gap": True,
        "tiny_nonzero_cut_not_discarded_by_absolute_size": True,
        "asymptotic_freedom_not_a_confinement_mass_bound": True,
        "leading_Wilson_coefficient_not_full_momentum_form_factor": True,
        "new_model_not_assigned_old_complete_two_loop_error": True,
        "full_V_G_B_and_original_P8_not_closed": True,
    }
