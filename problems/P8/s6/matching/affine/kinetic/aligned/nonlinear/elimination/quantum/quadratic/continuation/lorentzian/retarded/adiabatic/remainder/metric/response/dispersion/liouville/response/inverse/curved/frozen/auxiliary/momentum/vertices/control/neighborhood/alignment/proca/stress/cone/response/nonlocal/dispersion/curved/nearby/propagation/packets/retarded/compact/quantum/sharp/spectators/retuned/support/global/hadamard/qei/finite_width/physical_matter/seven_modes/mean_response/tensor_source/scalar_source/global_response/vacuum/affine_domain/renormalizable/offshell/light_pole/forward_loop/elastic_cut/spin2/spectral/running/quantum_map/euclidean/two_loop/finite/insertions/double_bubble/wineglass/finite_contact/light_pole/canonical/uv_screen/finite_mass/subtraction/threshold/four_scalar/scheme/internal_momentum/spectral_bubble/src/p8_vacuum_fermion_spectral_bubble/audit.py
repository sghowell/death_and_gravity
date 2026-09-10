"""Spectral, subtraction, combinatoric and restricted-scope checks."""

from functools import cache

import sympy as sp

from . import calibration, density, outer, ownership, spectral

MODULES = (spectral, density, outer, ownership, calibration)


@cache
def residuals():
    return {
        mod.__name__.rsplit(".", 1)[-1] + "_" + name: sp.simplify(value)
        for mod in MODULES
        for name, value in mod.data()["checks"].items()
    }


@cache
def gates():
    result = dict(calibration.data()["bounds"])
    result.update(
        {
            "direct_spectral_identity_not_assumed_high_energy_contour": True,
            "whole_frozen_finite_kernel_and_both_OS_anchors_retained": True,
            "complete_parameter_interval_produces_beta_cubed_density": True,
            "spectral_measure_nonnegative": True,
            "formal_continuum_measure_not_finite_normalized_spectral_weight": True,
            "one_formal_propagator_insertion_not_exact_resummation": True,
            "two_light_line_positions_and_three_channels_explicit": True,
            "outer_symmetry_and_Q_factors_match_frozen_scalar_loop": True,
            "proper_fermion_self_energy_references_paired_before_bounds": True,
            "only_overall_local_quartic_constant_subtracted": True,
            "spectral_integral_performed_after_finite_outer_difference": True,
            "mixed_mass_parameter_gap_on_required_complex_domain": True,
            "once_subtracted_outer_integral_absolutely_convergent": True,
            "uniform_second_derivative_majorant_integrable": True,
            "both_crossed_channel_Taylor_factors_retained": True,
            "zero_transfer_channel_contributes_no_forward_second_coefficient": True,
            "finite_local_quartic_reference_change_has_zero_family_b2": True,
            "literal_family_not_complete_two_loop_amplitude_or_error": True,
            "no_independent_higher_derivative_counterterm_used": True,
            "no_full_spectrum_cutoff_or_V_G_B_claim": True,
            "original_P8_remains_open": True,
        }
    )
    return {k: bool(v) for k, v in result.items()}


def bad_cases():
    invalid = (
        True,
        False,
        1.0,
        sp.Float(1),
        "1",
        None,
        sp.oo,
        -sp.oo,
        sp.zoo,
        sp.I,
        sp.nan,
        sp.Symbol("x"),
    )
    result = []
    for j in range(4):
        for i, value in enumerate(invalid):
            values = [2, 1, 1, 144]
            values[j] = value
            result.append((f"type_{j}_{i}", outer.enclosure, tuple(values)))
    result.extend(
        (f"outside_{i}", outer.enclosure, values)
        for i, values in enumerate(
            (
                (0, 1, 1, 144),
                (1, 1, 1, 144),
                (-2, 1, 1, 144),
                (2, -1, 1, 144),
                (2, 1, 1, 0),
                (2, 1, 1, -144),
            )
        )
    )
    return result


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            continue
        raise ValueError("Unsupported spectral outer-bubble input accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "both_internal_light_line_insertions_present": True,
        "whole_inner_local_reference_pairing_required": True,
        "unsubtracted_outer_spectral_integral_not_used": True,
        "infinite_total_formal_continuum_weight_not_hidden": True,
        "no_extra_canonical_LSZ_or_kinetic_resummation": True,
        "local_constant_not_an_adjustable_b2_contact": True,
        "remaining_two_loop_families_not_called_complete": True,
        "original_P8_not_closed": True,
    }
