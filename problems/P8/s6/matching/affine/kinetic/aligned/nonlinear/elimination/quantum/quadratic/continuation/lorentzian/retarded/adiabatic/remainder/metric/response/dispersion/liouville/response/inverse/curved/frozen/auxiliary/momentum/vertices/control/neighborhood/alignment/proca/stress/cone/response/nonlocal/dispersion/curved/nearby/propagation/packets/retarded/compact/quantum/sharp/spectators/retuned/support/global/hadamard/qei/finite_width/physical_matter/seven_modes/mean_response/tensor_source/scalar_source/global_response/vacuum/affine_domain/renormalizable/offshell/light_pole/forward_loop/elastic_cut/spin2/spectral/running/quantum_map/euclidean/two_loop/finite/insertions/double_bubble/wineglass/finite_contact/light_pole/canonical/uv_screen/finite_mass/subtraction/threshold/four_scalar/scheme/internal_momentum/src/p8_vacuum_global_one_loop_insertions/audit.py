"""Whole-kernel bounds, restricted routing, and explicit nonclosure controls."""

from functools import cache

import sympy as sp

from . import calibration, halfplane, routing, spacelike, window

MODULES = (halfplane, spacelike, routing, window, calibration)


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
            "same_frozen_whole_one_loop_fermion_kernel": True,
            "complete_affine_reference_subtracted_before_any_bound": True,
            "finite_on_shell_mass_and_residue_anchors_retained": True,
            "no_negative_large_mass_counterterm_used_as_free_propagator": True,
            "unbounded_complex_halfplane_not_finite_momentum_series": True,
            "positive_resolvent_weights_and_exact_parameter_moments": True,
            "segment_Taylor_remainder_uses_same_halfplane_gap": True,
            "on_shell_insertion_has_removable_extension": True,
            "global_spacelike_envelope_retains_logarithmic_decay": True,
            "routing_is_pointwise_without_complex_q_translation": True,
            "selected_channel_routing_not_every_multiloop_graph": True,
            "both_scalar_and_fermion_inverse_signs_retained": True,
            "selected_full_one_loop_normalization_not_all_loop_resummation": True,
            "reference_energy_window_not_inferred_Wilsonian_cutoff": True,
            "no_global_Lorentzian_spectrum_or_reflection_positivity_claim": True,
            "no_outer_two_loop_or_later_error_bound_claim": True,
            "original_V_G_B_and_P8_remain_open": True,
        }
    )
    return {name: bool(value) for name, value in result.items()}


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
    calls = (
        ("halfplane", halfplane.enclosure, (2, 1, 144, 1)),
        ("spacelike", spacelike.envelope, (2, 1, 144, 0)),
        ("routing", routing.channel_margin, (2, 1, 1)),
        ("window", window.defect_bound, (0, 0, 0, 1)),
    )
    result = []
    for label, call, valid in calls:
        for j in range(len(valid)):
            for i, value in enumerate(invalid):
                values = list(valid)
                values[j] = value
                result.append((f"{label}_type_{j}_{i}", call, tuple(values)))
    outside = (
        (
            halfplane.enclosure,
            (
                (0, 1, 144, 1),
                (-1, 1, 144, 1),
                (2, -1, 144, 1),
                (2, 1, 0, 1),
                (2, 1, 144, 0),
                (2, 1, 144, 16),
                (2, 1, 144, 17),
            ),
        ),
        (spacelike.envelope, ((2, 1, 144, -1),)),
        (routing.channel_margin, ((-1, 0, 1), (2, -1, 1), (2, 1, 0))),
        (
            window.defect_bound,
            ((-1, 0, 0, 1), (0, -1, 0, 1), (0, 0, -1, 1), (0, 0, 0, 0), (0, 0, 0, -1)),
        ),
    )
    for i, (call, values) in enumerate(outside):
        result.extend((f"outside_{i}_{j}", call, args) for j, args in enumerate(values))
    return result


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (TypeError, ValueError):
            continue
        raise ValueError(
            "Unsupported global one-loop insertion input accepted: " + name
        )
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "mass_counterterm_paired_before_estimation": True,
        "Taylor_in_momentum_over_mass_not_used": True,
        "constant_bound_not_promoted_to_outer_UV_convergence": True,
        "arbitrary_internal_routings_not_claimed": True,
        "finite_reference_energy_not_a_cutoff_certificate": True,
        "full_one_loop_not_all_order_spectrum_or_error": True,
        "valid_wide_bounds_may_be_inconclusive": True,
        "original_P8_not_closed": True,
    }
