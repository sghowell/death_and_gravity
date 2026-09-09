"""Exact amplitude controls for the new tensor state difference only."""

from functools import cache

from p8_proca_mean_response import band as inherited

from . import bounds


def calibrated(value=inherited.DEFAULT_AMPLITUDE):
    eta = inherited.exact_amplitude(value)
    d = bounds.data()
    return {
        "eta": eta,
        "global_proxy_density_upper": 4 * eta,
        "global_physical_tensor_density_absolute_upper": 12 * eta,
        "global_tensor_pressure_absolute_upper": 4 * eta,
        "global_weighted_mean_phase_upper": d[
            "global_weighted_mean_phase_upper_per_eta"
        ]
        * eta,
        "global_lapse_upper": d["global_lapse_upper_per_eta"] * eta,
        "global_linear_physical_log_scale_upper": d[
            "global_linear_physical_log_scale_upper_per_eta"
        ]
        * eta,
        "global_exact_frame_log_scale_upper": d[
            "global_actual_frame_log_scale_upper_per_eta"
        ]
        * eta,
        "global_matter_field_upper": d["global_matter_field_upper_per_eta"] * eta,
        "global_matter_density_absolute_upper": d["global_scalar_density_upper_per_eta"]
        * eta,
        "global_relative_matter_density_upper": d[
            "global_relative_scalar_density_upper_per_eta"
        ]
        * eta,
        "absolute_weighted_clock_integral_upper": d[
            "rounded_absolute_volume_weighted_clock_source_integral_upper_per_eta"
        ]
        * eta,
        "center_proper_acceleration_absolute_upper": 15 * eta,
    }


@cache
def bad_cases():
    return [(name, calibrated, args) for name, _, args in inherited.bad_cases()]
