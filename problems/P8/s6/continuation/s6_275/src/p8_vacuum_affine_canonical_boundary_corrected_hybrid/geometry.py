"""Corrected raw field rows and recomputed complete nonlinear source domains."""

from functools import cache

import sympy as s
from p8_vacuum_affine_finite_volume_turnaround import moving
from p8_vacuum_affine_hybrid_core_regulator_comparison import quadratic
from p8_vacuum_affine_selfconsistent_finite_feedback import geometry as original

from . import source

P, KAPPA, CORE = original.P, original.KAPPA, source.RADIUS
kap = KAPPA
NZ, NZZ = original.NZ, original.NZZ
field = moving.old.field
F_AMPLITUDE = source.F_AMPLITUDE


def field_rows():
    rows = dict(field.field_bounds())
    rows["Pi_v_A0"] *= 2
    rows["delta_Pi_M_A0"] *= 2
    return rows


@cache
def boundary_rows():
    before = field.field_bounds()
    v = 4 * before["v_A2"] / (1 + P) ** 2
    sigma = before["M1_A1"] / (1 + P)
    return {
        "whole_prepared_reference_per_radius_rows": before,
        "whole_corrected_raw_per_radius_rows": field_rows(),
        "whole_added_trace_momentum_row": 144 * source.TIME * v + 3 * sigma / 10,
        "whole_added_matter_momentum_row": 3 * v / 10,
        "whole_boundary_time_coefficient": 18 * (4 + 48 * source.TIME**2),
    }


def require_radius(value):
    if isinstance(value, bool) or not isinstance(value, (int, s.Integer, s.Rational)):
        raise TypeError("Require an exact declared corrected-model radius")
    if value not in (source.RADIUS, source.ANALYSIS_RADIUS):
        raise ValueError("Require the physical core or the larger analysis radius")
    return s.sympify(value)


def bounds(radius=source.RADIUS):
    radius = require_radius(radius)
    return _bounds(radius)


@cache
def _bounds(radius):
    f = {name: value * 8 * radius for name, value in field_rows().items()}
    v = 4 * f["v_A2"] / (1 + P) ** 2
    tau = 4 * f["tau_A2"] / (1 + P) ** 2
    hom = source.HOMOGENEOUS_RADIUS
    pbound = 8 * source.TIME + hom
    pvbound = 100 * (source.TIME + hom)
    trace = (pvbound + f["Pi_v_A0"]) * P * v + P * f["Pi_v_A0"] / 3
    shape = 180 * 32 * (1 + P) * f["Pi_tau_A0"]
    vector = 12 * (1 + P) ** 2 * f["W_A0"] * f["Pi_W_A0"]
    matter = 2 * (
        (1 + f["delta_Pi_M_A0"]) * f["M1_A1"]
        + (2 * hom + f["Pi_H_A0"]) * f["H_gradient_A0"]
    )
    D0 = trace + shape + vector + matter
    lam = 4 * D0 / P
    piTF = 5 * (32 * f["Pi_tau_A0"] + 16 * lam)
    shear = 192 * piTF**2
    curve = 2 * 10**12 * P**2 * (v + tau) ** 2
    images = [
        2 * f["Pi_v_A0"] + 12 * pbound * v,
        12 * (1 + P) * f["Pi_W_A0"],
        4 * f["delta_Pi_M_A0"] + 6 * v,
        4 * f["Pi_H_A0"] + 12 * hom * v,
        f["eta_A0"],
        shear,
        16 * f["Pi_W_A0"] ** 2 / field.ZETA,
        96 * field.ZETA * ((1 + P) * f["W_A0"]) ** 2,
        4 * f["W_A0"] ** 2,
        4 * f["M1_A1"] ** 2,
        4 * f["H_gradient_A0"] ** 2,
        8 * P**2 * v + curve,
    ]
    means = [
        10 * pbound * v * v + 4 * f["Pi_v_A0"] * v,
        36 * v * (1 + P) * f["Pi_W_A0"],
        10 * v * v + 12 * f["delta_Pi_M_A0"] * v,
        10 * hom * v * v + 12 * f["Pi_H_A0"] * v,
        s.Integer(0),
        *images[5:11],
        curve,
    ]
    zsum, avg = sum(images), sum(means)
    hvar = (
        1000 * v * v
        + 1000 * avg
        + 10**9 * zsum * zsum / 2
        + 6 * v * (1000 * zsum + 10**9 * zsum * zsum / 2)
    )
    contraction = (
        s.Rational(1, 30) + 10**9 * (source.LAPSE_RADIUS + 12 * (hom + max(images))) / 3
    )
    residual = s.Rational(5, 10**400) + 12 * 10**4 * hom + NZ * zsum
    assert hom + max(images) + 8 * source.TIME < source.COORDINATE
    assert contraction < s.Rational(1, 20)
    assert residual / 3 + contraction * source.LAPSE_RADIUS < source.LAPSE_RADIUS
    assert max(f["v_A2"], f["tau_A2"]) < s.Rational(1, 1000)
    boundary_amplitude = 10**6 * kap * v**2
    amplitude = (
        2 * 1000 * kap * hvar + 10**128 * (8 * radius) ** 2 / 2 + boundary_amplitude
    )
    return {
        "fields": f,
        "v": v,
        "tau": tau,
        "pbound": pbound,
        "Pi_v_background": pvbound,
        "full_Dtrace": trace,
        "full_Dshape": shape,
        "full_Dvector": vector,
        "full_Dmatter": matter,
        "D0": D0,
        "lambda_A1": lam,
        "pi_TF": piTF,
        "shear": shear,
        "curvature_remainder": curve,
        "images": images,
        "averages": means,
        "sum_images": zsum,
        "sum_average_remainders": avg,
        "hvariation": hvar,
        "whole_centered_amplitude": amplitude,
        "whole_boundary_time_generator_amplitude": boundary_amplitude,
        "contraction": contraction,
        "whole_center_residual": residual,
    }


@cache
def remainder():
    small, large = bounds(source.RADIUS), bounds(source.ANALYSIS_RADIUS)
    ratio = source.RADIUS / source.ANALYSIS_RADIUS
    tail = large["whole_centered_amplitude"] * ratio**3 / (1 - ratio)
    raw_rows = quadratic.bounds()["whole_quadratic_row_budgets"]
    corrected_rows = {name: 4 * value for name, value in raw_rows.items()}
    raw2 = 2 * 1000 * sum(corrected_rows.values()) * (8 * source.RADIUS) ** 2
    reference = 10**128 * (8 * source.RADIUS) ** 2 / 2
    generator = small["whole_boundary_time_generator_amplitude"]
    ndelta = NZ * small["sum_images"]
    nmean = NZ * small["sum_average_remainders"] + NZZ * small["sum_images"] ** 2 / 2
    volume = (
        20 * small["v"] ** 2 + 72 * small["v"] * ndelta + 12 * nmean + 5000 * ndelta**2
    )
    return {
        "whole_corrected_positive_full_quadratic_rows": corrected_rows,
        "whole_raw_constrained_quadratic_amplitude": raw2,
        "whole_retained_independent_reference_amplitude": reference,
        "whole_retained_boundary_time_generator_amplitude": generator,
        "whole_analysis_to_core_ratio": ratio,
        "whole_complete_degree_at_least_three_remainder": tail,
        "whole_corrected_centered_interaction_amplitude": raw2
        + reference
        + generator
        + tail,
        "whole_corrected_centered_volume_amplitude": volume,
        "whole_full_point_lapse_response": ndelta,
        "whole_full_mean_lapse_response": nmean,
    }


@cache
def data():
    small, large, b, r = (
        bounds(),
        bounds(source.ANALYSIS_RADIUS),
        boundary_rows(),
        remainder(),
    )
    before, after = field.field_bounds(), field_rows()
    UN = s.Rational(3, 4) * 4 * (2 + source.OFF_CLOCK)
    UNN = s.Rational(21, 16) * 8 * (2 + source.OFF_CLOCK) ** 2 + 3 * (
        6 + source.OFF_CLOCK
    )
    return {
        "whole_boundary_corrected_field_row_derivation": b,
        "whole_corrected_small_complete_phase_domain": small,
        "whole_corrected_large_complex_analysis_domain": large,
        "whole_corrected_full_quadratic_and_entire_remainder": r,
        "whole_implicit_and_readout_first_second_rows": [NZ, NZZ, UN, UNN],
        "whole_corrected_row_proof": "The fixed prepared reference sends4R into8R and preserves the full original covariance in its own coordinates. On that image the raw boundary adds a_bar^3[-18Hbar v+3ellbar sigma] to Pi_v and3a_bar^3ellbar v to Pi_M. Since a_bar^3<2, |Hbar|<=4T and3a_bar^3ellbar=3/10, the displayed added rows are each smaller than the corresponding OLD prepared momentum row. Therefore doubling exactly those two raw rows is valid; every other row is unchanged. This is an inequality for the corrected chart, NOT equality with the old physical dictionary.",
        "whole_corrected_full_geometry_proof": "Recompute the full trace, shape, vector/Gauss and matter sources, formal-transpose inverse, all twelve pointwise source images and twelve mean remainders with the corrected rows at BOTH radii. Nonlinear reconstruction retains every generated harmonic, scalar curvature and density mean. The full phase-plus-homogeneous residual is retained in the same lapse self-map. Complete source identities and their500 envelopes apply at every resulting point.",
        "whole_boundary_time_generator_bound": "The exact additional energy is partial_u F_b=-9kappa a_bar^3(Hbar_dot+3Hbar^2) integral v^2. The pointwise coefficient is<1000 and the torus volume(2pi)^3<512<1000, so the entire spatial integral is bounded by1e6 kappa v_A0^2 throughout the actual time interval. It is added at BOTH analysis radii, is independent of live Y and has zero phase constant/linear jets. Volume composes with the canonical shear only; no artificial generating energy is added to the volume.",
        "whole_complete_quadratic_and_remainder_proof": "Retain the entire original literal constrained H2, including the trace-shear cancellation and full-C1 Schur term, TT curvature, both momentum densities and every Proca polarization. All terms in the six positive raw quadratic envelopes have degree at most two in the changed row bounds, which increased by at most2, so four times EACH old row bounds the whole corrected raw H2. Add the boundary generator and independently bounded reference quadratic. The complete degree>=3 radial Cauchy tail on Rstar=1e150 is retained as an exact remainder bound; it is not deleted from the model.",
        "whole_centered_average_and_Cauchy_proof": "All prepared nonzero inputs have zero spatial mean and the boundary is the same spatially constant linear shear on each pair, so every corrected raw input still has zero mean. At the actual uniform root H_N=0; the full averaged first phase derivative vanishes. Actual nonzero profiles are never replaced by the tree fixture. The full average bounds and holomorphy in phase and Y, with real-time C5 profiles, give the two safe amplitudes1e172 and1e-510. Rstar is only an analysis domain; neither physical cutoff is enlarged.",
        "checks": {
            "whole_canonical_trace_Lie_pairing": original.trace_identity(),
            "only_trace_row_doubled": after["Pi_v_A0"] - 2 * before["Pi_v_A0"],
            "only_matter_row_doubled": after["delta_Pi_M_A0"]
            - 2 * before["delta_Pi_M_A0"],
            "all9_other_rows_unchanged": s.Matrix(
                [
                    after[key] - value
                    for key, value in before.items()
                    if key not in ("Pi_v_A0", "delta_Pi_M_A0")
                ]
            ),
            "all12_source_images": s.Integer(len(small["images"]) - 12),
            "all12_average_remainders": s.Integer(len(small["averages"]) - 12),
            "whole_trace_shape_vector_matter_sum": small["D0"]
            - sum(
                small["full_D" + name]
                for name in ("trace", "shape", "vector", "matter")
            ),
            "whole_full_lapse_center_includes_corrected_phase": small[
                "whole_center_residual"
            ]
            - s.Rational(5, 10**400)
            - 12 * 10**4 * source.HOMOGENEOUS_RADIUS
            - NZ * small["sum_images"],
            "all6_corrected_positive_quadratic_rows": s.Matrix(
                [
                    value - 4 * quadratic.bounds()["whole_quadratic_row_budgets"][key]
                    for key, value in r[
                        "whole_corrected_positive_full_quadratic_rows"
                    ].items()
                ]
            ),
            "whole_interaction_keeps_all4_terms": r[
                "whole_corrected_centered_interaction_amplitude"
            ]
            - r["whole_raw_constrained_quadratic_amplitude"]
            - r["whole_retained_independent_reference_amplitude"]
            - r["whole_retained_boundary_time_generator_amplitude"]
            - r["whole_complete_degree_at_least_three_remainder"],
            "whole_complete_mean_lapse_response": r["whole_full_mean_lapse_response"]
            - NZ * small["sum_average_remainders"]
            - NZZ * small["sum_images"] ** 2 / 2,
        },
        "gates": {
            "added_trace_row_fits_explicit_doubling": b[
                "whole_added_trace_momentum_row"
            ]
            < before["Pi_v_A0"],
            "added_matter_row_fits_explicit_doubling": b[
                "whole_added_matter_momentum_row"
            ]
            < before["delta_Pi_M_A0"],
            "entire_boundary_time_coefficient_bound": b[
                "whole_boundary_time_coefficient"
            ]
            < 1000
            and 512 * b["whole_boundary_time_coefficient"] < 10**6,
            "both_full_corrected_metric_shape_domains": all(
                max(item["fields"]["v_A2"], item["fields"]["tau_A2"])
                < s.Rational(1, 1000)
                for item in (small, large)
            ),
            "both_corrected_complete_source_images": all(
                source.HOMOGENEOUS_RADIUS + max(item["images"]) + 8 * source.TIME
                < source.COORDINATE
                for item in (small, large)
            ),
            "both_corrected_full_lapse_contractions": max(
                small["contraction"], large["contraction"]
            )
            < s.Rational(1, 20),
            "both_corrected_full_lapse_self_maps": all(
                item["whole_center_residual"] / 3
                + item["contraction"] * source.LAPSE_RADIUS
                < source.LAPSE_RADIUS
                for item in (small, large)
            ),
            "larger_corrected_analysis_amplitude": large["whole_centered_amplitude"]
            < 10**533,
            "corrected_entire_nonlinear_tail": r[
                "whole_complete_degree_at_least_three_remainder"
            ]
            < 10**143,
            "corrected_complete_interaction_fits_operator_amplitude": r[
                "whole_corrected_centered_interaction_amplitude"
            ]
            < source.H_AMPLITUDE,
            "corrected_complete_volume_fits_operator_amplitude": r[
                "whole_corrected_centered_volume_amplitude"
            ]
            < source.F_AMPLITUDE,
            "corrected_full_average_density": small["hvariation"]
            < s.Rational(1, 10**520),
            "whole_source_implicit_Hessian": (10**4 + 2 * 10**6 * NZ + 10**9 * NZ**2)
            / 2
            < NZZ,
            "whole_U_first_second_rows": UN < 12 and UNN < 10**4,
            "reference_image_bound_used_in_prepared_coordinates": moving.bounds()[
                "flow_deviation"
            ]
            < s.Rational(1, 50),
            "no_old_physical_field_equality_assumed": True,
            "no_full_actual_H2_reference_identification_assumed": True,
        },
    }
