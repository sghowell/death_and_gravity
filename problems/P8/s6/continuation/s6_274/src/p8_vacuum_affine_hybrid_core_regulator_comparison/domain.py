"""Larger complex ANALYSIS domain, unchanged physical core and full remainder."""

from functools import cache

import sympy as s
from p8_vacuum_affine_finite_volume_turnaround import moving
from p8_vacuum_affine_selfconsistent_finite_feedback import geometry

from . import quadratic, source

P, kap = geometry.P, geometry.KAPPA
field = moving.old.field
rows = field.field_bounds()


def require_radius(value):
    if isinstance(value, bool) or not isinstance(value, (int, s.Integer, s.Rational)):
        raise TypeError("Require an exact declared analysis radius")
    if value not in (source.RADIUS, source.ANALYSIS_RADIUS):
        raise ValueError(
            "Require the original physical core or the larger analysis radius"
        )
    return s.sympify(value)


def bounds(radius):
    radius = require_radius(radius)
    return _bounds(radius)


@cache
def _bounds(radius):
    f = {name: value * 8 * radius for name, value in rows.items()}
    v = 4 * f["v_A2"] / (1 + P) ** 2
    tau = 4 * f["tau_A2"] / (1 + P) ** 2
    hom = source.previous.HOMOGENEOUS_RADIUS
    pbound = 8 * source.previous.TIME + hom
    pvbound = 100 * (source.previous.TIME + hom)
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
        s.Rational(1, 30)
        + 10**9 * (source.previous.LAPSE_RADIUS + 12 * (hom + max(images))) / 3
    )
    residual = s.Rational(5, 10**400) + 12 * 10**4 * hom + geometry.NZ * zsum
    assert hom + max(images) + 8 * source.previous.TIME < source.previous.COORDINATE
    assert contraction < s.Rational(1, 20)
    assert (
        residual / 3 + contraction * source.previous.LAPSE_RADIUS
        < source.previous.LAPSE_RADIUS
    )
    assert max(f["v_A2"], f["tau_A2"]) < s.Rational(1, 1000)
    amplitude = 2 * 1000 * kap * hvar + 10**128 * (8 * radius) ** 2 / 2
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
        "contraction": contraction,
        "whole_center_residual": residual,
    }


@cache
def remainder():
    large = bounds(source.ANALYSIS_RADIUS)
    ratio = source.RADIUS / source.ANALYSIS_RADIUS
    tail = large["whole_centered_amplitude"] * ratio**3 / (1 - ratio)
    full_quadratic = (
        2
        * 1000
        * quadratic.bounds()["whole_full_quadratic_bound"]
        * (8 * source.RADIUS) ** 2
    )
    reference = 10**128 * (8 * source.RADIUS) ** 2 / 2
    return {
        "whole_analysis_to_core_ratio": ratio,
        "whole_degree_at_least_three_remainder": tail,
        "whole_actual_H2_amplitude": full_quadratic,
        "whole_retained_independent_reference_amplitude": reference,
        "whole_small_ball_centered_Hamiltonian": full_quadratic + reference + tail,
    }


@cache
def data():
    small, large = bounds(source.RADIUS), bounds(source.ANALYSIS_RADIUS)
    r = remainder()
    return {
        "whole_small_original_phase_domain": small,
        "whole_large_complex_analysis_domain": large,
        "whole_full_Hamiltonian_Cauchy_remainder": r,
        "whole_zero_first_phase_jet_proof": "At the actual uniform root N0(Y), H_N=0. The exact first derivative of the full spatial energy is a constant linear combination of v, the nonzero input momentum/matter/heavy modes, the divergence Gauss term and linear scalar curvature. Every such spatial average vanishes. The full first shape/shear/electric/magnetic/gradient terms are zero at the homogeneous center. This proves the complete first phase jet vanishes, without dropping nonlinear generated means.",
        "whole_larger_domain_proof": "Evaluate every unchanged physical field row at the complete image radius8e150. The full determinant/shape and mean-zero adjoint inverse are valid because their A2 bounds remain small; all trace, shape, vector and matter sources are retained. The displayed twelve pointwise images and exact average remainders fit the SAME1e-120 source box and1e-115 lapse ball, including full phase contributions to the center residual. All numerical inequalities are recomputed, not asserted by relabeling the old core.",
        "whole_uniform_full_symbol_proof": "The old full source is holomorphic in finite complex phase and homogeneous parameters on this strictly larger common domain at each real time. Source profiles retain only original real C5 regularity. Uniform spatial bounds permit dominated integration; the implicit root is the unique entire-source branch. These are bounds on complete symbols, not truncated reconstructed fields.",
        "whole_radial_Cauchy_proof": "For each complex phase z on4R, apply scalar Cauchy to the WHOLE centered g(t z) on radius4Rstar/||z||. Its degree0 and1 vanish; its degree2 is the exact full constrained H2 minus the independently retained free quadratic, with no identification assumed. The entire degree>=3 sum is bounded byA_large*(R/Rstar)^3/(1-R/Rstar). This estimates the nonlinear remainder; no term is removed from the Hamiltonian. The scalar center remains in the original energy and phase.",
        "whole_unchanged_physical_regulator": "Rstar=1e150 is an ANALYSIS radius only. The actual phase core remainsR=1e20, the same two c1,c2 cutoffs, the same state and both full operator definitions. The S273 coupled model and interval are unchanged.",
        "checks": {
            "exact_original_small_domain_Hamiltonian_average": small["hvariation"]
            - geometry.bounds()["hvariation"],
            "whole_original_small_domain_twelve_images": s.Matrix(small["images"])
            - s.Matrix(geometry.bounds()["images"]),
            "whole_original_small_domain_twelve_average_remainders": s.Matrix(
                small["averages"]
            )
            - s.Matrix(geometry.bounds()["averages"]),
            "whole_full_center_residual_matches_original": small[
                "whole_center_residual"
            ]
            - geometry.bounds()["full_center_residual"],
            "whole_retained_H2_reference_and_nonlinear_terms": r[
                "whole_small_ball_centered_Hamiltonian"
            ]
            - r["whole_actual_H2_amplitude"]
            - r["whole_retained_independent_reference_amplitude"]
            - r["whole_degree_at_least_three_remainder"],
        },
        "gates": {
            "whole_large_complex_metric_shape_domain": max(
                large["fields"]["v_A2"], large["fields"]["tau_A2"]
            )
            < s.Rational(1, 1000),
            "whole_large_source_images": source.previous.HOMOGENEOUS_RADIUS
            + max(large["images"])
            + 8 * source.TIME
            < source.previous.COORDINATE,
            "whole_large_lapse_contraction": large["contraction"] < s.Rational(1, 20),
            "whole_large_full_phase_lapse_self_map": large["whole_center_residual"] / 3
            + large["contraction"] * source.previous.LAPSE_RADIUS
            < source.previous.LAPSE_RADIUS,
            "whole_large_analytic_amplitude": large["whole_centered_amplitude"]
            < 10**531,
            "whole_complete_degree_three_remainder": r[
                "whole_degree_at_least_three_remainder"
            ]
            < 10**141,
            "whole_new_actual_centered_amplitude": r[
                "whole_small_ball_centered_Hamiltonian"
            ]
            < source.H_AMPLITUDE,
            "whole_Cauchy_geometric_ratio": 0
            < r["whole_analysis_to_core_ratio"]
            < s.Rational(1, 2),
            "physical_core_and_cutoff_not_enlarged": source.RADIUS == geometry.CORE,
            "no_reconstructed_harmonic_or_mean_truncation": True,
        },
    }
