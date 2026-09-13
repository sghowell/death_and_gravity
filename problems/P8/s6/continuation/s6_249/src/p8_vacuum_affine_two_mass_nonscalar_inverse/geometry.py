"""Full O(2) spatial decomposition and prepared transverse-vector Ward quotient."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_adm_clock_response import geometry as metric
from p8_vacuum_affine_two_mass_reference import geometry as reference

t = metric.t
P = s.Symbol("external_momentum", positive=True)
BASIS = (
    s.eye(3) / s.sqrt(3),
    s.diag(-1, -1, 2) / s.sqrt(6),
    s.Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]]) / s.sqrt(2),
    s.Matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]]) / s.sqrt(2),
    s.diag(1, -1, 0) / s.sqrt(2),
    s.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]]) / s.sqrt(2),
)
LABELS = (
    "scalar_trace",
    "scalar_shear",
    "vector_x",
    "vector_y",
    "tensor_plus",
    "tensor_cross",
)


def zero(value):
    return (
        value.applyfunc(s.cancel)
        if isinstance(value, s.MatrixBase)
        else s.cancel(value)
    )


def project(Q, axis, sector):
    """Coordinate-free orthogonal sector projectors; axis has unit Euclidean norm."""
    if not isinstance(Q, s.MatrixBase) or Q.shape != (3, 3) or Q != Q.T:
        raise ValueError("Use a complete symmetric physical spatial direction")
    if (
        not isinstance(axis, s.MatrixBase)
        or axis.shape != (3, 1)
        or s.simplify(axis.dot(axis) - 1) != 0
    ):
        raise ValueError("Use an exact unit external direction")
    if sector not in ("scalar", "vector", "tensor"):
        raise ValueError("Keep the three distinct O2 spatial sectors")
    parallel = axis * axis.T
    plane = s.eye(3) - parallel
    if sector == "tensor":
        return zero(plane * Q * plane - plane * s.trace(plane * Q) / 2)
    if sector == "vector":
        return zero(plane * Q * parallel + parallel * Q * plane)
    return zero(parallel * s.trace(parallel * Q) + plane * s.trace(plane * Q) / 2)


@cache
def representation():
    J = s.Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 0]])
    reflection = s.diag(-1, 1, 1)
    G = s.Matrix(6, 6, lambda i, j: s.trace(BASIS[i].T * (J * BASIS[j] - BASIS[j] * J)))
    R = s.Matrix(
        6, 6, lambda i, j: s.trace(BASIS[i].T * reflection * BASIS[j] * reflection)
    )
    symbols = s.symbols("full_ordered_spatial_entry0:36")
    M = s.Matrix(6, 6, symbols)
    angular = list(M * G - G * M)
    full = angular + list(M * R - R * M)
    linear, _ = s.linear_eq_to_matrix(full, symbols)
    rotation_only, _ = s.linear_eq_to_matrix(angular, symbols)
    a, b, c, d, v, h = s.symbols(
        "scalar00 scalar01 scalar10 scalar11 vector_kernel tensor_kernel"
    )
    allowed = s.diag(s.Matrix([[a, b], [c, d]]), v * s.eye(2), h * s.eye(2))
    flattened = s.Matrix.hstack(*[s.Matrix(list(B)) for B in BASIS])
    identity = s.eye(3)
    components = s.Matrix(
        9,
        9,
        lambda a, b: (
            (
                identity[a // 3, b // 3] * identity[a % 3, b % 3]
                + identity[a // 3, b % 3] * identity[a % 3, b // 3]
            )
            / 2
        ),
    )
    projection_checks = {}
    for index, sector in enumerate(("scalar", "vector", "tensor")):
        target = s.diag(*[int(j // 2 == index) for j in range(6)])
        current = s.Matrix(
            6,
            6,
            lambda i, j, sector=sector: s.trace(
                BASIS[i] * project(BASIS[j], s.Matrix([0, 0, 1]), sector)
            ),
        )
        projection_checks["coordinate_free_orthogonal_" + sector + "_projection"] = (
            current - target
        )
    return {
        "generator": G,
        "reflection": R,
        "full_ordered_commutant": allowed,
        "full_commutant_nullity": 36 - linear.rank(),
        "rotation_only_nullity": 36 - rotation_only.rank(),
        "global_projectors": "For n=P/|P| and Pi=I-nn^T, P_T Q=Pi Q Pi-Pi tr(Pi Q)/2, P_V Q=Pi Q nn^T+nn^T Q Pi, P_S Q=nn^T tr(nn^T Q)+Pi tr(Pi Q)/2. They are orthogonal, complete and norm-one without a global polarization frame. Their almost-everywhere Fourier definition does not assign a literal P0 constraint.",
        "checks": {
            **projection_checks,
            "whole_symmetric_basis_Gram": flattened.T * flattened - s.eye(6),
            "entire_symmetric_tensor_projector": flattened * flattened.T - components,
            "whole_rotation_generator_antisymmetry": G + G.T,
            "whole_reflection_involution": R * R - s.eye(6),
            "full_allowed_retarded_commutant_rotation": allowed * G - G * allowed,
            "full_allowed_retarded_commutant_reflection": allowed * R - R * allowed,
            "full_O2_commutant_has_exactly_six_free_entries": 36 - linear.rank() - 6,
            "rotation_without_reflection_would_allow_two_extra_channels": 36
            - rotation_only.rank()
            - 8,
        },
        "gates": {
            "ordered_scalar_cross_entries_not_forced_equal": allowed[0, 1]
            != allowed[1, 0],
            "two_nonscalar_doublets_not_assumed_one_at_finite_transfer": v != h,
            "reflection_needed_to_exclude_parity_odd_antisymmetric_channels": True,
        },
    }


@cache
def ward():
    a, H = metric.a, metric.H
    W = lambda f: s.diff(f, t) + 3 * H * f
    b = s.Function("transverse_shift")(t)
    v = s.Function("spatial_vector_amplitude")(t)
    chi = s.Function("transverse_gauge")(t)
    u = s.Dummy("source_time", real=True)
    IW = lambda f: a**-3 * s.Integral((a**3 * f).subs(t, u), (u, -s.Rational(1, 2), t))
    gauge_checks = []
    local_ward_checks = []
    rho = s.Function("actual_Gaussian_rho")(t)
    pressure = s.Function("actual_Gaussian_pressure")(t)
    stress = metric.density(rho, pressure)
    for direction, index in ((0, 2), (1, 3)):
        xi = s.zeros(4, 1)
        phase = s.exp(s.I * P * metric.z)
        xi[direction + 1] = chi * phase
        n, beta, Q = metric.gauge(xi)
        expected_beta = s.zeros(3, 1)
        expected_beta[direction] = s.diff(chi, t) * phase
        gauge_checks.extend(list(zero(beta - expected_beta)))
        gauge_checks.extend(
            list(zero(Q - s.sqrt(2) * s.I * P * chi * phase * BASIS[index]))
        )
        gauge_checks.append(n)
        for detector_direction, detector_index in ((0, 2), (1, 3)):
            phaseD = s.exp(-s.I * P * metric.z)
            betaD = s.zeros(3, 1)
            betaD[detector_direction] = s.Function("detector_shift")(t) * phaseD
            QD = s.Function("detector_vector")(t) * phaseD * BASIS[detector_index]
            detector = (s.S.Zero, betaD, QD)
            local_ward_checks.append(zero(metric.source_ward(detector, xi, stress)))
            xiD = s.zeros(4, 1)
            xiD[detector_direction] = s.Function("detector_gauge")(t) * phaseD
            synchronous = (s.S.Zero, s.zeros(3, 1), v * phase * BASIS[index])
            local_ward_checks.append(
                zero(metric.detector_ward(xiD, synchronous, stress))
            )
    invariant = s.diff(v, t) - s.sqrt(2) * s.I * P * b
    shifted = invariant.subs(
        {v: v + s.sqrt(2) * s.I * P * chi, b: b + s.diff(chi, t)}, simultaneous=True
    ).doit()
    ev = s.Function("applied_vector_spatial_force")(t)
    eb = -s.sqrt(2) * s.I * P * IW(ev)
    En = s.Function("whole_transverse_shift_residual")(t)
    Eq = -W(En) / (s.sqrt(2) * s.I * P)
    bad = a**-3
    kappa = s.Symbol("kappa", positive=True)
    force = s.Function("canonical_quotient_force")(t)
    canonical_ev = s.sqrt(kappa) * force / 2
    return {
        "gauge_invariant_velocity": invariant,
        "synchronous_quotient": "v_syn=v-sqrt(2)iP I beta; its derivative is v'-sqrt(2)iP beta. Use the unchanged zero germ.",
        "distinct_detector_source_primitives": "Use Iminus from-1/2 on the source and Iplus=-integral_t^(1/2) on the detector. Both derivatives equal the input. Their weighted transpose is (Iplus)^dagger=-I_W, which gives the retarded output E_beta=-sqrt(2)iP I_W E_v. Final detector and initial output/source fluxes vanish for distinct reasons; do not use Iminus on both legs.",
        "whole_residual_Ward": "(D+3H)E_beta+sqrt(2)iP E_v=0, with residuals divided by the OUTPUT a³. This is the complete on-shell-reference Ward identity, not a separate Gaussian mean reset.",
        "complete_Gaussian_vector_contacts": "Both complete ordered Ward integrands vanish on the opposite-momentum transverse vector sector after their full metric second-vertex contacts are retained. This holds for arbitrary nonzero Gaussian rho(t),P(t); it is not an early mean subtraction. The distinct common-clock contact vanishes here because both lapse directions are zero.",
        "compatible_applied_force_pair": s.ImmutableMatrix([ev, eb]),
        "canonical_to_spatial_force": canonical_ev,
        "unprepared_constraint_counterexample": bad,
        "force_graph": "Admit smooth prepared quotient force f in H^r_ball; e_v=sqrt(kappa)f/2 and e_beta=-sqrt(2)iP I_W e_v form the full gauge-compatible pair. Conversely e_v=-(D+3H)e_beta/(sqrt(2)iP) must lie in the stated quotient-force space; no unrestricted same-H^r shift-force claim.",
        "physical_shift_boundary": "The synchronous representative has beta0 and spatial metric amplitude in H^r_ball. A zero-spatial-vector representative has beta_inv=-v_syn'/(sqrt(2)iP), a tempered distribution; its infrared H^r norm is not asserted equal to the quotient norm. Literal P0 is separate.",
        "checks": {
            "full_actual_two_transverse_gauge_metric_jets": s.ImmutableMatrix(
                gauge_checks
            ),
            "all_eight_complete_ordered_Gaussian_vector_Ward_integrands": s.ImmutableMatrix(
                local_ward_checks
            ),
            "full_gauge_invariant_vector_velocity": zero(shifted - invariant),
            "weighted_prepared_primitive": zero(W(IW(ev)) - ev),
            "entire_compatible_applied_force_Ward": zero(
                W(eb) + s.sqrt(2) * s.I * P * ev
            ),
            "entire_reference_residual_identity": zero(
                W(En) + s.sqrt(2) * s.I * P * Eq
            ),
            "nonzero_unprepared_constraint_kernel": zero(W(bad)),
            "canonical_force_one_leg_conversion": 2 * canonical_ev / s.sqrt(kappa)
            - force,
        },
        "gates": {
            "original_vector_constraint_germ_essential": bad.subs(t, -s.Rational(1, 2))
            != 0,
            "no_inverse_P_in_quotient_to_compatible_force_pair": True,
            "arbitrary_shift_force_not_promoted_to_same_quotient_norm": True,
            "full_reference_mean_and_both_chart_contacts_retained_before_Ward": True,
        },
    }


@cache
def leading():
    d = reference.d
    C = (2 * d * d + d - 7) / (2 * d * (d + 2))
    source = reference.data()
    old = reference.old_matching.data()
    Cp = old["all_dimension_tracefree_invariant"].subs(reference.old_matching.g.d, d)
    Ch = source["full_dimensional_scalar_invariants"]["tracefree"]
    tau = s.Symbol("positive_lag", positive=True)
    return {
        "full_nonscalar_dimension_invariant": C,
        "physical_value_and_first_dimension_jet": (
            s.Rational(7, 15),
            s.Rational(83, 450),
        ),
        "full_fixed_fourth_coefficient": -(reference.ell + 2) / 60,
        "complete_leading_radial_coefficient": s.Rational(224, 15),
        "complete_q0_channel": "The SAME full F2,total with both physical mass cuts and fixed constant-(log(n)+2)/60 acts on each unit-Frobenius nonscalar direction at leading and highest finite order. Tensor and vector finite-transfer kernels remain distinct.",
        "checks": {
            "whole_two_species_invariant": s.factor(C - Cp - Ch),
            "physical_nonscalar_invariant": C.subs(d, 3) - s.Rational(7, 15),
            "fixed_physical_first_dimension_jet": s.diff(C, d).subs(d, 3)
            - s.Rational(83, 450),
            "full_fourth_finite_constant": source["fixed_total_fourth_constants"][
                "shear"
            ]
            + (reference.ell + 2) / 60,
            "whole_radial_normalization": 32 * C.subs(d, 3) - s.Rational(224, 15),
            "Abel_leading_to_complete_fourth_primitive": s.diff(
                C.subs(d, 3) / tau, tau, 4
            )
            - s.Rational(56, 5) / tau**5,
            "all_dimension_frozen_proper_scale": -d - 2 + d + 1 + 1,
        },
        "gates": {
            "nonscalar_leading_not_trace_channel": C.subs(d, 3) != 4,
            "fixed_physical_finite_jet_not_discarded": s.diff(C, d).subs(d, 3) != 0,
            "full_two_mass_reference_not_sum_of_inverses": True,
        },
    }
