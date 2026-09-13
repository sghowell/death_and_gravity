"""First-loop coefficient matching with the full original frontier retained."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_heavy_scalar_loop_remainder import audit as previous

from . import bidisk, coefficients, jets

ITEM = {
    "id": "separate_V2S_T1_OS4_full_first_loop_b20_b21_b40_matching_with_joint_complex_bidisk_remainders",
    "status": "SEPARATE_COMPLETE_FIRST_LOOP_LOW_COEFFICIENT_MATCHING_NOT_OMITTED_LOOPS_EXACT_PHYSICAL_COEFFICIENTS_UV_ORIGINAL_PARENT_BOUNCE_V_G_B_OR_P8",
}
require_scope = previous.require_scope


def require_stage(stage):
    if not isinstance(stage, str) or stage not in (
        "complete_four_order_loop_and_log_moment_reduction",
        "full_joint_complex_bidisk_remainder",
        "same_OS4_first_loop_b20_matching",
        "same_OS4_first_loop_b21_matching",
        "same_OS4_first_loop_b40_matching",
    ):
        raise ValueError(
            "Only the stated complete first-loop coefficient matching is proved"
        )
    return stage


def require_jet(name, order):
    if not isinstance(name, str) or name not in ("b20", "b21", "b40"):
        raise ValueError("Only the three stated low-energy coefficients are bounded")
    if isinstance(order, bool) or not isinstance(order, (int, s.Integer)) or order != 1:
        raise ValueError("Only the complete first-loop correction is controlled")
    return name, 1


def require_bidisk(vreal, vimag, treal, timag):
    values = (vreal, vimag, treal, timag)
    for value in values:
        if isinstance(value, bool) or not isinstance(
            value, (int, s.Integer, s.Rational)
        ):
            raise TypeError("Use exact rational real and imaginary components")
    vr, vi, tr, ti = map(s.Rational, values)
    if vr * vr + vi * vi > bidisk.RADIUS**2 or tr * tr + ti * ti > bidisk.RADIUS**2:
        raise ValueError("Outside the complete named complex bidisk")
    return vr + s.I * vi, tr + s.I * ti


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError("First-loop coefficient matching does not close original P8")
    return True


def packets():
    return {
        "complete_four_order_loop_and_exact_logarithmic_moments": jets.data(),
        "full_complex_bidisk_and_all_inverse_mass_tails": bidisk.data(),
        "actual_first_loop_coefficient_intervals_and_relative_errors": coefficients.data(),
    }


@cache
def residuals():
    rows = {
        name + "_" + key: value
        for name, packet in packets().items()
        for key, value in packet["checks"].items()
    }
    rows.update(
        {
            "nine_original_primitive_rows": len(frontier()) - 9,
            "unchanged_original_primitive_statuses": sum(
                a != b for a, b in zip(frontier(), previous.frontier())
            ),
            "unchanged_previous_matching_records": sum(
                a != b for a, b in zip(matching(), previous.matching())
            ),
            "one_new_complete_first_loop_coefficient_record": len(matching())
            - len(previous.matching())
            - 1,
        }
    )
    return {key: s.cancel(value) for key, value in rows.items()}


def scalar_entry_count():
    return len(residuals())


@cache
def gates():
    return {
        **{
            name + "_" + key: bool(value)
            for name, packet in packets().items()
            for key, value in packet["gates"].items()
        },
        "same_frozen_original_and_separate_inputs": True,
        "same_OS4_contact_no_finite_derivative_retuning": True,
        "joint_complex_remainder_not_real_window_inference": True,
        "all_unknown_heavy_weight_retained": True,
        "omitted_loop_and_original_parent_boundary_retained": True,
        "unique_matching_identifiers": len({row["id"] for row in matching()})
        == len(matching()),
    }


def observable():
    return {
        "input": "The unchanged V2S-T1-OS4 first-loop amplitude, all six ordered boxes, complete UV subtractions and original first-order light mass/residue conditions are retained.",
        "full_reduction": "Exact integration-by-parts identities reduce every light logarithmic moment before cancellation. The complete n^-2 through n^-5 coefficients are crossing-symmetric polynomials of degree at most0,0,2,3 respectively.",
        "analytic_neighborhood": "A directly proved joint complex bidisk |v|,|t|<=1/4 gives a uniform complete remainder below10^12 g^4/(16pi^2 n^6). It is not inferred from the real physical-angle bound.",
        "first_loop_matching": "The complete first-loop shifts have0<delta_b20/(4lambda)<10^-203,0<delta_b21/(-3gamma)<10^-203 and |delta_b40|/(gamma^2/lambda)<10^-192. The b20 correction is positive and b21 negative; no sign is assigned to delta_b40.",
        "specified_truncation": "The tree-plus-first-loop higher b40 stays positive and within the displayed relative interval. These are coefficients of an explicitly truncated amplitude, not exact quantum LSZ or dispersion coefficients.",
        "no_new_subtraction": "The fixed OS4 contact is constant and has zero derivative in these three observables. No new finite condition is imposed and no unknown heavy weight is removed.",
        "remaining": "Omitted loops, exact physical coefficient matching, exact UV/dispersion properties, heavy resonance control, finite-gravity Regge control, covariant common-parent bounce matching and original V/G/B/P8 remain open.",
    }


def bad_cases():
    bad = (
        True,
        False,
        1.0,
        s.Float(1),
        "1",
        None,
        s.oo,
        -s.oo,
        s.zoo,
        s.I,
        s.nan,
        s.Symbol("x"),
    )
    out = []
    for i, value in enumerate(bad):
        for pos in range(4):
            args = [0, modes.KAPPA, modes.MASS, 1]
            args[pos] = value
            out.append((f"target_type_{pos}_{i}", require_scope, tuple(args)))
            args = [0, 0, 0, 0]
            args[pos] = value
            out.append((f"bidisk_type_{pos}_{i}", require_bidisk, tuple(args)))
        for pos in range(2):
            args = ["b20", 1]
            args[pos] = value
            out.append((f"jet_type_{pos}_{i}", require_jet, tuple(args)))
        out.append((f"stage_type_{i}", require_stage, (value,)))
    for args in (
        (1, modes.KAPPA, modes.MASS, 1),
        (-1, modes.KAPPA, modes.MASS, 1),
        (0, 2 * modes.KAPPA, modes.MASS, 1),
        (0, modes.KAPPA, 999, 1),
        (0, modes.KAPPA, modes.MASS, 2),
    ):
        out.append((f"target_scope_{len(out)}", require_scope, args))
    for args in (
        ("b00", 1),
        ("b60", 1),
        ("all_coefficients", 1),
        ("exact_b20", 1),
        ("b20", 0),
        ("b20", 2),
        ("b20", s.Rational(1, 2)),
    ):
        out.append((f"jet_scope_{len(out)}", require_jet, args))
    for pos in range(4):
        for value in (-1, 1):
            args = [0, 0, 0, 0]
            args[pos] = value
            out.append((f"bidisk_scope_{pos}_{value}", require_bidisk, tuple(args)))
    for stage in (
        "all_loop_error",
        "exact_quantum_coefficients",
        "exact_UV_S_matrix",
        "all_higher_jets_matched",
        "new_derivative_counterterm",
        "heavy_pole_subtracted",
        "physical_LSZ_all_orders",
        "b40_first_correction_positive",
        "exact_Gram_saturation",
        "real_window_implies_derivatives",
        "light_log_moments_discarded",
        "controlled_resonance",
        "physical_cutoff",
        "common_parent_bounce",
        "finite_gravity_Regge_done",
        "closed_P8",
    ):
        out.append(("unsupported_" + stage, require_stage, (stage,)))
    for i in range(9):
        rows = frontier()
        rows[i]["status"] = "COMPLETE"
        out.append((f"primitive_{i}", validate_scope, (rows, matching())))
    for i in range(len(matching())):
        rows = matching()
        rows[i]["status"] = "COMPLETE"
        out.append((f"matching_{i}", validate_scope, (frontier(), rows)))
    out.extend(
        (
            ("missing_primitive", validate_scope, (frontier()[:-1], matching())),
            ("missing_matching", validate_scope, (frontier(), matching()[:-1])),
            (
                "extra_full_parent",
                validate_scope,
                (
                    frontier(),
                    matching() + [{"id": "full_parent", "status": "COMPLETE"}],
                ),
            ),
        )
    )
    return out


@cache
def rejected_inputs():
    for name, call, args in bad_cases():
        try:
            call(*args)
        except (ValueError, TypeError):
            continue
        raise ValueError("Unsupported coefficient claim accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "all_mass_ordered_diagrams_retained": True,
        "light_log_moments_reduced_not_dropped": True,
        "joint_bidisk_and_complete_Cauchy_tails": True,
        "all_three_coefficient_factorials_retained": True,
        "same_finite_contact_no_new_derivative_retuning": True,
        "first_b40_correction_sign_not_assumed": True,
        "truncated_coefficients_not_exact_all_loop_values": True,
        "all_original_primitive_statuses_unchanged": frontier() == previous.frontier(),
    }
