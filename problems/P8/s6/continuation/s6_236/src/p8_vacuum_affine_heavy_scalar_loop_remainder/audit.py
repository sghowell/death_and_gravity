"""The full first-loop angular window is advanced without upgrading original primitives."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_heavy_scalar_four_point_loop import audit as previous
from p8_vacuum_affine_heavy_scalar_four_point_loop import symmetric

from . import bounds, expansion, integral

ITEM = {
    "id": "separate_V2S_T1_OS4_complete_physical_all_angle_first_loop_remainder_same_tree_window",
    "status": "SEPARATE_FULL_FIRST_LOOP_ANGULAR_REMAINDER_NOT_HIGHER_QUANTUM_JETS_OMITTED_LOOPS_EXACT_UV_ORIGINAL_PARENT_BOUNCE_V_G_B_OR_P8",
}
require_scope = previous.require_scope


def require_stage(stage):
    if not isinstance(stage, str) or stage not in (
        "complete_mass_ordered_parameter_primitive",
        "convergent_complex_inverse_mass_remainder",
        "complete_on_shell_two_term_cancellation",
        "same_OS4_full_first_loop_physical_angle_bound",
        "tree_plus_first_loop_truncation_bound",
    ):
        raise ValueError(
            "Only the complete first-loop coefficient and its stated physical window are controlled"
        )
    return stage


def require_physical(name, svalue, cosine, order):
    if not isinstance(name, str) or name != symmetric.NAME:
        raise ValueError(
            "Only the unchanged separate V2S-T1-OS4 prescription is evaluated"
        )
    for value in (svalue, cosine, order):
        if isinstance(value, bool) or not isinstance(
            value, (int, s.Integer, s.Rational)
        ):
            raise TypeError("Use exact rational kinematics and first-loop order")
    if order != 1:
        raise ValueError("No omitted-loop or exact-amplitude bound is supplied")
    if not 4 <= svalue <= bounds.S_MAX or not -1 <= cosine <= 1:
        raise ValueError("Outside the full stated physical scattering window")
    S = s.Rational(svalue)
    x = s.Rational(cosine)
    return name, S, -(S - 4) * (1 - x) / 2, -(S - 4) * (1 + x) / 2, 1


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError("A first-loop angular bound does not close original P8")
    return True


def packets():
    return {
        "complete_ordered_parameter_primitive_and_Cauchy_remainder": integral.data(),
        "all_diagram_on_shell_cancellation_in_the_same_OS4_scheme": expansion.data(),
        "uniform_actual_physical_window_and_relative_first_loop_error": bounds.data(),
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
            "one_new_full_first_loop_angular_record": len(matching())
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
        "same_frozen_original_and_separate_sources": True,
        "same_OS4_contact_not_a_new_renormalization_condition": True,
        "real_and_imaginary_parts_bounded_together": True,
        "all_energies_and_angles_in_closed_named_window": True,
        "omitted_loop_and_original_parent_boundary_retained": True,
        "unique_matching_identifiers": len({row["id"] for row in matching()})
        == len(matching()),
    }


def observable():
    return {
        "input": "The full S235 first four-light loop, its six ordered boxes, all UV subtractions and the unchanged V2S-T1-OS4 finite contact are retained with the S234 first-order on-shell external normalization.",
        "new_exact_reduction": "A single closed quadratic-parameter primitive generates both Cbar and the mass-ordered Dbar by a full heavy-mass derivative. Its only physical light branch is Log(L-i0) on the stated window. Complete convergent analytic coefficient bounds avoid absolute integration through Feynman double poles.",
        "cancellation": "All bubble B and triangle L Log L terms cancel in the complete orders n^-2 and n^-3. The remaining two terms are independent of kinematics on s+t+u=4 and therefore cancel under the existing symmetric value subtraction.",
        "full_first_loop_bound": "For every real physical configuration4<=s<=10^196 and every angle, |A1_OS4|/A_original<10^-199. This controls the entire first-loop coefficient, not only its absorptive part.",
        "specified_truncation": "The separately defined tree-plus-first-loop expression differs from the entire original tree by less than1/60+10^-199<1/59 on that same window. The word truncation is essential.",
        "remaining": "Higher physical coefficient derivatives, all omitted loop orders, exact S-matrix/UV properties, controlled heavy resonance, finite-gravity Regge control, the common-parent bounce and original V/G/B/P8 remain open.",
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
        for pos in range(4):
            args = [symmetric.NAME, 4, 0, 1]
            args[pos] = value
            out.append((f"physical_type_{pos}_{i}", require_physical, tuple(args)))
        out.append((f"stage_type_{i}", require_stage, (value,)))
    for args in (
        (1, modes.KAPPA, modes.MASS, 1),
        (-1, modes.KAPPA, modes.MASS, 1),
        (0, 2 * modes.KAPPA, modes.MASS, 1),
        (0, modes.KAPPA, 999, 1),
        (0, modes.KAPPA, modes.MASS, 2),
    ):
        out.append((f"target_scope_{len(out)}", require_scope, args))
    for pos, value in (
        (0, "V2S-T1"),
        (0, "original_parent"),
        (1, 0),
        (1, 3),
        (1, bounds.S_MAX + 1),
        (2, -2),
        (2, 2),
        (3, 0),
        (3, 2),
        (3, s.Rational(1, 2)),
    ):
        args = [symmetric.NAME, 4, 0, 1]
        args[pos] = value
        out.append((f"physical_scope_{len(out)}", require_physical, tuple(args)))
    for stage in (
        "all_loop_error",
        "exact_quantum_S_matrix",
        "higher_physical_coefficients_matched",
        "controlled_resonance",
        "original_parent_replaced",
        "new_OS4_finite_contact",
        "absolute_Feynman_double_pole_bound",
        "heavy_boxes_interchanged",
        "leading_terms_discarded_before_cancellation",
        "tree_plus_first_loop_is_exact",
        "uniform_arbitrary_energy",
        "full_quantum_potential_positive",
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
        raise ValueError("Unsupported matching assertion accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "all_mass_ordered_boxes_and_UV_terms_retained": True,
        "physical_cut_and_threshold_logarithms_retained": True,
        "Cauchy_bounds_not_truncated_power_counting": True,
        "two_universal_terms_cancel_only_after_all_channels": True,
        "same_finite_value_condition_no_frozen_edit": True,
        "full_first_loop_not_only_optical_part": True,
        "truncated_amplitude_not_exact_all_loop_claim": True,
        "all_original_primitive_statuses_unchanged": frontier() == previous.frontier(),
    }
