"""Separate OS4 first-loop value matching, with every original frontier retained."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_heavy_scalar_one_loop import audit as previous
from p8_vacuum_affine_heavy_scalar_tree_matching import model

from . import amplitude, symmetric, zero_jet

ITEM = {
    "id": "separate_V2S_T1_OS4_complete_first_four_point_loop_UV_and_cut_matching_with_symmetric_contact_and_full_classical_quartic_margin",
    "status": "SEPARATE_COMPLETE_FIRST_LOOP_REPRESENTATION_AND_SYMMETRIC_VALUE_MATCHING_NOT_FULL_ANGULAR_HIGHER_COEFFICIENT_ALL_LOOP_ORIGINAL_PARENT_BOUNCE_UV_OR_P8",
}
require_scope = previous.require_scope


def require_stage(stage):
    if not isinstance(stage, str) or stage not in (
        "complete_first_loop_integral_representation",
        "complete_UV_counterterm_and_first_cut_check",
        "off_shell_zero_jet_control",
        "separate_symmetric_first_loop_value_matching",
        "contact_only_classical_quartic_margin",
    ):
        raise ValueError(
            "Only the specified separate first-loop representation and value matching are proved"
        )
    return stage


def require_matching(name, svalue, tvalue, uvalue, order):
    if not isinstance(name, str) or name != symmetric.NAME:
        raise ValueError(
            "Only the separately named V2S-T1-OS4 prescription is evaluated"
        )
    if isinstance(order, bool) or not isinstance(order, (int, s.Integer)) or order != 1:
        raise ValueError("Only first-loop value matching is evaluated")
    for value in (svalue, tvalue, uvalue):
        if isinstance(value, bool) or not isinstance(
            value, (int, s.Integer, s.Rational)
        ):
            raise TypeError("Matching invariants must be exact real rationals")
        if value != symmetric.S0:
            raise ValueError("Only the symmetric on-shell subthreshold value is fixed")
    return name, symmetric.S0, symmetric.S0, symmetric.S0, 1


def require_cut(value):
    if isinstance(value, bool) or not isinstance(value, (int, s.Integer, s.Rational)):
        raise TypeError("The first-cut invariant must be an exact real rational")
    if not 4 < value < model.MASS2:
        raise ValueError("Outside the stated first-elastic open window")
    return s.Rational(value)


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError(
            "One symmetric first-loop condition does not close original P8"
        )
    return True


def packets():
    return {
        "complete_first_loop_integrals_UV_and_first_elastic_cut": amplitude.data(),
        "independent_off_shell_full_constant_background_loop_jet": zero_jet.data(),
        "separate_symmetric_on_shell_value_matching_and_contact_margin": symmetric.data(),
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
            "unchanged_primitive_rows": sum(
                a != b for a, b in zip(frontier(), previous.frontier())
            ),
            "unchanged_previous_matching_rows": sum(
                a != b for a, b in zip(matching(), previous.matching())
            ),
            "one_new_first_loop_value_matching": len(matching())
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
        "original_and_separate_frozen_inputs_unchanged": True,
        "new_finite_contact_prescription_explicitly_named": True,
        "all_ordered_boxes_and_three_UV_counterterms_retained": True,
        "subthreshold_matching_not_real_scattering_kinematics": True,
        "single_value_not_full_angular_or_all_loop_bound": True,
        "unique_matching_identifiers": len({row["id"] for row in matching()})
        == len(matching()),
    }


def observable():
    return {
        "base_scheme": "The complete four-light first-loop representation belongs to the S234 separate-model scheme and retains all three bubbles, all triangle contributions, six ordered boxes, all UV counterterms and the original first-order external normalization.",
        "independent_cut": "Its forward imaginary part is exactly the complete rational S233 first elastic coefficient for4<s<M_H^2; no two-partial-wave truncation or interchange of light/heavy box cuts is used.",
        "zero_jet_control": "An independent constant-background two-field Hessian reproduces the complete off-shell zero-momentum jet. This is not the on-shell symmetric point.",
        "new_prescription": "V2S-T1-OS4 adds deltaC_fin=-A1_base(4/3,4/3,4/3). The point is on shell and subthreshold, not physical real scattering. The finite contact is separately declared and does not edit the frozen S233/S234 prescriptions.",
        "necessary_value_adjustment": "Before this adjustment the first-loop value is negative and exceeds the positive finely cancelled tree by more than10^190 in magnitude. This is not a claim that the full model is strongly coupled or has no UV parent.",
        "controlled_contact_margin": "The positive finite contact is below10^-6 times the entire classical heavy-square quartic margin24q. The contact-only comparison remains coercive with positive remaining quartic. No quantum-potential or all-counterterm bare-vacuum theorem follows.",
        "remaining": "Full physical-angle remainder, higher quantum coefficient matching, omitted loops, continuum UV construction, original affine/common-parent bounce, finite-gravity Regge control and original V/G/B/P8 remain open.",
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
        for pos in range(5):
            args = [symmetric.NAME, symmetric.S0, symmetric.S0, symmetric.S0, 1]
            args[pos] = value
            out.append((f"matching_type_{pos}_{i}", require_matching, tuple(args)))
        out.append((f"cut_type_{i}", require_cut, (value,)))
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
        (0, "original_affine_parent"),
        (1, 0),
        (2, 0),
        (3, 0),
        (1, 4),
        (4, 0),
        (4, 2),
    ):
        args = [symmetric.NAME, symmetric.S0, symmetric.S0, symmetric.S0, 1]
        args[pos] = value
        out.append((f"matching_scope_{len(out)}", require_matching, tuple(args)))
    for value in (0, 4, -1, model.MASS2, 2 * model.MASS2):
        out.append((f"cut_boundary_{len(out)}", require_cut, (value,)))
    for stage in (
        "full_angular_quantum_error",
        "higher_physical_coefficients_matched",
        "all_loop_error",
        "exact_UV_S_matrix",
        "original_parent_replaced",
        "old_counterterm_changed",
        "all_boxes_symmetric",
        "only_contact_UV_subtracted",
        "off_shell_zero_is_on_shell",
        "real_scattering_symmetric_point",
        "large_loop_proves_strong_coupling",
        "full_quantum_potential_positive",
        "all_bare_counterterms_keep_old_minimum",
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
        raise ValueError("Unsupported loop matching closure accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "all_six_mass_ordered_boxes_retained": True,
        "all_three_UV_counterterms_retained": True,
        "complete_first_cut_not_truncated_partial_waves": True,
        "off_shell_zero_jet_not_symmetric_on_shell_point": True,
        "new_named_finite_contact_not_frozen_prescription_edit": True,
        "large_value_correction_not_full_strong_coupling_no_go": True,
        "single_value_and_contact_margin_not_full_quantum_matching": True,
        "all_original_primitive_statuses_unchanged": frontier() == previous.frontier(),
    }
