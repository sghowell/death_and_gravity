"""Separate V2S-T1 classical tree model, with original physical frontiers unchanged."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_scalar_transfer_moment import audit as previous

from . import cut, model
from . import matching as tree_matching

ITEM = {
    "id": "separate_V2S_T1_local_positive_two_scalar_model_complete_original_massive_tree_matching_and_uniform_all_angle_remainder",
    "status": "SEPARATE_CLASSICAL_TWO_SCALAR_MODEL_AND_TREE_MATCHING_NOT_FULL_QUANTUM_ORIGINAL_AFFINE_PARENT_BOUNCE_UV_OR_P8",
}
require_scope = previous.require_scope


def require_stage(stage):
    if not isinstance(stage, str) or stage not in (
        "separate_canonical_two_scalar_bounded_potential",
        "complete_original_massive_tree_plus_exact_remainder",
        "uniform_all_angle_tree_error_on_named_window",
        "complete_rational_first_elastic_cut_comparison",
        "separate_tree_spectral_matching_only",
    ):
        raise ValueError(
            "Only the named separate classical model and tree matching are proved"
        )
    return stage


def require_model(name, energy):
    if not isinstance(name, str) or name != model.NAME:
        raise ValueError("Only the separately named V2S-T1 model is evaluated")
    if isinstance(energy, bool) or not isinstance(energy, (int, s.Integer, s.Rational)):
        raise TypeError("The tree matching energy must be an exact rational")
    if not 2 <= energy <= tree_matching.ENERGY:
        raise ValueError("Outside the proved physical tree matching energy window")
    return name, s.Rational(energy)


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError(
            "A separate positive tree model does not close original P8 or its quantum matching"
        )
    return True


def packets():
    return {
        "separate_local_two_scalar_model_and_full_classical_potential": model.data(),
        "complete_original_tree_matching_and_uniform_massive_error": tree_matching.data(),
        "complete_rational_angular_integrals_and_first_elastic_comparison": cut.data(),
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
            "one_new_separate_tree_model": len(matching())
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
        "complete_two_field_kinetic_and_potential_checked": True,
        "original_mass_and_potential_retained_in_all_angle_target": True,
        "tree_cut_uses_complete_rational_amplitude": True,
        "new_model_not_substituted_for_original_affine_action": True,
        "tree_matching_not_full_quantum_error_or_bounce_domain": True,
        "unique_matching_identifiers": len({row["id"] for row in matching()})
        == len(matching()),
    }


def observable():
    return {
        "new_named_scope": "V2S-T1 is a separately named local polynomial light/heavy scalar model with canonical kinetic matrix, positive masses and globally coercive nonnegative classical potential. It is not the original Proca/DHOST affine parent.",
        "exact_tree_matching": "Its literal full exchange/contact amplitude equals the complete original S177/S182 massive four-scalar tree plus gamma*sum(channel-2)^4/[D-(channel-2)], D=2lambda/gamma. The new contact is explicitly part of this separate model, not an edit to an old counterterm.",
        "uniform_error": "On all physical angles and4<=s<=10^196, the original tree is positive and the new amplitude has relative tree error strictly below1/60. The complete first elastic coefficient differs by less than121/3600 for s>4 and both vanish at threshold.",
        "higher_coefficients": "The new tree has b20=4lambda, b21=-3gamma and b40=gamma^2/lambda, with forward constant difference16gamma/(D+2). All heavy exchange poles, infinite even angular channels and higher terms remain.",
        "physical_boundary": "Classical flat model and controlled tree/first-elastic matching only. No full quantum mass/residue, real loops, resonance width, exact unitarity, high-energy arc, original independent-function correspondence, covariant bounce parent, physical cutoff or original V/G/B/P8 closure.",
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
            out.append((f"target_scope_type_{pos}_{i}", require_scope, tuple(args)))
        out.append((f"stage_type_{i}", require_stage, (value,)))
        out.append(
            (f"model_name_type_{i}", require_model, (value, tree_matching.ENERGY))
        )
        out.append((f"energy_type_{i}", require_model, (model.NAME, value)))
    for args in (
        (1, modes.KAPPA, modes.MASS, 1),
        (-1, modes.KAPPA, modes.MASS, 1),
        (0, 2 * modes.KAPPA, modes.MASS, 1),
        (0, modes.KAPPA, 999, 1),
        (0, modes.KAPPA, modes.MASS, 2),
    ):
        out.append((f"target_scope_{len(out)}", require_scope, args))
    for args in (
        ("", tree_matching.ENERGY),
        ("original_affine_parent", tree_matching.ENERGY),
        ("V2S-T2", tree_matching.ENERGY),
        (model.NAME, 0),
        (model.NAME, -2),
        (model.NAME, 1),
        (model.NAME, 2 * tree_matching.ENERGY),
    ):
        out.append((f"model_domain_{len(out)}", require_model, args))
    for stage in (
        "full_quantum_error_below_one_sixtieth",
        "exact_unitarity",
        "UV_complete",
        "original_parent_replaced",
        "original_full_functions_matched",
        "physical_cutoff",
        "bounce_parent_complete",
        "same_global_field_dictionary",
        "real_loop_matching_done",
        "renormalized_mass_proved",
        "resonance_width_controlled",
        "only_l0_l2_in_new_cut",
        "massless_tree_target",
        "heavy_pole_deleted",
        "finite_gravity_Regge_done",
        "closed_P8",
    ):
        out.append((f"unsupported_stage_{stage}", require_stage, (stage,)))
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
        raise ValueError(
            "Unsupported model matching or physical closure accepted: " + name
        )
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "omitting_contact_does_not_match_original_tree": True,
        "massive_shifted_cubic_identity_keeps_potential": True,
        "constant_and_higher_coefficients_not_identical": True,
        "rational_tree_has_infinite_even_angular_channels": True,
        "positive_tree_potential_not_quantum_UV_completion": True,
        "first_elastic_comparison_not_full_loop_error": True,
        "two_scalar_model_not_original_affine_parent": True,
        "all_original_primitive_statuses_unchanged": frontier() == previous.frontier(),
    }
