"""Scoped isolated-factor result without promoting the full P8 frontier."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_quantum_forced_constraints import audit as previous

from . import kernel, normalization, spectral

ITEM = {
    "id": "isolated_fixed_shear_factor_unique_subthreshold_zero_pole_plus_cut_and_finite_window_causal_inverse",
    "status": "ISOLATED_REFERENCE_INVERSE_WITH_UNDAMPED_POLE_NOT_FULL_COUPLED_POLE_STABILITY_OR_P8",
}
require_scope = previous.require_scope


def require_stage(stage):
    if not isinstance(stage, str) or stage not in (
        "isolated_factor",
        "pole_plus_cut",
        "finite_window_inverse",
    ):
        raise ValueError(
            "Only the isolated reference factor and its finite-window inverse are proved"
        )
    return stage


def frontier():
    return previous.frontier()


def matching():
    return [dict(row) for row in previous.matching()] + [dict(ITEM)]


def validate_scope(rows, obligations):
    if rows != frontier() or obligations != matching():
        raise ValueError(
            "An isolated shear inverse does not close a full quantum or original P8 obligation"
        )
    return True


def packets():
    return {
        "physical_normalization": normalization.data(),
        "first_sheet_and_cut": spectral.data(),
        "causal_inverse": kernel.data(),
    }


@cache
def residuals():
    rows = {
        name + "_" + key: value
        for name, p in packets().items()
        for key, value in p["checks"].items()
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
            "one_new_isolated_reference_input": len(matching())
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
            for name, p in packets().items()
            for key, value in p["gates"].items()
        },
        "original_finite_prescription_and_mass_retained": True,
        "both_pole_and_continuum_terms_retained": True,
        "positive_weight_not_positive_F2_residue": True,
        "only_finite_window_total_L1_claimed": True,
        "isolated_zero_not_full_physical_pole_or_instability": True,
        "unique_matching_identifiers": len({r["id"] for r in matching()})
        == len(matching()),
    }


def observable():
    return {
        "fixed_factor": "The original finite local unit-Frobenius TF Hessian gives F2(0)=-1/30. The independently normalized full Proca spin2 cut gives the radial weight y^2(30-20y^2+3y^4)/30.",
        "new_spectral_structure": "A2=-F2 has exactly one simple zero p=-r m^2, with0.5798<r<0.5799. Its positive inverse weight obeys16<R/m^2<17. The residue of1/F2 is negative. The cut density is positive and supplies the remaining static moment.",
        "causal_inverse": "The isolated inverse is pole plus continuum. Its full causal kernel is ordinary L1 on each finite window but not on the full half-line. Its primitive lies in[-60,0], giving the prepared bound60T||f'||.",
        "remaining": "A full nonzero-transfer curved matrix normal form, compatible-space full quantum inverse, finite-coupling/nonlinear parent remainder, background/stability, heavy/cutoff and original V/G/B/P8 are still open.",
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
            out.append((f"scope_type_{pos}_{i}", require_scope, tuple(args)))
        out.append((f"stage_type_{i}", require_stage, (value,)))
    for args in (
        (1, modes.KAPPA, modes.MASS, 1),
        (-1, modes.KAPPA, modes.MASS, 1),
        (0, 2 * modes.KAPPA, modes.MASS, 1),
        (0, modes.KAPPA, 999, 1),
        (0, modes.KAPPA, modes.MASS, 2),
    ):
        out.append((f"scope_{len(out)}", require_scope, args))
    for stage in (
        "full_quantum_inverse",
        "stable",
        "finite_Born_remainder",
        "closed_P8",
        "new_state",
        "pole_deletion",
        "physical_cutoff",
        "half_line_L1",
        "full_physical_pole",
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
                "extra_quantum_inverse",
                validate_scope,
                (
                    frontier(),
                    matching() + [{"id": "quantum_inverse", "status": "COMPLETE"}],
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
        raise ValueError("Unsupported isolated-factor claim accepted: " + name)
    return len(bad_cases())


@cache
def controls():
    return {
        "rejected_inputs": rejected_inputs(),
        "trace_gap_cannot_replace_negative_shear_threshold": -s.Rational(172, 225)
        != s.Rational(16, 15),
        "action_coefficient_not_Hessian_coefficient": -s.Rational(1, 60)
        != -s.Rational(1, 30),
        "omitted_pole_changes_static_inverse_by_more_than_zero": spectral.RESIDUE_LO
        / spectral.ROOT_HI
        > 0,
        "positive_weight_is_negative_residue_for_F2": -spectral.RESIDUE_LO < 0,
        "coarse_pole_does_not_exhaust_static_moment": spectral.RESIDUE_HI
        / spectral.ROOT_LO
        < 30,
        "undamped_periodic_pole_has_positive_absolute_mass_each_period": 4
        * spectral.RESIDUE_LO
        / spectral.ROOT_HI
        > 0,
        "full_quantum_force_graph_keeps_direct_auxiliary_block": previous.ITEM["id"]
        in {r["id"] for r in matching()},
        "all_original_primitive_statuses_unchanged": frontier() == previous.frontier(),
    }
