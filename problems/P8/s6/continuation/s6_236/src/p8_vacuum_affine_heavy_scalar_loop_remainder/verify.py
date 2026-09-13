"""Read-only complete first-loop physical angular matching remainder."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_heavy_scalar_four_point_loop import verify as parent
from p8_vacuum_affine_heavy_scalar_tree_matching import verify as model_input

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates" / "polynomial-vacuum-affine-heavy-scalar-loop-remainder.json"
)
PARENT_SHA = "1277412d19d9b022eb81c3728251d1836d88e4b047c8ad80794107c2aff7c911"
MODEL_SHA = "31716f2b701c0120c328dd8a2506f677eb30edeca956d531c57975e38d8d25cb"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_heavy_scalar_loop_remainder/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    for prior, expected in ((parent, PARENT_SHA), (model_input, MODEL_SHA)):
        if sha(prior.REPORT) != expected:
            raise ValueError(
                "A frozen full loop, OS4 condition or tree matching input changed"
            )
        prior.validate_report(
            json.loads(prior.REPORT.read_text()), prior.build_report()
        )
    return {
        "S6_235_complete_mass_ordered_loop_UV_cut_and_unchanged_OS4_value_condition_fully_rebuilt": PARENT_SHA,
        "S6_233_separate_model_complete_original_tree_and_full_window_bound_fully_rebuilt": MODEL_SHA,
        "original_target_and_all_first_loop_normalization_inputs_transitively_rebuilt": True,
        "all_frozen_original_and_separate_model_scientific_bytes_unchanged": True,
        "no_new_finite_counterterm_or_change_to_the_existing_OS4_prescription": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError(
            "A full branch, convergent remainder, cancellation or physical-window gate failed"
        )
    return {
        "schema": 1,
        "claim": "P8-S6.236.SEPARATE_V2S_T1_OS4_COMPLETE_MASS_ORDERED_FIRST_LOOP_PHYSICAL_ALL_ANGLE_REMAINDER_WITH_CONVERGENT_PARAMETER_BOUNDS_AND_UNCHANGED_SYMMETRIC_SUBTRACTION",
        "date": "2026-09-12",
        "status": "SEPARATE_FULL_FIRST_LOOP_PHYSICAL_ALL_ANGLE_MATCHING; NOT_HIGHER_QUANTUM_COEFFICIENT_JETS_OMITTED_LOOPS_EXACT_UV_ORIGINAL_AFFINE_PARENT_BOUNCE_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/integral.md",
            "notes/expansion.md",
            "notes/remainder.md",
            "notes/window.md",
            "notes/scope.md",
            "notes/literature.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "complete_parameter_primitive_and_convergent_remainder": serialize(
            payload(audit.integral.data())
        ),
        "full_on_shell_cancellation_and_actual_physical_matching_bound": serialize(
            {
                "expansion": payload(audit.expansion.data()),
                "bounds": payload(audit.bounds.data()),
            }
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The complete S235 mass-ordered first four-light loop is reduced to an exact parameter primitive and its heavy-mass derivative with the physical Log(L-i0) branch retained. Analytic coefficient bounds on a full complex inverse-mass disk control every term beyond the first two orders; the light threshold logarithm has a uniform integrable majorant. Only after combining all bubbles, triangles and six ordered boxes do all B and L Log L terms cancel from the orders n^-2 and n^-3. The remaining terms depend on n but not on the on-shell kinematics, so the already frozen V2S-T1-OS4 symmetric value subtraction removes both without a new finite condition. For every physical angle and4<=s<=10^196, the entire first-loop coefficient satisfies |A1_OS4|/A_original<10^-199. Thus the explicitly tree-plus-first-loop-truncated expression differs from the full original tree by less than1/60+10^-199<1/59 on the same energy window2<=E<=10^98. This is not an error bound on the omitted loops, an exact quantum S matrix, higher physical coefficient jets, a continuum UV theory, a controlled heavy resonance, finite-gravity Regge control or the original common-parent bounce. Original V/G/B/P8 remain open.",
        "not_established": [
            "Controlled higher physical coefficient derivatives b20,b21,b40 from this real-window bound alone",
            "An omitted-loop remainder, exact quantum amplitude, controlled resonance, high-energy arc or full UV construction",
            "A full quantum effective-potential or all-counterterm physical-vacuum theorem",
            "Identification with the original affine/DHOST/Proca parent, connected covariant bounce/state matching, finite-gravity Regge control or nonlinear stability",
            "A physical cutoff, exclusion of all UV parents or original V/G/B/P8 closure",
        ],
        "verification_boundary": "Exact complete primitive, derivative, routed loop cancellation and rational parameter checks accompany written complex-disk Cauchy and full real-window logarithmic bounds. These arguments are not FORMALIZED. Independent quadratures and boundary samples are diagnostics, not continuum certificates. All first-loop diagrams and the same finite contact are retained; tree-plus-first-loop is explicitly a truncation. Native/direct/ordinary/CLI retain original SymPy; only full regression uses the audited exact-GCD adapter.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The full first-loop angular remainder report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.236 full first-loop angular matching replay passed; omitted loops and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
