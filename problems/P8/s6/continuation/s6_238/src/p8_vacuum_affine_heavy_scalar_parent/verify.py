"""Read-only separate covariant parent and scoped classical matching report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_heavy_scalar_loop_coefficients import verify as parent
from p8_vacuum_affine_quantum_retuning import verify as original_input
from p8_vacuum_analytic_affine_parent import verify as affine_input

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-heavy-scalar-parent.json"
PARENT_SHA = "6d76badc63d988c9fae045c5a3628b20c7fb476ea9f45fdf98fd52a7a1b0e9fb"
ORIGINAL_SHA = "f754f6d62795d5a1c0ade4c576395f1587fc35f9c9895bc8b7b15b2b2dec2e16"
AFFINE_SHA = "e0cf546e77cf12ff4e0f2e6bef706ff6de1eddc40cacff3b0b151fcb7009a4f4"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_heavy_scalar_parent/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    for prior, expected in (
        (parent, PARENT_SHA),
        (original_input, ORIGINAL_SHA),
        (affine_input, AFFINE_SHA),
    ):
        if sha(prior.REPORT) != expected:
            raise ValueError(
                "A frozen original action, affine construction or separate quantum input changed"
            )
        prior.validate_report(
            json.loads(prior.REPORT.read_text()), prior.build_report()
        )
    return {
        "S6_237_original_frontier_and_separate_first_loop_input_fully_rebuilt": PARENT_SHA,
        "S6_182_literal_fixed_smooth_QG1_profile_and_vacuum_constant_fully_rebuilt": ORIGINAL_SHA,
        "S6_174_generic_full_connection_quotient_source_centering_and_chart_fully_rebuilt": AFFINE_SHA,
        "all_frozen_original_and_separate_model_scientific_bytes_unchanged": True,
        "new_full_candidate_not_automatic_transfer_of_old_quantum_results": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError(
            "A full parent, domain, clock-jet, potential or source-stress gate failed"
        )
    return {
        "schema": 1,
        "claim": "P8-S6.238.SEPARATE_SMOOTH_H8A420_COVARIANT_AFFINE_HEAVY_SCALAR_PARENT_WITH_FULL_VACUUM_TREE_MATCHING_CONNECTED_REGULAR_DOMAIN_COMPLETE_CLOCK_FOUR_JET_BOUND_AND_SCOPED_SOURCE_STRESS_CONTROL",
        "date": "2026-09-12",
        "status": "SEPARATE_SMOOTH_CLASSICAL_PARENT_AND_SCOPED_MATCHING; NOT_FULL_QUANTUM_PARENT_SAME_STATE_BOUNCE_UV_REGGE_ORIGINAL_V_G_B_OR_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/family.md",
            "notes/domain.md",
            "notes/clock.md",
            "notes/heavy.md",
            "notes/scope.md",
            "notes/literature.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "literal_new_parent_vacuum_and_full_clock_coefficient_tube": serialize(
            {
                "family": payload(audit.family.data()),
                "clock": payload(audit.clock.data()),
            }
        ),
        "complete_classical_potential_prescribed_stress_and_counterexample": serialize(
            payload(audit.heavy.data())
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "A separately named smooth CD-REG-AFFINE-ISO-QG1-H8A420 covariant action combines the literal fixed original f,r and QG1 profile with an explicit analytic clock-jet-null correction and a minimally coupled heavy scalar. All dependent Ia coefficients and generic full affine maps are recomputed. The new R keeps the complete original connected range, including its upper-edge margin, so the full quotient, source-centered mass and regular metric chart construct a new algebraic parent. Its decoupled flat degree-four tree is V2S-T1, while its full constant-field potential relative to the retained origin is coercive. Eighth-order independent clock zeros and seventh-order dependent zeros preserve the named classical jets; a directly proved complex neighborhood bounds all mixed four-jet differences of the full normalized F,R,A3,A4,A5,J by10^-2700 on the stated tube. Exact prescribed linear-light trajectories admit a closed heavy particular response. Including the full source metric derivative gives stress error below10^-8 relative to Y only on an explicit connected nonnegative-gradient two-chart corridor. An unrestricted-strip counterexample instead exceeds10^395 despite the small offset proxy. None of these classical/algebraic results transfers the old first-loop amplitudes or QG1 quantum mean background to the new model. Heavy-state stress, new higher-vertex loops, same-state nonlinear bounce, physical UV, finite-gravity Regge and original V/G/B/P8 remain open.",
        "not_established": [
            "A globally real-analytic complete QG1 scalar function; only the new additions and R are analytic",
            "Small heavy stress on the unrestricted field strip, an all-background physical kinetic gap or an interpolating coupled solution",
            "A new specified heavy quantum state and controlled curved stress, counterprofile or same-state nonlinear bounce",
            "Automatic S235-S237 first-loop transfer, a new-parent omitted-loop bound, a quantum vacuum theorem or physical UV and Regge completion",
            "A physical cutoff, all-parent exclusion or original V/G/B/P8 closure",
        ],
        "verification_boundary": "Literal original inputs and full new formulas accompany exact finite-jet, quotient-factor, source and stress identities. Written global real-domain, complex-neighborhood and elementary exponential inequalities give the continuum bounds. These are not FORMALIZED proofs, and numerical fixtures do not replace them. The failed unrestricted stress proposal is replaced by an explicit counterexample and a scoped corridor; no frozen claim is edited. Native/direct/ordinary/CLI retain original SymPy; only full regression uses the audited exact-GCD adapter.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete separate covariant parent report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.238 separate classical parent replay passed; new quantum matching and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
