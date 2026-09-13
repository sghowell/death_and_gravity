"""Read-only first-loop input for the separate V2S-T1 model, not original P8 closure."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_heavy_scalar_tree_matching import verify as parent
from p8_vacuum_affine_scalar_matching_scale import verify as normalization

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-heavy-scalar-one-loop.json"
PARENT_SHA = "31716f2b701c0120c328dd8a2506f677eb30edeca956d531c57975e38d8d25cb"
NORMALIZATION_SHA = "34da90ca48a00689b3b9d7db12b367db703c970236caabcecd293b1061559cd5"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_heavy_scalar_one_loop/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    for prior, expected in ((parent, PARENT_SHA), (normalization, NORMALIZATION_SHA)):
        if sha(prior.REPORT) != expected:
            raise ValueError(
                "A frozen separate model, original normalization or physical frontier changed"
            )
        prior.validate_report(
            json.loads(prior.REPORT.read_text()), prior.build_report()
        )
    return {
        "S6_233_separate_V2S_T1_full_classical_model_and_original_tree_matching_fully_rebuilt": PARENT_SHA,
        "S6_231_original_identical_normalization_and_full_error_boundary_fully_rebuilt": NORMALIZATION_SHA,
        "original_S6_177_S6_182_target_and_S6_232_moment_boundary_transitively_rebuilt": True,
        "all_frozen_original_and_separate_model_scientific_bytes_unchanged": True,
        "first_loop_prescription_explicitly_separate_not_original_affine_counterterm_change": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError(
            "A Gaussian, first-loop self-energy, analytic disk or decay gate failed"
        )
    return {
        "schema": 1,
        "claim": "P8-S6.234.SEPARATE_V2S_T1_EXACT_HEAVY_GAUSSIAN_INTEGRATION_COMPLETE_FIRST_LIGHT_SELF_ENERGY_COMPLEX_DISK_BOUND_AND_LEADING_HEAVY_DECAY_COEFFICIENT",
        "date": "2026-09-12",
        "status": "SEPARATE_FINITE_REGULATOR_AND_FIRST_LOOP_INPUT; NOT_ALL_LOOP_FOUR_POINT_MATCHING_ORIGINAL_AFFINE_PARENT_BOUNCE_UV_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/gaussian.md",
            "notes/loops.md",
            "notes/analytic.md",
            "notes/width.md",
            "notes/scope.md",
            "notes/literature.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "exact_heavy_integration_and_complete_one_light_loop_generator": serialize(
            payload(audit.gaussian.data())
        ),
        "complete_first_light_two_point_and_heavy_absorptive_coefficients": serialize(
            {
                "light": payload(audit.self_energy.data()),
                "heavy": payload(audit.width.data()),
            }
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The separate frozen V2S-T1 model permits exact finite-regulator Euclidean heavy Gaussian integration. The remaining nonlocal light action has a positive quartic lower bound, and its complete Hessian retains both local and bilocal terms. The full one-light-loop trace generates mixed heavy-light contributions despite the light-independent flat heavy determinant. An explicit heavy-one-point-zero condition cancels its associated source tadpole, while the full first light self-energy retains the quartic tadpole and mixed bubble. In the stated mu1 MSbar plus light on-shell scheme, abs(Pi_MS(1))<10^-7 and0<Pi_MS'(1)<10^-207. The subtracted first coefficient has a complete complex logarithm remainder and gives a propagator ratio error below10^-209 for the reciprocal of the one-loop-truncated inverse on abs(s-1)<=10^196. That truncated inverse has only its renormalization-imposed simple pole with unit residue in this disk, not a certified exact quantum LSZ pole. The mixed first bubble retains its physical threshold and does not acquire a pseudothreshold first-sheet cut. Independently normalized heavy and mixed cuts give the first heavy decay coefficient g^2 beta/(32pi M_H), with10^-208<Gamma_first/M_H<10^-207 and the outgoing formal second-sheet sign. No stable exact heavy atom, controlled full resonance, evaluated real four-point matching, all-loop remainder, continuum quantum construction, original affine/common-parent bounce or original V/G/B/P8 closure is established.",
        "not_established": [
            "A continuum interacting quantum UV measure, exact S matrix, high-energy arc or all-loop error bound",
            "Exact physical light mass/LSZ, complete quantum cut thresholds, or a controlled heavy resonance pole and width",
            "Evaluated complete renormalized real four-point matching or the full quantum angular error premise of S231/S232",
            "Identification with the original affine/DHOST/Proca action, a covariant common-parent bounce/state domain or finite-gravity Regge control",
            "A physical cutoff, full UV exclusion or original V/G/B/P8 closure",
        ],
        "verification_boundary": "Exact finite matrix Schur/Hessian identities, independently counted Wick factors, Feynman parameter and rational margin checks accompany written all-domain quadratic-form and complex-log proofs. These written arguments are not FORMALIZED. Finite-regulator integration, formal loop coefficients and algebraic inversion of a truncated inverse are distinguished from the unproved full quantum theory. Native/direct/ordinary/CLI retain original SymPy; only full regression uses the audited exact-GCD adapter.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The separate first-loop report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.234 separate first-loop input replay passed; full quantum matching and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
