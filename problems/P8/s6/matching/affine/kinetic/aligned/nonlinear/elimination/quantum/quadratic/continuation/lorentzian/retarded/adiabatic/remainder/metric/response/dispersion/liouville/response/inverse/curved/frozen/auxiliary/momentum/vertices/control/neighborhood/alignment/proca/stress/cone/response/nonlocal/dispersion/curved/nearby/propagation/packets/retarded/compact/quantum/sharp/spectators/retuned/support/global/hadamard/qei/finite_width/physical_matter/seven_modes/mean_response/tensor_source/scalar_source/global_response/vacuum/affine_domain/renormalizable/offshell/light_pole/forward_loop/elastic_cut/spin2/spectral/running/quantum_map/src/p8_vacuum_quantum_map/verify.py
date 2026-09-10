"""Read-only finite-order quantum map of the full matched polynomial vacuum."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_polynomial_vacuum_running import verify as parent

from . import audit, calibration, interactions, ward, wick

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-one-loop-quantum-field-map.json"
PARENT_SHA = "5d8e48a98de7013a1cd9c342ea5673cbd2a47269df4a6e678144e6be22cf428b"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_quantum_map/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen full-model reference-flow report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_117_fully_rebuilt": PARENT_SHA,
        "same_actual_S711_cubic_derivative_field_map": True,
        "same_parent_one_loop_mass_residue_counterterms_and_elastic_amplitude": True,
        "full_transformed_action_not_bare_quartic_derivative_truncation": True,
        "new_explicit_composite_prescription_no_old_rolling_transfer": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows = audit.residuals()
    gates = {k: bool(v) for k, v in audit.gates().items()}
    if not all(v is True for v in gates.values()):
        raise ValueError("An actual finite-order quantum-map gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.118.POLYNOMIAL_VACUUM_ONE_LOOP_QUANTUM_FIELD_MAP",
        "date": "2026-09-09",
        "status": "COMPLETE_FULL_MATCHED_ACTION_ONE_LOOP_SOURCE_AWARE_FIELD_MAP_ACTUAL_DERIVATIVE_WICK_MIXING_EVANESCENT_SUBTRACTION_POSITIVE_LOCAL_LIGHT_RESIDUE_AND_PROPER_LSZ_TRANSFER; NOT_BARE_TRUNCATED_TARGET_QUANTUM_EQUIVALENCE_GLOBAL_INVERSE_HIGHER_LOOP_FINITE_GRAVITY_ROLLING_PARENT_V_G_B_OR_ORIGINAL_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/interactions.md",
            "notes/wick.md",
            "notes/ward.md",
            "notes/bounds.md",
            "notes/literature.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": len(rows),
        "proof_checks": gates,
        "actual_derivative_map_Wick_mixing": serialize(payload(wick.data())),
        "all_generated_interactions_and_counterterms": serialize(
            payload(interactions.data())
        ),
        "source_identity_controls_and_LSZ": serialize(payload(ward.data())),
        "actual_composite_and_propagator_bounds": serialize(
            {
                "whole_disc": payload(calibration.data()),
                "origin": calibration.point(0),
                "light_mass_shell": calibration.point(1),
                "upper_real_disc_endpoint": calibration.point(2),
            }
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The actual cubic derivative map has an explicitly computed one-loop composite normalization, including its finite evanescent subtraction term. Full generated interactions, transformed fixed counterterms, regulated Jacobian and source identity give proper stable-light LSZ equivalence to the complete matched action. Its one-loop residue remains positive and its truncated inverse has no additional zeros in the parent unit disc. The parent one-loop b2 and elastic-cut bounds transfer to this full action, not to the bare derivative truncation. This is a formal finite-order flat-vacuum result, not a global or cosmological quantum-equivalence or original P8 closure theorem.",
        "not_established": [
            "Quantum equivalence of the bare earlier quartic derivative truncation without matched higher interactions",
            "Global differential invertibility, a nonperturbative functional measure or an all-order propagator spectrum",
            "Higher-loop scattering or field-map remainders and an all-energy UV completion",
            "Finite-gravity covariant field map, stress/curvature counterterms or IR/Regge contour and Delta_grav",
            "Old rolling-state transfer, absolute coupled renormalization or a common controlled propagating bounce parent",
            "Exclusion or completion of an entire original P8 row, full V/G/B or original P8 closure",
        ],
        "verification_boundary": "Native exact contraction of the actual four-dimensional 210-jet cubic map, covariant block assembly, Gamma-function pole and finite expansion, complete generated interaction degrees, topology, Gaussian source/Jacobian/sextic controls, off-shell counterexample, LSZ factors and actual rational local bounds. Continuous functional change-of-variable and analytic estimates are explicit source-pinned written perturbative proofs, not proof-assistant formalized, nonperturbative or peer reviewed. Every ancestor, own source and report field is read-only rebuilt with unmodified scientific SymPy.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The full matched vacuum one-loop quantum-map report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.118.POLYNOMIAL_VACUUM_ONE_LOOP_QUANTUM_FIELD_MAP replay passed; full V/G/B and original P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
