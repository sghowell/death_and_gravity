"""Read-only actual global seven-mode physical quadratic state certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as certificate
from p8_proca_physical_matter import verify as parent

from . import audit, bridge, covariance, deformation, model, projectors, state

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "actual-global-seven-mode-state.json"
PARENT_SHA = "0ee647fed9d5710696ad3015c8894dcd3d0e9e735a224fda79e43402b6976b48"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_proca_seven_modes/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen actual physical matter certificate changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_102_fully_rebuilt": PARENT_SHA,
        "same_original_global_retuned_constant_Proca_action": True,
        "same_chosen_preexisting_scalar_state_factor": True,
        "no_old_tadpole_or_subtraction_reference_reassigned": True,
    }


def payload(data):
    return {key: value for key, value in data.items() if key != "checks"}


@cache
def build_report():
    prior = prior_checks()
    identities = audit.residuals()
    gates = {key: bool(value) for key, value in audit.gates().items()}
    if not all(value is True for value in gates.values()):
        raise ValueError("An actual joint seven-mode state proof gate failed")
    exact = certificate.certify_residuals(identities)
    return {
        "schema": 1,
        "claim": "P8-S6.103.ACTUAL_GLOBAL_SEVEN_MODE_STATE",
        "date": "2026-09-09",
        "status": "COMPLETE_POSITIVE_GLOBAL_QUASIFREE_STATE_FOR_ALL_SEVEN_ACTUAL_PHYSICAL_QUADRATIC_MODES_WITH_UNCHANGED_SCALAR_MARGINAL_CORRECTED_PROCA_PROPAGATION_FULL_DENSITY_CCR_AND_TWO_SIGNED_CONES; RENORMALIZED_WARD_COMPLETED_SOURCE_BACKREACTION_V_G_B_AND_ORIGINAL_P8_OPEN",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/model.md",
            "notes/covariance.md",
            "notes/preparation.md",
            "notes/infrared.md",
            "notes/hadamard.md",
            "notes/canonical.md",
            "notes/scope.md",
        ],
        "exact_residuals": exact,
        "named_exact_check_count": len(identities),
        "checked_scalar_entries": sum(
            v.rows * v.cols if isinstance(v, sp.MatrixBase) else 1
            for v in identities.values()
        ),
        "proof_checks": gates,
        "actual_fourteen_phase_hamiltonian": serialize(
            {**payload(model.data()), **payload(model.energy())}
        ),
        "actual_cartesian_polarization_projectors": serialize(
            payload(projectors.data())
        ),
        "positive_reference_covariances": serialize(payload(covariance.data())),
        "smooth_reference_preparation": serialize(
            {**payload(deformation.data()), "whole_compact_strip": deformation.strip()}
        ),
        "actual_global_product_state": serialize(payload(state.data())),
        "canonical_physical_stress_bridge": serialize(payload(bridge.data())),
        "controls": serialize(audit.controls()),
        "verdict": "The unchanged scalar state and explicitly prepared two-tensor/three-Proca factors define a positive global quasifree state for the actual seven-mode quadratic action on its original global clock. All fourteen density commutators, literal action normalizations, the constrained positive three-column Proca covariance, Cartesian TT projector and exact canonical physical-observable bridge are checked. Whole-strip energy and infrared bounds establish distributional existence; corrected Proca Cauchy propagation and a blockwise wavefront argument retain the clock and physical matter signed cones. This fills the free tensor/vector state gap but not renormalized stress, second-order mean response, Ward completion, backreaction or V/G/B UV matching. Original P8 remains open.",
        "not_established": [
            "A complete covariant renormalized physical stress, completed quantum Ward identity, second-order mean solution or self-consistent semiclassical background",
            "A physical-matter QEI transferred from the positive test energy, positivity of a stress second variation, or any arbitrary-amplitude perturbative statement",
            "A nonlinear interacting quantum gravity or covariant BRST state, or spacelike locality of raw spatial TT gauge components",
            "Interaction loop, omitted-operator, heavy-threshold or finite-gravity contour error control, actual V/G/B matching, or a UV completion",
            "Original P8(b) classification or original P8 closure",
        ],
        "verification_boundary": "Native exact actual tensor/Proca Hamiltonians and temporal constraint; complete scalar-plus-spectator fourteen-phase CCR; positive covariances and physical Proca polarization sum on the exact mass shell; global Cartesian projectors and minimal-scalar tensor normalization; literal all-time old-action coefficient and canonical density bridges; explicit smooth preparation and exact analytic-bound inputs. The corrected external Hadamard theorem and written global, energy, infrared and microlocal arguments are pinned but not proof-assistant formalized. Every ancestor, own source/proof/test byte and report field is replayed with native scientific SymPy.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The actual global seven-mode state report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.103.ACTUAL_GLOBAL_SEVEN_MODE_STATE replay passed; source, V/G/B and original P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
