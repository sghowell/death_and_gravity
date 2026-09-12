"""Read-only complete reference spatial and prepared Gaussian metric response."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_homogeneous_trace_anchor import verify as parent
from p8_vacuum_affine_proca_gaussian import verify as gaussian

from . import assembly, audit, dimension, remainder, shapes, vertices

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-full-spatial-remainder.json"
PARENT_SHA = "0c12575dca9e1ba1361343e7832e93e377e6198d7923ec49fe7d7f51e54b1153"
GAUSSIAN_SHA = "1175702d6980db090396863b807af0e767241daeb7d8dabf138ed11c67e1ae75"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_full_spatial_remainder/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    for previous, expected in ((parent, PARENT_SHA), (gaussian, GAUSSIAN_SHA)):
        if sha(previous.REPORT) != expected:
            raise ValueError("A source-pinned historical Gaussian input changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_218_full_homogeneous_trace_anchor_fully_rebuilt": PARENT_SHA,
        "S6_176_original_Gaussian_state_transitively_rebuilt": GAUSSIAN_SHA,
        "same_actual_reference_parent_Hessian_not_full_nonlinear_source_free_family": True,
        "no_frozen_scientific_test_or_report_bytes_changed": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A full reference Gaussian metric response gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.219.FULL_SPATIAL_SCALAR_REMAINDERS_ORIGINAL_REGULATOR_AND_COMPLETE_PREPARED_GAUSSIAN_METRIC_RESPONSE",
        "date": "2026-09-12",
        "status": "COMPLETE_REFERENCE_GAUSSIAN_METRIC_RESPONSE_AND_WEAK_BOUND; NOT_REDUCED_INVERSE_FULL_SOURCED_PARENT_BACKGROUND_STABILITY_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/vertices.md",
            "notes/remainder.md",
            "notes/shapes.md",
            "notes/dimension.md",
            "notes/assembly.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "full_constraint_vertices_pairs_and_original_scalar_shapes": serialize(
            {
                "features": payload(vertices.feature_data()),
                "pairs": payload(vertices.pair_data()),
                "shapes": payload(shapes.data()),
            }
        ),
        "full_remainders_dimension_current_and_prepared_response": serialize(
            {
                "remainder": payload(remainder.data()),
                "dimension": payload(dimension.data()),
                "spatial": payload(assembly.data()),
                "prepared": payload(assembly.ward_data()),
            }
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The full spatial trace-shift and temporal-constraint feature is exactly minus the independent constrained ADM first vertex, and the complete noncommuting second contact has the correct positive Hamiltonian signs. Their Frobenius/operator norm contractions allow all ten-field actual-state, occupation, time-IBP and UV-subtracted endpoint bounds and original one-/two-leg tails to transfer without enlargement. The full nine-pair product retains ten geometries, all ordered longitudinal cross products and the constraint square. Pair exchange and corrected Schwarz parity concern centered total grades, not endpoint labels. Direct complete scalar endpoint extraction and independent centered integration give the exact three ordered cubic/linear original-cutoff shapes; all quadratic/finite shapes cancel while the q4 log shell remains. The scalar odd nonlog coefficients are nonzero and are not treated as the tracefree special case. Explicit96 complex-dimensional geometry/inverse-basis majorants,35 grouped transverse moments and all-four-order phase/frequency controls give the full scalar dominated dimension limit at fixed physical source invariants. The actual full homogeneous trace/tracefree anchor and original six-invariant finite prescription assemble the complete spatial reference Gaussian current with bound2e95 M Z136 and homogeneous-anchored original-regulator error3e54 M Z136/K. The normalized trace/trace kernel has bound2e95 and both independently ordered cross kernels have bounds3e49; this does not equate the cross kernels. The original prepared Ward maps and complete density/chart corrections therefore supply the full reference Gaussian metric response, bounded by1e118 V04 U138. Independent literal metric determinant/field-strength variations, nondegenerate six-phase-space covariance comparisons, all-source-jet original-band integrals including the log shell, grazing-strip and transpose controls supplement the proof. This establishes a reference Gaussian weak operator, not the genuinely reduced inverse, full finite-amplitude sourced-parent remainder, corrected quantum background/stability, physical cutoff or original V/G/B/P8 closure.",
        "not_established": [
            "A genuinely constrained/reduced scalar-mixed same-space inverse or contraction",
            "A finite-amplitude remainder for the full nonlinear affine sourced parent",
            "A quantum-corrected background, bounce or stability from small canonical weak bounds",
            "Physical EFT cutoff, coupled-loop/heavy-sector and original V/G/B obligations",
            "Original P8 closure or exclusion of an original primitive by this conditional reference result",
        ],
        "verification_boundary": "Exact full ADM/physical-feature bridges, complete pairs, scalar shapes, safe dimension denominators, norm arithmetic, original-regulator assembly and prepared response are checked. Independent literal metric/covariance/original-band tests do not substitute for continuous proofs. Infinite-mode, angular, dimensional, Sobolev and Ward arguments are written proofs, not FORMALIZED. Native/direct/ordinary/CLI use original SymPy; only full regression uses the audited exact-GCD adapter. Frozen predecessors are unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The full reference Gaussian metric report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.219 full reference Gaussian metric response replay passed; reduced inverse and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
