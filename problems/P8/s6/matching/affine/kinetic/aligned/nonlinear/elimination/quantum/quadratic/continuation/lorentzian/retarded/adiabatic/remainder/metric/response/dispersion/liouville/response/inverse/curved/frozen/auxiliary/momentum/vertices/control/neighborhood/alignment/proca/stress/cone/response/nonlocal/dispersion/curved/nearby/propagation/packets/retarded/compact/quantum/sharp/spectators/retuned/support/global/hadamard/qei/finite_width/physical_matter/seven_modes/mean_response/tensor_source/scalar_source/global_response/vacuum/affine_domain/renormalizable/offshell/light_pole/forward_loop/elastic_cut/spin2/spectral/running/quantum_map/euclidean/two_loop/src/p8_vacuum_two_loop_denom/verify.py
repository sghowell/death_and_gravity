"""Read-only complete bare two-loop topology and compact denominator verifier."""

import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as certificate
from p8_vacuum_euclidean_control import verify as parent

from . import audit, calibration, graphs, positive, subgraphs

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-two-loop-denominators.json"
PARENT_SHA = "a24747e99e5002aed275317aa68da4aa0bf41bd90f590e5ab8c47970c76f0e61"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_two_loop_denom/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen restricted Euclidean-control report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_119_fully_rebuilt": PARENT_SHA,
        "same_canonical_full_Gaussian_heavy_quartic_action": True,
        "same_actual_mass_and_fixed_light_on_shell_scheme": True,
        "known_kernel_iterations_not_full_two_loop_error": True,
        "no_rolling_state_or_finite_gravity_transfer": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows = audit.residuals()
    gates = {k: bool(v) for k, v in audit.gates().items()}
    if not all(v is True for v in gates.values()):
        raise ValueError("A complete two-loop denominator or graph gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.120.POLYNOMIAL_VACUUM_TWO_LOOP_DENOMINATORS",
        "date": "2026-09-09",
        "status": "COMPLETE_BARE_TWO_LOOP_LIGHT_1PI_FOUR_POINT_TOPOLOGY_FULL_HEAVY_REFINEMENTS_EXACT_COMPACT_FORWARD_DISC_INTERIOR_DENOMINATOR_GAP_AND_UV_FOREST_INVENTORY; NOT_INTEGRATED_SUBTRACTED_TWO_LOOP_AMPLITUDE_LSZ_HIGHER_LOOP_REGGE_FINITE_GRAVITY_V_G_B_OR_ORIGINAL_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/graphs.md",
            "notes/denominators.md",
            "notes/positive.md",
            "notes/subgraphs.md",
            "notes/literature.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": len(rows),
        "proof_checks": gates,
        "complete_bare_graph_inventory": serialize(payload(graphs.data())),
        "all_exact_positive_polynomial_witnesses": serialize(positive.summaries()),
        "all_local_UV_subgraphs_and_restricted_forests": serialize(
            {
                "summaries": subgraphs.summary(),
                "bare_families": {
                    kind: subgraphs.data(kind, (0, 0, 0)) for kind in graphs.KINDS
                },
            }
        ),
        "actual_mass_and_compact_disc_enclosures": serialize(
            {
                "general": payload(calibration.data()),
                "actual_bare_family_centers": {
                    kind: calibration.point(kind, (0, 0, 0), (sp.Rational(1, 4),) * 4)
                    for kind in graphs.KINDS
                },
                "actual_fully_refined_tadpole_imaginary_endpoint": calibration.point(
                    "tadpole_insertion", (3, 3, 3), (sp.Rational(1, 7),) * 7, 0, 1
                ),
            }
        ),
        "controls": serialize(audit.controls()),
        "verdict": "All bare two-loop light-1PI four-point skeletons of the actual full-heavy quartic parent are enumerated. Every one of 192 exact heavy refinements has an independently checked first graph polynomial and an exact rational positive-polynomial witness giving Re F at least U/4 on the interior simplex throughout the unit forward disc, with a positive actual-mass lift and crossing-derivative ratio at most four. Divergent local subgraphs and vertex-aware compatible forests are recorded. This eliminates an interior-denominator obstacle on this compact domain, but does not evaluate UV subtractions, finite scheme matching, two-loop LSZ or an integrated forward-coefficient error.",
        "not_established": [
            "Ultraviolet boundary convergence or evaluated forest subtraction integrals",
            "Finite fixed-scheme counterterm coefficients and complete two-loop mass residue LSZ or forward coefficient",
            "All higher-loop errors or an all-energy complex scattering contour",
            "Finite-gravity IR/Regge contour or Delta_grav",
            "Old rolling-state transfer absolute coupled renormalization or a common controlled bounce parent",
            "Exclusion or completion of an entire original P8 row full V/G/B or original P8 closure",
        ],
        "verification_boundary": "Native exhaustive finite graph and refinement enumeration, independent weighted-Laplacian and spanning-tree identities, full two-forest invariants, exact rational positive-polynomial witnesses, UV subgraph/forest enumeration and actual parameter enclosures. The continuous simplex/complex-disc argument and diagram interpretation are explicit source-pinned written proofs, not proof-assistant formalized, nonperturbative or peer reviewed. All ancestors, own sources and report fields rebuild read-only with unmodified scientific SymPy.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete two-loop denominator report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.120.POLYNOMIAL_VACUUM_TWO_LOOP_DENOMINATORS replay passed; full V/G/B and original P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
