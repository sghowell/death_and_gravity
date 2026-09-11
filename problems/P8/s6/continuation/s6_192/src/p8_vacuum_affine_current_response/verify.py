"""Read-only actual physical current comparison response and finite-amplitude remainder."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_matrix_response_tail import verify as parent
from p8_vacuum_affine_proca_gaussian import verify as gaussian

from . import adiabatic, audit, low, tail, vertices

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-current-response.json"
PARENT_SHA = "fddfd18474a054b1ceadfe5d35c16355b010ee83487c2ae204ab345453570746"
GAUSSIAN_SHA = "1175702d6980db090396863b807af0e767241daeb7d8dabf138ed11c67e1ae75"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_current_response/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    for previous, expected in ((parent, PARENT_SHA), (gaussian, GAUSSIAN_SHA)):
        if sha(previous.REPORT) != expected:
            raise ValueError(
                "The source-pinned actual affine Gaussian state chain changed"
            )
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_191_fully_rebuilt": PARENT_SHA,
        "S6_176_actual_CD_Gaussian_sector_transitively_rebuilt": GAUSSIAN_SHA,
        "fixed_stress_prescription_and_full_source_retained": True,
        "all_frozen_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("An actual physical current comparison response gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.192.ACTUAL_PHYSICAL_SHEAR_CURRENT_C2_MODE_COMPARISON_AND_FULL_FINITE_AMPLITUDE_REMAINDER",
        "date": "2026-09-11",
        "status": "ACTUAL_UNCHANGED_STATE_PHYSICAL_SHEAR_CURRENT_C2_MODE_COMPARISON; NOT_FIXED_COVARIANT_MATCHING_FULL_RESPONSE_INVERSE_BACKGROUND_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/vertices.md",
            "notes/low.md",
            "notes/parity.md",
            "notes/contour.md",
            "notes/tail.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "physical_metric_vertices_and_exact_finite_band": serialize(
            {"vertices": payload(vertices.data()), "low": payload(low.data())}
        ),
        "ordered_adiabatic_comparison_and_complete_C2_current": serialize(
            {"adiabatic": payload(adiabatic.data()), "tail": payload(tail.data())}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The unchanged actual three-mode state admits a complete physical shear-current C2 fourth-order mode comparison on the admitted homogeneous history family. Full metric vertices, frame and frequency derivatives are retained. A momentum-independent covariance growth estimate controls the finite proof-band complement; ordered matrix parity and complex Cauchy bounds control the infinite reference-minus-fourth-order tail. Combining them with the actual state tail gives full comparison bounds 1e77,1e94,1e111 through derivative order2 and Taylor remainder at most epsilon^2*1e111/2, zero at epsilon0 and strict otherwise. J_ad4 is not yet matched to the original fixed covariant mu=m subtraction and finite contacts. This is not a complete physical self-energy, feedback inverse, finite-coupling background, physical cutoff or original V/G/B and P8 closure.",
        "not_established": [
            "Finite local/contact matching of the fourth-order comparison to the original fixed covariant mu=m prescription",
            "The full spatial/mixed physical response, feedback inverse or finite-coupling quantum background",
            "Full interacting parent state/measure/loops, physical heavy/cutoff/threshold and omitted-order matching",
            "Complete vacuum cuts/contour/truncation, finite-gravity IR/Regge or original V/G/B and P8 closure",
        ],
        "verification_boundary": "Exact ordered vertex derivatives, balanced symplectic and full covariance response equations, reference/current parity, complete finite/infinite integrals and zero-safe Taylor factor. Independent actual matrix-root/metric-vertex, propagator-response and holomorphic coefficient fixtures supplement the written proofs. Native, direct science, ordinary and CLI use original SymPy; only full regression uses the audited exact GCD adapter. Analytic and functional arguments are not FORMALIZED.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The actual physical current comparison report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.192.ACTUAL_PHYSICAL_SHEAR_CURRENT_C2_MODE_COMPARISON_AND_FULL_FINITE_AMPLITUDE_REMAINDER replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
