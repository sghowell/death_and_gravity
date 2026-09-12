"""Read-only actual spatial matching difference and fixed curved pole identity."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_curved_linear_conversion import verify as parent
from p8_vacuum_affine_proca_gaussian import verify as gaussian

from . import audit, contact, density, jets, matching

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates" / "polynomial-vacuum-affine-spatial-matching-difference.json"
)
PARENT_SHA = "7aa73d597f35e68b485e1aed2f0e6ed8c3a72648d992f1f555a8d5fed02d3612"
GAUSSIAN_SHA = "1175702d6980db090396863b807af0e767241daeb7d8dabf138ed11c67e1ae75"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_spatial_matching_difference/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {key: value for key, value in d.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    for previous, expected in ((parent, PARENT_SHA), (gaussian, GAUSSIAN_SHA)):
        if sha(previous.REPORT) != expected:
            raise ValueError(
                "The source-pinned actual affine Gaussian state chain changed"
            )
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_210_fully_rebuilt": PARENT_SHA,
        "S6_176_actual_CD_Gaussian_sector_transitively_rebuilt": GAUSSIAN_SHA,
        "fixed_stress_prescription_and_full_source_retained": True,
        "all_frozen_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("An actual spatial matching difference gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.211.ACTUAL_SPATIAL_UV_DIFFERENCE_COMPLETE_CONTACT_CANCELLATION_AND_FIXED_CURVED_POLE_IDENTITY",
        "date": "2026-09-11",
        "status": "ACTUAL_SPATIAL_UV_DIFFERENCE_COMPLETE_CONTACT_CANCELLATION_AND_FIXED_CURVED_POLE_IDENTITY; NOT_FINITE_DIMENSIONAL_MATCHING_FULL_RESPONSE_INVERSE_BACKGROUND_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/jets.md",
            "notes/density.md",
            "notes/contact.md",
            "notes/matching.md",
            "notes/norms.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "actual_complete_WKB_spatial_difference_and_contact_cancellation": serialize(
            {
                "jets": payload(jets.data()),
                "density": payload(density.data()),
                "contact": payload(contact.data()),
            }
        ),
        "fixed_spatial_pole_identity_finite_matching_boundary_and_norms": serialize(
            {"matching": payload(matching.data())}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The full actual unit-W8 spatial UV symbol is now differenced relative to its unchanged homogeneous anchor. Seven original physical-pair contractions and all105 angular channel/source/degree entries are computed with detector time fixed. Exactly11 spatial differences are nonzero. The full original one-leg contact cancels before internal integration only in the P-minus-P0 difference; its generally nonzero homogeneous anchor remains. The actual curved logarithmic symbol agrees exactly with the existing fixed Proca pole Hessian, including its Einstein sign, longitudinal curvature and first/second source-time derivatives. The original one-ball power coefficient is noncovariant and the exact finite lower-band polynomial is retained. Complete original pair-band conversion and its uniform tail remain unchanged. Explicit power/log coefficient norms are1/200 and1e4 in Z24, with log canonical display4e-796, not cutoff-uniform bounds. Independent full four-order W8 coefficients, all source jets, literal constrained ten-field/nine-pair azimuths, complete angular integration at two Cauchy resolutions, noncommuting contacts and pole/finite-weight controls check the result. Physical dimension3 pole agreement does not fix the evanescent finite dimensional match or supply a full response, reduced mixed inverse, nonlinear background/stability, remaining parent loops, cutoff, finite-gravity IR/Regge or original V/G/B and P8 closure.",
        "not_established": [
            "Original finite spatial dimensional matching with all evanescent Hamiltonian, polarization and angular contributions",
            "Permission to change the reference/state, original cutoff or finite coefficients, delete homogeneous contacts or odd endpoints, or add the fixed local target twice",
            "Full matched response, reduced mixed inverse, finite-coupling background/stability or cutoff-uniform smallness from coefficient normalization",
            "Remaining parent loops/cutoff/scattering/IR/Regge or original V/G/B and P8 closure",
        ],
        "verification_boundary": "Exact full actual WKB UV jets and source differentiation, all physical angular coefficients and general tracefree reconstruction, original contact cancellation only in the spatial difference, exact radial lower endpoint and full proper-clock fixed pole identity. Independent all-four-W8 finite coefficient extraction, full ten-field/nine-pair azimuth checks, all angular coefficients and two inverse-radius resolutions, noncommuting contact and norm/Green controls supplement the written proof. Native, direct science, ordinary and CLI retain original SymPy; only full regression uses the audited exact GCD adapter. Continuous Fourier, locality, dimensional-boundary and regulated-current arguments are not FORMALIZED.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The actual spatial matching difference report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.211.ACTUAL_SPATIAL_UV_DIFFERENCE_COMPLETE_CONTACT_CANCELLATION_AND_FIXED_CURVED_POLE_IDENTITY replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
