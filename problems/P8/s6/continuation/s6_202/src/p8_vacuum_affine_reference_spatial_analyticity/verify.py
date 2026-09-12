"""Read-only joint spatial analyticity and far-endpoint Taylor remainder."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_proca_gaussian import verify as gaussian
from p8_vacuum_affine_reference_state_prefactor import verify as parent

from . import audit, boundary, domain, frame, readouts

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT
    / "certificates"
    / "polynomial-vacuum-affine-reference-spatial-analyticity.json"
)
PARENT_SHA = "b75148fac1694e6a689f955d76cb30fc94734f063f85c4ec952215baa2d35eda"
GAUSSIAN_SHA = "1175702d6980db090396863b807af0e767241daeb7d8dabf138ed11c67e1ae75"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_reference_spatial_analyticity/*.py"))
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
        "S6_201_fully_rebuilt": PARENT_SHA,
        "S6_176_actual_CD_Gaussian_sector_transitively_rebuilt": GAUSSIAN_SHA,
        "fixed_stress_prescription_and_full_source_retained": True,
        "all_frozen_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A joint spatial reference analyticity gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.202.ACTUAL_CURVED_UNIT_W8_JOINT_SPATIAL_ANALYTICITY_FULL_READOUTS_AND_FAR_ENDPOINT_TAYLOR_REMAINDER",
        "date": "2026-09-11",
        "status": "ACTUAL_CURVED_UNIT_W8_JOINT_SPATIAL_ANALYTICITY_FULL_READOUTS_AND_FAR_ENDPOINT_TAYLOR_REMAINDER; NOT_NEAR_REGION_LOCAL_COVARIANT_MATCHING_FULL_INVERSE_BACKGROUND_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/domain.md",
            "notes/frame.md",
            "notes/readouts.md",
            "notes/endpoints.md",
            "notes/far-tail.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "joint_complex_domain_and_complete_analytic_frame": serialize(
            {"domain": payload(domain.data()), "frame": payload(frame.data())}
        ),
        "full_constrained_readouts_and_far_endpoint_Taylor_tail": serialize(
            {"readouts": payload(readouts.data()), "boundary": payload(boundary.data())}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "For the same CD sector and preparation, the unit-W8 reference has a common complex time and external-momentum domain ||p||<=1e-6nu(k) for real |k|>=1000. Both frequency branches, a local complex-orthogonal polarization frame with genuine Euclidean norm bounds, all ten constrained field readouts, all nine physical pairs and the inverse summed phase are controlled. Analytic Schwarz continuation is applied to the complete pre-current coefficients, not to an imaginary-part operation or to the Borel state. The complete first-five-endpoint row from time jets0..4 is below2e27nu; its fourth spatial Taylor remainder is below4e57|P|^5nu^-4 when |P|<=1e-6nu/2. Only this far-region remainder is integrated. It is below1e54||D||L2 X46[Gamma], with the original both-created-mode regulator tail below1e58||D||L2 X46[Gamma]/K. The two canonical factors give4e-746 and4e-742/K. The complementary near-momentum region, full finite and divergent Taylor coefficients, complete reference contact and original fixed covariant spatial matching remain open. No full inverse, interacting background, stability, physical cutoff, remaining parent loops or finite-gravity IR/Regge follows. Original V/G/B and P8 remain open.",
        "not_established": [
            "The complementary near-momentum endpoint region or a full all-momentum Taylor-subtracted reference response",
            "Finite and divergent Taylor coefficients, full reference contact and original fixed covariant spatial matching",
            "Borel-state momentum analyticity, a new exact reference state, a full mixed inverse, interacting background or stability",
            "Full parent/cutoff/scattering/IR/Regge matching or original V/G/B and P8 closure",
        ],
        "verification_boundary": "Exact joint-domain constants, complex frame identities, complete constrained-readout bridges, all endpoint recursion rows, Cauchy powers, radial bounds and canonical factors. Independent actual full W8 complex evaluations, genuine Euclidean frames, both Schwarz signs, all nine stress pairs, nonzero fifth Taylor remainders and removed two-leg integrals supplement written continuum estimates. Native, direct science, ordinary and CLI use original SymPy; only full regression uses the audited exact GCD adapter. Joint holomorphy, Cauchy estimates and the far continuum argument are not FORMALIZED.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The joint spatial reference analyticity report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.202.ACTUAL_CURVED_UNIT_W8_JOINT_SPATIAL_ANALYTICITY_FULL_READOUTS_AND_FAR_ENDPOINT_TAYLOR_REMAINDER replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
