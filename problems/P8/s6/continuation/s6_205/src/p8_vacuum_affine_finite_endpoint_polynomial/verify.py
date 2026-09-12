"""Read-only finite endpoint Taylor restoration with original regulator tail."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_local_spatial_hessian import verify as parent
from p8_vacuum_affine_proca_gaussian import verify as gaussian

from . import audit, degrees, limit, radials, tail

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates" / "polynomial-vacuum-affine-finite-endpoint-polynomial.json"
)
PARENT_SHA = "2c7edc8da05b175c0df786284f1730989be0860859f08ea98416ca2229345414"
GAUSSIAN_SHA = "1175702d6980db090396863b807af0e767241daeb7d8dabf138ed11c67e1ae75"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_finite_endpoint_polynomial/*.py"))
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
        "S6_204_fully_rebuilt": PARENT_SHA,
        "S6_176_actual_CD_Gaussian_sector_transitively_rebuilt": GAUSSIAN_SHA,
        "fixed_stress_prescription_and_full_source_retained": True,
        "all_frozen_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A finite endpoint Taylor restoration gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.205.TEN_FINITE_CURVED_ENDPOINT_SPATIAL_TAYLOR_CELLS_RESTORED_WITH_ORIGINAL_REGULATOR_TAIL",
        "date": "2026-09-11",
        "status": "TEN_FINITE_CURVED_ENDPOINT_SPATIAL_TAYLOR_CELLS_RESTORED_WITH_ORIGINAL_REGULATOR_TAIL; NOT_REMAINING_UV_CONTACT_MATCHING_FULL_RESPONSE_INVERSE_BACKGROUND_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/degrees.md",
            "notes/moments.md",
            "notes/coefficients.md",
            "notes/tail.md",
            "notes/restoration.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "complete_endpoint_degree_partition_and_finite_radial_bounds": serialize(
            {"degrees": payload(degrees.data()), "radials": payload(radials.data())}
        ),
        "original_two_leg_tail_and_restored_actual_finite_piece": serialize(
            {"tail": payload(tail.data()), "limit": payload(limit.data())}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The complete first-five-endpoint spatial Taylor sector splits into ten absolutely integrable cells j+n>=5 and fifteen candidate UV cells j+n<=4, retaining all75 source time-jet entries and all nine physical pairs. Exact massive radial moments and the complete source-jet Cauchy bounds control the finite continuum polynomial by1e50||D||L2 X46[Gamma]. Its original both-created-mode regulator error is below1e55||D||L2 X46[Gamma]/K, after separate low/high external-transfer estimates in the same norm. Only after dominated removal is this convergent sector a genuine finite spatial polynomial. Every retained term has positive spatial degree, leaving zero transfer unchanged. These finite oversubtractions are restored to the S203 known actual-current piece, giving5e54 M[D]Y[Gamma] and tail7e60 M[D]Y[Gamma]/K, with canonical displays2e-745 and28e-740/K. The S204 fixed local target is not double-counted as an independently matched quantum contribution. The complete contact and fifteen candidate cells still require their actual UV expansion, finite/divergent coefficients, regulator artifacts and original fixed covariant matching; candidate is not a claim that each cell diverges. No full matched response, mixed inverse, interacting background, stability, cutoff, remaining parent loops or finite-gravity IR/Regge is established. Original V/G/B and P8 remain open.",
        "not_established": [
            "Actual UV coefficients, finite/divergent and regulator-artifact matching of the remaining fifteen candidate cells with the complete contact",
            "Polynomiality or locality of the remaining moving sharp-band integral, or actual divergence of every candidate cell",
            "A full matched response, derivative-compatible inverse, interacting background or stability",
            "Full parent/cutoff/scattering/IR/Regge matching or original V/G/B and P8 closure",
        ],
        "verification_boundary": "Exact complete degree/time-jet partition, finite massive radial moments, removed-union power bounds, restoration identity and canonical factors. Independent actual mixed spatial/time coefficients at multiple momentum scales and two Cauchy resolutions retain all nine physical pairs. Original two-leg angular integrals, half-band tails, exact order/prewarmed-cache guards and a borderline-majorant control supplement written continuum estimates. Native, direct science, ordinary and CLI use original SymPy; only full regression uses the audited exact GCD adapter. Cauchy estimates, dominated removal and the finite-polynomial restoration are not FORMALIZED.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The finite endpoint Taylor restoration report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.205.TEN_FINITE_CURVED_ENDPOINT_SPATIAL_TAYLOR_CELLS_RESTORED_WITH_ORIGINAL_REGULATOR_TAIL replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
