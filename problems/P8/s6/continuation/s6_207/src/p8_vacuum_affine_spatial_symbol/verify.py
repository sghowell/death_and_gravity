"""Read-only full spatial UV symbol with independent source and detector times."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_proca_gaussian import verify as gaussian
from p8_vacuum_affine_sharp_band_artifact import verify as parent

from . import audit, benchmark, extraction, sectors, symbol

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-spatial-symbol.json"
PARENT_SHA = "baef3e1fbb6564045eecc11b245bb7ab984bbc11cfd0de62a115e0edb2e90c31"
GAUSSIAN_SHA = "1175702d6980db090396863b807af0e767241daeb7d8dabf138ed11c67e1ae75"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_spatial_symbol/*.py"))
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
        "S6_206_fully_rebuilt": PARENT_SHA,
        "S6_176_actual_CD_Gaussian_sector_transitively_rebuilt": GAUSSIAN_SHA,
        "fixed_stress_prescription_and_full_source_retained": True,
        "all_frozen_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A complete spatial UV symbol gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.207.COMPLETE_FOUR_SECTOR_SPATIAL_UV_SYMBOL_AND_FINITE_ENDPOINT_COEFFICIENT_EXTRACTION",
        "date": "2026-09-11",
        "status": "COMPLETE_FOUR_SECTOR_SPATIAL_UV_SYMBOL_AND_FINITE_ENDPOINT_COEFFICIENT_EXTRACTION; NOT_FULL_CONTACT_REGULATOR_COVARIANT_MATCHING_ALL_P_NORM_INVERSE_BACKGROUND_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/sectors.md",
            "notes/normalization.md",
            "notes/extraction.md",
            "notes/benchmark.md",
            "notes/curved-tests.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "complete_two_time_polarization_factorization_and_normalized_symbol": serialize(
            {"sectors": payload(sectors.data()), "symbol": payload(symbol.data())}
        ),
        "finite_UV_extraction_and_nonzero_massive_benchmark": serialize(
            {
                "extraction": payload(extraction.data()),
                "benchmark": payload(benchmark.data()),
            }
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The complete tracefree spatial Proca pair kernel has an exact four-sector projector factorization retaining all nine physical pairs, the longitudinal constraint, independent detector/source times and distinct inverse phases. Setting x=1/r and using scaled momenta n,-n+xP gives an exact analytic normalization without moving individual polarization frames. Each endpoint is x^(j-1) times its normalized analytic source-jet row. Degrees0,...,4-j give all fifteen potentially nonintegrable endpoint/power slots and all35 source time-jet entries; each degree-d coefficient is polynomial in external momentum of degree at most d. The actual sixth/eighth W8 frequency terms begin at x^6/x^8 and do not affect these retained Taylor jets, without replacing the finite-momentum reference or prepared state. The subtracted high-momentum remainder is O(k^-4) at fixed P, with a fixed-P integrable tail; an explicit uniform all-P norm is not supplied. Exact massive flat tensor/vector/scalar coefficients include nonzero mass and momentum corrections, and independent full-field actual curved W8 calculations validate the factorization, source-time jets, all retained coefficients at two Cauchy resolutions, higher-W8 independence and nonzero fourth-power remainders. A third-odd-endpoint UV cancellation in the chosen fixture does not remove its nonzero finite endpoint. The original pair/contact bands and S206 leading nonpolynomial artifact remain. Full subleading regulator conversion, actual contact and finite/divergent fixed covariant matching, full response/inverse/background/stability, physical cutoff and original V/G/B are not established. Original P8 remains open.",
        "not_established": [
            "Complete subleading sharp-band conversion and quantum contact/candidate-cell finite and divergent coefficients in the original covariant prescription",
            "An explicit all-external-momentum remainder norm, renormalized comparator, changed reference/state or permission to alter matching",
            "Deletion of all odd endpoints from the curved response or a full physical current divergence or model exclusion",
            "Full matched response, derivative-compatible inverse, finite-coupling background/stability, parent loops/cutoff/IR/Regge or original V/G/B and P8 closure",
        ],
        "verification_boundary": "Exact general two-time four-sector identities, inverse-radius normalization, finite slot/time-jet counts, source-only derivative order and massive flat coefficient recurrences. Independent original physical W8 modes, full field vs projector calculations, all source time jets, two actual curved coefficient resolutions, higher-W8 jet independence, massive flat field-level Cauchy coefficients, fixed-P remainder scaling and wrong-coincidence/common-frequency/odd-endpoint controls supplement written analyticity and power-counting arguments. Native, direct science, ordinary and CLI use original SymPy; only full regression uses the audited exact GCD adapter. The fixed-P analytic remainder and finite UV extraction are not FORMALIZED.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete spatial UV symbol report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.207.COMPLETE_FOUR_SECTOR_SPATIAL_UV_SYMBOL_AND_FINITE_ENDPOINT_COEFFICIENT_EXTRACTION replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
