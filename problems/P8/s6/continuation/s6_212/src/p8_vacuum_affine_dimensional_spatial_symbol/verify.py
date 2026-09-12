"""Read-only actual dimensional spatial UV finite part."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_proca_gaussian import verify as gaussian
from p8_vacuum_affine_spatial_matching_difference import verify as parent

from . import audit, density, geometry, jets, matching

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates" / "polynomial-vacuum-affine-dimensional-spatial-symbol.json"
)
PARENT_SHA = "94514af725bfc242cf045064e48742c22e5282eb966daa110483c46928ad3f65"
GAUSSIAN_SHA = "1175702d6980db090396863b807af0e767241daeb7d8dabf138ed11c67e1ae75"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_dimensional_spatial_symbol/*.py"))
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
        "S6_211_fully_rebuilt": PARENT_SHA,
        "S6_176_actual_CD_Gaussian_sector_transitively_rebuilt": GAUSSIAN_SHA,
        "fixed_stress_prescription_and_full_source_retained": True,
        "all_frozen_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("An actual dimensional spatial UV finite part gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.212.ACTUAL_DIMENSIONAL_SPATIAL_UV_FINITE_PART_WITH_EVANESCENT_MODE_AND_COUNTERTERM_VARIATIONS",
        "date": "2026-09-11",
        "status": "ACTUAL_DIMENSIONAL_SPATIAL_UV_FINITE_PART_WITH_EVANESCENT_MODE_AND_COUNTERTERM_VARIATIONS; NOT_FULL_DIMENSION_LIMIT_ASSEMBLED_RESPONSE_INVERSE_BACKGROUND_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/geometry.md",
            "notes/modes.md",
            "notes/dimension.md",
            "notes/curvature.md",
            "notes/finite.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "actual_dimensional_geometry_modes_and_fixed_source_UV_coefficients": serialize(
            {
                "jets": payload(jets.data()),
                "density": payload(density.data()),
                "geometry": payload(geometry.data()),
            }
        ),
        "continued_fixed_counterterm_original_MSbar_finite_part_and_norms": serialize(
            {"matching": payload(matching.data())}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The actual dimension-dependent spatial ultraviolet finite difference is evaluated in the original MSbar prescription. The full field-strength vertex and active/passive sphere contractions retain the actual dimensional mode multiplicity. Actual continued canonical rates and all WKB/source jets reproduce all105 physical angular coefficients. Invariant reconstruction fixes the physical source before taking dimension derivatives; the evanescent odd endpoint is nonzero and retained. Full continued R_old,R_squared,Ricci_squared,Riemann_squared and volume Hessians are varied at the original fixed four-dimensional scalar pole weights, retaining the Euler and log(a) contributions. The exact finite expression keeps the original comoving lower-band polynomial and scale. Its first-time coefficient is the proper-time derivative of the second one, and an explicit rational log-polynomial norm is below1e5 Z24, canonically4e-795. Independent literal field-strength/sphere and noncommuting full metric second variations in dimensions3,4,5,6, all-four-order WKB Cauchy coefficients, complex-dimension radial finite-part extraction at two resolutions and omission controls check the result. No full assembled-current dimension limit, reduced mixed inverse, finite-coupling background/stability, remaining parent loops/cutoff, finite-gravity IR/Regge or original V/G/B and P8 closure is supplied.",
        "not_established": [
            "Full subtracted spatial reference dimension-limit interchange and assembly with the fixed homogeneous current and actual known remainders",
            "Permission to reset the state/reference, change the original regulator or finite coefficients, remove evanescent endpoints/Euler/volume terms, or count the local target twice",
            "A full matched spatial/scalar/metric response or reduced inverse, controlled finite-coupling background/stability or physical cutoff",
            "Remaining parent loops/scattering/IR/Regge or original V/G/B and P8 closure",
        ],
        "verification_boundary": "Exact genuine general-dimensional field-strength geometry, actual dimensional canonical rates and UV jets, all physical coefficient limits and fixed-source invariant derivatives, continued fixed pole Hessians and original MSbar radial normalization, finite proper-time Green identity and explicit coefficient norm. Independent integer-dimensional matrix/sphere checks, literal mixed metric curvature variations, full four-order WKB coefficient tests, complex-epsilon finite parts and two resolutions supplement the written analytic derivation. Native, direct science, ordinary and CLI retain original SymPy; only full regression uses the audited exact GCD adapter. Continuous Fourier, local-variation and dimensional arguments are not FORMALIZED; full current assembly remains outside this checkpoint.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The actual dimensional spatial UV finite part report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.212.ACTUAL_DIMENSIONAL_SPATIAL_UV_FINITE_PART_WITH_EVANESCENT_MODE_AND_COUNTERTERM_VARIATIONS replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
