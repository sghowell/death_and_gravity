"""Read-only full continuum response on strictly separated supports."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_proca_gaussian import verify as gaussian
from p8_vacuum_affine_spatial_current import verify as parent

from . import audit, projection, response, support, tail

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-separated-response.json"
PARENT_SHA = "ac027ba606e11566837d2186e3cd8ad28888ff2669a585bd2085c3f0bdc12a64"
GAUSSIAN_SHA = "1175702d6980db090396863b807af0e767241daeb7d8dabf138ed11c67e1ae75"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_separated_response/*.py"))
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
        "S6_195_fully_rebuilt": PARENT_SHA,
        "S6_176_actual_CD_Gaussian_sector_transitively_rebuilt": GAUSSIAN_SHA,
        "fixed_stress_prescription_and_full_source_retained": True,
        "all_frozen_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A separated-support continuum weak response gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.196.ACTUAL_SPATIAL_CONTINUUM_SEPARATED_SUPPORT_WEAK_PROCA_CURRENT_RESPONSE_WITH_QUANTITATIVE_REGULATOR_TAIL",
        "date": "2026-09-11",
        "status": "ACTUAL_SPATIAL_CONTINUUM_SEPARATED_SUPPORT_WEAK_PROCA_CURRENT_RESPONSE_WITH_QUANTITATIVE_REGULATOR_TAIL; NOT_COINCIDENT_RETARDED_EXTENSION_OPERATOR_NORM_INVERSE_BACKGROUND_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/tail.md",
            "notes/projection.md",
            "notes/response.md",
            "notes/support.md",
            "notes/limit.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "complete_stress_vector_and_projection_tail": serialize(
            {"tail": payload(tail.data()), "projection": payload(projection.data())}
        ),
        "ordered_support_continuum_weak_response": serialize(
            {"support": payload(support.data()), "response": payload(response.data())}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "For the unchanged actual conditional Gaussian Proca sector, the full centered smeared stress creates a two-particle vector. Its both-leg momentum projection has a uniform quantitative tail below1e52 N[f]^2/K, including high external momentum via the existing spatial H1 norm. The same orthogonal projection for source and readout gives an exact discarded-tail inner product, hence a weak response cutoff error below5e51 N[D]N[Gamma]/K. On strictly time-ordered compact supports all metric and unchanged finite-curvature contacts vanish, and the full first-order response is the convergent i/4 stress commutator. Its complete bound is5e49, its two-factor canonical display2e-750, and the canonical tail error2e-748/K, in the stated test norms. Full Proca locality gives zero on mutually spacelike supports; a finite spatial projection need not preserve that property. No coincident-time retarded extension, general operator-norm differentiability, finite-amplitude remainder, full response inverse, interacting background, physical cutoff or original V/G/B and P8 closure follows.",
        "not_established": [
            "The fixed renormalized spatial/mixed retarded response on overlapping supports or the diagonal",
            "Operator-norm metric differentiability, a finite-amplitude remainder, full inverse or interacting background",
            "Full parent state/measure/loops, physical heavy/cutoff/threshold and omitted-order matching",
            "Complete vacuum cuts/contour/truncation, finite-gravity IR/Regge or original V/G/B and P8 closure",
        ],
        "verification_boundary": "Exact all-momentum majorant constants, common-projection algebra, full physical stress normalization and support conditions. Independent massive radial and nonzero-transfer integrals, literal multimode Fock operators and an overlapping-time counterexample supplement continuous written arguments. Native, direct science, ordinary and CLI use original SymPy; only full regression uses the audited exact GCD adapter. The analytic Hilbert-space and causal-support arguments are not FORMALIZED.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The separated-support weak current report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.196.ACTUAL_SPATIAL_CONTINUUM_SEPARATED_SUPPORT_WEAK_PROCA_CURRENT_RESPONSE_WITH_QUANTITATIVE_REGULATOR_TAIL replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
