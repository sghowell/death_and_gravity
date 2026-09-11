"""Read-only actual C2 covariance-tail response and finite-amplitude remainder."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_matrix_adiabatic import verify as parent
from p8_vacuum_affine_proca_gaussian import verify as gaussian

from . import audit, covariance, mixed, reference, variation

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-matrix-response-tail.json"
PARENT_SHA = "861466f6f14c73ba732ee81265a21fed3e237adcada37e1a78749f9b6c21d015"
GAUSSIAN_SHA = "1175702d6980db090396863b807af0e767241daeb7d8dabf138ed11c67e1ae75"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_matrix_response_tail/*.py"))
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
        "S6_190_fully_rebuilt": PARENT_SHA,
        "S6_176_actual_CD_Gaussian_sector_transitively_rebuilt": GAUSSIAN_SHA,
        "fixed_stress_prescription_and_full_source_retained": True,
        "all_frozen_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("An actual C2 covariance-tail response gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.191.ACTUAL_MATRIX_COVARIANCE_TAIL_C2_SHEAR_RESPONSE_AND_FINITE_AMPLITUDE_TAYLOR_REMAINDER",
        "date": "2026-09-11",
        "status": "ACTUAL_UNCHANGED_STATE_MATRIX_COVARIANCE_TAIL_C2_PARAMETER_RESPONSE_AND_TAYLOR_REMAINDER; NOT_COMPLETE_COVARIANT_STRESS_RESPONSE_INVERSE_FINITE_COUPLING_BACKGROUND_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/mixed.md",
            "notes/reference.md",
            "notes/variation.md",
            "notes/covariance.md",
            "notes/tail.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "full_mixed_jets_and_reference": serialize(
            {"mixed": payload(mixed.data()), "reference": payload(reference.data())}
        ),
        "actual_C2_graph_covariance_and_weighted_tail": serialize(
            {
                "variation": payload(variation.data()),
                "covariance": payload(covariance.data()),
            }
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The unchanged actual three-mode state has complete first and second shear-amplitude derivative bounds for its finite-reference covariance error. Mixed positive-root and ordered reference derivatives retain polarization mixing without a condition-number shortcut. Exact nonlinear error equations retain the quadratic first-response contribution, and the full covariance and frequency-weight derivatives are included. Common infinite majorants give a C2 high-band weighted covariance-tail integral with norm bounds1e-65,1e-47,1e-28 through derivative order2 and finite-amplitude Taylor remainder at most epsilon^2*1e-28/2 (zero at epsilon0 and strict otherwise). This is not the full fixed-covariant stress response, its remaining contacts, low band, full feedback inverse, interacting background or physical cutoff. Original V/G/B and P8 remain OPEN.",
        "not_established": [
            "The complete fixed-covariant UV-subtracted stress response and all finite contacts",
            "The low-band and full physical observable response or a full gravitational feedback remainder",
            "Full spatial/mixed/nonlinear inverse, interacting quantum background and parent state/measure/loops",
            "Physical heavy/cutoff/threshold, omitted-order matching, finite-gravity IR/Regge or original V/G/B and P8 closure",
        ],
        "verification_boundary": "Exact mixed positive-root and reference recurrences, complete noncommuting first/second graph error equations, full covariance derivatives and continuous infinite-tail/Taylor estimates. Independent matrix, parameter-reference, covariance and weighted-integral fixtures supplement the written proofs. Native, direct science, ordinary and CLI use original SymPy; only full regression uses the audited exact GCD adapter. Analytic and functional arguments are not FORMALIZED.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The actual C2 covariance-tail report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.191.ACTUAL_MATRIX_COVARIANCE_TAIL_C2_SHEAR_RESPONSE_AND_FINITE_AMPLITUDE_TAYLOR_REMAINDER replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
