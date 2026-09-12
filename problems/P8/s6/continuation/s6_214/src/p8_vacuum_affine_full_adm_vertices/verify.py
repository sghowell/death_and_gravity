"""Read-only full constrained ADM Proca vertices and form bounds."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_matched_spatial_current import verify as parent
from p8_vacuum_affine_proca_gaussian import verify as gaussian

from . import adm, audit, bounds, chart, vertices

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-full-adm-vertices.json"
PARENT_SHA = "2da5ec3dbafea291c83232eb475e619080efe7315853deaa12d79b11d854809a"
GAUSSIAN_SHA = "1175702d6980db090396863b807af0e767241daeb7d8dabf138ed11c67e1ae75"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_full_adm_vertices/*.py"))
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
        "S6_213_fully_rebuilt": PARENT_SHA,
        "S6_176_actual_CD_Gaussian_sector_transitively_rebuilt": GAUSSIAN_SHA,
        "fixed_stress_prescription_and_full_source_retained": True,
        "all_frozen_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A full constrained ADM Proca vertex gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.214.FULL_CONSTRAINED_ADM_PROCA_VERTICES_CHART_CONTACTS_AND_ALL_MOMENTUM_FORM_BOUNDS",
        "date": "2026-09-12",
        "status": "FULL_CONSTRAINED_ADM_PROCA_VERTICES_CHART_CONTACTS_AND_ALL_MOMENTUM_FORM_BOUNDS; NOT_FULL_NEW_SECTOR_CONTINUUM_MATCHING_REDUCED_INVERSE_BACKGROUND_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/adm.md",
            "notes/vertices.md",
            "notes/shift.md",
            "notes/chart.md",
            "notes/bounds.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "actual_full_ADM_constraint_and_complete_metric_vertices": serialize(
            {"adm": payload(adm.data()), "vertices": payload(vertices.data())}
        ),
        "four_metric_contact_chain_and_all_momentum_form_bounds": serialize(
            {"chart": payload(chart.data()), "bounds": payload(bounds.data())}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The full constrained Gaussian Proca ADM Hamiltonian and all first/second lapse, shift and symmetric spatial metric vertices are derived in the unchanged actual CD parent. Literal four-metric contraction and Legendre reduction retain the lapse/trace-dependent positive temporal-constraint energy and the longitudinal shift transport term. Full noncommuting spatial contacts recover the frozen tracefree sector. The complete shift is the covector Lie generator with both momenta and reverse-pair reality retained. Its ten-field feature form has exact norm a||beta|| for arbitrary complex Fourier beta, established by the singular-value projector identity, not only real Hermiticity. The complete unreduced direction norm sigma=|n|+(5/2)||Q||op+a||beta|| bounds first and second energy-relative vertices uniformly in all momenta. The same actual-state finite-band covariance response and contact bounds then apply. Nonlinear ADM-to-four-metric second derivatives keep the actual one-point contact chain; zero ADM shift/lapse contacts are not deleted covariant contacts. The full nonlinear affine source is not set to zero: the existing clock-source vanishing identity supplies only the stated Gaussian linear-response bridge. S213's matched tracefree continuum response and historical display corrections are preserved. Independent arbitrary-SPD action/Legendre, literal matrix-exponential, Fourier, complex reality, singular-value and omission checks support the formulas. No full new-sector ultraviolet matching, reduced scalar/mixed inverse, finite-amplitude spatial C2 response, controlled finite-coupling background/stability, other parent loops/cutoff, finite-gravity IR/Regge or original V/G/B and P8 closure is supplied.",
        "not_established": [
            "Full lapse/shift/spatial-trace Gaussian continuum state/time/UV/contact matching and associated Ward structure",
            "Genuinely reduced scalar/mixed inverse, canonical reduced-mode normalization or stable nonlinear quantum solution",
            "Finite-amplitude spatial C2 response, interacting finite-coupling background/stability, remaining parent measure/sectors/loops or physical cutoff",
            "Permission to delete longitudinal constraint/shift terms, four-metric chart contacts, full nonlinear source or unchanged preparation",
            "Canonical vacuum cuts/contour/truncation, finite-gravity IR/Regge or original V/G/B and P8 closure",
        ],
        "verification_boundary": "Exact ADM Legendre/constraint/boundary algebra, full two-momentum first/second vertices, sharp complex shift norm, noncommuting four-metric chain contacts and all-momentum energy forms. Independent four-metric action contractions, arbitrary-SPD magnetic and Legendre checks, literal full spatial expansions, complex reverse-pair reality and energy-normalized singular values through widely separated momenta supplement the written proofs. These continuous field/Fourier/form-norm arguments are not FORMALIZED. Native/direct/ordinary/CLI retain original SymPy; full regression alone uses the audited exact-GCD adapter. The full new-sector continuum and genuinely reduced inverse remain open.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The full constrained ADM Proca vertices report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.214.FULL_CONSTRAINED_ADM_PROCA_VERTICES_CHART_CONTACTS_AND_ALL_MOMENTUM_FORM_BOUNDS replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
