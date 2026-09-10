"""Read-only first spectral covariance insertion in the scalar quadratic sector."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_fermion_spectral_bubble import verify as parent

from . import audit, calibration, enclosure, kernel, quadratic, selected

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-fermion-spectral-pole.json"
PARENT_SHA = "be4188f2bcd3ff3b0d68a20fd2b385c0c9e936c2ab4df52626cefabf8e21e5c2"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_fermion_spectral_pole/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen first spectral outer-bubble report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_135_fully_rebuilt": PARENT_SHA,
        "same_complete_first_fermion_spectral_insertion": True,
        "prior_scientific_sources_and_reports_unchanged": True,
        "selected_quadratic_family_not_complete_two_loop_or_P8_closure": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows = audit.residuals()
    gates = audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A spectral quadratic-family gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.136.FIRST_FERMION_SPECTRAL_QUADRATIC_SECTOR",
        "date": "2026-09-10",
        "status": "FIRST_SPECTRAL_COVARIANCE_QUADRATIC_SECTOR_AND_SELECTED_POLE_CONTROL; NOT_COMPLETE_TWO_LOOP_HIGHER_LOOP_ERROR_SPECTRUM_CUTOFF_V_G_B_OR_ORIGINAL_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/quadratic.md",
            "notes/kernel.md",
            "notes/bounds.md",
            "notes/selected.md",
            "notes/literature.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": len(rows),
        "proof_checks": gates,
        "literal_quadratic_covariance_insertion": serialize(payload(quadratic.data())),
        "spectral_outer_kernel_and_enclosure": serialize(
            {"kernel": payload(kernel.data()), "enclosure": payload(enclosure.data())}
        ),
        "selected_functional_order_ledger": serialize(payload(selected.data())),
        "actual_selected_pole_bounds": serialize(payload(calibration.data())),
        "controls": serialize(audit.controls()),
        "verdict": "Differentiating the frozen complete scalar quadratic half-trace once in its light covariance includes all local tadpole terms and exactly one spectral light line in the mixed H-Phi bubble. The fixed H one-point reference and complete inner/outer on-shell subtractions are retained before estimation. Integrable spectral derivative bounds give a positive finite family slope below 10^-616 and a uniform quadratic remainder coefficient below 10^-1017. The complete one-loop quadratic functional plus precisely this sector, with its explicit mass/field references and formal-order normalization products, retains a unit-residue mass-one local pole, positive local Phi curvature, and the stated finite Euclidean reference-window control. It is not the complete enlarged-model two-loop inverse, a higher-loop/truncation error budget, exact spectrum, reflection positivity, derived cutoff, justified V contour, finite-gravity G, common-parent B or original P8 closure.",
        "not_established": [
            "Other enlarged-model two-loop quadratic primitives and interaction conversions",
            "Complete two-loop inverse and higher-loop or truncation error budget",
            "Global potential, exact stable-heavy/gauge spectrum or reflection positivity",
            "Derived Wilsonian cutoff or high-energy V contour and strict full V verdict",
            "Finite-gravity regulated Regge/IR Delta for G",
            "Controlled common bounce-parent field/state/cutoff dictionary for B",
            "Original P8 closure or exclusion of a complete old ladder row",
        ],
        "verification_boundary": "Exact full-half-trace differentiation, frozen mixed-bubble normalization, local forest/reference identities, spectral parameter derivatives, formal-order field bookkeeping and independent rational calibration. Analytic continuation and dominated continuum bounds are written proofs, not proof-assistant formalization or all-orders QFT. Native and ordinary scientific replay uses unmodified SymPy with interpreter-only runtime allowances.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The spectral quadratic-sector report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.136.FIRST_FERMION_SPECTRAL_QUADRATIC_SECTOR replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
