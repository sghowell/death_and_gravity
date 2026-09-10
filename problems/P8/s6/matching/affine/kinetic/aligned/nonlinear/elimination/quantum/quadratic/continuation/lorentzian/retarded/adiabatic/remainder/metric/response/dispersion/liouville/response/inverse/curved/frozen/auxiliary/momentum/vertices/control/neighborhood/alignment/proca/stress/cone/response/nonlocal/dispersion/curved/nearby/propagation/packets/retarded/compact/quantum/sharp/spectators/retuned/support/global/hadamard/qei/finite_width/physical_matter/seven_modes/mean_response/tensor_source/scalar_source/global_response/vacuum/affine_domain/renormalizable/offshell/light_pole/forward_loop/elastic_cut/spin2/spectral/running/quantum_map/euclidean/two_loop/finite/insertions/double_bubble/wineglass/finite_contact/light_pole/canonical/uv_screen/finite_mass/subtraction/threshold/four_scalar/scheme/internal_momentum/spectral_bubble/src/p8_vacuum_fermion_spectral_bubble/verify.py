"""Read-only first fermion spectral insertion in the local-quartic bubble."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_global_one_loop_insertions import verify as parent

from . import audit, calibration, density, outer, ownership, spectral

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-fermion-spectral-bubble.json"
PARENT_SHA = "a6d91c84be131cfd36d877a870213f15c61c16e487ec00816c1eff98cfb98404"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_fermion_spectral_bubble/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen global one-loop insertion report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_134_fully_rebuilt": PARENT_SHA,
        "same_whole_on_shell_fermion_kernel_and_reference_boundary": True,
        "prior_scientific_sources_and_reports_unchanged": True,
        "literal_family_not_complete_two_loop_or_original_P8_closure": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows = audit.residuals()
    gates = audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A first spectral outer-bubble gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.135.FIRST_FERMION_SPECTRAL_OUTER_BUBBLE",
        "date": "2026-09-10",
        "status": "POSITIVE_FINITE_FIRST_FERMION_SPECTRAL_QUARTIC_BUBBLE_FAMILY; NOT_COMPLETE_TWO_LOOP_HIGHER_LOOP_ERROR_EXACT_SPECTRUM_CUTOFF_V_G_B_OR_ORIGINAL_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/spectral.md",
            "notes/outer.md",
            "notes/ownership.md",
            "notes/bounds.md",
            "notes/literature.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": len(rows),
        "proof_checks": gates,
        "complete_spectral_insertion": serialize(payload(spectral.data())),
        "positive_measure_and_outer_subtraction": serialize(
            {"density": payload(density.data()), "outer": payload(outer.data())}
        ),
        "six_placement_reference_ownership": serialize(payload(ownership.data())),
        "actual_family_coefficient_bounds": serialize(payload(calibration.data())),
        "controls": serialize(audit.controls()),
        "verdict": "A direct parameter-integral identity gives a nonnegative spectral representation of the entire first on-shell fermion propagator insertion. Its total formal continuum weight is infinite; no normalized exact spectral theorem is claimed. The local-quartic outer bubble includes both internal insertion positions and all three channels, with the inherited symmetry and loop factors. After pairing the complete inner mass/residue references and subtracting one overall local quartic constant, the spectral integral and its derivatives have explicit integrable complex-domain majorants. The literal six-contribution family's forward second coefficient is strictly positive and below 10^-1418, or 10^-819 times the positive tree reference. A finite constant conversion cannot tune this family's b2. Other two-loop graph and parameter-conversion families, higher-loop/truncation errors, full spectrum, justified cutoff, V contour, finite-gravity G, common-parent B and original P8 closure remain unestablished.",
        "not_established": [
            "Complete enlarged-model two-loop amplitude and reference conversion",
            "Other heavy-exchange, fermion-box and primitive or insertion families",
            "Higher-loop and truncation error budgets or strict full V verdict",
            "Normalized exact spectral measure, reflection positivity or full spectrum",
            "Derived cutoff and high-energy V contour",
            "Finite-gravity IR/Regge Delta for G and common bounce-parent B",
            "Original P8 closure or exclusion of any full old ladder row",
        ],
        "verification_boundary": "Exact spectral differentiation, frozen-kernel agreement, parameter interval algebra, proper/overall local-reference identities, literal insertion counts and independent rational bounds. Analytic continuation, dominated integral operations and continuum positivity are written proofs, not proof-assistant formalized or all-orders quantum-field-theory claims. Native and ordinary scientific replay uses unmodified SymPy with interpreter-only runtime allowances.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The spectral outer-bubble report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.135.FIRST_FERMION_SPECTRAL_OUTER_BUBBLE replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
