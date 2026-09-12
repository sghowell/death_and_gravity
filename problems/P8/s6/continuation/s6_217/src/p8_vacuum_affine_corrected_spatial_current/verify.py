"""Read-only corrected complete spatial current and prepared Ward input."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_ordered_scalar_symbol import verify as parent
from p8_vacuum_affine_proca_gaussian import verify as gaussian

from . import assembly, audit, response, subtraction

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates" / "polynomial-vacuum-affine-corrected-spatial-current.json"
)
PARENT_SHA = "70e8538425cfab3035cea3c4c9d3a5da63797182625f79811bbc9ad2cd9ed815"
GAUSSIAN_SHA = "1175702d6980db090396863b807af0e767241daeb7d8dabf138ed11c67e1ae75"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_corrected_spatial_current/*.py"))
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
            raise ValueError("A source-pinned historical Gaussian input changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_216_corrected_phase_and_full_scalar_UV_fully_rebuilt": PARENT_SHA,
        "S6_176_original_Gaussian_state_transitively_rebuilt": GAUSSIAN_SHA,
        "successor_reassembly_not_endorsement_of_withdrawn_historical_formulas": True,
        "no_frozen_scientific_test_or_report_bytes_changed": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A corrected full spatial current gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.217.CORRECTED_COMPLETE_RETARDED_TRACEFREE_SPATIAL_CURRENT_ORIGINAL_REGULATOR_AND_PREPARED_WARD_INPUT",
        "date": "2026-09-12",
        "status": "CORRECTED_COMPLETE_TRACEFREE_SPATIAL_WEAK_CURRENT_ORIGINAL_ANCHORED_REGULATOR_AND_KNOWN_WARD_INPUT; NOT_FULL_SCALAR_ANCHOR_REDUCED_INVERSE_BACKGROUND_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/phase.md",
            "notes/decomposition.md",
            "notes/subtraction.md",
            "notes/dimension.md",
            "notes/assembly.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "complete_corrected_retarded_decomposition_and_subtraction": serialize(
            {
                "canonical": payload(response.canonical_data()),
                "decomposition": payload(response.decomposition_data()),
                "physical_UV": payload(subtraction.physical_data()),
                "dimension": payload(subtraction.dimension_data()),
            }
        ),
        "full_current_original_regulator_and_prepared_Ward_input": serialize(
            {
                "assembly": payload(assembly.data()),
                "Ward": payload(assembly.ward_data()),
            }
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The independent unequal-mode covariance/CCR bridge fixes the annihilation-source orientation throughout the complete current, not only its local UV part. Every first-five endpoint and its UV subtraction acquires the same (-1)^j factor; the fifth endpoint reverses and the sixth bulk retains its full positive-detector-phase form. Actual-minus-reference comparison, endpoint and original-mask tails keep their absolute bounds. All unaveraged odd nonlog physical UV coefficients vanish; seven transfer-independent odd log entries, including source-jet1 magnetic terms, leave a decaying shell error and cannot change A3,A1,K2,K0. The corrected evanescent finite term remains mandatory. Both complex-dimensional phase branches are bounded explicitly and continued at the same parameter. The unchanged actual homogeneous one-point derivative supplies the anchor, giving Jcorrect=H0+Known_correct(P)-Known_correct(0)+Ffinite_correct(P), a2e95 M Z136 weak bound and3e54 M Z136/K original homogeneous-anchored regulator error. Both canonical external factors give8e-705 and12e-746/K. The independent prepared Ward maps and full chart terms now have a corrected known tracefree input below1e118 V04 U138. Finite time-dependent unequal-mode Fock/covariance fixtures independently test the entire retarded integral, all six endpoints, nonzero sixth bulk and wrong-phase controls. Two successor status changes explicitly restore only this new corrected tracefree current and known Ward input, without re-endorsing S212/S213's withdrawn formulas. Three ordered scalar kernels and the homogeneous trace anchor, full reduced inverse, nonlinear quantum background and remaining original V/G/B obligations remain open.",
        "not_established": [
            "The complete three ordered scalar currents, homogeneous trace anchor and their full state/time/contact/remainder norms",
            "Validity of the withdrawn S212 finite value or S213 uncorrected physical current formula",
            "An unanchored original-cutoff convergence rate, a physical EFT cutoff identification or a new finite counterterm",
            "A same-space reduced scalar inverse or nonlinear stability from small weak tracefree canonical bounds",
            "Full reduced response, finite-amplitude quantum background, remaining parent/loop/heavy/cutoff or original V/G/B and P8 closure",
        ],
        "verification_boundary": "Exact finite CCR/covariance algebra and independent unequal time-dependent oscillator Fock calculations check the full phase/endpoint/bulk bridge. Unaveraged UV parity, all physical tracefree spatial differences, evanescent finite change, dimension coefficient majorants, norm arithmetic and scope are checked. Continuous infinite-mode, dimensional, Sobolev and Ward arguments are written proofs, not FORMALIZED. Native/direct/ordinary/CLI retain original SymPy; only full regression uses the audited exact-GCD adapter. Frozen historical scientific, proof, test and report bytes remain unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The corrected complete spatial current report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.217 corrected complete tracefree current and known Ward input replay passed; three scalar kernels and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
