"""Read-only MS fermion self-energy-chord certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_fermion_proper_references import verify as parent

from . import audit, calibration, catalog, domain, kernel, tail

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-fermion-self-energy-chord.json"
PARENT_SHA = "f7ae3eeb624a0d733eb0a0d27fd7b17c56aec0e393c3cebad9021a0947f351a8"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_fermion_self_energy_chord/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen MS fermion self-energy-chord report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_139_fully_rebuilt": PARENT_SHA,
        "same_GY14_nongravitational_reference_model": True,
        "all_prior_scientific_sources_and_reports_unchanged": True,
        "self_energy_subsets_not_complete_primitives_or_original_P8": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows = audit.residuals()
    gates = audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A MS fermion self-energy-chord gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.140.MS_FERMION_SELF_ENERGY_CHORD_BOUNDS",
        "date": "2026-09-10",
        "status": "RENORMALIZED_SCALAR_AND_GAUGE_SELF_ENERGY_CHORD_BOUNDS; NOT_COMPLETE_PRIMITIVE_ROWS_TWO_LOOP_MATCHING_HIGHER_LOOP_V_G_B_OR_ORIGINAL_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/catalog.md",
            "notes/kernel.md",
            "notes/tail.md",
            "notes/domain.md",
            "notes/literature.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_MS_kernel_and_global_domain": serialize(
            {"kernel": payload(kernel.data()), "domain": payload(domain.data())}
        ),
        "primitive_word_ownership": serialize(payload(catalog.data())),
        "soft_tail_and_radial_bound": serialize(payload(tail.data())),
        "actual_partial_primitive_frontier": serialize(
            {
                "actual": payload(calibration.data()),
                "primitive_frontier": audit.frontier(),
            }
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The complete renormalized one-loop MS fermion self-energy, including scalar and gauge exchange and its proper mass/kinetic counterterms, is inserted into every marked propagator of the first four-scalar fermion loop. The six cyclic boxes and four marks give 24 self-energy-chord words in each primitive sector. A global complex parameter gap and Dirac norm bound yield a momentum-dependent soft Cauchy radius. The degree-four-and-higher soft tail is projected before outer integration; exact integration gives the radial moment 41/36. Complete-subset Lorentz/permutation symmetry makes the lower soft degrees irrelevant to b2. The combined scalar/gauge bound is below 1e-1400, relative to the tree below 1e-800, in the shared proper MS reference. Only two partial primitive rows advance: the vertex and other words, remaining conversion/insertion tasks, complete two-loop error, higher loops, V, G, B and original P8 remain open.",
        "not_established": [
            "Twenty-four vertex and twelve other words in each primitive scalar/gauge quartic sector",
            "Five other wholly unevaluated primitive fermion-sector rows",
            "Other counterterm, canonical-field and parameter-conversion contributions",
            "S6.138 outer-reference finite interaction conversion",
            "Complete enlarged-model two-loop amplitude, pole or error",
            "Higher loops, exact spectrum, cutoff, V contours, G, B or original P8 closure",
        ],
        "verification_boundary": "Exact cyclic-word/marked-propagator bijection and determinant variation, complete dimensional-MS subgraph algebra, global parameter-gap and Dirac identities, momentum-dependent Cauchy-tail constants and exact radial integration with a mutation-checked partial frontier. The continuum analytic and symmetry arguments are written proofs, not proof-assistant formalization or independent peer review. Native, ordinary, CLI and direct science use unmodified SymPy with interpreter-only allowances.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The MS fermion self-energy-chord differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.140.MS_FERMION_SELF_ENERGY_CHORD_BOUNDS replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
