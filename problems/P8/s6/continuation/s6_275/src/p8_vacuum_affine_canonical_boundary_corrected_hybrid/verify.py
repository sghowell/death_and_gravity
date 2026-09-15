"""Read-only original canonical-boundary correction and new finite hybrid report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_hybrid_core_regulator_comparison import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT
    / "certificates/polynomial-vacuum-affine-canonical-boundary-corrected-hybrid.json"
)
PARENT_SHA = "4854ad5b56788d8f90c01d24055b3f0337f5285b9f80926e71ffc30286099ecd"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(
            ROOT.glob("src/p8_vacuum_affine_canonical_boundary_corrected_hybrid/*.py")
        )
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen archived finite hybrid comparison report changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    return {
        "S6_274_archived_comparison_and_full_ancestry_rebuilt": PARENT_SHA,
        "archived_numeric_replay_does_not_reendorse_refuted_physical_dictionary": True,
        "S269_to_S274_physical_identifications_explicitly_qualified": True,
        "new_corrected_solutions_not_identical_to_archived_solutions": True,
        "original_raw_source_and_prepared_seed_retained_with_complete_transport": True,
        "S261_refutation_and_S265_integrability_boundary_unchanged": True,
        "all_original_primitive_and_matching_frontiers_retained": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A corrected canonical-boundary hybrid proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.275.ORIGINAL_ADM_CANONICAL_BOUNDARY_DICTIONARY_ERRATUM_AND_NEW_CORRECTED_FINITE_HYBRID_TURNAROUND_COMPARISONS",
        "date": "2026-09-15",
        "status": "CANONICAL_BOUNDARY_CORRECTION_AND_NEW_FINITE_HYBRID; NOT_ARCHIVED_SAME_PHYSICAL_STATE_OR_ORIGINAL_V_G_B_P8",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/boundary.md",
            "notes/source.md",
            "notes/geometry.md",
            "notes/quantum.md",
            "notes/dynamics.md",
            "notes/validation.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_canonical_boundary_erratum_source_and_corrected_geometry": serialize(
            {name: payload(packets[name]) for name in names[:3]}
        ),
        "whole_corrected_quantum_hybrid_and_regulator_comparisons": serialize(
            {name: payload(packets[name]) for name in names[3:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The original raw ADM versus prepared momenta differ by the S220 boundary generating function f=-9Hv^2+3ell v sigma. Both independent raw-velocity and full tree-H2 derivations confirm it. Its two bounce momentum shifts change every positive configuration covariance and cannot be a global scalar phase. The archived S269-S274 raw/prepared physical identification is explicitly qualified without rewriting frozen records. After exact fixed-reference canonical transport, both original cutoffs and full orderings define NEW corrected finite hybrid solutions on |u|<=1e-180. Homogeneous deviation is<1e-390, volume relative error<1e-380, endpoint gap>5e-360 and every minimum lies in |u|<1e-188. Positive core leakage is<1e-6; coupled cutoff(Y,phase-factored state,volume) differences are<1e-580,1e-12,1e-514 and ordering differences<1e-630,1e-65,1e-566. The full original source, primitive, heavy, Gauss, nonlinear reconstruction, scalar-center forces and complete positive-heat remainders are retained. These claims concern the corrected finite hybrid, not archived-solution equality, a full nonlinear boundary derivation, homogeneous quantization, sharp support, unique or strict minimum, regulator removal, unlocalized dynamics, physical UV matching, omitted-loop/Regge control or global completion. Original V/G/B/P8 remain OPEN.",
        "not_established": [
            "Equality of the corrected physical state or solutions with the archived unshifted model",
            "An unproved full nonlinear ADM boundary generating function or nonlinear metaplectic covariance",
            "An irrelevant global boundary phase or actual live-profile H2 identification with the reference",
            "Exact joint support, unique or strict volume minimum, or homogeneous quantization",
            "Unlocalized dynamics, uniform mode/volume/cutoff removal, original fully quantum means or physical UV matching",
            "Omitted-loop/Regge/global completion or formal verification of the written analytic proofs",
        ],
        "verification_boundary": "Exact original-source identities, independent boundary/flow/covariance checks, corrected outward full-domain bounds and independent operator/dynamical diagnostics support explicit written proofs for a NEW corrected finite hybrid. Full ancestry replays its frozen contents, not its refuted physical interpretation. Native/direct/ordinary/CLI use original SymPy; only an audited captured complete regression may use the exact-GCD adapter. No historical source or report is rewritten and original P8 is not closed.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete canonical-boundary corrected hybrid report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.275 canonical-boundary corrected finite HYBRID replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
