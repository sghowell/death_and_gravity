"""Read-only same-model finite hybrid core and regulator certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_selfconsistent_finite_feedback import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates/polynomial-vacuum-affine-hybrid-core-regulator-comparison.json"
)
PARENT_SHA = "1ef7c3a7165db3688c4141313ed8fb7ee5a2c43ee0613ddd787144a8b201cdce"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_hybrid_core_regulator_comparison/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen self-consistent finite hybrid report changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    return {
        "S6_273_full_selfconsistent_finite_hybrid_and_ancestry_rebuilt": PARENT_SHA,
        "same_four_finite_hybrid_solutions_not_replacement_dynamics": True,
        "original_source_seed_parameters_reference_and_two_cutoffs_unchanged": True,
        "S261_refutation_and_S265_integrability_boundary_unchanged": True,
        "all_original_primitive_and_matching_frontiers_retained": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A full finite hybrid comparison gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.274.SAME_FULL_SOURCE_FINITE_HYBRID_POSITIVE_CORE_PROBABILITY_AND_COUPLED_CUTOFF_ORDERING_PHYSICAL_VOLUME_COMPARISONS",
        "date": "2026-09-15",
        "status": "SAME_FINITE_HYBRID_CORE_POVM_AND_REGULATOR_COMPARISONS; NOT_EXACT_SUPPORT_REGULATOR_REMOVAL_OR_ORIGINAL_V_G_B_P8",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/quadratic.md",
            "notes/domain.md",
            "notes/operators.md",
            "notes/dynamics.md",
            "notes/validation.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_source_full_quadratic_and_complete_nonlinear_remainder": serialize(
            {name: payload(packets[name]) for name in names[:3]}
        ),
        "whole_exact_operators_and_coupled_core_regulator_comparisons": serialize(
            {name: payload(packets[name]) for name in names[3:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "For the same four S273 finite hybrid solutions, the complete constrained quadratic Hamiltonian and its whole nonlinear Cauchy remainder give a sharper full centered amplitude1e172. The exact sixth antiheat identity retains its entire positive-heat remainder and both original orderings. On the same |u|<=1e-180 interval the evolved positive coherent outside-core probability is<1e-6. Fully coupled cutoff comparisons give five-coordinate, phase-factored state and physical-volume differences<1e-580,1e-12,1e-514; ordering comparisons give<1e-630,1e-65,1e-566. Every original scalar-center force, full trace/shear constraint, tensor curvature, heavy source, longitudinal Gauss term and regulator remainder remains. The complex analysis radius is not a changed physical cutoff. No exact support, un-factored global-state distance, homogeneous quantization, regulator removal, unlocalized Hamiltonian, unique or strict minimum, physical matching, omitted-loop/Regge/UV control or nonlinear global completion follows. Original V/G/B/P8 remain OPEN.",
        "not_established": [
            "Exact joint phase support or zero coherent core leakage",
            "A vector-norm comparison of un-factored states with different scalar phases",
            "Quantized homogeneous variables or the original fully quantum mean",
            "Unlocalized dynamics or a uniform mode, volume or cutoff-removal limit",
            "A unique or strict volume minimum, matching, omitted-loop/Regge/UV or global completeness",
            "Formal verification of the written analytic and operator proofs",
        ],
        "verification_boundary": "Exact full-source and curvature identities, outward rational envelopes, complete finite derivative evaluations and independent analytic/operator diagnostics support written proofs for the unchanged finite hybrid. Full ancestry is rebuilt. Native/direct/ordinary/CLI use original SymPy; only a captured complete regression may use the audited exact-GCD adapter. No successful historical replay changes archived physical errata or closes original frontiers.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete same-model hybrid comparison report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.274 same finite HYBRID core and regulator comparison replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
