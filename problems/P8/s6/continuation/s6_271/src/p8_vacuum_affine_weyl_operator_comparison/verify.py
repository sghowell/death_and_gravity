"""Read-only finite Weyl operator ordering, positivity and unitary comparison."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_quantitative_local_time import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates/polynomial-vacuum-affine-weyl-operator-comparison.json"
PARENT_SHA = "5b6adbb4be256ec2df0051d781d43a2e94b8ab12d099c79295e87c821fb0fc01"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_weyl_operator_comparison/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen full finite-time regulator report changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    return {
        "S6_270_entire_finite_time_source_domain_and_ancestry_rebuilt": PARENT_SHA,
        "original_source_parameters_and_pure_nonzero_preparations_unchanged": True,
        "S261_refutation_and_S265_Gaussian_integrability_boundary_unchanged": True,
        "S269_Omega_star_clarification_and_complex_formal_adjoint_retained": True,
        "finite_Weyl_ordering_comparison_not_original_physical_matching": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A complete Weyl operator or unitary comparison gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.271.COMPLETE_EXPLICIT_FULL_WEYL_OPERATOR_NORM_ORDERING_COMPARISON_POSITIVE_PHYSICAL_VOLUME_AND_ENTIRE_UNITARY_ORDERING_AND_TWO_CUTOFF_REGULATOR_COMPARISON_WITH_UNCHANGED_PARENT",
        "date": "2026-09-14",
        "status": "EXPLICIT_FINITE_WEYL_OPERATOR_ORDERING_POSITIVE_VOLUME_AND_UNITARY_COMPARISON; NOT_ORIGINAL_PHYSICAL_MATCHING_CONTINUUM_LOOP_OR_V_G_B_P8",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/kernel.md",
            "notes/derivatives.md",
            "notes/comparison.md",
            "notes/dynamics.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_full_source_and_explicit_Weyl_operator_bound": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_high_derivatives_and_evaluated_ordering_unitary_comparison": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The same finite48-pair full-source regulator admits an explicit Gaussian-frame Schur bound with all196 mixed phase/cutoff derivatives retained. The complete calibrated-coherent to Weyl volume OPERATOR error is<1e-200, the actual Weyl volume is>1/2 and within1e-199 ofI. The full real bounded Weyl Hamiltonian generates an exact unitary; its OPERATOR difference from the calibrated-coherent evolution is<1e-944 on |u|<=1e-2000. Two defined Weyl cutoffs have state difference<1e-943 and readout mean difference<1e-199. Every original pure-state, source, scalar phase, time, primitive, nonlinear implicit, spatial, matter/vector and residual-translation contact remains. No original interacting mean, unregularized Hamiltonian, physical Wilsonian matching, uniform continuum limit, omitted loops or global completeness follows. Original V/G/B/P8 remain OPEN.",
        "not_established": [
            "Positivity-preserving Weyl quantization for arbitrary positive symbols",
            "An original interacting volume mean or unregularized singular Hamiltonian",
            "Bounds for arbitrary cutoffs outside the two explicitly differentiated choices",
            "Uniform mode-count, torus-volume or physical regulator-removal limits",
            "A homogeneous quantum state or literal identity with the R3 state",
            "Original physical matching, omitted loops, all-energy UV positivity or global completeness",
            "The stronger old coherent-only cutoff bounds unchanged for Weyl orderings",
            "Formal verification of the written analytic and operator proofs",
        ],
        "verification_boundary": "Direct normalized Gaussian integration, every product-Schur weight, exact high-derivative induction/partition/Cauchy inequalities and all196 growth ratios support an identified written OPERATOR theorem. This is not a symbol-sup shortcut or an unknown-dimensional constant. The high derivatives are in PHASE, not unproved real-time profiles. Full exact source bindings and ancestor reports remain. Numerical kernel, quadrature, Schur and unitary fixtures are independent diagnostics, not the full P8 quantum evolution. Native/direct/ordinary/CLI use original SymPy; only the captured full regression uses the audited exact-GCD adapter. Both historical errata and the original open frontiers are unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete finite Weyl operator comparison report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.271 complete finite Weyl operator and unitary comparison replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
