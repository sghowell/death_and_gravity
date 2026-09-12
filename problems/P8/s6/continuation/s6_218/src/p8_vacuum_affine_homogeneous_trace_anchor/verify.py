"""Read-only full homogeneous trace anchor and original finite-volume matching."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_corrected_spatial_current import verify as parent
from p8_vacuum_affine_proca_gaussian import verify as gaussian

from . import anchor, audit, comparison, homogeneous, local

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates" / "polynomial-vacuum-affine-homogeneous-trace-anchor.json"
)
PARENT_SHA = "db4e6bbce7c09d54ef1705581f898bc10fecc5fe860e28c510e587b106fe0a0a"
GAUSSIAN_SHA = "1175702d6980db090396863b807af0e767241daeb7d8dabf138ed11c67e1ae75"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_homogeneous_trace_anchor/*.py"))
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
        "S6_217_corrected_full_tracefree_current_fully_rebuilt": PARENT_SHA,
        "S6_176_original_Gaussian_state_transitively_rebuilt": GAUSSIAN_SHA,
        "same_actual_reference_Hessian_not_full_nonlinear_source_free_parent": True,
        "no_frozen_scientific_test_or_report_bytes_changed": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A full homogeneous trace anchor gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.218.FULL_HOMOGENEOUS_TRACE_GAUSSIAN_ANCHOR_ORIGINAL_VOLUME_MATCHING_AND_ORDERED_CROSS_ZEROS",
        "date": "2026-09-12",
        "status": "COMPLETE_HOMOGENEOUS_TRACE_ANCHOR_AND_FULL_HOMOGENEOUS_HILBERT_LIFT; NOT_NONZERO_TRANSFER_SCALAR_REMAINDERS_REDUCED_INVERSE_FULL_SOURCED_PARENT_BACKGROUND_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/hamiltonian.md",
            "notes/domination.md",
            "notes/comparison.md",
            "notes/matching.md",
            "notes/anchor.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "full_trace_Hamiltonian_domination_and_original_comparison": serialize(
            {
                "Hamiltonian": payload(homogeneous.data()),
                "domination": payload(homogeneous.domination_data()),
                "comparison": payload(comparison.data()),
            }
        ),
        "general_dimension_volume_current_and_actual_homogeneous_anchor": serialize(
            {
                "local": payload(local.data()),
                "anchor": payload(anchor.data()),
            }
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The full constrained homogeneous trace family Q=(2/3)phi I has A=a exp(phi/3), varying volume and longitudinal constraint. Its exact scalar T/T/L normalized current vertices and all mixed frequency/squeeze jets satisfy the hypotheses of the existing full finite-reference and actual-state comparison estimates. These hypotheses, the unchanged initial graph and complete second contact are checked explicitly; no tracefree Hamiltonian is simply relabeled. Independent arbitrary-dimensional scalar K/omega differentiation and radial continuation reproduce the original fixed local invariant action with the full A^d Jacobian. The physical finite local trace Euler current includes the nonzero m^4 volume term, has all five formal density Green identities, and its first three raw amplitude bounds are below1e11. The complete Gaussian determinant trace current has coordinate-density C0,C1,C2 displays(1e78,1e95,2e111). Its first derivative is the actual reference P0 trace/trace metric Hessian by finite-mode tangent uniqueness, full contact and S176's vanishing reference source/first source variation. Both ordered trace/tracefree cross anchors vanish by rotation covariance. Normalization to I/sqrt3 and the orthogonal initial-germ Hilbert lift give the full homogeneous anchor bound1e95||D||L2 Z130, with canonical displays3e-705 for unit trace and4e-705 for the full anchor. Independent literal ADM matrices, full covariance-flow derivatives, scale-factor Euler variations, full-chart integrated Hessians and deleted-constraint/volume controls supplement the written proof. Finite-amplitude C2 is only a determinant-component result, not a remainder for the full nonlinear sourced parent. The three nonzero-transfer scalar currents still need their full state/time/endpoint/contact/all-transfer assembly; reduced inverse, nonlinear quantum background and original V/G/B remain open.",
        "not_established": [
            "Complete three nonzero-transfer scalar currents with state/time/endpoint/contact and all-transfer regulator/dimension remainders",
            "A finite-amplitude C2 remainder for the full affine sourced parent rather than its Gaussian determinant component",
            "Equality of the two nonzero-transfer ordered cross retarded kernels from their zero-transfer anchors",
            "A same-space reduced inverse, nonlinear quantum background or stability from derivative-losing canonical anchor bounds",
            "Remaining parent/loop/heavy/cutoff, original V/G/B or original P8 closure",
        ],
        "verification_boundary": "Exact full trace Hamiltonian and current vertices, mixed-jet dominance, scalar K/omega local action, general-dimensional matching, full volume Euler current, local majorants, tangent contacts, normalization and scope are checked. Independent literal ADM/covariance and scale-factor variational tests do not import the formulas they test. Continuous all-momentum, dimensional, C2 dominated differentiation, rotation and Hilbert arguments are written proofs, not FORMALIZED. Native/direct/ordinary/CLI retain original SymPy; only full regression uses the audited exact-GCD adapter. Frozen predecessors remain unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The homogeneous trace anchor report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.218 full homogeneous trace anchor replay passed; nonzero-transfer scalar remainders and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
