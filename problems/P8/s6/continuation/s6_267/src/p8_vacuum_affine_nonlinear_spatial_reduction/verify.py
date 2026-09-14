"""Read-only full local spatial reduction with unchanged original ancestry."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_nonlinear_lapse_branch import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates/polynomial-vacuum-affine-nonlinear-spatial-reduction.json"
PARENT_SHA = "2736c2e1e5be41840a14973e3ec7fa70a20580072961a3f0663f14cd2499e663"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_nonlinear_spatial_reduction/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError(
            "The frozen complete nonlinear auxiliary-branch report changed"
        )
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    return {
        "S6_266_entire_nonlinear_auxiliary_report_and_ancestry_rebuilt": PARENT_SHA,
        "whole_original_sources_H_Proca_and_Gaussian_preparations_unchanged": True,
        "S265_Gaussian_integrability_boundary_and_S261_refutation_unchanged": True,
        "spatial_cotangent_and_translation_factor_not_full_quantum_closure": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A full spatial reduction or original-reference gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.267.COMPLETE_FULL_CONVOLUTION_NONLINEAR_SPATIAL_SHAPE_AND_SOURCE_PINNED_COTANGENT_MOMENTUM_REDUCTION_WITH_RESIDUAL_TRANSLATIONS_AND_UNCHANGED_PURE_REFERENCE_ZERO_CHARGE_FACTOR_WITH_UNCHANGED_PARENT",
        "date": "2026-09-14",
        "status": "FULL_LOCAL_NONLINEAR_SPATIAL_COTANGENT_REDUCTION_AND_PURE_NONZERO_TRANSLATION_FACTOR; NOT_NONLINEAR_QUANTUM_STATE_DOMAIN_EVOLUTION_CONTINUUM_REGULATOR_CUTOFF_OR_ORIGINAL_V_G_B_P8",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/shape.md",
            "notes/cotangent.md",
            "notes/translations.md",
            "notes/reference.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_full_nonlinear_shape_source_and_cotangent_reduction": serialize(
            {name: payload(packets[name]) for name in names[:3]}
        ),
        "whole_residual_translation_and_same_pure_reference_factor": serialize(
            payload(packets[names[3]])
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The explicit full-convolution weighted Fourier contraction constructs a positive unit-determinant divergence-free nonlinear spatial shape, including generated homogeneous corrections and its complete cubic remainder. The full original matter/vector/metric generator has a unique mean-zero cotangent momentum lift with every inverse, adjoint and boundary contact retained. The remaining three global translation constraints commute and are retained. The same finite nonzero pure Gaussian mode prescriptions give a common zero-charge reference factor, without conditioning or selecting a homogeneous state. This is a classical local reduction and an original linear-CCR reference symmetry result, not a nonlinear Gaussian pushforward, support theorem or interacting quantum construction. Every original source, state prescription, auxiliary-domain qualification and historical refutation remains. Original V/G/B/P8 remain OPEN.",
        "not_established": [
            "A globally unique diffeomorphism quotient or inverse uniform in torus radius",
            "A new homogeneous quantum state or two-polarization prescription at zero momentum",
            "Nonlinear transport or support of the original Gaussian inside the classical chart",
            "A self-adjoint interacting physical Hamiltonian, complete ordering or quantum measure",
            "A continuum gauge/BRST regulator, interacting physical mean or omitted-loop estimate",
            "A physical Wilsonian matching cutoff or original vacuum and finite-gravity gates",
            "A time-dependent nonlinear inhomogeneous bounce or global completeness theorem",
            "Formal verification of the written Banach, elliptic, cotangent and operator arguments",
        ],
        "verification_boundary": "Exact full source bindings, complete determinant/projector identities, all cotangent inverse/adjoint contacts, literal matter/vector translation generators and unrestricted coupled pure-Gaussian symmetry support the written local proof. Independent finite Fourier and matrix diagnostics are not the infinite-convolution theorem or a regulator for the full diffeomorphism algebra. Native/direct/ordinary/CLI use original SymPy; only captured full regression uses the audited exact-GCD adapter. Every frozen ancestor and both historical errata remain unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete nonlinear spatial reduction report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.267 full nonlinear spatial reduction and residual translations replay passed; nonlinear quantum domain, regulator and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
