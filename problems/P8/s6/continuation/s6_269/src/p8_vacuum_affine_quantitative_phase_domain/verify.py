"""Read-only evaluated finite bounce domain and same-seed volume comparison."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_local_quantum_regulator import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates/polynomial-vacuum-affine-quantitative-phase-domain.json"
PARENT_SHA = "410af7db450afb132e77a4e2c696ddd3a027a78c324becb3209c78509a4b1df0"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_quantitative_phase_domain/*.py"))
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
            "The frozen explicit finite local quantum regulator report changed"
        )
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    return {
        "S6_268_entire_explicit_finite_quantum_regulator_and_ancestry_rebuilt": PARENT_SHA,
        "original_source_parameters_and_pure_nonzero_preparations_unchanged": True,
        "S261_refutation_and_S265_Gaussian_integrability_boundary_unchanged": True,
        "evaluated_finite_bounce_comparison_not_original_continuum_matching": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError(
            "A full quantitative phase-domain or same-seed volume gate failed"
        )
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.269.COMPLETE_EVALUATED_FULL_NONLINEAR_BOUNCE_PHASE_DOMAIN_ALL_TWELVE_INVARIANT_IMAGES_SAME_PURE_REFERENCE_COHERENT_TAIL_AND_POSITIVE_FINITE_VOLUME_REGULATOR_COMPARISON_WITH_UNCHANGED_PARENT",
        "date": "2026-09-14",
        "status": "EVALUATED_FINITE_BOUNCE_PHASE_DOMAIN_SAME_SEED_COHERENT_TAIL_AND_POSITIVE_VOLUME_COMPARISON; NOT_TIME_EVOLUTION_ORDERING_MATCHING_REGULATOR_REMOVAL_OR_ORIGINAL_V_G_B_P8",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/reference.md",
            "notes/geometry.md",
            "notes/invariants.md",
            "notes/quantum.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_source_reference_and_evaluated_canonical_field_bounds": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_nonlinear_geometry_invariant_images_and_quantum_comparison": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "For the explicitly declared L=1, P=1e64, all-eight-channel finite torus family at the bounce, the entire full-covariance whitened ball of radius2e20 reconstructs complete nonlinear spatial fields whose twelve original auxiliary invariants are all below1e-260. Every source, density, adjoint ghost, curvature and generated Fourier contact remains. The SAME pure reference has coherent core leakage below exp(-1e39), strictly positive rather than zero. The positive finite coherent physical volume has operator distance fromI below1e-255. Two admissible same-core volume cutoffs have evaluated exponentially small seed vector and mean differences. This is a finite bounce-slice regulated observable, not evaluated interaction/ordering dynamics, a new original action or state, an original quantum mean or a continuum limit. Original V/G/B/P8 remain OPEN.",
        "not_established": [
            "An evaluated positive common time interval or later-time state leakage",
            "The complete Hamiltonian interaction norm or first-Weyl ordering error",
            "Positivity or original matching of the separately calibrated volume",
            "An original interacting physical quantum mean or unregularized Hamiltonian",
            "A homogeneous quantum state or identity with the original R3 state",
            "Uniform mode-count/volume estimates, regulator removal or loop bounds",
            "Original Wilsonian matching, finite-gravity gates or nonlinear completeness",
            "Formal verification of the written Banach and operator arguments",
        ],
        "verification_boundary": "Exact source/canonical bindings, complete correlated-state row normalization, generic full adjoint and non-diagonal trace contacts, rational outward field/invariant/volume budgets and full gamma-tail identities support the independent written proofs. Galerkin inverse, covariance, radial integration and arbitrary-precision tail checks are diagnostics, not a full diffeomorphism regulator or original time simulation. Native/direct/ordinary/CLI use original SymPy; only the captured full regression uses its audited exact-GCD adapter. Frozen ancestors and both historical errata remain unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete quantitative phase-domain report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.269 evaluated finite bounce phase-domain and volume comparison replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
