"""Read-only complete finite-regulator volume turnaround certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_weyl_operator_comparison import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates/polynomial-vacuum-affine-finite-volume-turnaround.json"
PARENT_SHA = "970319d0dc460d46e9474e2ef9c17c1624aa7244b584461ac6cab1ff46771c7e"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_finite_volume_turnaround/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen finite Weyl operator comparison report changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    return {
        "S6_271_full_Weyl_operator_source_and_ancestry_rebuilt": PARENT_SHA,
        "original_source_parameters_state_and_two_cutoffs_unchanged": True,
        "S261_refutation_and_S265_integrability_boundary_unchanged": True,
        "S269_Omega_star_and_complex_formal_adjoint_clarifications_retained": True,
        "finite_turnaround_not_self_consistent_background_or_original_matching": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A whole finite volume turnaround gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.272.COMPLETE_FULL_SOURCE_LONGER_TRANSLATED_DOMAIN_AND_STATE_INDEPENDENT_FINITE_VOLUME_TURNAROUND_FOR_BOTH_OPERATOR_ORDERINGS_AND_TWO_ORIGINAL_CUTOFFS",
        "date": "2026-09-14",
        "status": "FINITE_REGULATOR_VOLUME_TURNAROUND_WITH_EXTERNAL_MEANS; NOT_UNIQUE_STRICT_OR_SELF_CONSISTENT_BOUNCE_OR_ORIGINAL_V_G_B_P8",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/reference.md",
            "notes/geometry.md",
            "notes/auxiliary.md",
            "notes/quantum.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_clock_adapted_source_and_full_phase_geometry": serialize(
            {name: payload(packets[name]) for name in names[:3]}
        ),
        "whole_translated_auxiliary_and_operator_volume_turnaround": serialize(
            {name: payload(packets[name]) for name in names[3:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The same finite48-pair full-source regulator has an evaluated real interval |u|<=1e-130. A fixed-balance reference bound, full scaled Fourier reconstruction, exact moving homogeneous clock center and complete translated auxiliary contraction give a full complex volume-symbol error<1e-320. The actual calibrated-coherent volume is within2e-320 ofI; the full Weyl operator is within1e-264 ofI and aboveI/2. For both full unitary orderings and both original cutoffs, the physical-volume mean endpoints exceed its center by>5e-260, forcing every global minimum strictly inside |u|<1e-132. The proof uses exact unitarity, not small state displacement. Every original source, state, scalar phase, primitive, implicit, spatial, matter/Gauss, covariance and translation contact remains. The homogeneous variables are external; no unique or strict bounce, self-consistent feedback, original interacting mean, continuum limit, matching, omitted loops or global completion follows. Original V/G/B/P8 remain OPEN.",
        "not_established": [
            "A unique volume minimum or strict positive second derivative",
            "Self-consistent homogeneous quantum mean or background equations",
            "Small longer-time state, cutoff, ordering-unitary or leakage errors",
            "An unregularized singular Hamiltonian or its interacting volume mean",
            "Uniform mode-count, torus-volume or regulator-removal limits",
            "Original physical matching, omitted loops, UV or global completion",
            "Formal verification of the written analytic and operator proofs",
        ],
        "verification_boundary": "Exact full clock identities, complete source derivatives, outward rational estimates and independently tested geometric/operator identities support the written theorem. The same parent and all frozen ancestor reports are rebuilt. Independent numerical and algebraic fixtures are diagnostics, not a finite simulation defining the full quantum system. Native/direct/ordinary/CLI use original SymPy; only a captured full regression may use the audited exact-GCD adapter. Both historical errata and every original open frontier remain.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete finite volume turnaround report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.272 complete finite volume turnaround replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
