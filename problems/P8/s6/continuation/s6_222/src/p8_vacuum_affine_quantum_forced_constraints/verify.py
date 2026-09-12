"""Read-only full quantum-forced scalar interface and ordered graph-domain certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_scalar_tame_propagator import verify as parent

from . import audit, estimates, feedback, forces

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates" / "polynomial-vacuum-affine-quantum-forced-constraints.json"
)
PARENT_SHA = "7560a2ba6fcb470f784662046e2561d641bf4f6f0034104a158d977093075255"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_quantum_forced_constraints/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("A source-pinned classical propagator input changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_221_uniform_tame_classical_propagator_fully_rebuilt": PARENT_SHA,
        "S6_220_complete_scalar_Gaussian_response_and_both_contacts_transitively_rebuilt": True,
        "S6_182_current_fixed_QG1_and_original_quantum_preparation_unchanged": True,
        "no_frozen_scientific_proof_test_or_report_bytes_changed": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A full forced-constraint or domain gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.222.FULL_QUANTUM_FORCED_SCALAR_CONSTRAINT_INTERFACE_ORDERED_SCHUR_GRAPH_EQUIVALENCE_AND_FIRST_RESPONSE_DOMAIN",
        "date": "2026-09-12",
        "status": "EXACT_FORCED_CONSTRAINTS_AND_ORDERED_GRAPH_EQUIVALENCE; NOT_QUANTUM_INVERSE_CONTRACTION_FINITE_COUPLING_REMAINDER_STABILITY_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/forces.md",
            "notes/contacts.md",
            "notes/feedback.md",
            "notes/domains.md",
            "notes/bounds.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "complete_force_dependent_constraints_metric_interface_and_contacts": serialize(
            {
                "forces": payload(forces.data()),
                "contacts": payload(forces.contact_data()),
            }
        ),
        "ordered_causal_Schur_graph_and_precise_first_response_spaces": serialize(
            {"feedback": payload(feedback.data()), "bounds": payload(estimates.data())}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "Solving the full scalar coefficient action with independent normalized metric forces gives b=(pv-3gb)/2 and n=(Lc-gn-delta gzeta)/(2Jc), without Theta or inverse-transfer division. Both scalar momenta and weighted equations are retained. The exact phase and metric interfaces are Z'=KZ+Fg and s=CZ+Dg, with F=-Jcan C^T and a nonzero symmetric rank2 direct auxiliary block D. The full S220 Gaussian scalar response includes the ADM and additional clock contacts exactly once. Its normalization Qbar=(kappa a^3)^-1 Rhat is an output-density factor and cannot be commuted through the retarded kernel. Classical variation of constants yields [I-Qbar(D+C G0F)]g=e+Qbar C G0h, with a reversible reconstruction on the explicitly stated common graph domain. This does not establish that domain's density, invariance, surjectivity or a bounded quantum inverse. The force map has bound400(1+|P|^2), the direct block below100, and the normalized scalar response maps prepared H13_t H^(r+8) toL2_t H^(r-3) with bound1e-683. A prescribed first response has a classical phase output inC H^(r-17), with the exp(1e29) propagator constant; a comparison-phase source uses the unchanged finite C13 and loses27 spatial derivatives. The resulting scalar feedback bound gives onlyL2H^(r-19), not its required H13H^(r+8) source domain. Independent literal force saddles, complete local clock compositions, ordered two-time block solves and nonzero deleted-auxiliary/normalization controls supplement the written proof. A comparison multiplier shows why derivative-losing norm data alone neither establish a Neumann contraction nor exclude a different bounded inverse. No quantum state, finite counterterm, higher-derivative mode or physical cutoff is changed. Full quantum inversion, finite-coupling/nonlinear parent remainder, background/stability, heavy-sector control and original V/G/B/P8 remain open.",
        "not_established": [
            "Existence, uniqueness, boundedness or a contraction for the full quantum Schur inverse",
            "Density or feedback invariance of the stated graph domain from the available weak estimates",
            "A finite-coupling Born remainder, convergence, order reduction or full quantum canonical Hamiltonian",
            "Numerical high-time-jet stress bounds, independent quantum-state source maps or a literal global homogeneous solution",
            "Finite nonlinear parent/background stability, heavy/physical cutoff or original V/G/B/P8 closure",
        ],
        "verification_boundary": "Exact complete forced Legendre/constraint algebra, symplectic adjoint, rank2 direct response, clock composition, ordered block factorization and norm arithmetic are checked independently. Riesz extension, spatial Sobolev shifts, classical Duhamel and graph-domain equivalence are written proofs, not FORMALIZED. A finite-dimensional block fixture is not a continuum inverse certificate. Native/direct/ordinary/CLI retain original SymPy; only separately audited full regression uses the exact-GCD adapter. Frozen predecessors are unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The quantum-forced scalar graph report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.222 quantum-forced scalar interface replay passed; full quantum inverse and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
