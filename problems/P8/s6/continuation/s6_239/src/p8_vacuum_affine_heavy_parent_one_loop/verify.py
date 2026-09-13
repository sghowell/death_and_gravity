"""Read-only complete formal first-loop matching of the new classical limit."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_heavy_scalar_loop_coefficients import verify as polynomial_input
from p8_vacuum_affine_heavy_scalar_one_loop import verify as self_input
from p8_vacuum_affine_heavy_scalar_parent import verify as parent

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-heavy-parent-one-loop.json"
PARENT_SHA = "934b6c035fb5de0eb89ee1e1a03c2624c78a33985c0503bce883c41240b05905"
POLYNOMIAL_SHA = "6d76badc63d988c9fae045c5a3628b20c7fb476ea9f45fdf98fd52a7a1b0e9fb"
SELF_SHA = "d6646d25d158bcfeb32bf07895e8fcc2ca08e05988919f803ce50aa5e45f58a5"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_heavy_parent_one_loop/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    for prior, expected in (
        (parent, PARENT_SHA),
        (polynomial_input, POLYNOMIAL_SHA),
        (self_input, SELF_SHA),
    ):
        if sha(prior.REPORT) != expected:
            raise ValueError(
                "A frozen full parent or separate polynomial loop input changed"
            )
        prior.validate_report(
            json.loads(prior.REPORT.read_text()), prior.build_report()
        )
    return {
        "S6_238_literal_full_new_classical_parent_fully_rebuilt": PARENT_SHA,
        "S6_237_separate_polynomial_complete_loop_coefficients_fully_rebuilt": POLYNOMIAL_SHA,
        "S6_234_full_source_onepoint_and_on_shell_prescription_fully_rebuilt": SELF_SHA,
        "all_frozen_original_and_separate_model_scientific_bytes_unchanged": True,
        "new_loop_of_classical_limit_not_quantum_limit_interchange": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A complete-loop, matching, counterterm or scope gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.239.COMPLETE_FORMAL_FIRST_LOOP_OF_H8A420_CLASSICAL_LIMIT_WITH_ALL_HIGHER_VERTEX_CONTRACTIONS_NEW_OS4_MATCHING_FULL_ANGLE_AND_FORWARD_COEFFICIENT_BOUNDS_AND_SEPARATE_FINITE_CLOCK_EXTENSION_BUDGET",
        "date": "2026-09-12",
        "status": "FORMAL_COMPLETE_FIRST_LOOP_OF_CLASSICAL_LIMIT_WITH_SCOPED_ERRORS; NOT_QUANTUM_DECOUPLING_ALL_LOOPS_CURVED_BOUNCE_UV_REGGE_ORIGINAL_V_G_B_OR_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/germs.md",
            "notes/contractions.md",
            "notes/renormalization.md",
            "notes/matching.md",
            "notes/scope.md",
            "notes/literature.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "literal_complete_action_germs_graph_counting_and_one_loop": serialize(
            {"germs": payload(audit.germs.data()), "loops": payload(audit.loops.data())}
        ),
        "complete_matching_bounds_and_finite_clock_extensions": serialize(
            payload(audit.bounds.data())
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The complete formal one-loop four-Phi coefficient of S238's fixed-canonical CLASSICAL limiting Phi/H action is derived, with its free tensor, M1 and Proca spectators retained. Exact connected graph counting includes all six local degree-six tadpoles and all cubic/five-field contractions: the full-source onepoint cancellation, the tadpole-dressed heavy exchange and the mixed3+1 bubble. Dimensional tensor averaging is performed before pole-only MSbar, retaining its finite evanescent terms. A NEW full-parent OS4 value prescription leaves an explicit full rational exchange correction and two local symmetric polynomials. The complete first loop is below10^-199 of the original full tree throughout2<=E<=10^98 and all physical angles, and the combined tree plus first loop stays within1/59. The complete forward shifts retain0<delta_b20/(4lambda)<10^-203 and0<delta_b21/(-3gamma)<10^-203, while abs(delta_b40)/(gamma²/lambda)<10^-192 with no sign assigned. The nonzero new b40 exchange correction is included. The new positive finite contact fits the separate contact-only global potential margin. Explicit finite contact, quadratic and onepoint extensions by1-T preserve the relevant vacuum jets and have a separately proved full mixed four-jet bound10^-400 on the clock tube. None of this interchanges quantization with gravitational decoupling, controls omitted loops or the full curved quantum action/state/response, supplies a nonlinear same-state bounce, proves physical UV or finite-gravity Regge completion, or closes original V/G/B/P8.",
        "not_established": [
            "A regulator-removed quantum decoupling limit of the finite-gravity parent or a proof that quantization and the classical limit commute",
            "Omitted-loop control, an exact heavy resonance, quantum vacuum theorem or physical UV completion",
            "A full curved UV counterfunctional, specified new heavy state, renormalized curved stress/response or same-state nonlinear bounce",
            "The classical10^-2700 tube bound for the separately specified finite counterterms; their new budget is10^-400",
            "Finite-gravity Regge estimates, a physical cutoff, all-parent exclusion or original V/G/B/P8 closure",
        ],
        "verification_boundary": "Literal full-action jets, complete graph patterns, exhaustive720-permutation dimensional contractions and exact rational matching bounds are accompanied by written continuum proofs. Independent multilinear plane-wave coefficient extraction with exact-moment signed sphere cubature checks all local vertices in dimensions3,4,5. Its numerical evaluation is a diagnostic, not an interval proof. Native/direct/ordinary/CLI retain original SymPy; only full regression uses the audited exact-GCD adapter. No frozen scientific or prescription bytes are changed.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete new-parent formal first-loop report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.239 complete classical-limit first-loop replay passed; quantum parent and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
