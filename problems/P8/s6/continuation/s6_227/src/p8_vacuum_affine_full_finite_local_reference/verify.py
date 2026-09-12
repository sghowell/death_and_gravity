"""Read-only actual full finite local reference with ordered causal inverse."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_curved_scalar_reference import verify as local_input
from p8_vacuum_affine_weighted_channel_resolvent import verify as parent

from . import audit, coefficients, factorization, resolvent

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates" / "polynomial-vacuum-affine-full-finite-local-reference.json"
)
PARENT_SHA = "98ac59825e67545f5afd2ce5a430d0c381f49811d85ee5bf805a3939b8c9de3d"
LOCAL_SHA = "5c4a3649bc435475d6bfd65901af47b63da2543faa9bc9bcd6364e00fac9a02d"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_full_finite_local_reference/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    for prior, expected in ((parent, PARENT_SHA), (local_input, LOCAL_SHA)):
        if sha(prior.REPORT) != expected:
            raise ValueError(
                "A source-pinned original local Hessian or weighted factor input changed"
            )
        prior.validate_report(
            json.loads(prior.REPORT.read_text()), prior.build_report()
        )
    return {
        "S6_226_original_channel_log_coercivity_and_weighted_reference_fully_rebuilt": PARENT_SHA,
        "S6_225_actual_complete_original_finite_local_Hessian_fully_rebuilt": LOCAL_SHA,
        "S6_223_original_shear_pole_and_S6_86_trace_factor_transitively_rebuilt": True,
        "S6_222_actual_full_quantum_force_graph_transitively_rebuilt_not_inverted_here": True,
        "original_state_prescription_and_all_frozen_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError(
            "A complete local factorization, coefficient or inverse gate failed"
        )
    return {
        "schema": 1,
        "claim": "P8-S6.227.ACTUAL_COMPLETE_FINITE_LOCAL_HESSIAN_CURVATURE_CONJUGATION_UNIFORM_BOUND_AND_FULL_FINITE_REFERENCE_INVERSE",
        "date": "2026-09-12",
        "status": "ACTUAL_FULL_FINITE_LOCAL_REFERENCE_CAUSAL_INVERSE; NOT_FULL_NONLOCAL_CURVED_QUANTUM_SYSTEM_UNWEIGHTED_SMALLNESS_STABILITY_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/factorization.md",
            "notes/coefficients.md",
            "notes/reciprocity.md",
            "notes/conjugation.md",
            "notes/inverse.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "actual_complete_local_Hessian_factorization_and_coefficients": serialize(
            {
                "factorization": payload(factorization.data()),
                "coefficients": payload(coefficients.data()),
            }
        ),
        "source_time_reciprocity_conjugated_remainder_and_full_finite_inverse": serialize(
            payload(resolvent.data())
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The complete original S225 finite local remainder has the exact Euler-Hessian identity Rloc=Bc*A+A^T Bc+D*HD+D*J+J^T D+V0, with the coefficient-derivative matrix J and the complete density boundary retained. The source-pinned alpha/gamma coefficients include all original background-curvature, second-metric-variation, Einstein and volume terms. Actual whole-slab bounds give ||A||<110m^2, ||H||<185m^2, ||J||<286m^2 and ||V0||<160m^4 at the unchanged m=1000. Retarded/advanced reciprocity and reversed-time primal Volterra bounds give ||Z(t,s)||<=(5/2)elapsed exp(20elapsed) and the SOURCE-time derivative bound ||partial_s Z||<=(5/2)exp(20elapsed). The ordinary kernel of ZD* is partial_s Z; its initial-delta cancellation is explicit. The six correctly ordered terms of Omega=Z Rloc Y therefore have weighted same-H^r norm at most6825m^2/[4(sigma-20)^2]+3575m^2/(sigma-20)^3+1000m^4/(sigma-20)^4. At sigma_M=2m exp(32+8M), this is below1/4 uniformly in all spatial momenta and every realr. Including the optional S226 channel matrix V of norm at mostM gives a convergent causal middle inverse of Fdiag+V+Omega with norm at most20/(123+32M). ComposingY andZ gives the full finite local reference inverse with norm at most125/[(123+32M)(sigma_M-20)^4]. Both ordered distribution graph inverse identities retain the full original finite local Hessian, original shear pole and initial boundary. The physical conformal density stays on the right and the bound retains kappa; removing the time weight costs exp(sigma_M T). Independent arbitrary-coefficient Euler/boundary identities, original coefficient bounds, exact constant-coefficient and nonconstant sixth-jet reciprocity, initial-atom omission controls and full noncommuting finite block solves supplement the written functional proof. Actual full nonlocal curved mass/state/contact/tree/matter matching and the S222 coupled inverse remain open, as do physical stability and original V/G/B/P8.",
        "not_established": [
            "Identification or compatible bounds for the complete actual nonlocal curved mass/state/contact/classical/matter response",
            "A full actual S222 auxiliary/clock graph inverse or finite-coupling/nonlinear common-parent remainder",
            "Unweighted small response, quantum stability, absence of growing poles or a physical cutoff",
            "Original V/G/B or P8 closure; no original local finite term, pole, state or preparation was changed",
        ],
        "verification_boundary": "Exact full local density/Euler factorization, original finite coefficient bridge, all interval and convolution constants, source-time initial-boundary identity and ordered inverse algebra are checked. Retarded/advanced reciprocity, Volterra bounds, causal distribution graph composition and weighted functional analysis are written proofs, not FORMALIZED. Finite causal matrices and finite Taylor jets are independent diagnostics, not continuum discretization certificates. Native/direct/ordinary/CLI retain original SymPy; only full regression uses the audited exact-GCD adapter. Frozen predecessors are unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The full finite local reference report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.227 full finite local reference replay passed; actual full quantum inverse and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
