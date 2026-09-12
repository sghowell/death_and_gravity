"""Read-only original transfer derivative and mandatory spectral-weight requirement."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_scalar_matching_scale import verify as parent
from p8_vacuum_canonical_affine_decoupling import verify as tree_input

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-scalar-transfer-moment.json"
PARENT_SHA = "34da90ca48a00689b3b9d7db12b367db703c970236caabcecd293b1061559cd5"
TREE_SHA = "b85c821b59bc41489418c0c190a2cc573cd9852cb19254627c3e2150d1dd6446"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_scalar_transfer_moment/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    for prior, expected in ((parent, PARENT_SHA), (tree_input, TREE_SHA)):
        if sha(prior.REPORT) != expected:
            raise ValueError(
                "A source-pinned original scalar amplitude, normalization or physical frontier changed"
            )
        prior.validate_report(
            json.loads(prior.REPORT.read_text()), prior.build_report()
        )
    return {
        "S6_231_original_identical_scalar_normalization_unitarity_and_conditional_forward_matching_fully_rebuilt": PARENT_SHA,
        "S6_177_complete_fixed_canonical_classical_decoupling_and_original_scalar_tree_fully_rebuilt": TREE_SHA,
        "S6_182_same_retained_vacuum_and_unchanged_nonconstant_low_jets_transitively_rebuilt": True,
        "original_state_prescription_and_all_frozen_scientific_bytes_unchanged": True,
        "positive_rational_atom_not_substituted_for_original_parent": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError(
            "An original transfer normalization, positive moment or diagnostic-scope gate failed"
        )
    return {
        "schema": 1,
        "claim": "P8-S6.232.ORIGINAL_SCALAR_TRANSFER_MATCHING_MANDATORY_LOW_ABSORPTIVE_WEIGHT_AND_HIGHER_FORWARD_COEFFICIENT_WITH_POSITIVE_SPECTRAL_CONTROL",
        "date": "2026-09-12",
        "status": "EXACT_ORIGINAL_TRANSFER_DATA_AND_CONDITIONAL_FULL_SPECTRAL_WEIGHT_AND_HIGHER_COEFFICIENT_REQUIREMENT; NOT_FULL_LOOP_BOUND_PHYSICAL_CUTOFF_EXACT_UNITARY_COMPLETION_UV_NO_GO_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/transfer.md",
            "notes/split.md",
            "notes/first-elastic.md",
            "notes/positive-pole.md",
            "notes/scope.md",
            "notes/literature.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "original_crossing_transfer_data_and_first_elastic_derivative": serialize(
            payload(audit.coefficients.data())
        ),
        "mandatory_full_moment_first_elastic_bound_and_separate_positive_atom": serialize(
            {
                "moment": payload(audit.moment.data()),
                "atomic_control": payload(audit.atomic.data()),
            }
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "For the same original S182 vacuum and complete S177 massive scalar tree, the crossing coefficient is exactly b2_tree(t)=4lambda-3gamma*t. The complete first elastic nonforward coefficient from the original normalized identical channels has transfer slope beta*b(S)^2/[60pi(S-4)], with all massive threshold and potential factors retained. For any actual scalar S matrix satisfying the explicitly unproved physical-mass, spectral-gap, positive endpoint derivative and differentiated fixed-t dispersion/arc premises, an arbitrary split K=M^2>4 gives b21+3b20/[2(K-2)]>=-3J4_low(K)/2. All low spectral weight and any retained positive heavy atoms remain. Under physical coefficient errors |b20-4lambda|<=4lambda*delta0 and |b21+3gamma|<=3gamma*delta1, this requires J4_low>=2gamma(1-delta1)-4lambda(1+delta0)/(K-2). At M10^99, delta0<=1 and delta1<=1/2 require FULL J4_low>gamma/5. Independently, the complete original first elastic contribution obeys0<J4_first<[48lambda^2K+(9/16)gamma^2K^3+gamma^2]/pi^2<gamma*10^-200. Thus the full low moment must exceed that contribution by more than2*10^199, although the original tree amplitude is below10^-202 on that energy window. No small full-loop error is assumed to derive this required enhancement. Independently, the positive full spectral Gram matrix gives J4_total^2<=b20*b40, while positive transfer weight gives J4_total>=-2*b21/3 for b21<0. Therefore b40>=gamma^2(1-delta1)^2/[lambda(1+delta0)], including the named lower gamma^2/(8lambda)>0 and the exact-match lower gamma^2/lambda. The original tree has b40=0, so a full matching model must supply this higher coefficient. A separately named positive rational pole with M_H^2=2+2lambda/gamma and g^2=1/2^26 matches these two low coefficients and supplies J4=2gamma, while its v4 coefficient gamma^2/lambda is nonzero. That rational tree is not an exact unitary amplitude or a matched parent; it only demonstrates compatibility and Gram-bound saturation at the level of positive spectral data, not physical unitary optimality. The split is not a physical cutoff, the first elastic upper bound is not a full absorptive bound, and no pole, state, regulator or original parent action is changed. Full quantum matching, finite-gravity control, nonlinear common-parent bounce estimates and original V/G/B/P8 remain OPEN.",
        "not_established": [
            "The original full quantum amplitude, physical mass/LSZ matching, absence of additional low cuts, positive derivative convergence or the differentiated fixed-t dispersion and arc conditions",
            "Physical b20 or b21 error bounds, a full absorptive calculation, control of higher loops or heavy thresholds, or dominance of the first elastic contribution",
            "An exact unitary UV completion or a matched common parent from the separate positive rational pole; its full higher coefficients differ from the original tree",
            "A physical cutoff, UV-parent exclusion, finite-gravity forward/Regge bound, full curved nonlinear or bounce control, original V/G/B or P8 closure",
        ],
        "verification_boundary": "Original crossing and first elastic coefficients, positive partial-wave endpoint factors, exact low/high moment decomposition, rational inequalities, conservative massive integral bounds and the separate atom identities are replayed independently. The continuum endpoint convergence and differentiated dispersion implication have explicit unproved physical premises and written conditional proofs, not FORMALIZED. Finite positive-measure fixtures and the rational atom are diagnostics, not S matrices or numerical continuum certificates. Native/direct/ordinary/CLI retain original SymPy; only full regression uses the audited exact-GCD adapter. No frozen predecessor changes.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The original transfer spectral-moment report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.232 original transfer and conditional spectral-weight replay passed; full matching and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
