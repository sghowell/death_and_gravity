"""Read-only separate local two-scalar tree matching, not original-parent completion."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_scalar_transfer_moment import verify as parent
from p8_vacuum_canonical_affine_decoupling import verify as tree_input

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates" / "polynomial-vacuum-affine-heavy-scalar-tree-matching.json"
)
PARENT_SHA = "6a6c17275c0ed6c2dec6e65deb596fc264945168d2ebd058a53277a13d0057ff"
TREE_SHA = "b85c821b59bc41489418c0c190a2cc573cd9852cb19254627c3e2150d1dd6446"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_heavy_scalar_tree_matching/*.py"))
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
                "A source-pinned original tree, required spectral coefficient or physical frontier changed"
            )
        prior.validate_report(
            json.loads(prior.REPORT.read_text()), prior.build_report()
        )
    return {
        "S6_232_required_full_spectral_weight_higher_coefficient_and_separate_positive_atom_fully_rebuilt": PARENT_SHA,
        "S6_177_original_complete_canonical_tree_and_classical_decoupling_fully_rebuilt": TREE_SHA,
        "S6_231_original_identical_normalization_and_S6_182_retained_target_vacuum_transitively_rebuilt": True,
        "original_state_prescription_and_all_frozen_scientific_bytes_unchanged": True,
        "V2S_T1_is_separately_named_not_the_original_affine_parent": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError(
            "A separate-model potential, original tree remainder or first-elastic gate failed"
        )
    return {
        "schema": 1,
        "claim": "P8-S6.233.SEPARATE_LOCAL_CANONICAL_TWO_SCALAR_MODEL_WITH_COERCIVE_POTENTIAL_COMPLETE_ORIGINAL_MASSIVE_TREE_MATCHING_AND_UNIFORM_ALL_ANGLE_REMAINDER",
        "date": "2026-09-12",
        "status": "SEPARATE_V2S_T1_CLASSICAL_MODEL_AND_COMPLETE_TREE_MATCHING; NOT_FULL_QUANTUM_ORIGINAL_AFFINE_PARENT_BOUNCE_UV_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/model.md",
            "notes/amplitude.md",
            "notes/error.md",
            "notes/cut.md",
            "notes/scope.md",
            "notes/literature.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "separate_local_two_scalar_classical_model_and_potential": serialize(
            payload(audit.model.data())
        ),
        "complete_original_tree_remainder_and_full_first_elastic_comparison": serialize(
            {
                "matching": payload(audit.tree_matching.data()),
                "cut": payload(audit.cut.data()),
            }
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The separately named V2S-T1 model has canonical light/heavy scalar kinetic terms, light mass1, M_H^2=D+2 with D=2lambda/gamma, g=1/8192 and contact C=-g^2(3/D-2/D^2). Its full potential completes into a positive heavy square, phi^2/2 and a strictly positive remaining phi^4 term g^2(D-1)/[6D^2(D+2)], giving a coercive nonnegative classical potential and unique flat vacuum with positive masses. A second exact positive square gives the global bound V>=(phi^2+H^2)/6 at the actual D>2. Its literal exchange/contact tree amplitude equals the COMPLETE original S177/S182 massive scalar tree plus the exact rational remainder gamma*sum(channel-2)^4/[D-(channel-2)]. The massive shifted-cubic identity retains the potential term and the original coefficients; no original action is edited. Over all physical angles and4<=s<=10^196, the original tree is positive and the relative TREE error lies strictly between0 and1/60, using6r^2/(1-r) with r=(s-2)/D<32/625. The new tree matches b20=4lambda and b21=-3gamma, has b40=gamma^2/lambda and forward constant difference16gamma/(D+2); all higher terms and the heavy pole remain. Complete rational angular primitives, including the threshold limits and infinitely many even channels, give the first elastic comparison rho_first_original<=rho_first_V2S_T1<(61/60)^2 rho_first_original above threshold in the named window. Heavy single and pair production lie above that window, and the mixed two-body channel is forbidden by phi parity. This is a local classical model and a controlled complete-tree/first-elastic comparison, not a full quantum amplitude error, exact unitary S matrix, resonance-width calculation, high-energy arc theorem, original independent-function matching, covariant common-parent bounce construction, physical cutoff or original V/G/B/P8 closure.",
        "not_established": [
            "A mathematically constructed full quantum UV theory, physical mass/LSZ and real-loop matching, all-loop error or exact resonance width",
            "Identification of V2S-T1 with the original affine/Proca/DHOST parent or a match of its complete independent canonical functions and higher-field vertices",
            "A covariant common-parent field/state/derivative domain, controlled heavy stress and nonlinear bounce, finite-gravity forward or Regge control",
            "A physical cutoff, UV-parent exclusion, original V/G/B or P8 closure",
        ],
        "verification_boundary": "The complete potential square, canonical kinetic and mass entries, literal vertex factors, full original massive amplitude identity, rational all-angle margins and complete angular primitives are checked independently. Coercivity, uniform physical error, convergence of the even angular expansion and the first elastic comparison are written classical/perturbative-coefficient arguments, not FORMALIZED or full quantum theorems. All target predecessors remain frozen. Native/direct/ordinary/CLI retain original SymPy; only full regression uses the audited exact-GCD adapter.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The separate heavy-scalar tree matching report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.233 separate two-scalar tree matching replay passed; full quantum and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
