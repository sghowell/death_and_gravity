"""Read-only leading gauge-Yukawa matching-screen checkpoint."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_two_loop_canonical import verify as parent

from . import audit, calibration, cuts, flow, model, threshold

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-gauge-yukawa-matching-screen.json"
PARENT_SHA = "362ca51c5f1b51733b0c12b5fce5e003834a3f2c0a748ec4926e4966b0d662cb"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_gauge_yukawa_screen/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen complete canonical two-loop report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_127_fully_rebuilt": PARENT_SHA,
        "old_scalar_model_and_scientific_reports_unchanged": True,
        "new_gauge_candidate_separately_named_not_an_old_amplitude_bound": True,
        "original_finite_EFT_positivity_contract_not_replaced_by_full_UV_construction": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows = audit.residuals()
    gates = {k: bool(v) for k, v in audit.gates().items()}
    if not all(v is True for v in gates.values()):
        raise ValueError("A gauge-Yukawa matching-screen gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.128.GAUGE_YUKAWA_LEADING_MATCHING_SCREEN",
        "date": "2026-09-10",
        "status": "POSITIVE_ONE_LOOP_MARGINAL_FIXED_FLOW_WITH_NONZERO_LEADING_HEAVY_FERMION_GAUGE_THRESHOLD_AND_OVERLAPPING_MASSLESS_CUT_SCREEN; NOT_FULL_FINITE_MASS_MATCHING_UV_NO_GO_V_G_B_OR_ORIGINAL_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/flow.md",
            "notes/threshold.md",
            "notes/cuts.md",
            "notes/bounds.md",
            "notes/literature.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entries(),
        "proof_checks": gates,
        "new_candidate_and_one_loop_marginal_flow": serialize(
            {"model": payload(model.data()), "flow": payload(flow.data())}
        ),
        "leading_heavy_fermion_threshold": serialize(payload(threshold.data())),
        "leading_matched_operator_gauge_cut": serialize(
            {
                "data": payload(cuts.data()),
                "interior": cuts.point(),
                "center": cuts.point(2),
            }
        ),
        "prospective_boundary_calibration_not_full_matching": serialize(
            payload(calibration.data())
        ),
        "controls": serialize(audit.controls()),
        "verdict": "A separately named unbroken SU(3) extension with fourteen vectorlike Dirac fundamentals has a positive explicit one-loop marginal fixed-flow ray. The two opposite scalar Yukawas cancel the odd gauge threshold but generate a nonzero leading Phi^2 F^2 coefficient. Independent physical-polarization and Lorentz contractions give the leading matched-operator massless gauge cut, first present at three parent loops. Its s and u cuts overlap the old massive scalar disc; the boundary mismatch vanishes at the center but not in its neighborhood. The tiny prospective coefficient cannot justify inheriting the original massive-gap dispersion hypotheses. This is a leading matching screen, not a full finite-mass on-shell amplitude, confinement statement, exclusion of the candidate or row, proof of UV consistency or P8 closure.",
        "not_established": [
            "A full finite-mass momentum or higher-loop threshold remainder for the new candidate",
            "Relevant-coupling running global quantum potential stability or nonperturbative ultraviolet completion",
            "A Yang-Mills confinement gap or its precise light spectrum",
            "A new-model complete forward-amplitude error budget or gapped dispersion theorem",
            "Finite-gravity IR/Regge Delta controlled common bounce parent or original P8 closure",
        ],
        "verification_boundary": "Exact SU(3) generator algebra, matrix/tensor one-loop beta-function reduction and explicit ODE ray, heavy flavor determinant and leading threshold dictionary, independent cut contractions and rational prospective scales. Renormalization formulae are sourced inputs; threshold/phase-space/cut arguments and scope limits are written proofs, not proof-assistant formalized or nonperturbative. All ancestor, own source and report fields are rebuilt read-only using unmodified scientific SymPy.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The gauge-Yukawa matching-screen report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.128.GAUGE_YUKAWA_LEADING_MATCHING_SCREEN replay passed; full V/G/B and original P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
