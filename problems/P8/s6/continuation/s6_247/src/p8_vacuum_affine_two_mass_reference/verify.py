"""Read-only complete unchanged-scheme two-mass Gaussian flat reference inverse."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_band_coupled_inverse import verify as matching_input
from p8_vacuum_affine_flat_scalar_quotient_inverse import verify as quotient_input
from p8_vacuum_affine_heavy_adm_clock_response import verify as parent
from p8_vacuum_affine_heavy_curved_state import verify as state_input

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-two-mass-reference.json"
PARENT_SHA = "c1b11dd8afe198c5a176e0e66b64065fd1a3a4644f6f8054cf80e9ec38d20123"
STATE_SHA = "d4a6081d0482b2168255299cdd0080c90a1716eeffbc70f2645e6e2452e78eda"
QUOTIENT_SHA = "a2287a5083d5dceab6886f4f7de3304ce3a7af56c772e76dc869115ad8ac3fb5"
MATCHING_SHA = "076648eee2112f5b86302a072e11265350e155a92877f5b68b5007ab0cf09432"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_two_mass_reference/*.py"))
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
        (state_input, STATE_SHA),
        (quotient_input, QUOTIENT_SHA),
        (matching_input, MATCHING_SHA),
    ):
        if sha(prior.REPORT) != expected:
            raise ValueError(
                "A frozen actual Gaussian response, prescription or reference geometry changed"
            )
        prior.validate_report(
            json.loads(prior.REPORT.read_text()), prior.build_report()
        )
    return {
        "S6_246_complete_actual_heavy_ADM_common_clock_response_fully_rebuilt": PARENT_SHA,
        "S6_240_actual_heavy_mass_and_full_covariant_finite_prescription_fully_rebuilt": STATE_SHA,
        "S6_224_original_complete_projector_geometry_and_causal_quotient_fully_rebuilt": QUOTIENT_SHA,
        "S6_229_original_full_dimensional_Proca_leading_match_fully_rebuilt": MATCHING_SHA,
        "all_frozen_original_and_separate_model_scientific_bytes_unchanged": True,
        "old_Proca_inverse_not_transferred_to_the_sum_or_new_curved_model": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A complete two-mass reference inverse gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.247.COMPLETE_UNCHANGED_SCHEME_TWO_MASS_GAUSSIAN_FLAT_REFERENCE_WITH_FULL_SCALAR_CUT_FIXED_FINITE_CURVATURE_TWO_THRESHOLD_CUT_ONLY_RECIPROCAL_AND_UNIFORM_CAUSAL_QUOTIENT_INVERSE",
        "date": "2026-09-13",
        "status": "COMPLETE_TWO_MASS_GAUSSIAN_FLAT_QUOTIENT_REFERENCE_INVERSE_WITH_NO_FIRST_SHEET_ISOLATED_POLE; NOT_ACTUAL_CURVED_COUPLED_INVERSE_SMALLNESS_STABILITY_LIGHT_LOOPS_NONLINEAR_UV_REGGE_ORIGINAL_V_G_B_OR_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/cut.md",
            "notes/finite.md",
            "notes/measure.md",
            "notes/kernels.md",
            "notes/inverse.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "complete_two_mass_input": serialize(
            {
                "full_minimal_scalar_cut_and_fixed_finite_prescription": payload(
                    audit.geometry.data()
                ),
                "complete_two_threshold_reciprocal": payload(audit.spectral.data()),
            }
        ),
        "complete_reference_inverse": serialize(payload(audit.inverse.data())),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The actual conditional Proca-plus-heavy Gaussian fourth-order flat reference retains both full physical cuts and their independently fixed finite prescriptions. The minimal-scalar stress pair, full spin0/spin2 densities and general-dimensional fixed-physical leading/finite matching are derived without a conformal replacement or new counterterm. The combined factors have no first-sheet zero: the already specified heavy finite term changes the actual total threshold sign. Complete cut-only reciprocal measures retain both physical thresholds and the two-sided interior trace cusp. No old Proca pole is deleted from its frozen result or inserted into the inverse of the sum. Ordinary factor inverses belong to half-line L1 with finite unevaluated constants. Their absolutely convergent shifted primitives give the full flat scalar quotient a two-sided causal inverse on its complete initial-boundary graph, uniformly over all external momenta with no spatial derivative loss. The normalized operator bound375*T^4/[32*(log(n)+2)] is below125*T^4/4224. The physical force inverse restores64*pi^2*kappa, and the amplitude-to-metric Gram factor remains. This does not establish a controlled actual curved coupled inverse, interacting light/mixed loops, nonlinear same-state bounce, quantum gravitational limit, physical UV/Regge or original V/G/B/P8 closure.",
        "not_established": [
            "An actual curved tree/state/profile/matter normal form or compatible coupled inverse on the new model",
            "Small physical forcing errors, stability, no full-gravity growing roots or finite inhomogeneous nonlinear feedback",
            "An inverse of the gauge-singular full three-source block or literal homogeneous lapse/shift constraint solution",
            "Interacting light/mixed state, omitted-loop control or quantization commuting with the gravitational limit",
            "A physical cutoff, exact UV S matrix, finite-gravity Regge estimate or all-parent exclusion",
            "A changed finite prescription, suppressed heavy cut, omitted threshold cusp or original V/G/B/P8 closure",
        ],
        "verification_boundary": "Complete exact scalar stress, curvature, dimensional, spectral and ordered inverse identities support the written contour, endpoint and causal graph arguments. Independent full angular/boosted stress checks, literal curvature Hessians, complex radial quadrature, complete two-mass reciprocal reconstructions including the actual huge mass, threshold diagnostics and both matrix inverse products test the retained formulas. Primary scalar cut context is Martin-Verdaguer gr-qc/0001098v1 Eq5.8, with an explicit convention conversion; their subtraction and full gravity propagator prescriptions are not transferred. Native/direct/ordinary/CLI use original SymPy; only full regression uses the separately audited exact-GCD adapter. Frozen scientific bytes remain unchanged. The continuum proofs are not FORMALIZED.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete two-mass flat reference report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.247 complete two-mass flat reference replay passed; curved inverse and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
