"""Read-only complete first-loop low-energy coefficient matching."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_heavy_scalar_loop_remainder import verify as parent
from p8_vacuum_affine_heavy_scalar_tree_matching import verify as model_input

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT
    / "certificates"
    / "polynomial-vacuum-affine-heavy-scalar-loop-coefficients.json"
)
PARENT_SHA = "7a5e56979d2a755b7df1fc26f8e6d3d6a4e1bf760d4aa20911b139c250d28cb8"
MODEL_SHA = "31716f2b701c0120c328dd8a2506f677eb30edeca956d531c57975e38d8d25cb"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_heavy_scalar_loop_coefficients/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    for prior, expected in ((parent, PARENT_SHA), (model_input, MODEL_SHA)):
        if sha(prior.REPORT) != expected:
            raise ValueError(
                "A frozen loop, finite prescription or original tree target changed"
            )
        prior.validate_report(
            json.loads(prior.REPORT.read_text()), prior.build_report()
        )
    return {
        "S6_236_complete_primitive_and_physical_first_loop_remainder_fully_rebuilt": PARENT_SHA,
        "S6_233_exact_separate_tree_and_original_low_coefficient_targets_fully_rebuilt": MODEL_SHA,
        "S6_234_S6_235_full_loop_normalization_counterterms_and_OS4_condition_transitively_rebuilt": True,
        "all_frozen_original_and_separate_model_scientific_bytes_unchanged": True,
        "no_finite_derivative_retuning_or_unknown_heavy_weight_subtraction": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError(
            "A complete logarithmic moment, complex-bidisk or coefficient-error gate failed"
        )
    return {
        "schema": 1,
        "claim": "P8-S6.237.SEPARATE_V2S_T1_OS4_COMPLETE_FIRST_LOOP_LOW_ENERGY_B20_B21_B40_MATCHING_WITH_JOINT_COMPLEX_BIDISK_REMAINDERS_AND_UNCHANGED_FINITE_PRESCRIPTION",
        "date": "2026-09-12",
        "status": "SEPARATE_COMPLETE_FIRST_LOOP_LOW_COEFFICIENT_MATCHING; NOT_OMITTED_LOOPS_EXACT_PHYSICAL_COEFFICIENTS_FULL_UV_ORIGINAL_AFFINE_PARENT_BOUNCE_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/jets.md",
            "notes/bidisk.md",
            "notes/products.md",
            "notes/coefficients.md",
            "notes/scope.md",
            "notes/literature.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "complete_four_order_loop_and_logarithmic_moment_reduction": serialize(
            payload(audit.jets.data())
        ),
        "full_joint_bidisk_and_actual_coefficient_matching_bounds": serialize(
            {
                "bidisk": payload(audit.bidisk.data()),
                "coefficients": payload(audit.coefficients.data()),
            }
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The complete mass-ordered first-loop amplitude of the unchanged V2S-T1-OS4 prescription is expanded through n^-5 with every light logarithmic moment reduced by an exact full integration-by-parts identity. The complete orders are crossing-symmetric polynomials of degree at most0,0,2,3. A separately proved joint complex bidisk |v|,|t|<=1/4 gives a full remainder below10^12 g^4/(16pi^2 n^6), including every higher inverse-mass term and all vertex/triangle/box cross products. Cauchy bounds with the correct coefficient factors16,64,256 then control the complete first-loop corrections:0<delta_b20/(4lambda)<10^-203,0<delta_b21/(-3gamma)<10^-203 and |delta_b40|/(gamma^2/lambda)<10^-192. The first b20 shift is positive and b21 negative; the sign of the first b40 shift is not determined. The explicitly tree-plus-first-loop higher coefficient remains positive and near its separate tree value. The existing constant OS4 contact has zero derivative in these observables, so no new finite condition or heavy-weight subtraction is used. These are complete first-loop matching statements, not omitted-loop bounds, exact physical LSZ or dispersion coefficients, a controlled resonance, full UV/Regge theory or original common-parent bounce and V/G/B/P8 closure.",
        "not_established": [
            "An omitted-loop remainder or exact all-order physical coefficient matching",
            "A sign for the first b40 correction, exact Gram saturation, exact quantum S matrix, controlled resonance or full high-energy dispersion arc",
            "A full quantum potential or all-counterterm physical-vacuum theorem",
            "The original affine/DHOST/Proca common-parent bounce/state matching, finite-gravity Regge control or nonlinear stability",
            "A physical cutoff, all-parent exclusion or original V/G/B/P8 closure",
        ],
        "verification_boundary": "Exact complete coefficient jets, full logarithmic moment identities, ordered diagram cancellations and rational margins accompany written joint-holomorphy and complete Cauchy-tail proofs. These arguments are not FORMALIZED. Independent quadratures and high-precision diagnostics do not replace the continuum estimates. The real-window result is not used as an analytic derivative bound. Native/direct/ordinary/CLI retain original SymPy; only full regression uses the audited exact-GCD adapter.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The full first-loop coefficient report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.237 full first-loop coefficient matching replay passed; omitted loops and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
