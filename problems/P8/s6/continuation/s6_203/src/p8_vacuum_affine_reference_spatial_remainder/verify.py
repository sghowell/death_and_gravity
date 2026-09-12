"""Read-only all-momentum endpoint remainder and original regulator tail."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_proca_gaussian import verify as gaussian
from p8_vacuum_affine_reference_spatial_analyticity import verify as parent

from . import audit, limit, near, real, tail

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates" / "polynomial-vacuum-affine-reference-spatial-remainder.json"
)
PARENT_SHA = "2536af1c794aa5e5f23bf6c02aa7a7bcc5a51391431d4c02a7ee3e7214a1ca69"
GAUSSIAN_SHA = "1175702d6980db090396863b807af0e767241daeb7d8dabf138ed11c67e1ae75"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_reference_spatial_remainder/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    for previous, expected in ((parent, PARENT_SHA), (gaussian, GAUSSIAN_SHA)):
        if sha(previous.REPORT) != expected:
            raise ValueError(
                "The source-pinned actual affine Gaussian state chain changed"
            )
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_202_fully_rebuilt": PARENT_SHA,
        "S6_176_actual_CD_Gaussian_sector_transitively_rebuilt": GAUSSIAN_SHA,
        "fixed_stress_prescription_and_full_source_retained": True,
        "all_frozen_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("An all-momentum endpoint remainder gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.203.ACTUAL_CURVED_UNIT_W8_ALL_MOMENTUM_SPATIAL_TAYLOR_SUBTRACTED_ENDPOINT_REMAINDER_AND_ORIGINAL_REGULATOR_TAIL",
        "date": "2026-09-11",
        "status": "ACTUAL_CURVED_UNIT_W8_ALL_MOMENTUM_SPATIAL_TAYLOR_SUBTRACTED_ENDPOINT_REMAINDER_AND_ORIGINAL_REGULATOR_TAIL; NOT_FIXED_LOCAL_COVARIANT_MATCHING_FULL_INVERSE_BACKGROUND_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/real.md",
            "notes/near.md",
            "notes/tail.md",
            "notes/decomposition.md",
            "notes/regulator.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "all_real_endpoint_rows_and_complete_near_radial_bounds": serialize(
            {"real": payload(real.data()), "near": payload(near.data())}
        ),
        "all_momentum_regulator_limit_and_known_actual_finite_piece": serialize(
            {"tail": payload(tail.data()), "limit": payload(limit.data())}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The full first-five-endpoint unit-W8 remainder is controlled at all real internal and external momenta. The original |k|<m band is retained unexpanded; for |k|>=m only the degree-four spatial Taylor integrand is subtracted. The complete near region is bounded byC4|P|^4+C5|P|^5, including all five Taylor orders and the logarithmic radial term. Combining the low, near and previously controlled far pieces gives3e54||D||L2 X46[Gamma], with original both-created-mode regulator error5e60||D||L2 X46[Gamma]/K. The same six-spatial-derivative norm controls the tail. With the unchanged actual-to-unit correction and sixth-time-derivative bulk, the known actual finite piece is below4e54 M[D]Y[Gamma] and has tail6e60 M[D]Y[Gamma]/K. Both canonical factors give16e-746 and24e-740/K for that combined piece. The remaining Taylor-polynomial INTEGRAND sector need not be polynomial after integration over a moving sharp two-leg intersection; a constant-integrand geometric counterexample makes this explicit without claiming an actual current divergence coefficient. Its full finite/divergent content, regulator artifacts and fixed covariant matching with the distinct complete one-leg contact remain unresolved. No physical state, prescription or cutoff is changed, and no full matched response, inverse, interacting background, stability, remaining parent loops or finite-gravity IR/Regge is established. Original V/G/B and P8 remain open.",
        "not_established": [
            "Finite and divergent Taylor-integrand coefficients, regulator artifacts and fixed covariant matching with the complete distinct contact",
            "Polynomiality after a moving sharp band is integrated, or an actual current divergence coefficient from the geometric fixture alone",
            "A full matched response, derivative-compatible mixed inverse, interacting background or stability",
            "Full parent/cutoff/scattering/IR/Regge matching or original V/G/B and P8 closure",
        ],
        "verification_boundary": "Exact complete real row and radial constants, near power and logarithm bounds, all spatial weights, removed-union constants, sharp-intersection geometric counterexample and canonical factors. Independent actual all-five endpoint time jets, large external-momentum readouts, all five near radial integrals, original low and near removed unions and complete sphere intersections supplement written continuum estimates. Native, direct science, ordinary and CLI use original SymPy; only full regression uses the audited exact GCD adapter. The all-momentum estimates and weak-limit argument are not FORMALIZED.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The all-momentum endpoint remainder report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.203.ACTUAL_CURVED_UNIT_W8_ALL_MOMENTUM_SPATIAL_TAYLOR_SUBTRACTED_ENDPOINT_REMAINDER_AND_ORIGINAL_REGULATOR_TAIL replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
