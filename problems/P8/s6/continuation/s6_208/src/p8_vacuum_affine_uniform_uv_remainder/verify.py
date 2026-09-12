"""Read-only uniform actual UV-subtracted endpoint remainder and original tail."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_proca_gaussian import verify as gaussian
from p8_vacuum_affine_spatial_symbol import verify as parent

from . import audit, domain, estimates, limit, tail

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-uniform-uv-remainder.json"
PARENT_SHA = "d3eed45830d5b1af1a26bfe969bdff93a3aafefddebc61e7ce9c9d375cfaf9dc"
GAUSSIAN_SHA = "1175702d6980db090396863b807af0e767241daeb7d8dabf138ed11c67e1ae75"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_uniform_uv_remainder/*.py"))
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
        "S6_207_fully_rebuilt": PARENT_SHA,
        "S6_176_actual_CD_Gaussian_sector_transitively_rebuilt": GAUSSIAN_SHA,
        "fixed_stress_prescription_and_full_source_retained": True,
        "all_frozen_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A uniform actual UV remainder gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.208.UNIFORM_ACTUAL_CURVED_UV_SUBTRACTED_ENDPOINT_REMAINDER_AND_ORIGINAL_TWO_LEG_TAIL",
        "date": "2026-09-11",
        "status": "UNIFORM_ACTUAL_CURVED_UV_SUBTRACTED_ENDPOINT_REMAINDER_AND_ORIGINAL_TWO_LEG_TAIL; NOT_FULL_UV_CONTACT_REGULATOR_COVARIANT_MATCHING_RESPONSE_INVERSE_BACKGROUND_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/domain.md",
            "notes/readouts.md",
            "notes/remainder.md",
            "notes/tail.md",
            "notes/repartition.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "complete_uniform_domain_and_all_momentum_endpoint_bounds": serialize(
            {"domain": payload(domain.data()), "estimates": payload(estimates.data())}
        ),
        "original_regulator_tail_and_exact_actual_current_repartition": serialize(
            {"tail": payload(tail.data()), "limit": payload(limit.data())}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "A new explicit joint inverse-radius/time domain rho(P)=1/[100(m+|P|)] controls the complete actual unit-W8 normalized frequencies, genuine complex polarization frames and all ten constrained readouts. All four W8 coefficients and both Schwarz signs remain; normalized field norm20 and full pair norm1000 are proved only in that domain. Complete Cauchy source-jet estimates retain all fifteen UV coefficient slots and all35 time-jet entries. Far, complementary near and unexpanded low regions give a uniform actual UV-subtracted endpoint bound2e48||D||L2 X46[Gamma], with original two-created-mode regulator error4e53||D||L2 X46[Gamma]/K. Every near power and logarithmic term is retained. The removed union adds at most one momentum power, so the same six-spatial-derivative X46 norm suffices. An exact finite-regulator repartition rebases from S201's actual-state correction and finite time remainder, without double-counting earlier finite endpoint pieces or adding the unmatched local target. The new known actual piece has5e48 M[D]Y[Gamma] and tail5e53 M[D]Y[Gamma]/K, with canonical displays2e-751 and2e-746/K. The complete actual UV-symbol integral and full one-leg contact still require all finite/divergent coefficients, subleading sharp-band artifacts and original fixed covariant matching. No full matched response, reduced mixed inverse, finite-coupling background/stability, physical cutoff, remaining parent loops or finite-gravity IR/Regge result is established. Original V/G/B and P8 remain open.",
        "not_established": [
            "Complete quantum UV/contact finite and divergent coefficients, subleading sharp-band conversion and original fixed covariant matching",
            "Permission to replace the original masks, alter the full W8 reference or actual state, change finite counterterms or double-count earlier finite pieces",
            "A full matched response, reduced lapse/shift or derivative-compatible mixed inverse, finite-coupling background or stability",
            "Remaining parent loops/cutoff/scattering/IR/Regge or original V/G/B and P8 closure",
        ],
        "verification_boundary": "Exact scaled domain constants, full source-jet row bounds, all near power/log and far integrals, six-degree weight identities, original removed-union tail and exact finite-regulator repartition/canonical factors. Independent actual full W8 domains/readouts, all UV coefficients at separated external scales and two Cauchy resolutions, actual far and large-transfer near remainders, exact original far/near/low removed-region quadratures and exact-order/repartition controls supplement written continuous estimates. Native, direct science, ordinary and CLI use original SymPy; only full regression uses the audited exact GCD adapter. The Cauchy, global Fourier-norm and regulator-removal arguments are not FORMALIZED.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The uniform actual UV remainder report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.208.UNIFORM_ACTUAL_CURVED_UV_SUBTRACTED_ENDPOINT_REMAINDER_AND_ORIGINAL_TWO_LEG_TAIL replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
