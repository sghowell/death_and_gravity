"""Read-only full massive matter F1 and scoped gravity-vertex matching report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_gaussian_metric_pole_matching import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates/polynomial-vacuum-affine-massive-matter-graviton-vertex.json"
)
PARENT_SHA = "608a9394ac4d06f0f2e5c22a52e830f023364a906ea8e811bd2b1c7415e0c019"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_massive_matter_graviton_vertex/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError(
            "The frozen fixed H/Proca and unmatched light metric report changed"
        )
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    return {
        "S6_285_fixed_Gaussian_metric_matching_and_entire_ancestry_rebuilt": PARENT_SHA,
        "S279_remains_rejected_and_archived": True,
        "S275_S276_S277_and_scoped_P8a_unchanged": True,
        "all6_historical_physical_qualifications_unchanged": True,
        "all_original_primitive_and_matching_frontiers_retained": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A generated matter vertex or Ward proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.286.ORIGINAL_MASSIVE_MATTER_F1_POSITIVE_SPECTRAL_MEASURE_BOUNDED_T_CHANNEL_GRAVITY_PIECE_AND_CONDITIONAL_BACKGROUND_WARD_NORMALIZATION",
        "date": "2026-09-15",
        "status": "SCOPED_GENERATED_MATTER_F1_AND_CONDITIONAL_WARD_MATCHING; NOT_FULL_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/ward.md",
            "notes/vertex.md",
            "notes/spectral.md",
            "notes/matching.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_conditional_Ward": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_massive_matter_F1_and_t_channel_matching": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "Both original massive scalar stress triangles and cubic metric contacts have the complete mixed Ward normalization, including the inherited scalar OS counterterm. The entire generated matter F1 has a positive gapped spectral measure, F1prime(0)=Pi_second(mu)/6 and a controlled complex forward limit. Its finite t-channel spin2 coefficient is-Pi_second(mu)/(3kappa); its magnitude is below10^-1205 on abs(t)<=2 at the actual hierarchy. The general background identity cancels the entire regulated raw GR vertex/LSZ residue conditionally, not as a proof of the ordinary massless physical vertex. Explicit transverse Ricci-derivative freedom disproves a universal Ward-only slope determination. Full F2, H-metric mixing, curved local matching, the crossed finite amplitude and improved b20, massless IR/Regge and original-state/bounce obligations remain. Original V/G/B/P8 remain OPEN.",
        "not_established": [
            "Complete physical scalar-graviton F2, heavy-metric mixing and finite curvature-improvement/Ricci-derivative matching",
            "S285 light pure-curvature polynomial or full physical Newton pole sign; complete four-point local/tadpole matching",
            "Full finite crossed amplitude, improved b20, omitted loops or the finite Regge contour remainder",
            "Ordinary internal-gravity vertex, regulated-to-physical massless IR/detector/LSZ limit or exact interacting vacuum",
            "Physical EFT cutoff, original quantum state/domain/measure/bounce or V/G/B/P8 closure",
        ],
        "verification_boundary": "Whole original source, inherited scalar OS conditions, literal two massive stress triangles and both cubic metric contacts, all-component Ward identities, arbitrary symmetric constant-metric variation, positive spectral/complex-domain proofs and complete harmonic tensor contraction. The written analytic and covariance arguments are not FORMALIZED. Native/direct/ordinary/CLI retain original SymPy; the exact-GCD adapter is FULL-only. No frozen ancestor, raw report, finite counterterm or state is edited.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The generated massive matter vertex report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.286 generated massive matter F1 and conditional Ward replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
