"""Read-only full unequal-mass Regge-error certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_complete_matter_graviton_endpoint import verify as metric_input
from p8_vacuum_affine_one_newton_inclusive_assembly import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates/polynomial-vacuum-affine-unequal-mass-regge-residue.json"
PARENT_SHA = "6206f5d56959f01a8467e2c2b74ff43a2fad3cc75df0cee67e9bb8ce82f39564"
METRIC_SHA = "191329017c78af9fcb017a3bf77df87bebb2fb47c231948f37dc1e743200c9b7"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_unequal_mass_regge_residue/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {k: v for k, v in data.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen inclusive-assembly parent report changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    if sha(metric_input.REPORT) != METRIC_SHA:
        raise ValueError("The frozen complete matter metric endpoint changed")
    metric_input.validate_report(
        json.loads(metric_input.REPORT.read_text()), metric_input.build_report()
    )
    return {
        "S6_297_inclusive_assembly_entire_ancestry_rebuilt": PARENT_SHA,
        "S6_290_complete_matter_metric_endpoint_rebuilt": METRIC_SHA,
        "S279_remains_rejected_and_archived": True,
        "S275_S276_S277_and_scoped_P8a_unchanged": True,
        "all6_historical_physical_qualifications_unchanged": True,
        "all_original_primitive_and_matching_frontiers_retained": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A complete unequal-mass Regge proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.298.UNEQUAL_MASS_SPIN2_SPECTRAL_MOMENT_AND_CONDITIONAL_REGGE_ERROR",
        "date": "2026-09-15",
        "status": "SCOPED_COMPLETE_UNEQUAL_MASS_MOMENT_AND_CONDITIONAL_REGGE_ERROR_MAP; NOT_FULL_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/waves.md",
            "notes/moments.md",
            "notes/regge.md",
            "notes/matching.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_unequal_mass_spin2": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_positive_spectral_moment_and_conditional_Regge_budget": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "Both original unequal-mass cuts reproduce the complete generated matter spin2 slope. The formal HH cut contributes a quarter within2*10^-98, whereas the low light window4<T<=16 captures less than10^-195 of that slope. Coupled-pole unitarity at spin2 reproduces both endpoint factors. A positive Legendre-Q representation yields a whole-function order error and explicit finite-transfer contour budget retaining the unknown arc, residue-ratio and trajectory inputs. These inputs and the high-energy contour are not matched. Original V/G/B/P8 remain OPEN.",
        "not_established": [
            "Actual Regge trajectory or species residue-ratio matching",
            "A numerical transfer arc norm or absence of other Regge singularities",
            "High-energy complex-energy remainder and analytic trajectory terms",
            "Full physical massless IR observable and omitted-loop bounds",
            "Original state/domain/measure/bounce or V/G/B/P8 closure",
        ],
        "verification_boundary": "Exact tree projections, source-preserving massive graph moments, positive integral inequalities and written finite-contour arguments prove the stated scoped identities and conditional bound. Independent quadrature checks the full functions but is not the proof of a uniform bound. Integer-spin unitarity does not establish the assumed complex-J continuation or a UV completion. Native/direct/ordinary/CLI retain original SymPy; the exact-GCD adapter remains FULL-only. All frozen ancestors, historical qualifications and original closure statuses remain unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete unequal-mass Regge report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.298 full unequal-mass Regge error-map replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
