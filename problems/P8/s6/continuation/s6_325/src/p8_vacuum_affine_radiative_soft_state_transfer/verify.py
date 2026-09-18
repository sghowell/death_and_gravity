"""Read-only uniform radiative finite-conversion and paired soft-state report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_radiative_angular_finite import verify as angular_input
from p8_vacuum_affine_single_residual_soft_dressing import verify as soft_input
from p8_vacuum_affine_three_real_probability_overlap import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates/polynomial-vacuum-affine-radiative-soft-state-transfer.json"
)
PARENT_SHA = "d1e3fdffd82b802be0fd0a8721f527f14697ea4a4543c300055c1f4ada321509"
ANGULAR_SHA = "5e18135d1a49c2846047f2a71c3de00f7c82b84a131d8a6f924890308f4ccd99"
SOFT_SHA = "22e89710b44e3df51ab35cf4ef541f70c42c68ee740a29f261f9f9d16694b59c"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_radiative_soft_state_transfer/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {k: v for k, v in data.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    for module, digest in (
        (previous, PARENT_SHA),
        (angular_input, ANGULAR_SHA),
        (soft_input, SOFT_SHA),
    ):
        if sha(module.REPORT) != digest:
            raise ValueError("A frozen radiative-state soft input changed")
        module.validate_report(
            json.loads(module.REPORT.read_text()), module.build_report()
        )
    return {
        "S6_324_probability_overlap_and_state_change_boundary_rebuilt": PARENT_SHA,
        "S6_301_full_dimensional_radiative_angular_factor_rebuilt": ANGULAR_SHA,
        "S6_309_same_state_leading_soft_pairing_convention_rebuilt": SOFT_SHA,
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
        raise ValueError("A uniform radiative-state soft proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.325.UNIFORM_RADIATIVE_FINITE_CONVERSION_AND_PAIRED_SOFT_STATE_TRANSFER",
        "date": "2026-09-17",
        "status": "SCOPED_UNIFORM_RADIATIVE_CONTINUITY_AND_SOFT_TRANSFER; NOT_FULL_LOOP_ALL_N_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/current.md",
            "notes/angular.md",
            "notes/regulator.md",
            "notes/matching.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_and_uniform_state_continuity": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_paired_soft_transfer": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The full S300/S301 radiative current has uniform finite-conversion continuity controlled by total radiated energy, independent of finite null multiplicity and collinear splitting. An energy-dependent regulator quotient bound supplies a genuine integrable dominator for the named Born-seeded same-state real/virtual soft difference. Its total variation is below200000*x*(1+ln(1/x))/kappa^2 after regulator removal. The physical marked state and full additional-soft dimensional projector/phase are distinguished. This is not a complete radiative hard loop, full NNLO inclusive rate, all-N hard matching, quantum state, Regge, bounce or original P8 closure.",
        "not_established": [
            "Complete radiative hard-loop reconstruction and whole NNLO real-virtual assembly",
            "Evanescent outer-state and finite hard matching",
            "All-N hard subtraction/summation or positive normalized detector measure",
            "Interacting quantum state, unitarity or absolute complex Regge",
            "Original common-parent bounce or V/G/B/P8 closure",
        ],
        "verification_boundary": "Exact projected-current derivative and minimum-majorant integral;18 rational original recoil/radial/angular calibrations with exact collinear splitting; full fixed-ball trace/radial and phase estimates including a regulator-uniform state difference; same-state soft probability pole with retained amplitude phase, endpoint integrability and calorimetric Born Gamma coefficient. The written proof establishes uniformity, not samples or integer gates alone. Only the named additional-soft insertion is evaluated, not all hard-loop graphs. Not kernel-formalized; original SymPy outside FULL.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The radiative soft-state report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.325 uniform radiative soft-state transfer replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
