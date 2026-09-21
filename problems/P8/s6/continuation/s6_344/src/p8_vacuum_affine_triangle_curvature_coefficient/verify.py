"""Read-only finite selected-triangle curvature matching in an explicit lift."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_curvature_contact_basis import verify as input0
from p8_vacuum_affine_heavy_parent_one_loop import verify as input3
from p8_vacuum_affine_heavy_scalar_four_point_loop import verify as input4
from p8_vacuum_affine_radiative_curvature_matching import verify as input2
from p8_vacuum_affine_triangle_box_radiation import verify as input1

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates/polynomial-vacuum-affine-triangle-curvature-coefficient.json"
)
PARENT_SHA = "dd49e3ecf06b03e905779e22625ad57147f3ab1d2b69f5d49896cea55338a89e"
RADIATION_SHA = "b73d47b8bb15941ebb82540312d1a681cea47295c0cef7f993aae7174aad1d04"
CURVATURE_SHA = "f8f9cd57c0240fe6008efd743d399c27dfc393a14d50501e8362289aae204db0"
SOURCE_SHA = "a09852126d79ca223af06555e10283328cb144a14945c127abee73b7ce8b9c05"
FLAT_SHA = "1277412d19d9b022eb81c3728251d1836d88e4b047c8ad80794107c2aff7c911"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_triangle_curvature_coefficient/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    for module, digest in (
        (input0, PARENT_SHA),
        (input1, RADIATION_SHA),
        (input2, CURVATURE_SHA),
        (input3, SOURCE_SHA),
        (input4, FLAT_SHA),
    ):
        if sha(module.REPORT) != digest:
            raise ValueError("A frozen finite-triangle matching input changed")
        module.validate_report(
            json.loads(module.REPORT.read_text()), module.build_report()
        )
    return {
        "S6_343_complete_first_curvature_basis_and_selected_poles_rebuilt": PARENT_SHA,
        "S6_342_complete_selected_triangle_and_box_radiation_rebuilt": RADIATION_SHA,
        "S6_336_literal_curvature_normalization_and_extra_parent_chi_rebuilt": CURVATURE_SHA,
        "S6_239_complete_original_limiting_source_inventory_rebuilt": SOURCE_SHA,
        "S6_235_original_flat_triangle_integral_and_mass_assignment_rebuilt": FLAT_SHA,
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
        raise ValueError("A finite triangle matching proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.344.FINITE_SELECTED_TRIANGLE_CURVATURE_COEFFICIENT_IN_EXPLICIT_LIFT",
        "date": "2026-09-18",
        "status": "SCOPED_FINITE_SELECTED_TRIANGLE_DEGREE6_CURVATURE_MATCHING; NOT_FULL_HARD_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/lift.md",
            "notes/triangle.md",
            "notes/moment.md",
            "notes/dressing.md",
            "notes/calibration.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_lift_and_triangle_coefficient": serialize(
            {name: payload(packets[name]) for name in names[:4]}
        ),
        "whole_complete_dressing_and_original_calibrations": serialize(
            {name: payload(packets[name]) for name in names[4:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "For the original selected matter triangle and the fully specified labeled scalar-box lift, the complete local radiative difference vanishes below six derivatives. Its first nonzero term is fixed generically over the entire Feynman simplex, including all three literal loop-line insertions. The outer-heavy-line product identity and all24 Bose labels give chi_triangle=g^2*(C+g^2/n)*c_triangle(n)/(16pi^2), with an exact convergent integral and rational-log formula for c_triangle. The original known coefficient is positive and its magnitude divided by the positive Born coefficient is below10^-206. This local derivative coefficient is already contained in the full S342 loop, is not a truncation bound at physical timelike kinematics, and does not assign or bound independent extra parent chi. Original V/G/B/P8 remains open.",
        "not_established": [
            "The independent added parent-theory curvature coefficient chi or full physical matching",
            "A complete aggregate curvature coefficient including boxes and every other matter class in this convention",
            "A physical above-threshold amplitude approximation or local Taylor truncation-error bound",
            "Internal-graviton loops, the full curved-background counterfunctional or finite-gravity quantum decoupling",
            "Full virtual or inclusive matching, exact LSZ, quantum unitarity, complex Regge, same-parent bounce, UV or original V/G/B/P8 closure",
        ],
        "verification_boundary": "Exact symbolic off-shell Gram, soft-dot and polarization algebra with all three cyclic line positions, coefficient-wise gamma integration and full simplex polynomial equality. Independently exact mass-integral primitive/endpoints, equal-mass limit, heavy-mass asymptotics and all24 Bose normalization; plus original massive four-vector calibrations and54 actual frozen S342 line kernels. Written uniform analytic-origin argument justifies the local expansion; no sampled fit substitutes for the whole polynomial. Not kernel-formalized, not the full physical nonlocal matching problem, and all frozen ancestors unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete finite-triangle curvature report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.344 finite selected triangle curvature coefficient replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
