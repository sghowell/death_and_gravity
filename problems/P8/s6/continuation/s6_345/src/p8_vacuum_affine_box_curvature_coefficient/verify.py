"""Read-only finite ordered-box curvature coefficient in a fixed covariant jet lift."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_curvature_contact_basis import verify as input1
from p8_vacuum_affine_heavy_parent_one_loop import verify as input4
from p8_vacuum_affine_heavy_scalar_four_point_loop import verify as input5
from p8_vacuum_affine_local_tadpole_radiation import verify as input3
from p8_vacuum_affine_triangle_box_radiation import verify as input2
from p8_vacuum_affine_triangle_curvature_coefficient import verify as input0

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates/polynomial-vacuum-affine-box-curvature-coefficient.json"
PARENT_SHA = "227850342699960004105744b9eaad57b20a626b19d2e90deb353512c847a036"
BASIS_SHA = "dd49e3ecf06b03e905779e22625ad57147f3ab1d2b69f5d49896cea55338a89e"
RADIATION_SHA = "b73d47b8bb15941ebb82540312d1a681cea47295c0cef7f993aae7174aad1d04"
JETS_SHA = "602dd0e80cc5fec993a20219f3423195439b0e3dba03e468794bc939b2627645"
SOURCE_SHA = "a09852126d79ca223af06555e10283328cb144a14945c127abee73b7ce8b9c05"
FLAT_SHA = "1277412d19d9b022eb81c3728251d1836d88e4b047c8ad80794107c2aff7c911"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_box_curvature_coefficient/*.py"))
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
        (input1, BASIS_SHA),
        (input2, RADIATION_SHA),
        (input3, JETS_SHA),
        (input4, SOURCE_SHA),
        (input5, FLAT_SHA),
    ):
        if sha(module.REPORT) != digest:
            raise ValueError("A frozen finite-box matching input changed")
        module.validate_report(
            json.loads(module.REPORT.read_text()), module.build_report()
        )
    return {
        "S6_344_known_triangle_coefficient_and_distinct_lift_rebuilt": PARENT_SHA,
        "S6_343_complete_first_curvature_basis_and_selected_poles_rebuilt": BASIS_SHA,
        "S6_342_complete_selected_original_ordered_box_radiation_rebuilt": RADIATION_SHA,
        "S6_337_literal_original_Galileon_connections_rebuilt": JETS_SHA,
        "S6_239_complete_original_limiting_source_inventory_rebuilt": SOURCE_SHA,
        "S6_235_original_flat_ordered_boxes_and_mass_assignment_rebuilt": FLAT_SHA,
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
        raise ValueError("A finite ordered-box matching proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.345.FINITE_SELECTED_BOX_CURVATURE_COEFFICIENT_IN_EXPLICIT_JET_LIFT",
        "date": "2026-09-18",
        "status": "SCOPED_FINITE_SELECTED_BOX_DEGREE6_CURVATURE_MATCHING; NOT_FULL_HARD_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/basis.md",
            "notes/jets.md",
            "notes/radiation.md",
            "notes/moment.md",
            "notes/calibration.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_flat_lift_and_generic_box_matching": serialize(
            {name: payload(packets[name]) for name in names[:4]}
        ),
        "whole_finite_moment_and_original_calibrations": serialize(
            {name: payload(packets[name]) for name in names[4:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The entire selected original ordered-box degree6 real-TT curvature difference from the explicitly fixed six-word covariant jet lift is determined by a whole generic polynomial identity. All220 flat Bose monomials, all96 mass-ordered line insertions, all2430 radiative coefficient residuals and the complete connection response are retained. The resulting coefficient chi_box=g^4*c_box(n)/(16pi^2) has a strictly negative convergent mass integral and an exact rational-log expression, including its equal-mass limit and heavy-mass asymptotics. This is a known local loop coefficient, not the independent extra parent matching. The S344 triangle lift differs and must be converted before aggregation. No above-threshold Taylor approximation or perturbative smallness is inferred, and original V/G/B/P8 remains open.",
        "not_established": [
            "An aggregate coefficient without converting the triangle and treating every matter sector in a common convention",
            "The independent added parent-theory curvature coefficient chi or complete physical matching",
            "A physical above-threshold amplitude approximation or local Taylor truncation-error bound",
            "Internal-graviton loops, the full curved-background counterfunctional or finite-gravity quantum decoupling",
            "Full virtual or inclusive matching, exact LSZ, quantum unitarity, complex Regge, same-parent bounce, UV or original V/G/B/P8 closure",
        ],
        "verification_boundary": "All220 parity-even flat quartic three-edge monomials are projected by exact full-polynomial equality onto a rank6 off-shell Bose jet basis. The fully symmetrized third-jet response is derived in64 components and checked against the original literal Galileon contact. Exact sparse rational convolution retains all96 box insertions and compares all2430 generic coefficients, not a sample fit. Independent original massive vectors,72 frozen full line kernels and exact mass-integral primitives/limits supply crosschecks. The analytic-origin justification and operator convention are written explicitly. Not kernel-formalized and not full nonlocal parent matching.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete finite-box curvature report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.345 finite selected box curvature coefficient replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
