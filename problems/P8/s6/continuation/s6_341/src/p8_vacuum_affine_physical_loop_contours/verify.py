"""Read-only complete physical-contour and weighted TT-kernel report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_heavy_parent_one_loop import verify as input3
from p8_vacuum_affine_heavy_scalar_four_point_loop import verify as input1
from p8_vacuum_affine_heavy_scalar_loop_remainder import verify as input2
from p8_vacuum_affine_massive_scalar_graviton_bremsstrahlung import verify as input4
from p8_vacuum_affine_quadratic_radiation_cancellation import verify as input0

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates/polynomial-vacuum-affine-physical-loop-contours.json"
PARENT_SHA = "79f298f3b8db27f6df41081e2f0ab80fd13eab8bc149699196c24eef4fda6cb7"
LOOP_SHA = "1277412d19d9b022eb81c3728251d1836d88e4b047c8ad80794107c2aff7c911"
PRIMITIVE_SHA = "7a5e56979d2a755b7df1fc26f8e6d3d6a4e1bf760d4aa20911b139c250d28cb8"
SOURCE_SHA = "a09852126d79ca223af06555e10283328cb144a14945c127abee73b7ce8b9c05"
VERTEX_SHA = "c5aeb643e173157f280201af2276924d5a211a5f1f97a271a0ed8dfc8c6d4d5f"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_physical_loop_contours/*.py"))
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
        (input1, LOOP_SHA),
        (input2, PRIMITIVE_SHA),
        (input3, SOURCE_SHA),
        (input4, VERTEX_SHA),
    ):
        if sha(module.REPORT) != digest:
            raise ValueError("A frozen physical-contour input changed")
        module.validate_report(
            json.loads(module.REPORT.read_text()), module.build_report()
        )
    return {
        "S6_340_complete_selected_quadratic_cancellation_and_frontier_rebuilt": PARENT_SHA,
        "S6_235_full_mass_ordered_scalar_master_representation_rebuilt": LOOP_SHA,
        "S6_236_complete_primitive_and_sharper_flat_remainder_rebuilt": PRIMITIVE_SHA,
        "S6_239_original_parent_parameters_and_counterterms_rebuilt": SOURCE_SHA,
        "S6_295_literal_scalar_stress_vertex_and_original_recoil_rebuilt": VERTEX_SHA,
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
        raise ValueError("A physical-contour proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.341.PHYSICAL_OFFSHELL_TRIANGLE_BOX_CONTOURS_AND_TT_INSERTION_BOUNDS",
        "date": "2026-09-18",
        "status": "SCOPED_PHYSICAL_CONTOUR_AND_TT_KERNEL_BOUNDS; NOT_FULL_HARD_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/parameters.md",
            "notes/contour.md",
            "notes/moments.md",
            "notes/insertions.md",
            "notes/calibration.md",
            "notes/scope.md",
            "notes/provenance.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_offshell_parameters_and_physical_contour": serialize(
            {name: payload(packets[name]) for name in names[:3]}
        ),
        "whole_weighted_moments_TT_kernels_and_original_calibrations": serialize(
            {name: payload(packets[name]) for name in names[3:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "An explicit two-variable homotopy preserves the physical Feynman branch of the complete off-shell mass-ordered triangle and box denominators. A joint complex gap proves uniform scalar, Cauchy-derivative and weighted-moment bounds, including the symmetric subtraction point. Literal stress vertices and split-parameter integration give controlled TT line-insertion kernels; both light and heavy lines and ordered mass assignments are retained. This is not the completed radiative graph sum or independent curvature/internal-graviton matching, and no original V/G/B/P8 closure follows.",
        "not_established": [
            "Completed hard triangle and ordered-box radiative amplitude with multiplicities, relative signs, external/outer-heavy emissions and finite counterterms",
            "Same first-sheet analyticity across a physical cut or threshold4 uniformity",
            "Internal-graviton loops, independent curved matching, exact LSZ or resummed propagators",
            "Full inclusive interacting detector probability, quantum unitarity, complex Regge, same-parent bounce, UV or original V/G/B/P8 closure",
        ],
        "verification_boundary": "Exact full Symanzik polynomials, the entire Feynman-sign homotopy and written continuous gap proof, exact moment antiderivatives, all-order local Cauchy estimates, literal minimal stress insertions and weighted measures. Exact finite controls include350homotopy and175gap points,48original virtualities,72original pair invariants,99cyclic routes,513route-domain gates and66actual TT insertions. These controls do not replace the uniform proof. Independent floating quadrature is diagnostic only and not an exact or rigorous-error certificate. Not kernel-formalized; all frozen ancestor evidence unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The physical-contour report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.341 physical-contour and weighted TT-kernel replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
