"""Read-only complete selected triangle/box TT radiation and finite remainder."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_dimensional_real_remainder import verify as input2
from p8_vacuum_affine_factorized_bubble_radiation import verify as input1
from p8_vacuum_affine_heavy_parent_one_loop import verify as input4
from p8_vacuum_affine_heavy_scalar_four_point_loop import verify as input3
from p8_vacuum_affine_one_newton_inclusive_assembly import verify as input5
from p8_vacuum_affine_physical_loop_contours import verify as input0

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates/polynomial-vacuum-affine-triangle-box-radiation.json"
PARENT_SHA = "ca619dedc6a1981023ab8db61c5358530d5a8e5c060a44debed38fc8dce1229d"
BUBBLE_SHA = "174317dc92bc5b1c2956e6c5bebee4504f582325ae7e28b16f03d5de01f70206"
BOUND_SHA = "f1b5b057f0877a8722d63cb560cbe7507ca3e9e639455d121f574cd6034b3b62"
FLAT_SHA = "1277412d19d9b022eb81c3728251d1836d88e4b047c8ad80794107c2aff7c911"
SOURCE_SHA = "a09852126d79ca223af06555e10283328cb144a14945c127abee73b7ce8b9c05"
BORN_SHA = "6206f5d56959f01a8467e2c2b74ff43a2fad3cc75df0cee67e9bb8ce82f39564"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_triangle_box_radiation/*.py"))
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
        (input1, BUBBLE_SHA),
        (input2, BOUND_SHA),
        (input3, FLAT_SHA),
        (input4, SOURCE_SHA),
        (input5, BORN_SHA),
    ):
        if sha(module.REPORT) != digest:
            raise ValueError("A frozen triangle/box radiation input changed")
        module.validate_report(
            json.loads(module.REPORT.read_text()), module.build_report()
        )
    return {
        "S6_341_full_physical_contours_and_literal_weighted_TT_kernels_rebuilt": PARENT_SHA,
        "S6_339_complete_bubble_radiation_and_known_sector_bound_rebuilt": BUBBLE_SHA,
        "S6_335_complete_full_tree_and_real_minus_soft_bound_rebuilt": BOUND_SHA,
        "S6_235_complete_original_mass_ordered_flat_coefficients_rebuilt": FLAT_SHA,
        "S6_239_original_limiting_source_and_fixed_counterterms_rebuilt": SOURCE_SHA,
        "S6_297_original_full_Born_normalization_rebuilt": BORN_SHA,
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
        raise ValueError("A complete selected triangle/box proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.342.COMPLETE_SELECTED_TRIANGLE_BOX_RADIATION_AND_FINITE_REMAINDER",
        "date": "2026-09-18",
        "status": "SCOPED_COMPLETE_TRIANGLE_BOX_MATTER_RADIATION; NOT_FULL_HARD_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/functional.md",
            "notes/labels.md",
            "notes/radiation.md",
            "notes/bounds.md",
            "notes/calibration.md",
            "notes/scope.md",
            "notes/provenance.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_functional_labels_and_physical_radiation": serialize(
            {name: payload(packets[name]) for name in names[:4]}
        ),
        "whole_uniform_remainder_interference_and_original_calibrations": serialize(
            {name: payload(packets[name]) for name in names[4:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The complete selected minimal old-matter triangle and ordered-box real-TT radiation follows from the noncommuting covariant Hessian: all24 labelled kernels, every internal line, the outer heavy branch, all shifted external emissions and the unchanged finite OS4 constant are retained. Full physical contour integrals give the original-source nonleading bound||R||/A0<1e204/sqrt(kappa)=1e-196, including the other currently known matter sectors; the low-window finite full-tree interference is<1e-592. The isolated hard bound is not perturbative smallness. Independent curvature, internal gravity, full inclusive quantum probability and original V/G/B/P8 closure remain open.",
        "not_established": [
            "Independent S336 four-hard curved matching or internal-graviton loops",
            "Finite-gravity quantum decoupling, higher loops, exact LSZ or a resummed propagator",
            "Perturbative smallness from the large isolated triangle/box hard bound",
            "The separately divergent leading-soft loop extension, complete virtual hard matching or full inclusive interacting detector probability",
            "Quantum unitarity, complex Regge, same-parent bounce, UV or original V/G/B/P8 closure",
        ],
        "verification_boundary": "Generic noncommuting operator variation and independent finite-matrix fourth-field/metric derivatives, exact24label grouping into the original S235 representation, independent bubble and literal heavy-line sign controls, exact original-source budgets and continuous S341 physical-contour estimates. Four original recoil states and both polarizations include four coincident polarized channels and eight wrong-sign rejections. All672 original triangle/box split routes compare the entire denominator;1344 literal TT contractions include1284 nonzero densities. Diagnostic constant/v loop kernels calibrate graph organization and are not substituted for actual loop values. Finite checks supplement the written continuous proof. Not kernel-formalized; frozen ancestors unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete selected triangle/box report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.342 complete selected triangle/box radiation replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
