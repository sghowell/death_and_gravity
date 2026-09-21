"""Read-only complete selected quadratic-radiation-cancellation radiation report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_factorized_bubble_radiation import verify as input0
from p8_vacuum_affine_heavy_parent_one_loop import verify as input3
from p8_vacuum_affine_heavy_scalar_one_loop import verify as input2
from p8_vacuum_affine_massive_scalar_graviton_bremsstrahlung import verify as input4
from p8_vacuum_affine_mixed_source_radiation import verify as input1

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT
    / "certificates/polynomial-vacuum-affine-quadratic-radiation-cancellation-radiation.json"
)
PARENT_SHA = "174317dc92bc5b1c2956e6c5bebee4504f582325ae7e28b16f03d5de01f70206"
KERNEL_SHA = "faf709f6d72af441b9ff14d6ceb57fe573cf22d6b37a45d80ccf58205d5c9bc6"
LOOP_SHA = "d6646d25d158bcfeb32bf07895e8fcc2ca08e05988919f803ce50aa5e45f58a5"
SOURCE_SHA = "a09852126d79ca223af06555e10283328cb144a14945c127abee73b7ce8b9c05"
VERTEX_SHA = "c5aeb643e173157f280201af2276924d5a211a5f1f97a271a0ed8dfc8c6d4d5f"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_quadratic_radiation_cancellation/*.py"))
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
        (input1, KERNEL_SHA),
        (input2, LOOP_SHA),
        (input3, SOURCE_SHA),
        (input4, VERTEX_SHA),
    ):
        if sha(module.REPORT) != digest:
            raise ValueError(
                "A frozen quadratic-radiation-cancellation radiation input changed"
            )
        module.validate_report(
            json.loads(module.REPORT.read_text()), module.build_report()
        )
    return {
        "S6_339_full_factorized_bubbles_and_unassigned_curvature_frontier_rebuilt": PARENT_SHA,
        "S6_338_both_literal_mixed_massive_triangles_rebuilt": KERNEL_SHA,
        "S6_234_complete_two_point_inverse_and_fixed_OS_conditions_rebuilt": LOOP_SHA,
        "S6_239_original_source_graph_count_and_counterterm_extension_rebuilt": SOURCE_SHA,
        "S6_295_literal_minimal_scalar_vertex_normalization_rebuilt": VERTEX_SHA,
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
        raise ValueError(
            "A selected quadratic-radiation-cancellation radiation proof gate failed"
        )
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.340.COMPLETE_SELECTED_QUADRATIC_RADIATION_CANCELLATION",
        "date": "2026-09-18",
        "status": "SCOPED_COMPLETE_QUADRATIC_RADIATION_CANCELLATION; NOT_FULL_HARD_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/ward.md",
            "notes/loops.md",
            "notes/cancellation.md",
            "notes/calibration.md",
            "notes/scope.md",
            "notes/provenance.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_inverse_and_actual_mixed_loops": serialize(
            {name: payload(packets[name]) for name in names[:3]}
        ),
        "whole_exact_quadratic_cancellation_and_original_calibrations": serialize(
            {name: payload(packets[name]) for name in names[3:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The complete selected old massive-matter first-loop external quadratic radiation is exactly zero in physical real TT projection. Both literal mixed loop insertions and all fixed counterterms give the inverse divided difference; its vertex contribution cancels the single self-energy propagator graph with the correct Dyson sign, legwise for arbitrary shifted hard subamplitudes. The quartic tadpole and heavy-onepoint source pair are retained and canceled under the existing conditions. Positive massive parameter gaps and a double on-shell zero justify the continuous soft limit. No exact LSZ, resummed propagator, remaining hard triangle/box or internal-graviton loop, independent curvature matching, inclusive probability or original P8 closure follows.",
        "not_established": [
            "Remaining hard triangle, ordered-box or internal-graviton radiation",
            "Exact LSZ, resummed propagator, unprojected/off-shell metric response or independent curved matching",
            "Full inclusive interacting detector probability or loop rate",
            "Quantum unitarity, complex Regge, same-parent bounce, UV or original V/G/B/P8 closure",
        ],
        "verification_boundary": "Full generic inverse Ward and literal free scalar tensor identities, explicit both-mass exact-D triangle insertions, all fixed mass/residue/source terms, complete noncommuting variations of eight quadratic powers, exact full logarithmic double zero and original positive parameter gaps. Twelve original recoil external legs, two polarizations and three rational Feynman parameters independently compare both finite triangle endpoints with the propagator insertion. Finite controls do not replace the written uniform analytic proof. Not kernel-formalized; all frozen ancestor evidence unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The quadratic-radiation-cancellation radiation report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.340 complete selected quadratic-radiation-cancellation radiation replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
