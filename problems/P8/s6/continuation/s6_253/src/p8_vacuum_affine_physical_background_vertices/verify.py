"""Read-only entire homogeneous physical-background Gaussian vertex report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_gaussian_measure_response import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates" / "polynomial-vacuum-affine-physical-background-vertices.json"
)
PARENT_SHA = "504045fc426a73d2de1d22d5532cdf0238b6659cbbe7c24e162abf8ba2a5e3da"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_physical_background_vertices/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen full reference Gaussian functional changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    return {
        "S6_252_whole_Gaussian_functional_and_entire_source_ancestry_fully_rebuilt": PARENT_SHA,
        "all_frozen_original_and_separate_model_scientific_bytes_unchanged": True,
        "same_full_product_preparations_QG2_fixed_profiles_and_M1_mean": True,
        "no_interacting_state_or_QG3_retuning_substituted": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A complete physical-background vertex or scope gate failed")
    packets = audit.packets()
    return {
        "schema": 1,
        "claim": "P8-S6.253.COMPLETE_CURRENT_HOMOGENEOUS_PHYSICAL_GAUSSIAN_BACKGROUND_VERTICES_WITH_JOINT_SCALAR_LONGITUDINAL_PROCA_CONSTRAINT_MEASURE_ALL_EIGHT_REFERENCE_MODES_AND_HELD_VECTOR_EMBEDDING_CONTACT",
        "date": "2026-09-13",
        "status": "EXACT_COMPLETE_HOMOGENEOUS_GAUGE_FIXED_GAUSSIAN_VERTICES_JOINT_MEASURE_AND_FIXED_PRODUCT_RESPONSE; NOT_FULL_COVARIANT_VERTEX_MEASURE_CURVED_SUBTRACTION_PHYSICAL_MEAN_LOOP_BOUND_NONLINEAR_UV_REGGE_ORIGINAL_V_G_B_OR_P8",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/adm.md",
            "notes/coupled.md",
            "notes/physical.md",
            "notes/ctp.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "entire_current_parent_germs_coupled_reduction_and_joint_measure": serialize(
            {
                name: payload(packets[name])
                for name in (
                    "entire_parent_source_germs",
                    "all_clock_lapse_jets_and_reference_coefficients",
                    "whole_aligned_three_mode_constraint_reduction",
                    "whole_joint_physical_Gaussian_measure",
                )
            }
        ),
        "whole_physical_vertices_embedding_contacts_and_fixed_product_response": serialize(
            {
                name: payload(packets[name])
                for name in (
                    "whole_held_vector_source_and_embedding_contact",
                    "whole_homogeneous_physical_background_vertices",
                    "all_remaining_modes_and_existing_state_transports",
                    "whole_fixed_product_Gaussian_response",
                )
            }
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The complete current parent licenses all required clock lapse jets while retaining both fixed profiles and all three vacuum constants. Its homogeneous aligned scalar and longitudinal-Proca Gaussian is reduced jointly, including the whole source square, nonzero lapse/vector auxiliary coupling, lower constraint bracket and physical-density measure. Every6x6 physical lapse, scale and matter first/second vertex is specified by its exact whole-matrix chain rule and full source-pinned substitution; all off-reference Ward contacts, momentum/scale chains and time boundaries remain. The full held-W0 square gives the nonzero second embedding contact and its normal-vector light-loop vertex. Both TT, both transverse Proca and H modes retain their entire physical vertices and unchanged state maps, completing the eight-mode finite Gaussian response. This does not construct arbitrary inhomogeneous covariant vertices, nonlinear gauge/connection/field-map measure, selected Ward/variational curved subtraction, renormalized physical means, a loop-size bound, a compatible nonlinear bounce, quantum UV/Regge or original V/G/B/P8 closure.",
        "not_established": [
            "The full arbitrary inhomogeneous covariant physical background Hessian and nonlinear fluctuation-map contacts",
            "A nonlinear BRST, ghost, affine connection or disformal quantum measure",
            "A selected covariant Ward/variational curved counterfunctional, physical metric/vector/light mean or quantitative omitted-loop norm",
            "A uniform perturbed principal-symbol neighborhood, interacting state, convergent loop expansion or compatible nonlinear bounce",
            "Quantum gravitational decoupling, physical UV scattering, finite-gravity infrared/Regge control or original V/G/B/P8 closure",
            "Formalization of the written gauge, functional or continuum arguments",
        ],
        "verification_boundary": "Exact complete symbolic identities and independently re-entered ADM, physical-path, Fourier, Dirac, Fock and evolution diagnostics support the explicitly scoped Gaussian result. A same-reference Wick functional is not a physical covariant stress prescription. The finite interface does not certify an unrestricted background or a loop norm. Native/direct/ordinary/CLI use original SymPy; only full regression uses the audited exact-GCD adapter. No frozen scientific file is modified.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete homogeneous physical-background report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.253 full homogeneous physical-background vertex replay passed; physical curved subtraction and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
