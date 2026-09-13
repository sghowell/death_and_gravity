"""Read-only full current nonscalar conditional inverse report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_two_mass_coupled_inverse import verify as parent
from p8_vacuum_affine_two_mass_reference import verify as reference_input

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates" / "polynomial-vacuum-affine-two-mass-nonscalar-inverse.json"
)
PARENT_SHA = "f34f33bc1568492f2b4ce375d57f9a416ec5e1560435012245cae4ca80d7273c"
REFERENCE_SHA = "dce8b1593e15778f9b5ecb2ea224e0cb9b834d2e23a6d1c6fad126c358a8afd0"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_two_mass_nonscalar_inverse/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    for prior, expected in ((parent, PARENT_SHA), (reference_input, REFERENCE_SHA)):
        if sha(prior.REPORT) != expected:
            raise ValueError(
                "A frozen actual conditional or complete two-mass reference input changed"
            )
        prior.validate_report(
            json.loads(prior.REPORT.read_text()), prior.build_report()
        )
    return {
        "S6_248_actual_current_scalar_clock_M1_inverse_and_full_response_ancestry_fully_rebuilt": PARENT_SHA,
        "S6_247_complete_same_prescription_two_mass_reference_fully_rebuilt": REFERENCE_SHA,
        "all_frozen_original_and_separate_model_scientific_bytes_unchanged": True,
        "old_QG1_tensor_inverse_or_scalar_force_domain_not_transferred": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("An actual nonscalar inverse gate failed")
    packets = audit.packets()
    return {
        "schema": 1,
        "claim": "P8-S6.249.ACTUAL_CURRENT_TWO_MASS_CONDITIONAL_TENSOR_AND_TRANSVERSE_VECTOR_CAUSAL_INVERSE_WITH_FULL_NONSCALAR_CURVATURE_AND_COMPLETE_WARD_COMPATIBLE_SHIFT_FORCE_GRAPH",
        "date": "2026-09-13",
        "status": "ACTUAL_CONDITIONAL_NONSCALAR_TWO_MASS_PREPARED_FOURIER_BALL_INVERSE; NOT_ARBITRARY_SHIFT_FORCE_UNRESTRICTED_GRAPH_SMALLNESS_STABILITY_INTERACTING_LOOPS_NONLINEAR_UV_REGGE_ORIGINAL_V_G_B_OR_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/sectors.md",
            "notes/local.md",
            "notes/ward.md",
            "notes/remainder.md",
            "notes/inverse.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_nonscalar_assembly": serialize(
            {
                name: payload(data)
                for name, data in packets.items()
                if name != "complete_nonscalar_inverse_and_constraint_recovery"
            }
        ),
        "actual_normal_form_and_inverse": serialize(
            payload(packets["complete_nonscalar_inverse_and_constraint_recovery"])
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The same actual conditional QG2-H8A420 reference admits two-sided tensor and transverse-vector causal inverses on each stated smooth prepared external Fourier ball. The complete O(2) representation, including reflection, and coordinate-free orthogonal projectors separate the full ordered scalar block and distinct tensor/vector doublets. Literal four-dimensional Weyl and spatial Ricci variations give different full finite-transfer nonscalar curvature actions; the whole current classical/profile action and both canonical chain factors remain. Complete ordered Gaussian vector Ward integrands cancel only after all metric second-vertex contacts, with nonzero one-current means kept. Distinct retarded-source and advanced-detector primitives give the retarded weighted shift residual. The same full two-mass shear factor, not a sum of inverses or old Proca pole, inverts the full-mode strong remainder by an ordered weighted Volterra argument with finite unevaluated physical constants. Both original vector constraints follow on the explicit compatible force graph and unchanged initial germ. With S248 this gives the direct sum of the stated conditional metric/M1 amplitude and compatible-force subgraphs. It does not give arbitrary same-norm shift-force surjectivity, evaluated smallness or stability, interacting light/mixed loops, nonlinear same-state bounce, quantum gravitational limit, physical UV/Regge or original V/G/B/P8 closure.",
        "not_established": [
            "Arbitrary same-Hr shift-force or zero-spatial-vector shift norm on a full Fourier ball",
            "Uniform all-external-momentum bounds, unrestricted completed-graph surjectivity or a literal homogeneous vector theorem",
            "Evaluated physical inverse smallness, stability or finite inhomogeneous nonlinear same-state bounce",
            "Interacting light/mixed quantum state, omitted-loop bounds or quantization commuting with the gravitational limit",
            "Physical cutoff, exact UV S matrix, finite-gravity Regge remainder or all-parent exclusion",
            "Original V/G/B/P8 closure or machine formalization of the continuum arguments",
        ],
        "verification_boundary": "Exact full representation, literal curvature, whole weighted local action, both ordered vector Ward/contact, fixed-dimensional matching, ordered inverse and physical constraint identities support the written full-mode and Volterra proofs. Independent Christoffel curvature, arbitrary-axis projector, unexpanded metric, unequal-time weighted-transpose, infrared force, and noncommuting causal-matrix tests check the formulas without replacing continuum proofs. Native/direct/ordinary/CLI use original SymPy; only full regression uses the separately audited exact-GCD adapter. Every frozen scientific byte remains unchanged. The continuum arguments are not FORMALIZED.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The actual nonscalar two-mass inverse report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.249 actual nonscalar conditional inverse replay passed; unrestricted physics and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
