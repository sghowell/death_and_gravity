"""Read-only prepared Ward reconstruction and remaining scalar response boundary."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_full_adm_vertices import verify as parent
from p8_vacuum_affine_proca_gaussian import verify as gaussian

from . import audit, kinematics, norms, scalar, ward

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates" / "polynomial-vacuum-affine-prepared-ward-reconstruction.json"
)
PARENT_SHA = "2d2600e003055aec80398fdafa69e8e01a54c7c52e915832a51a0d23aa12db01"
GAUSSIAN_SHA = "1175702d6980db090396863b807af0e767241daeb7d8dabf138ed11c67e1ae75"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_prepared_ward_reconstruction/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {key: value for key, value in d.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    for previous, expected in ((parent, PARENT_SHA), (gaussian, GAUSSIAN_SHA)):
        if sha(previous.REPORT) != expected:
            raise ValueError(
                "The source-pinned actual affine Gaussian state chain changed"
            )
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_214_fully_rebuilt": PARENT_SHA,
        "S6_176_actual_CD_Gaussian_sector_transitively_rebuilt": GAUSSIAN_SHA,
        "fixed_stress_prescription_and_full_source_retained": True,
        "all_frozen_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A prepared Ward reconstruction gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.215.PREPARED_ORDERED_WARD_RECONSTRUCTION_KNOWN_PIECE_BOUND_AND_THREE_REMAINING_SCALAR_KERNELS",
        "date": "2026-09-12",
        "status": "PREPARED_ORDERED_WARD_RECONSTRUCTION_AND_KNOWN_PIECE_BOUNDS; NOT_MISSING_SCALAR_MATCHING_FULL_RESPONSE_REDUCED_INVERSE_BACKGROUND_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/kinematics.md",
            "notes/ward.md",
            "notes/boundary.md",
            "notes/scalar.md",
            "notes/norms.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "prepared_synchronous_maps_ordered_Ward_and_chart_terms": serialize(
            {"maps": payload(kinematics.data()), "Ward": payload(ward.data())}
        ),
        "three_scalar_reduction_and_explicit_known_piece_norms": serialize(
            {"scalar": payload(scalar.data()), "norms": payload(norms.data())}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The full conditional Gaussian metric response problem is reduced by an exact prepared synchronous decomposition and ordered Ward reconstruction. Retarded source primitives preserve the initial state germ; advanced detector primitives kill the final flux, while the initial differentiated flux vanishes from the source and current tangent. Direct physical-metric Lie derivatives fix every lapse, shift and spatial component without inverse H or momentum. The actual contravariant weight-one one-point density and complete nonlinear ADM chart contacts remain. Source covariance and detector conservation are derived separately, without making a retarded response a symmetric ordinary Hessian or asserting diffeomorphism invariance of a sharp computational band. The original causal proofs extend the matched tracefree response to the needed initial-only source domain. Rotation covariance leaves exactly three additional ordered scalar kernels, not two; both cross entries vanish at zero transfer but the homogeneous trace anchor remains. Explicit Cauchy/primitive/projection estimates and original one-point bounds control the known tracefree-plus-Ward piece by1e118 V04 U138. The displayed full estimate is conditional on three missing scalar matching/norm inputs and does not assert them. Independent inverse-metric/determinant volume-current variation, complete action reconstruction, endpoint and contact counterexamples and normalized rotational checks support the algebra and written proofs. The actual nonlinear parent source, finite prescription, prepared state and corrected historical frontier remain unchanged. No full scalar continuum matching, genuinely reduced inverse, nonlinear quantum solution, remaining parent loop/cutoff matching, finite-gravity IR/Regge or original P8 closure is supplied.",
        "not_established": [
            "The three ordered scalar kernels, homogeneous trace anchor, their original state/time/UV/contact matching and numerical norms",
            "Strong stress-norm differentiability from a weak Weyl-observable relative-Cauchy-evolution theorem, or full operator existence from a formal Ward identity",
            "Genuinely reduced canonical scalar/mixed inverse, contraction, finite-amplitude spatial control or stable finite-coupling quantum background",
            "Remaining parent sectors/measure/loops, cutoff/heavy/threshold/omitted-order matching, canonical vacuum cuts/contour/truncation or finite-gravity IR/Regge",
            "Permission to reset the state, retune fixed coefficients, delete chart contacts or declare original V/G/B and P8 closed",
        ],
        "verification_boundary": "Exact synchronous ADM/metric algebra, ordered weight-one density Ward and full contact chain, differentiated endpoint flux, normalized rotation projection, exact Cauchy majorants and known-piece norm arithmetic. Independent covariant volume-action, inverse/determinant, boundary/contact negative controls, primitive norms and arbitrary-axis checks supplement the written causal/domain/covariance proofs. These continuous distributional and Sobolev arguments are not FORMALIZED. Native/direct/ordinary/CLI retain original SymPy; only full regression uses the audited exact-GCD adapter. Scalar construction/matching and full reduced feedback remain open.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The prepared Ward reconstruction report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.215.PREPARED_ORDERED_WARD_RECONSTRUCTION_KNOWN_PIECE_BOUND_AND_THREE_REMAINING_SCALAR_KERNELS replay passed; original V/G/B and P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
