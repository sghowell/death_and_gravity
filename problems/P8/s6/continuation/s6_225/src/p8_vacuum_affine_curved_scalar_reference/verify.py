"""Read-only actual curved scalar geometry and ordered reference inverse."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_flat_scalar_quotient_inverse import verify as parent
from p8_vacuum_affine_local_tensor_response import verify as original

from . import audit, geometry, local, resolvent

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-curved-scalar-reference.json"
PARENT_SHA = "a2287a5083d5dceab6886f4f7de3304ce3a7af56c772e76dc869115ad8ac3fb5"
ORIGINAL_SHA = "567a8a84e981c637b000ef671c350b63c54f4fc357ca60eaf1ce1de9c1293b02"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_curved_scalar_reference/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    for prior, expected in ((parent, PARENT_SHA), (original, ORIGINAL_SHA)):
        if sha(prior.REPORT) != expected:
            raise ValueError(
                "A source-pinned flat quotient or original finite local input changed"
            )
        prior.validate_report(
            json.loads(prior.REPORT.read_text()), prior.build_report()
        )
    return {
        "S6_224_full_flat_scalar_quotient_uniform_inverse_fully_rebuilt": PARENT_SHA,
        "S6_189_unchanged_original_finite_local_Hessian_fully_rebuilt": ORIGINAL_SHA,
        "S6_223_original_shear_pole_and_S6_86_trace_factor_transitively_rebuilt": True,
        "S6_222_full_quantum_force_graph_transitively_rebuilt_not_inverted_here": True,
        "original_state_prescription_and_all_frozen_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError(
            "An actual geometry, local Hessian or ordered inverse gate failed"
        )
    return {
        "schema": 1,
        "claim": "P8-S6.225.ACTUAL_CURVED_SCALAR_CURVATURE_COORDINATES_COMPLETE_FINITE_LOCAL_HESSIAN_AND_ORDERED_UNIFORM_REFERENCE_INVERSE",
        "date": "2026-09-12",
        "status": "ACTUAL_LOCAL_GEOMETRY_AND_CURVATURE_ADAPTED_REFERENCE_INVERSE; NOT_ACTUAL_FULL_CURVED_QUANTUM_INVERSE_KAPPA_SMALLNESS_STABILITY_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/geometry.md",
            "notes/local.md",
            "notes/coordinates.md",
            "notes/adjoint.md",
            "notes/inverse.md",
            "notes/spaces.md",
            "notes/scope-validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "actual_curved_geometry_and_complete_original_local_Hessian": serialize(
            {
                "geometry": payload(geometry.data()),
                "local": payload(local.data()),
            }
        ),
        "ordered_curvature_coordinate_reference_inverse_and_physical_density": serialize(
            payload(resolvent.data())
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The literal four-dimensional exponential spatial scalar chart on the unchanged a=(1+t^2)^2 background gives S=a^2 deltaR_old/6 and W=3C0101/(a^2 phase), with Bc=B0+C D in fixed conformal time. Continuous whole-slab bounds give ||C||<8, ||Cprime||<35 and conformal length at most1. The complete original finite local scalar Hessian is the curvature-square factor -4SD SG-(4/45)WD WG plus the explicitly retained background-curvature, second-variation, Einstein and volume terms of total differential order at most2. Its homogeneous unit-Frobenius tracefree restriction agrees with the source-pinned S189 fourth and kinetic terms. Both causal curvature-coordinate inverses exist by separate Volterra velocity equations; the formal conformal-density adjoint retains -D(C^T), including its coefficient derivative. Their first-derivative kernels are bounded by (5/2)exp(20 elapsed) and (5/2)exp(108 elapsed). The ordered reference Aref=Bc* diag(Ftrace(D^2+q),(8/3)F2(D^2+q)) Bc has inverse Y Kdiag Z, represented by an ordinary Y Jdiag partial Z kernel. The retained shifted trace and pole-plus-cut shear primitives give kernel bound375 elapsed^3 exp(108 elapsed)/16 and inverse norm375T^4 exp(108T)/64<6T^4 exp(108T) on C_tH^r and L2_tH^r, every realr and all comoving momenta without spatial derivative loss. Both inverse products include the full causal initial boundary. Conformal physical force normalization keeps right multiplication by64pi^2 kappa a^4, yielding a bound below2250pi^2 kappa T^4 exp(108T), not kappa smallness. Independent complete linearized Weyl components, an ADM/warped-metric local action reconstruction, actual interval bounds, integrated adjoints and noncommuting finite causal fixtures supplement the written functional proof. The scalar factors remain the specified flat reference factors. Actual curved mass/state/contact/tree/matter differences, the full S222 force graph inverse and original V/G/B/P8 closure are not established.",
        "not_established": [
            "Identification of the specified flat scalar factors with the full actual curved mass and state response",
            "Compatible-space boundedness or smallness of the complete local, state, contact, tree and matter remainder",
            "A complete S222 auxiliary/clock graph bridge, full quantum inverse or finite-coupling/nonlinear parent remainder",
            "Pole deletion, order reduction, small backreaction, stability, physical cutoff or original V/G/B/P8 closure",
        ],
        "verification_boundary": "Exact actual geometry, complete finite Hessian, original tensor bridge, adjoint boundary and all convolution constants are checked. Causal variable-coefficient Volterra inversion, distribution graph extension and uniform function-space estimates are written proofs, not FORMALIZED. Finite causal matrix controls establish ordering diagnostics, not a continuum discretization limit. Native/direct/ordinary/CLI retain original SymPy; only full regression uses the audited exact-GCD adapter. All frozen predecessors are unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The curved scalar reference report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.225 curved scalar reference replay passed; actual full quantum inverse and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
