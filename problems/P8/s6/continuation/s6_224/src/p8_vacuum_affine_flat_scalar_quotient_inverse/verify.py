"""Read-only uniform full flat scalar quotient reference inverse."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_flat_spatial_conversion import verify as spatial
from p8_vacuum_affine_isolated_shear_resolvent import verify as parent

from . import audit, estimates, geometry, resolvent

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates" / "polynomial-vacuum-affine-flat-scalar-quotient-inverse.json"
)
PARENT_SHA = "c7c0eb49287a4456328fe86b20a8362b70f540309439ba98b74ad1f191e0c1e6"
SPATIAL_SHA = "739ab8954fc755149312818d2c6db8e5d38d943a65cb085160a337891bf30507"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_flat_scalar_quotient_inverse/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {k: v for k, v in data.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    for prior, expected in ((parent, PARENT_SHA), (spatial, SPATIAL_SHA)):
        if sha(prior.REPORT) != expected:
            raise ValueError(
                "A source-pinned isolated factor or full spatial projector input changed"
            )
        prior.validate_report(
            json.loads(prior.REPORT.read_text()), prior.build_report()
        )
    return {
        "S6_223_fixed_shear_pole_plus_cut_and_finite_window_inverse_fully_rebuilt": PARENT_SHA,
        "S6_200_complete_nonzero_transfer_spatial_projectors_fully_rebuilt": SPATIAL_SHA,
        "S6_86_trace_and_S6_189_original_finite_Hessian_transitively_rebuilt": True,
        "S6_222_full_quantum_force_graph_transitively_rebuilt_not_inverted_here": True,
        "original_state_prescription_and_all_frozen_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError(
            "A full quotient, inverse or uniform normalization gate failed"
        )
    return {
        "schema": 1,
        "claim": "P8-S6.224.FULL_FLAT_SCALAR_GAUGE_QUOTIENT_TWO_CHANNEL_CAUSAL_REFERENCE_INVERSE_AND_UNIFORM_ALL_MOMENTUM_BOUND",
        "date": "2026-09-12",
        "status": "FULL_FLAT_QUOTIENT_REFERENCE_INVERSE_NO_SPATIAL_DERIVATIVE_LOSS; NOT_FULL_CURVED_COUPLED_INVERSE_KAPPA_SMALLNESS_STABILITY_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/geometry.md",
            "notes/coordinates.md",
            "notes/spectral.md",
            "notes/inverse.md",
            "notes/spaces.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "full_covariant_projectors_ordered_gauge_and_reference_factorization": serialize(
            payload(geometry.data())
        ),
        "complete_causal_matrix_inverse_uniform_spaces_and_physical_units": serialize(
            {
                "resolvent": payload(resolvent.data()),
                "bounds": payload(estimates.data()),
            }
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The full nonzero-transfer covariant projectors and unchanged finite Weyl/scalar-curvature Hessian give the complete two-channel flat spatial scalar quotient B^T diag(Ftrace(D^2+q),(8/3)F2(D^2+q)) B. Both opposite time/Fourier gauge legs are retained, and the original three-source block has a gauge kernel rather than an inverse. The causal curvature-coordinate matrix B has an explicit inverse using wave and double-primitive kernels, with a continuousq0 extension and no inverse-transfer singularity. Its kernel Rq vanishes at the initial boundary and has first derivative norm below5/2 uniformly inq. The actual shifted trace and shear spectral measures, including the full shear pole, give a diagonal primitive Jq bounded by45/2. The ordinary matrix convolution Eq=Rq'*Jq*Rq^T has both reference inverse identities and norm at most375t^3/16. Consequently the normalized isolated reference inverse is bounded by375T^4/64<6T^4 onC_tH^r andL2_tH^r for every realr, uniformly over all spatial momenta without derivative loss. Its first time derivative is also bounded. The causal graph includes the full initial boundary; no initial atoms or quantum pole are deleted. The physical Hessian inverse retains64pi^2, and the flat S222 force normalization also retainskappa, giving375pi^2 kappa T^4, not kappa smallness. Three independent rotated four-dimensional curvature/projector fixtures, ordered gauges, direct matrix solves, time-domain forced wave solutions and both shifted continuum densities supplement the written causal and functional proof. This is a full isolated flat quotient reference inverse, not the complete curved/state/contact/tree/matter quantum-force inverse or a finite-coupling/nonlinear parent remainder, stability, cutoff or original V/G/B/P8 closure.",
        "not_established": [
            "The full three-source inverse before removing its gauge kernel, or a literal homogeneous lapse/shift solution",
            "A complete curved/state/contact/tree/matter normal form or compatible-space inverse of the S222 quantum-force system",
            "A kappa-small feedback estimate, uniform half-line scalarL1 bound, pole deletion or order reduction",
            "Finite-coupling/nonlinear parent control, quantum background/stability, heavy/physical cutoff or original V/G/B/P8 closure",
        ],
        "verification_boundary": "Exact full tensor/gauge factorization, finite local normalization, both rational inverse products, source-domain coordinate limits and all convolution constants are checked. Causal distribution graph inversion, spectral primitive bounds, Plancherel/Minkowski and dominated continuity are written proofs, not FORMALIZED. High-precision shifted dispersion is an independent diagnostic, not a continuum existence proof. Native/direct/ordinary/CLI retain original SymPy; only full regression uses the audited exact-GCD adapter. All frozen predecessors are unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The flat scalar quotient inverse report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.224 full flat scalar quotient inverse replay passed; curved coupled inverse and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
