"""Read-only corrected retarded phase and full ordered scalar UV matching."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_prepared_ward_reconstruction import verify as parent
from p8_vacuum_affine_proca_gaussian import verify as gaussian

from . import audit, density, geometry, jets, matching

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-ordered-scalar-symbol.json"
PARENT_SHA = "bba0d9e34fd6645444a52940587dced5500746a8fb49124434cce1ad819059a5"
GAUSSIAN_SHA = "1175702d6980db090396863b807af0e767241daeb7d8dabf138ed11c67e1ae75"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_ordered_scalar_symbol/*.py"))
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
            raise ValueError("A source-pinned historical Gaussian input changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_215_historical_report_fully_rebuilt": PARENT_SHA,
        "S6_176_original_Gaussian_state_transitively_rebuilt": GAUSSIAN_SHA,
        "historical_replay_not_endorsement_of_withdrawn_physical_identifications": True,
        "no_frozen_scientific_test_or_report_bytes_changed": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A corrected phase/scalar UV gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.216.CORRECTED_RETARDED_PAIR_PHASE_FULL_ORDERED_SCALAR_UV_AND_FIXED_DIMENSIONAL_FINITE_PART",
        "date": "2026-09-12",
        "status": "CORRECTED_RETARDED_PAIR_PHASE_FULL_ORDERED_SCALAR_SPATIAL_UV_FINITE_PART; NOT_FULL_CURRENT_REASSEMBLY_SCALAR_ANCHOR_REDUCED_INVERSE_BACKGROUND_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/phase.md",
            "notes/vertices.md",
            "notes/geometry.md",
            "notes/curvature.md",
            "notes/matching.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "corrected_phase_full_trace_vertices_and_dimensional_geometry": serialize(
            {
                "phase": payload(jets.phase_data()),
                "vertices": payload(jets.data()),
                "geometry": payload(geometry.data()),
            }
        ),
        "ordered_scalar_UV_fixed_pole_and_finite_matching": serialize(
            {"density": payload(density.data()), "matching": payload(matching.data())}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "An independent finite Fock and real symplectic covariance audit fixes the actual annihilation-source retarded orientation: sharp detector, positive detector phase and +Im. The S198 creation-amplitude identity is exact, but S207/S212 combined it with annihilation amplitudes; every odd endpoint must reverse. The full arbitrary-trace spatial Proca pair includes the temporal-constraint amplitude and its two ordered cross terms and square. Ten field-strength geometries give350 radial coefficients and105 new scalar channel entries. Six invariant reconstruction precedes the fixed-source dimension derivative. Full exponential-volume curvature Hessians, continued Euler and the original prescribed pole reproduce all physical corrected logarithmic coefficients. The same original radial MSbar finite part includes lower-band, dimensional-mode and full volume terms. Its corrected tracefree difference is p^2(3t^2+1)(-31T+30V)Gamma0/(420pi^2), and the full local spatial finite operator obeys1e6||D||L2 Z24. Ordered off-diagonal coefficients need not coincide; they obey the proper-time formal Green relation. Independent literal ADM polarization pairs, integer-dimensional field strengths and metric curvature, full four-order WKB coefficient extraction and wrong-phase/volume controls supplement the exact symbolic proof. Four historical status rows explicitly record the branch error and withdrawn physical identifications. No full scalar current, homogeneous trace anchor, corrected S213 finite endpoint/bulk/regulator assembly or S215 known-input restoration is inferred. The original state, parent and finite prescription remain fixed; full reduced response/inverse, nonlinear quantum background, remaining parent/loop/heavy/cutoff and canonical vacuum/finite-gravity obligations remain open.",
        "not_established": [
            "The complete scalar current, homogeneous trace anchor and full actual state/time/UV/contact remainders with norms",
            "Restoration of S213 physical current/anchored-regulator identification or S215 known tracefree input before consistent finite endpoint and bulk reassembly",
            "A retarded kernel symmetric under source/readout exchange from a local finite operator's formal Green identity",
            "A new finite counterterm, different state, changed mass/parent, metric initial-state reset or deletion of trace/chart constraints",
            "Full reduced canonical inverse, finite-amplitude response, stable quantum background, remaining parent matching/cutoff or original V/G/B and P8 closure",
        ],
        "verification_boundary": "Finite CCR/covariance and retarded ODE algebra fix phase independently of UV matching. Full trace field-strength and ADM vertices, exact six-invariant general-dimensional extraction, literal curved-volume checks, independent four-order WKB, original finite prescription and local coefficient majorants are checked. Continuous dimensional, covariant and Sobolev arguments are written proofs, not FORMALIZED. Native/direct/ordinary/CLI retain original SymPy; only full regression uses the audited exact-GCD adapter. Historical report replay does not restore withdrawn physical claims.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The corrected phase/scalar UV report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.216 corrected pair phase and ordered scalar UV replay passed; withdrawn current identifications and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
