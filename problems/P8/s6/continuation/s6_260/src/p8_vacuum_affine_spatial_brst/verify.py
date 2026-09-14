"""Read-only complete spatial/projective BRST and conditional source report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_complement_measure import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-spatial-brst.json"
PARENT_SHA = "c259d101b5029c4b50dacc29cfc9597653885bfd2ef172dbe6a0849119658c65"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_spatial_brst/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen complete affine-complement parent changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    return {
        "S6_259_full_connection_complement_and_complete_ancestry_rebuilt": PARENT_SHA,
        "same_current_parent_functions_profiles_sources_and_boundaries": True,
        "classical_BRST_and_conditional_Ward_not_replacement_quantum_preparation": True,
        "finite_algebra_not_anomaly_free_continuum_measure_or_original_P8": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A complete classical BRST or conditional source gate failed")
    packets = audit.packets()
    classical = tuple(
        name
        for name in packets
        if name
        not in (
            "whole_convergent_two_source_gauge_orbit",
            "whole_differentiated_conditional_Nielsen_identity",
        )
    )
    return {
        "schema": 1,
        "claim": "P8-S6.260.COMPLETE_RETAINED_SPATIAL_PROJECTIVE_BRST_COMPENSATED_TRACE_GAUGE_AND_CONDITIONAL_FIXED_SOURCE_WARD_MEAN_TRANSPORT_WITH_UNCHANGED_PARENT",
        "date": "2026-09-14",
        "status": "EXACT_FULL_CLASSICAL_SPATIAL_PROJECTIVE_BRST_AND_CONDITIONAL_SOURCE_WARD_TRANSPORT; NOT_COMPLETE_QUANTUM_REGULATOR_STATE_FIXED_MEAN_CUTOFF_ORIGINAL_V_G_B_OR_P8",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/algebra.md",
            "notes/projective.md",
            "notes/zero-modes.md",
            "notes/ward.md",
            "notes/state.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_classical_spatial_projective_BRST_and_gauge_source_bridge": serialize(
            {name: payload(packets[name]) for name in classical}
        ),
        "whole_conditional_Ward_mean_and_finite_orbit_diagnostics": serialize(
            {name: payload(packets[name]) for name in packets if name not in classical}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The unchanged parent has an explicit complete local spatial/projective classical BRST algebra on the regular reduced branch. The time-dependent shift, density weights and all 64 connection components are retained. The original trace gauge requires a compensating projective transformation; its full ghost shift gives contractible trace pairs and the finite triangular seven-ghost matrix. Every original source contact and physical map is retained. The complete canonical spatial flux and off-gauge ghost action are BRST closed. Translation projection and finite Fourier closure have explicit nonzero defects. Conditional fixed-source Ward/Nielsen identities retain state, measure, regulator and endpoint defects and the off-shell onepoint Hessian contact. A convergent two-source orbit integral gives a gauge-dependent composite coordinate mean without changing the physical Gaussian, while a noninvariant boundary weight changes a physical expectation. These statements do not evaluate the actual interacting P8 state, quantum measure, Nielsen vector, fixed mean, cutoff, nonlinear bounce or UV/Regge gates.",
        "not_established": [
            "An off-constraint BFV prescription for all unreduced auxiliary and shift-primary variables",
            "An anomaly-free continuum/time regulator, complete quantum ordering or renormalized spatial/covariant gauge measure",
            "A BRST-invariant interacting extension of the fixed free preparations or vanishing of all endpoint and state defects",
            "The actual P8 Nielsen vector, physical curved subtraction, interacting fixed mean or a bound on its correction",
            "A closed mean-zero or finite-Fourier diffeomorphism regulator, global gauge slice or inverse on translation modes",
            "A controlled Wilsonian cutoff/matching error, nonlinear bounce, quantum gravitational limit, physical UV scattering or finite-gravity IR/Regge bound",
            "Formalization of the written continuum algebra, conditional Ward integration and state-boundary arguments",
        ],
        "verification_boundary": "Exact full Grassmann jet identities and independent finite diagnostics support the written classical and conditional statements. External Horava or Higgs models do not supply this parent, state or quantum measure. Native/direct/ordinary/CLI use original SymPy; only captured full regression uses the audited exact-GCD adapter. Every earlier source, profile, preparation and report remains unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete spatial/projective BRST report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.260 complete classical spatial/projective BRST and conditional source Ward replay passed; full quantum measure, fixed mean and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
