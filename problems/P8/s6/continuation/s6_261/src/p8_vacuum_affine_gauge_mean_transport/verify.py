"""Read-only full gauge-orbit and fixed-reference mean-transport certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_spatial_brst import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-gauge-mean-transport.json"
PARENT_SHA = "323cfcb726f73905cc22c63d7d935e7149cdf04a3dc9964ffdb23202ef19c85c"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_gauge_mean_transport/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen complete spatial/projective BRST parent changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    return {
        "S6_260_complete_BRST_Ward_parent_and_whole_ancestry_rebuilt": PARENT_SHA,
        "same_parent_profiles_sources_preparations_and_canonical_boundaries": True,
        "finite_Weyl_gauge_jets_not_interacting_Nielsen_mean": True,
        "original_quantum_ordering_regulator_cutoff_and_P8_still_open": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError(
            "A whole gauge-orbit or fixed-reference mean-transport gate failed"
        )
    packets = audit.packets()
    geometric = tuple(
        name
        for name in packets
        if name
        not in (
            "whole_unchanged_coupled_reference_rows_and_commutator",
            "whole_actual_bounce_regular_canonical_covariance_symbols",
        )
    )
    return {
        "schema": 1,
        "claim": "P8-S6.261.WHOLE_SPATIAL_GAUGE_WEIGHT_ORBIT_FIXED_COUPLED_REFERENCE_MEAN_TRANSPORT_PHYSICAL_VOLUME_CANCELLATION_AND_ACTUAL_BOUNCE_COVARIANCE_SYMBOLS",
        "date": "2026-09-14",
        "status": "EXACT_FULL_CLASSICAL_GAUGE_ORBIT_FIXED_REFERENCE_WEYL_JETS_AND_BOUNCE_SYMBOLS; NOT_INTERACTING_NIELSEN_MEAN_ORDERING_REGULATOR_CUTOFF_ORIGINAL_V_G_B_OR_P8",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/orbit.md",
            "notes/fourier.md",
            "notes/reference.md",
            "notes/crossing.md",
            "notes/volume.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_spatial_gauge_orbit_Fourier_and_physical_volume_transport": serialize(
            {name: payload(packets[name]) for name in geometric}
        ),
        "whole_fixed_coupled_reference_and_actual_bounce_symbols": serialize(
            {name: payload(packets[name]) for name in packets if name not in geometric}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The unchanged spatial density gauge has an explicit weight-varying orbit. Its exact periodic monotone map retains two tensor profiles, and its second-order three-direction displacement keeps all generated Fourier modes. The unchanged full coupled Gaussian covariance gives the finite Weyl-jet coordinate means -9 Cvv/4 and -9 Cvn/4; both are retained in the physical lapse-dependent volume identity and cancel its two-point variation. The generic reduced lapse/volume commutator is nonzero, vanishing only at the linear bounce in this calculation. The complete regular central chart and whole fixed-profile lapse constraint yield explicit actual-reference bounce ultraviolet symbols: metric variance k^-3, lapse cross covariance k^-1, and lapse variance k. The corresponding radial logarithmic, quadratic and quartic asymptotics do not determine a physical cutoff, finite interacting error, nonlinear quantum ordering, actual Nielsen vector, interacting mean or P8 closure.",
        "not_established": [
            "The actual interacting Nielsen vector, vanishing Ward/state/endpoint/regulator defects or a gauge-independent interacting physical source functional",
            "A nonlinear quantum ordering prescription, complete BRST regulator or Gaussian support inside the nonlinear auxiliary branch",
            "An interacting fixed quantum mean, physical curved subtraction or quantitative omitted-loop and nonlinear correction bounds",
            "A finite-Fourier Lie algebra, global general three-dimensional gauge slice, infinite-volume inverse gap or prescription for the homogeneous constraints",
            "A physical Wilsonian cutoff from the unrenormalized reference asymptotics, controlled matching or a cutoff identified with the heavy mass",
            "A compatible nonlinear bounce, quantum gravitational limit, vacuum UV scattering or finite-gravity infrared/Regge bound",
            "Formalization of the written local gauge, fixed-reference symbol and asymptotic integration arguments",
        ],
        "verification_boundary": "Exact full identities and independent quadrature, Gaussian cubature, Fourier, Fock and full-matrix comparisons support the scoped classical and fixed-reference results. Finite comparison matrices are diagnostics, not a reselected physical state. Every earlier source, preparation and report is unchanged. Native/direct/ordinary/CLI use original SymPy; only captured full regression uses the audited exact-GCD adapter.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The full gauge-mean-transport report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.261 complete finite gauge-orbit and fixed-reference mean transport replay passed; interacting mean, cutoff and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
