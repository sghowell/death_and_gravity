"""Read-only source-pinned physical-volume correction and explicit erratum."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_gauge_mean_transport import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates" / "polynomial-vacuum-affine-physical-volume-correction.json"
)
PARENT_SHA = "4ca9381026403ff0b759c68c56c4683a5abf2610ccb2a3d949dad6d2e63e557e"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_physical_volume_correction/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The immutable historical S261 evidence changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    return {
        "S6_261_immutable_historical_evidence_and_whole_ancestry_rebuilt": PARENT_SHA,
        "S261_R_minus_half_physical_identification_is_explicitly_refuted": True,
        "historical_replay_is_not_renewed_physical_binding_certification": True,
        "same_original_action_profiles_sources_preparations_and_boundaries": True,
        "corrected_classical_source_jets_not_interacting_mean_or_P8_closure": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError(
            "A source-pinned physical metric or explicit correction gate failed"
        )
    packets = audit.packets()
    actual = tuple(
        name
        for name in packets
        if name != "whole_explicit_erratum_and_retained_reference_symbols"
    )
    return {
        "schema": 1,
        "claim": "P8-S6.262.SOURCE_PINNED_WHOLE_PHYSICAL_ADM_VOLUME_AND_MATTER_BINDING_CORRECTED_FIXED_REFERENCE_CONTACTS_AND_EXPLICIT_S261_PHYSICAL_IDENTIFICATION_REFUTATION",
        "date": "2026-09-14",
        "status": "CORRECTED_ACTUAL_CLASSICAL_PHYSICAL_VOLUME_BINDING_AND_WEYL_CONTACTS_WITH_EXPLICIT_S261_REFUTATION; NOT_QUANTUM_ORDERING_REGULATOR_INTERACTING_MEAN_CUTOFF_ORIGINAL_V_G_B_OR_P8",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/adm.md",
            "notes/jets.md",
            "notes/contact.md",
            "notes/erratum.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_actual_physical_metric_action_volume_and_corrected_contacts": serialize(
            {name: payload(packets[name]) for name in actual}
        ),
        "whole_explicit_S261_erratum_and_retained_reference_symbols": serialize(
            payload(packets["whole_explicit_erratum_and_retained_reference_symbols"])
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The literal unchanged parent has C=R**(-1/2)=Cchi**2=M**(-2), not R-1/2. Its full shifted ADM metric, inverse, scalar and Maxwell contractions, mass and original source density, and extended point Jacobian agree with the independent frozen coefficients. Its spatial density is U=R**(-3/4), with bounce lapse slope +3/2. The corrected finite fixed-reference physical onepoint contact is -27 Cvv/4-27 Cvn/8 and the two-point contact is its negative. The immutable S261 identification C=R-1/2 and its physical slope -6 are explicitly REFUTED. Generic cancellations can pass with either factor and are not physical binding tests. S261 generic orbits, Weyl coordinate means and bounce covariance symbols are retained within their original restricted scope. No original action, state or frozen evidence is edited. Quantum ordering, interacting physical mean, regulator, subtraction, actual cutoff, nonlinear bounce and original P8 remain OPEN.",
        "not_established": [
            "Renewed physical certification of the false S261 factor, slope or derived physical contact",
            "An interacting Nielsen vector or fixed mean, vanishing Ward/state/endpoint defects, physical subtraction or full gauge-independent quantum source functional",
            "A nonlinear quantum ordering, a complete BRST regulator or Gaussian support inside the nonlinear auxiliary branch",
            "A physical Wilsonian cutoff, controlled all-sector matching, finite omitted-loop bounds or identification of the heavy mass with the cutoff",
            "A compatible nonlinear bounce, quantum gravity limit, vacuum UV scattering or finite-gravity IR/Regge bound",
            "Formal verification of the full written classical, source-binding and reference-symbol arguments",
        ],
        "verification_boundary": "Independent full-metric and full-two-form contractions, finite-difference coordinate Jacobians and source Hessians, clock-jet comparisons, periodic inverse-map quadrature, correlated Gaussian cubature and unaliased Fourier sampling test the corrected classical and finite Weyl identities. Numeric and polynomial fixtures do not replace the unchanged physical parent or preparation. S261 source/report bytes remain immutable historical evidence with an explicit physical-binding refutation. Native/direct/ordinary/CLI use original SymPy; only captured full regression uses the audited exact-GCD adapter.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The source-pinned physical-volume correction report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.262 corrected physical-parent volume replay passed; S261 false binding REFUTED; interacting mean, cutoff and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
