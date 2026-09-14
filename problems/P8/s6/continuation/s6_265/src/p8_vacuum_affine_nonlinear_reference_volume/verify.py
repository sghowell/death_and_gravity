"""Read-only full nonlinear reference composites and exact integrability boundary."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_quantitative_gaussian_window import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates/polynomial-vacuum-affine-nonlinear-reference-volume.json"
PARENT_SHA = "7460bf04a813c31f844f73be22ca1d5d44e6d6f2dd3ffe4fbd25fa35e21aef8d"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_nonlinear_reference_volume/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen unchanged-Gaussian covariance report changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    return {
        "S6_264_complete_fixed_Gaussian_and_whole_ancestry_rebuilt": PARENT_SHA,
        "same_original_R_both_switches_state_preparation_and_canonical_maps": True,
        "S261_false_physical_binding_is_not_renewed": True,
        "same_original_profiles_masses_H_Proca_states_and_determinants": True,
        "reference_composite_and_coefficient_obstruction_not_interacting_mean_or_action_no_go": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError(
            "A complete source, nonlinear remainder, covariance or integrability gate failed"
        )
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.265.COMPLETE_POSITIVE_LAPSE_NONLINEAR_REFERENCE_COEFFICIENT_AND_COMMUTING_BOUNCE_VOLUME_MOMENTS_WITH_EXACT_GAUSSIAN_INDIVIDUAL_COEFFICIENT_INTEGRABILITY_BOUNDARY_WITH_UNCHANGED_PARENT",
        "date": "2026-09-14",
        "status": "FULL_NONPOLYNOMIAL_REFERENCE_COMPOSITES_AND_SINGULAR_INDIVIDUAL_COEFFICIENT_OBSTRUCTION; NOT_NONLINEAR_AUXILIARY_RECONSTRUCTION_FULL_ACTION_NO_GO_INTERACTING_MEAN_REGULATOR_CUTOFF_OR_ORIGINAL_V_G_B_P8",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/derivatives.md",
            "notes/covariance.md",
            "notes/gaussian.md",
            "notes/integrability.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_coefficients_global_bounds_and_full_lapse_jets": serialize(
            {name: payload(packets[name]) for name in names[:3]}
        ),
        "whole_complete_reference_composites_and_integrability_boundary": serialize(
            {name: payload(packets[name]) for name in names[3:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The entire original R with both actual switches gives bounded positive-lapse inverse-power spectral factors. Under the SAME original band-limited Gaussian reference, their full nonlinear means differ from the exact quadratic center moments by less than1e-960 on the tiny slab. At the bounce the complete hat-log-scale and lapse Weyl smears strongly commute, so the full positive-lapse reference volume exp(3v_hat)U(1+n) has a positive integrable joint spectral meaning. Its mean differs from1 plus the complete S264 quadratic contact by less than1e-960, keeping the cross covariance. The actual lapse variance is above1e-576 throughout the slab and1e-572 at the bounce, not zero. Conversely unrestricted positive-lapse Gaussian functional substitution gives an infinite mean for the original INDIVIDUAL ADM kinetic coefficient M/N, and a finite first but infinite second moment for the complete canonical matter coefficient N/U. This is a precise integrability obstruction to naive substitution, not a divergence or no-go theorem for the full constrained Hamiltonian, action or candidate. No original source, state, H/Proca preparation, determinant, quantum ordering or physical cutoff is changed; original V/G/B/P8 remain OPEN.",
        "not_established": [
            "A nonlinear auxiliary reconstruction identifying its lapse with the linear Gaussian reference lapse",
            "A full constrained Hamiltonian/action divergence, candidate inconsistency or original-P8 no-go from individual coefficient integrability",
            "A stationary interacting physical mean, full quantum gauge ordering/regulator or vanishing state/endpoint/Ward defects",
            "A projected, conditioned, support-restricted or newly selected replacement of the original Gaussian state",
            "Global auxiliary/Legendre pivot regularity from global positivity of the single coefficient R",
            "Joint spectral probabilities away from the commuting bounce slice or a uniform spacetime/infinite-volume event",
            "A physical Wilsonian cutoff, physical curved subtraction, complete omitted-loop bounds or full H/Proca interacting onepoint",
            "Nonlinear inhomogeneous stability, original V/G/B/P8 closure or reopening the completed scoped P8(a) result",
            "Formal verification of the written functional-calculus, Gaussian remainder, energy/covariance or integrability arguments",
        ],
        "verification_boundary": "The complete source and lapse weights, exact outward derivative bounds, full pure-covariance normalization, finite Gaussian moments and rational tail/remainder budgets are checked. Independent finite spectral, Gaussian-integral and singular-cutoff diagnostics support but do not replace the written arguments. The positive-lapse localization belongs to the stated observable, not to a changed state. Native/direct/ordinary/CLI use original SymPy; only captured complete P8 regression uses the separately audited exact-GCD adapter. All frozen predecessor bytes and the explicit S261 physical-binding erratum are preserved.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete nonlinear reference-volume and integrability report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.265 full nonlinear reference composites and individual-coefficient integrability replay passed; constrained quantum mean, cutoff and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
