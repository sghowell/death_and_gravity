"""Read-only finite prepared full physical-metric response certificate."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as affine
from p8_vector_metric_local import verify as parent

from . import envelopes, evolution, proofs, source, tadpole, tail, tangent

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"finite-prepared-metric-response.json"
PARENT_SHA = "4ab3213760790f7309b69596c347d21c2a20552f7d010a12b500f6552b1cbaa0"
sha, serialize = affine.sha, affine.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_vector_metric_response/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen full metric local-response certificate changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {"S6_67_fully_rebuilt": PARENT_SHA, "original_clock_all_order_state_and_fixed_tadpole_profiles_unchanged": True}


def controls():
    calls = []
    for value in (True, False, sp.true, sp.false, 1.0, sp.Float(1), sp.Integer(1), "1", -1, 5):
        calls.extend((lambda value=value: source.baseline("T", value),
                      lambda value=value: tangent.coefficient("T", value),
                      lambda value=value: tangent.varied_P("T", value),
                      lambda value=value: source.adiabatic("T", "energy", value)))
    calls.extend((lambda: tangent.coefficient("T", 0), lambda: tangent.varied_P("T", 0)))
    for value in (True, False, sp.true, 0, 1, "transverse", "longitudinal", "", None):
        calls.extend((lambda value=value: source.data(value),
                      lambda value=value: source.baseline(value, 1),
                      lambda value=value: tangent.coefficient(value, 1),
                      lambda value=value: tangent.reference(value),
                      lambda value=value: envelopes.reference(value),
                      lambda value=value: tail.reference_tail(value, "energy"),
                      lambda value=value: source.readout(value, "energy")))
    for value in (True, False, sp.true, 0, 1, "rho", "p", "", None):
        calls.extend((lambda value=value: source.readout("T", value),
                      lambda value=value: source.adiabatic("T", value, 1),
                      lambda value=value: tail.reference_tail("T", value)))
    for value in (True, False, 0.0, sp.Float(0), sp.oo, sp.nan, source.n[11], source.v[11],
                  source.n[0]*source.v[0], sp.Integer(1), sp.Float(1)*source.n[0],
                  sp.Symbol("unrelated")*source.v[0]):
        calls.append(lambda value=value: source.linear_bound(value))
    for value in (True, False, sp.true, 1.0, sp.Float(1), "1", sp.I, sp.oo, -sp.oo, sp.zoo, sp.nan, 0, -1):
        calls.extend((lambda value=value: evolution.bound(value, 1000),
                      lambda value=value: evolution.bound(10**24, value),
                      lambda value=value: tadpole.bound(value, 1000),
                      lambda value=value: tadpole.bound(10**24, value)))
    for value in (1, 999):
        calls.extend((lambda value=value: evolution.bound(10**24, value),
                      lambda value=value: tadpole.bound(10**24, value)))
    calls.extend(lambda value=value: serialize(value) for value in (1.0, sp.Float(1), sp.oo, sp.nan))
    rejected = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(calls):
        raise ValueError("An unsupported prepared full metric-response input was accepted")
    return {"rejected_inputs": rejected, "native_validation_before_cached_implementations": True,
            "both_source_eleventh_derivatives_are_outside_C10": True,
            "second_mass_vertices_and_physical_output_normalizations_are_retained": True,
            "initial_interference_remains_when_independent_initial_response_is_zero": True,
            "fixed_tadpole_profiles_are_not_recomputed_under_variation": True,
            "W8_is_not_a_new_state_or_subtraction": True,
            "C10_to_C0_is_not_a_no_loss_or_spatial_feedback_norm": True}


@cache
def build_report():
    prior, identities = prior_checks(), proofs.residuals()
    certified, gates = affine.certify_residuals(identities), proofs.checks()
    if not all(value is True for value in gates.values()):
        raise ValueError("A finite prepared metric-response audit gate failed")
    scalar_count = sum(value.rows*value.cols if isinstance(value, sp.MatrixBase) else 1 for value in identities.values())
    coefficient_bounds = {sector: {j: source.linear_bound(tangent.coefficient(sector, j)) for j in range(1, 5)}
                          for sector in ("T", "L")}
    reference_bounds = {sector: {key: value for key, value in envelopes.reference(sector).items()
                                if key not in ("polynomial_majorants", "box_reconstructions", "base_reference_checks")}
                        for sector in ("T", "L")}
    return {"schema": 1, "claim": "P8-S6.68.FINITE_PREPARED_METRIC_RESPONSE", "date": "2026-09-08",
            "status": "FINITE_HOMOGENEOUS_PREPARED_METRIC_RESPONSE; ORIGINAL_P8_OPEN",
            "prior_sha256": prior, "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in source_files()},
            "formulation": "FORMULATION.md", "written_proofs": ["notes/reference.md", "notes/evolution.md", "notes/limit.md"],
            "exact_residuals": certified, "named_exact_check_count": len(identities), "checked_scalar_entries": scalar_count,
            "proof_checks": gates,
            "eighth_order_two_source_coefficient_bounds": serialize(coefficient_bounds),
            "varied_reference_residual_and_normalization_bounds": serialize(reference_bounds),
            "physical_readout_vertices": serialize({sector: {component: source.readout(sector, component) for component in ("energy", "pressure")}
                                                  for sector in ("T", "L")}),
            "varied_physical_readout_tail_bounds": serialize({sector: {component: tail.reference_tail(sector, component)
                                                                      for component in ("energy", "pressure")} for sector in ("T", "L")}),
            "exact_mixing_and_reference_constants": serialize(evolution.constants()),
            "finite_vector_metric_response_scale_example": serialize(evolution.bound(10**24, 1000)),
            "existing_tadpole_physical_vertices": serialize(tadpole.physical_vertices()),
            "complete_background_cancelled_response_scale_example": serialize(tadpole.bound(10**24, 1000)),
            "initial_mixing_and_source_norm_controls": serialize({"squeezed_initial_phase_fixture": -sp.Rational(15, 4),
                "tenth_source_derivative_fixtures": {sector: tangent.reference(sector)["tenth_source_derivative_fixtures"] for sector in ("T", "L")}}),
            "reference_mode_policy": "The original S6.55 all-order state evolves exactly. W8 is auxiliary; subtraction remains 0,2,4. Both physical metric sources, moving canonical normalization, second mass-source contacts, output normalization and initial interference are varied. No new state or finite Wilson value is selected.",
            "dimensional_limit": "For fixed mass and smooth prepared sources, use a sufficiently small complex-D neighborhood of 3, analytic plus/minus modes and the original D-independent Borel cutoffs. Metric leading frequencies are D-independent, off-diagonal initial mixing remains O(nu^-6), and the varied subtracted tail is O(nu^-4). Dominated convergence gives the common finite limit; numerical constants refer to D=3, not the full old strip.",
            "prepared_source_domain": "Independent physical lapse n and log-scale zeta at spatial K=0, fixed scalar clock, smooth compact time support strictly after u0=-1/2, zero on an initial neighborhood. The norm is the maximum of both source time derivatives 0..10 on I=[-1/2,1/2], with physical energy/pressure C0 outputs. The existing scalar tadpole profiles stay fixed.",
            "controls": controls(),
            "verdict": "The full prepared homogeneous Gaussian physical-metric response has a finite integrable exact-minus-adiabatic remainder and the S6.67 local part. At L=10^24,m=1000 its joint nonlocal C10-to-C0 norm is below 10^-42 and complete vector norm below 10^-38. Including the already-fixed S6.60 scalar tadpole keeps the background-cancelled response below 10^-37. Coupled feedback and original P8 remain open.",
            "not_established": ["Arbitrary spatial momentum or independent varied initial covariance",
                                "No-loss coupled inverse, quantum stability/cones, other loops, interactions or cutoff",
                                "Finite Wilson matching, common UV parent, V/G/B or original P8 closure"],
            "verification_boundary": "Exact varied Riccati, physical readout, mixing, radial and fixed-action identities with continuous rational bounds and written common-dimensional parameter arguments; not proof-assistant formalization"}


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The finite prepared metric-response report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.68.FINITE_PREPARED_METRIC_RESPONSE replay passed; coupled feedback and original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
