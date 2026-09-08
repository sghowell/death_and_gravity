"""Read-only finite homogeneous prepared mass-response certificate."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as affine
from p8_vector_mass_adiabatic import verify as parent

from . import envelopes, evolution, proofs, source, tail, tangent

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"finite-prepared-mass-response.json"
PARENT_SHA = "1a6a16fcfd54827aefdfbe6e5a42f5a90fd001e611ae926c22bca01e29a66c70"
sha, serialize = affine.sha, affine.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_vector_mass_remainder/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen finite local mass-response certificate changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {"S6_65_fully_rebuilt": PARENT_SHA, "selected_clock_state_and_fixed_profiles_unchanged": True}


def controls():
    calls = []
    for value in (True, False, sp.true, sp.false, 1.0, sp.Float(1), sp.Integer(1), "1", -1, 5):
        calls.extend((lambda value=value: source.baseline("T", value),
                      lambda value=value: tangent.coefficient("T", value),
                      lambda value=value: tangent.varied_P("T", value)))
    calls.extend((lambda: tangent.coefficient("T", 0), lambda: tangent.varied_P("T", 0)))
    for value in (True, False, sp.true, 0, 1, "transverse", "longitudinal", "", None):
        calls.extend((lambda value=value: source.data(value),
                      lambda value=value: source.baseline(value, 1),
                      lambda value=value: tangent.coefficient(value, 1),
                      lambda value=value: tangent.reference(value),
                      lambda value=value: envelopes.reference(value),
                      lambda value=value: tail.reference_tail(value)))
    for value in (True, False, 0.0, sp.Float(0), sp.oo, sp.nan, source.n[11],
                  source.n[0]*source.n[1], sp.Integer(1), sp.Float(1)*source.n[0],
                  sp.Symbol("unrelated")*source.n[0]):
        calls.append(lambda value=value: source.linear_bound(value))
    for value in (True, False, sp.true, 1.0, sp.Float(1), "1", sp.I, sp.oo, -sp.oo, sp.zoo, sp.nan, 0, -1):
        calls.extend((lambda value=value: evolution.bound(value, 1000),
                      lambda value=value: evolution.bound(10**24, value)))
    for value in (1, 999):
        calls.append(lambda value=value: evolution.bound(10**24, value))
    calls.extend(lambda value=value: serialize(value) for value in (1.0, sp.Float(1), sp.oo, sp.nan))
    rejected = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(calls):
        raise ValueError("An unsupported prepared mass-response input was accepted")
    return {"rejected_inputs": rejected, "native_reference_order_validation_before_cache": True,
            "C10_norm_does_not_accept_uncontrolled_eleventh_derivatives": True,
            "nonzero_initial_mixing_retained_despite_zero_initial_response": True,
            "eighth_order_reference_not_new_subtraction_or_state": True,
            "common_dimensional_limit_before_separate_finite_assignments": True,
            "homogeneous_C10_to_C0_bound_not_no_loss_feedback": True}


@cache
def build_report():
    prior, identities = prior_checks(), proofs.residuals()
    exact, gates = affine.certify_residuals(identities), proofs.checks()
    if not all(value is True for value in gates.values()):
        raise ValueError("A finite prepared mass-response audit gate failed")
    scalar_count = sum(value.rows*value.cols if isinstance(value, sp.MatrixBase) else 1 for value in identities.values())
    coefficient_bounds = {sector: {j: source.linear_bound(tangent.coefficient(sector, j)) for j in range(1, 5)}
                          for sector in ("T", "L")}
    reference_bounds = {sector: {key: value for key, value in envelopes.reference(sector).items()
                                if key not in ("polynomial_majorants", "box_reconstructions", "base_reference_checks")}
                        for sector in ("T", "L")}
    return {"schema": 1, "claim": "P8-S6.66.FINITE_PREPARED_MASS_RESPONSE", "date": "2026-09-08",
            "status": "FINITE_HOMOGENEOUS_PREPARED_MASS_RESPONSE; ORIGINAL_P8_OPEN",
            "prior_sha256": prior, "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in source_files()},
            "formulation": "FORMULATION.md", "written_proofs": ["notes/reference.md", "notes/evolution.md", "notes/limit.md"],
            "exact_residuals": exact, "named_exact_check_count": len(identities), "checked_scalar_entries": scalar_count,
            "proof_checks": gates,
            "eighth_order_source_coefficient_bounds": serialize(coefficient_bounds),
            "varied_reference_residual_and_normalization_bounds": serialize(reference_bounds),
            "varied_readout_tail_bounds": serialize({sector: tail.reference_tail(sector) for sector in ("T", "L")}),
            "exact_mixing_and_reference_constants": serialize(evolution.constants()),
            "finite_response_scale_example": serialize(evolution.bound(10**24, 1000)),
            "initial_mixing_and_source_norm_controls": serialize({"squeezed_initial_phase_fixture": -sp.Rational(15, 4),
                "tenth_source_derivative_fixtures": {sector: tangent.reference(sector)["tenth_source_derivative_fixture"] for sector in ("T", "L")}}),
            "reference_mode_policy": "The original S6.55 all-order prepared Gaussian state is evolved exactly. W8 is only a comparison reference. Subtraction remains orders 0,2,4. Its source variation includes moving canonical normalization, physical readout weights, temporal contact and nonzero initial interference.",
            "dimensional_limit": "The written proof uses a sufficiently small complex-D neighborhood of 3 inside the frozen preparation domain, analytic plus/minus modes and D-independent fixed Borel cutoffs. Initial mixing remains uniformly O(nu^-6), the varied finite tail is O(nu^-4), and the radial measure is integrable for |Re D-3|<1/4 after shrinking the neighborhood as needed. Numerical bounds are on D=3, not uniform numerical claims throughout the old strip.",
            "prepared_source_domain": "Spatial K=0, smooth compact time support strictly after u0=-1/2, zero on an initial neighborhood. The source norm is the maximum over time derivatives 0..10 on I=[-1/2,1/2]. Outputs are normalized physical C0 readouts. No independent varied initial covariance is added.",
            "controls": controls(),
            "verdict": "The finite selected-state prepared homogeneous mass-response block has an integrable exact-minus-adiabatic remainder and the frozen finite local component. At L=10^24,m=1000 the nonlocal C10-to-C0 bound is below 10^-43 and the complete mass-block bound is below 10^-39. Full spatial/metric response, coupled feedback and original P8 remain open.",
            "not_established": ["Arbitrary spatial-momentum finite response, full metric or second mass-source blocks, or arbitrary varied initial states",
                                "A no-loss coupled feedback norm, quantum stability/cones, interactions or cutoff",
                                "Finite Wilson matching, a common UV parent, V/G/B or original P8 closure"],
            "verification_boundary": "Exact varied Riccati, tail, Wronskian, mixing and radial identities with continuous rational majorants, plus written parameter-dependence and common-dimensional dominated-convergence arguments; not proof-assistant formalization"}


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The finite prepared mass-response report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.66.FINITE_PREPARED_MASS_RESPONSE replay passed; coupled response and original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
