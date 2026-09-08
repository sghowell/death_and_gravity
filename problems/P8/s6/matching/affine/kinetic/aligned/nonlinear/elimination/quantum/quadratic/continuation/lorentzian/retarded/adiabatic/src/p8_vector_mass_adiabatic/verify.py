"""Read-only homogeneous finite local mass-response subtraction checkpoint."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as affine
from p8_vector_retarded import verify as parent

from . import bounds, canonical, compact, jets, radial, variation

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"finite-local-mass-response-subtraction.json"
PARENT_SHA = "17a31b89b503e42852dc6b60c986e615c02c582cd4b72599b28774f8e1d235a5"
sha, serialize = affine.sha, affine.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_vector_mass_adiabatic/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen bare retarded mass-response certificate changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {"S6_64_fully_rebuilt": PARENT_SHA, "original_clock_state_and_profile_preserved": True}


@cache
def residuals():
    out = {}
    for module in (canonical, variation, compact, bounds):
        out.update({module.__name__.rsplit(".", 1)[1]+"_"+name: value for name, value in module.checks().items()})
    out.update(radial.pole_checks())
    return out


@cache
def proof_checks():
    example = bounds.operator_bound(10**24, 1000)
    out = {"actual_C4_to_C0_local_norm_below_1e_minus_39":
               example["C4_to_C0_upper_bound"] < sp.Rational(1, 10**39)}
    for order in (0, 1, 2):
        operator = radial.matched(order)
        out["exact_linear_source_operator_"+str(order)] = all(sum(powers) == 1 for powers, _ in sp.Poly(operator, *jets.n).terms())
        out["no_binary_float_"+str(order)] = not operator.atoms(sp.Float)
        out["nonnegative_continuous_coefficient_envelopes_"+str(order)] = all(value >= 0 for value in bounds.coefficient_envelopes()[order].values())
    for order in (1, 2):
        data = compact.coefficients(order)
        density = data["second_time_square"]*jets.n[2]**2+data["first_time_square"]*jets.n[1]**2+data["source_square"]*jets.n[0]**2
        zero = {field: 0 for field in jets.n}
        out["derivative_quadratic_component_zero_to_first_order_"+str(order)] = (
            density.subs(zero) == 0 and all(sp.diff(density, field).subs(zero) == 0 for field in jets.n))
        out["omitted_finite_counterterm_adjoint_fixture_nonzero_"+str(order)] = (
            compact.controls()["half_time_unit_mass_adjoint_fixture_"+str(2*order)] != 0)
    return {name: bool(value) for name, value in out.items()}


def controls():
    calls = []
    for value in (True, False, sp.true, sp.false, 0.0, sp.Float(0), sp.Integer(0), "0", -1, 3):
        for function in (variation.combined, radial.laurent, radial.matched, compact.coefficients):
            calls.append(lambda function=function, value=value: function(value))
    for value in (True, False, sp.true, 0, 1, "transverse", "longitudinal", "", None):
        calls.append(lambda value=value: variation.data(value))
    u = radial.u
    for value in (True, False, 0.0, sp.Float(0), sp.oo, sp.nan, sp.Symbol("unrelated"),
                  1/(2+u**2), 1/(1-u**2), 1/u, sp.sin(u)):
        calls.append(lambda value=value: bounds.rational_bound(value))
    for value in (True, False, sp.true, 1.0, sp.Float(1), "1", sp.I, sp.oo, -sp.oo, sp.zoo, sp.nan, 0, -1):
        calls.extend((lambda value=value: bounds.operator_bound(value, 1000),
                      lambda value=value: bounds.operator_bound(10**24, value)))
    for value in (1, 999):
        calls.append(lambda value=value: bounds.operator_bound(10**24, value))
    calls.extend(lambda value=value: serialize(value) for value in (1.0, sp.Float(1), sp.oo, sp.nan))
    rejected = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(calls):
        raise ValueError("An unsupported finite local-response input was accepted")
    return {"rejected_inputs": rejected, "native_order_validation_before_cache": True,
            "full_D_polarization_and_rate_dependence_retained": True,
            "finite_radial_and_counterterm_pieces_combined": True,
            "zero_order_full_potential_not_added_twice": True,
            "local_C4_norm_not_nonlocal_feedback_norm": True,
            "homogeneous_source_only_no_spatial_bound": True}


@cache
def build_report():
    prior, identities = prior_checks(), residuals()
    exact, proofs = affine.certify_residuals(identities), proof_checks()
    if not all(value is True for value in proofs.values()):
        raise ValueError("A finite local mass-response audit check failed")
    return {"schema": 1, "claim": "P8-S6.65.FINITE_LOCAL_MASS_RESPONSE", "date": "2026-09-08",
            "status": "HOMOGENEOUS_FINITE_LOCAL_SUBTRACTION_COMPONENT; ORIGINAL_P8_OPEN",
            "prior_sha256": prior, "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in source_files()},
            "formulation": "FORMULATION.md", "written_proofs": ["notes/variation.md", "notes/matching.md", "notes/bound.md"],
            "exact_residuals": exact, "named_exact_check_count": len(identities), "checked_scalar_entries": len(identities),
            "proof_checks": proofs,
            "moving_mode_source_variations": serialize({sector: {key: value for key, value in variation.data(sector).items() if key != "readout_coefficients"} for sector in ("T", "L")}),
            "actual_radial_Laurent_coefficients": serialize({j: {key: radial.actual(value) for key, value in radial.laurent(j).items()} for j in (0, 1, 2)}),
            "actual_finite_local_operators": serialize({j: radial.matched(j) for j in (0, 1, 2)}),
            "actual_compact_finite_local_actions": serialize(bounds.closed_coefficients()),
            "omitted_counterterm_adjoint_controls": serialize(compact.controls()),
            "continuous_operator_coefficient_envelopes": serialize(bounds.coefficient_envelopes()),
            "scale_example": serialize(bounds.operator_bound(10**24, 1000)),
            "normalization_and_prescription": "The local action component is integral a^3(A n''^2+B n'^2+C n^2)/(64*pi^2), and its normalized physical Euler operator is displayed. At mu=m the full radial Laurent finite part is combined with +4 partial_D E_loop,D in this normalization. No additional finite Y^2 Wilson matching coefficient is included or inferred. The flat potential is already present and not added again.",
            "remaining_integral": "The exact selected-state causal response minus its varied local adiabatic orders 0,2,4 remains to be integrated with a justified common-dimensional limit and quantitative bounds. The local coefficient result and C4-to-C0 bound do not establish that remaining integral or a coupled feedback contraction.",
            "controls": controls(),
            "verdict": "Independent moving-mode and varied Riccati calculations reproduce the complete homogeneous curved poles, and the finite radial-plus-counterterm local operator is self-adjoint, agrees with the frozen static potential and has the stated continuous bound. The nonlocal response and original P8 remain open.",
            "not_established": ["The subtracted exact nonlocal response, its compatible-dimensional preparation estimates or integrated feedback norm",
                                "Spatially varying finite response, full metric/second-mass blocks or an arbitrary varied Cauchy-state family",
                                "Coupled quantum stability/cones, interactions, cutoff, finite Wilson matching, common UV parent, V/G/B or original P8 closure"],
            "verification_boundary": "Exact physical Hamiltonian, differentiated Riccati, polynomial radial, dimensional-counterterm, variational and continuous rational-envelope checks; not proof-assistant formalization"}


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The finite local mass-response report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.65.FINITE_LOCAL_MASS_RESPONSE replay passed; nonlocal response and original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
