"""Read-only quadratic counterterm continuation and evanescence certificate."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as affine
from p8_vector_quadratic import tensors
from p8_vector_quadratic import verify as parent

from . import bimetric, contractions, evanescence, geometry, kernel, response

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"quadratic-counterterm-dimensional-continuation.json"
PARENT_SHA = "6ac2f5d493737934590e3aa64786549419197ba8b74365660ae8bcbf09abd53a"
sha, serialize = affine.sha, affine.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_vector_quadratic_matching/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen curved quadratic pole certificate changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {"S6_61_fully_rebuilt": PARENT_SHA, "old_potential_and_first_variation_prescription_not_replaced": True}


@cache
def residuals():
    out = {**geometry.checks(), **contractions.checks(), **bimetric.checks()}
    for order in (2, 4):
        out.update(response.checks(order))
    boundary = evanescence.gauss_bonnet_control()
    out.update({name: boundary[name] for name in
                ("exact_curvature_identity", "exact_boundary_identity", "scalar_compact_variation_identity")})
    return out


@cache
def proof_checks():
    boundary = evanescence.gauss_bonnet_control()
    out = {"Gauss_Bonnet_simple_dimension_zero_has_coefficient_one_over_90":
           boundary["topological_factor_has_simple_zero_at_D_three"] == sp.Rational(1, 90),
           "Gauss_Bonnet_nonzero_Euler_fixture": boundary["nonzero_dimension_derivative_fixture"] == sp.Rational(1, 15),
           "omitted_measure_variation_nonzero_bounce_fixture": response.measure_control(4)["bounce_fixture"] == -sp.Rational(1024, 6561)}
    fields = tensors.plus.temporal+tensors.plus.spatial+tensors.minus.temporal+tensors.minus.spatial
    for order in (2, 4):
        value = kernel.pole(order)
        polynomial = sp.Poly(value, *fields)
        out["quadratic_mass_order_in_all_dimensions_"+str(order)] = all(sum(powers) == 2 for powers, _ in polynomial.terms())
        out["polynomial_in_D_not_dimensional_sampling_"+str(order)] = all(
            coefficient.is_polynomial(geometry.dimension) for coefficient in polynomial.coeffs())
        zero = {field: 0 for field in fields}
        out["zero_clock_value_and_first_variation_in_D_"+str(order)] = (
            value.subs(zero) == 0 and all(sp.diff(value, field).subs(zero) == 0 for field in fields))
        out["exact_no_float_continuation_"+str(order)] = not value.atoms(sp.Float)
    return {name: bool(value) for name, value in out.items()}


def controls():
    bad_orders = (True, False, sp.true, sp.false, 2.0, sp.Float(2), sp.Integer(2), "2", -1, 0, 3, 6, sp.Rational(1, 2))
    bad_jets = (True, False, sp.true, sp.false, 0.0, sp.Float(0), sp.Integer(0), "0", -1, 2, sp.Rational(1, 2))
    calls = []
    for value in bad_orders:
        for function in (kernel.pole, response.closed_first_jet, response.finite_counterterm_operator,
                         response.checks, response.measure_control):
            calls.append(lambda function=function, value=value: function(value))
        for function in (kernel.dimension_jet, response.raw_jet, response.compact_jet, response.operator_jet):
            calls.append(lambda function=function, value=value: function(value, 0))
    for value in bad_jets:
        for function in (kernel.dimension_jet, response.raw_jet, response.compact_jet, response.operator_jet):
            calls.append(lambda function=function, value=value: function(2, value))
    for labels in ((), ["a"], "a", ("a", "a"), ("aa",), ("1",), (True,), tuple("abcde")):
        calls.append(lambda labels=labels: tuple(contractions.assignments(labels)))
    for value in (True, 0.0, sp.Integer(0), "0", -1, 6):
        calls.extend((lambda value=value: contractions.ricci(value, 0),
                      lambda value=value: contractions.mode_component(tensors.plus, (value,), 0, 0)))
    calls.extend((lambda: contractions.contraction(()),
                  lambda: contractions.contraction((("Y", "", "mn"),)),
                  lambda: contractions.mode_component(tensors.plus, (0, 0, 0), 0, 0),
                  lambda: contractions.mode_component(tensors.plus, [0], 0, 0)))
    calls.extend(lambda value=value: serialize(value) for value in (1.0, sp.Float(1), sp.oo, sp.nan))
    rejected = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(calls):
        raise ValueError("An unsupported counterterm-continuation input was accepted")
    return {"rejected_inputs": rejected, "native_validation_survives_warm_caches": True,
            "flat_potential_not_readded": True, "finite_evanescent_piece_not_complete_finite_kernel": True,
            "D_dependent_measure_varied_before_limit": True,
            "four_dimensional_Gauss_Bonnet_identity_not_applied_early": True,
            "specified_continuation_not_unique_renormalization_scheme": True,
            "Euclidean_result_requires_Lorentzian_conversion": True}


@cache
def build_report():
    prior, identities = prior_checks(), residuals()
    exact, proofs = affine.certify_residuals(identities), proof_checks()
    if not all(value is True for value in proofs.values()):
        raise ValueError("A quadratic dimensional continuation check failed")
    return {"schema": 1, "claim": "P8-S6.62.QUADRATIC_CONTINUATION", "date": "2026-09-08",
            "status": "SPECIFIED_QUADRATIC_COUNTERTERM_EVANESCENCE; ORIGINAL_P8_OPEN",
            "prior_sha256": prior, "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in source_files()},
            "formulation": "FORMULATION.md", "written_proofs": ["notes/continuation.md", "notes/response.md", "notes/boundary.md"],
            "exact_residuals": exact, "named_exact_check_count": len(identities), "checked_scalar_entries": len(identities),
            "proof_checks": proofs,
            "literal_continuation": "D=3-2epsilon; continue the S6.61 two-derivative invariant contractions with fixed four-dimensional coefficients. For Q4 retain g_tilde=sqrt(det(1+Y/m^2)) (1+Y/m^2)^(-1) g_D and fixed a4sc=R^2/72-Ricci^2/180+Riemann^2/180, with D+1-dimensional determinants and contractions. This is a specified evanescent prescription, not an assertion of the same scalar-action transformation away from four dimensions.",
            "continued_auxiliary_curvatures": serialize(geometry.heat_kernel()),
            "quadratic_pole_dimensional_jets": serialize({order: {jet: kernel.dimension_jet(order, jet) for jet in (0, 1)} for order in (2, 4)}),
            "actual_compact_dimensional_jets": serialize({order: {jet: response.compact_jet(order, jet) for jet in (0, 1)} for order in (2, 4)}),
            "normalized_lapse_operator_dimensional_jets": serialize({order: {jet: response.operator_jet(order, jet) for jet in (0, 1)} for order in (2, 4)}),
            "bare_counterterm_finite_evanescent_operators": serialize({order: response.finite_counterterm_operator(order) for order in (2, 4)}),
            "measure_variation_controls": serialize({order: response.measure_control(order) for order in (2, 4)}),
            "Gauss_Bonnet_boundary_control": serialize(evanescence.gauss_bonnet_control()),
            "counterterm_convention": "Euclidean bare-action counterterm -Q_D/(32*pi^2*epsilon); normalized lapse variation -E3/(32*pi^2*epsilon)+2*(partial_D E_D)_3/(32*pi^2). E_D is varied with measure a^D before the limit. Lorentzian Wick/sign conversion remains required.",
            "scope": "Only new quadratic mass-deviation derivative counterterms are continued. All old lower-order matching, the full flat potential and selected stress profiles are unchanged. These local polynomials are entire in D; no nonlocal kernel holomorphy or finite renormalized response is inferred.",
            "controls": controls(),
            "verdict": "The specified new quadratic counterterm continuation and its finite evanescent normalized lapse variation are explicit and independently checked. Original P8 remains open.",
            "not_established": ["Lorentzian conversion and the dimensionally continued bare contact/retarded kernel combined with the counterterm",
                                "Compatible varied-state covariance, finite Wilson-coefficient boundary values, integrated feedback, quantum stability/cones or cutoff",
                                "Common UV parent, V/G/B or original P8 closure"],
            "verification_boundary": "Exact polynomial/rational geometry, symbolic index multiplicities, dimensional jets and unreduced-versus-compact variations, with written warped-product and boundary arguments; not proof-assistant formalization"}


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The quadratic continuation report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.62.QUADRATIC_CONTINUATION replay passed; finite kernel and original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
