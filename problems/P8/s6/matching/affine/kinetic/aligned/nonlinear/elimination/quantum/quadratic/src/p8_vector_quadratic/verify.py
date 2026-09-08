"""Read-only curved quadratic vector mass-insertion UV-pole certificate."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as affine
from p8_clock_tadpole import verify as parent

from . import bimetric, checks, geometry, invariants, kernel, response, tensors

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"curved-quadratic-vector-pole.json"
PARENT_SHA = "dba6c1f44957adaa54e15557bbf14eac6eeee5f911524caea301315d91425241"
sha, serialize = affine.sha, affine.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_vector_quadratic/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen selected-state clock-profile certificate changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {"S6_60_fully_rebuilt": PARENT_SHA, "unchanged_selected_action_state_and_clock_profiles": True}


@cache
def residuals():
    out = {**tensors.checks(), **geometry.checks(), **checks.flat(),
           **checks.scalar_conformal()["compact_variation_checks"], **checks.general_second(),
           **checks.constant_anisotropic(), **response.actual_mass_checks()}
    for order in (0, 2, 4):
        out.update(response.checks(order))
        value = kernel.pole(order)
        out["real_Fourier_bilinear_order_"+str(order)] = sp.expand(sp.im(value))
        out["opposite_leg_symmetry_order_"+str(order)] = sp.expand(value-tensors.swap_legs(value))
    return out


@cache
def proof_checks():
    out = {"withheld_reference_has_nonzero_exact_fixture": checks.withheld_control()["nonzero_exact_fixture"] == sp.Rational(237, 20),
           "accepted_fourth_order_not_literal_control": kernel.pole(4) != invariants.evaluate(4)/(960*tensors.mass**4),
           "reference_basis_has_expected_term_counts": (len(invariants.ZERO), len(invariants.SECOND), len(invariants.FOURTH)) == (2, 9, 48)}
    mass_fields = tensors.plus.temporal+tensors.plus.spatial+tensors.minus.temporal+tensors.minus.spatial
    for order in (0, 2, 4):
        value = sp.Poly(kernel.pole(order), *mass_fields)
        out["homogeneous_quadratic_mass_order_"+str(order)] = all(sum(powers) == 2 for powers, _ in value.terms())
        zero = {field: 0 for field in mass_fields}
        out["zero_clock_value_order_"+str(order)] = value.as_expr().subs(zero) == 0
        out["zero_first_mass_variation_order_"+str(order)] = all(sp.diff(value.as_expr(), field).subs(zero) == 0 for field in mass_fields)
        out["exact_rational_not_float_order_"+str(order)] = not kernel.pole(order).atoms(sp.Float)
    out["geometric_density_has_two_jet_locality"] = not (bimetric.local_quadratic_density().free_symbols
                                                        -set(bimetric.temporal.values())-set(bimetric.spatial.values())-set(tensors.H))
    return {name: bool(value) for name, value in out.items()}


def controls():
    bad_orders = (True, False, sp.true, sp.false, 1.0, sp.Float(0), sp.Integer(0),
                  "0", -1, 1, 6, sp.Rational(1, 2))
    calls = []
    for value in bad_orders:
        for function in (kernel.pole, invariants.evaluate, response.raw, response.local_form, response.closed_form, response.checks):
            calls.append(lambda value=value, function=function: function(value))
    bad_indices = (True, False, sp.true, sp.false, 0.0, sp.Float(0), sp.Integer(0), "0", -1, 4, sp.Rational(1, 2))
    for value in bad_indices:
        calls.extend((lambda value=value: geometry.metric(value),
                      lambda value=value: geometry.connection(value, 0, 0),
                      lambda value=value: geometry.curvature(0, value, 0, 0),
                      lambda value=value: geometry.ricci(0, value),
                      lambda value=value: geometry.derivative(sp.Integer(0), value),
                      lambda value=value: tensors.plus.covariant((value,), 0, 0),
                      lambda value=value: tensors.plus.covariant((), value, 0)))
    calls.extend(lambda name=name, sign=sign: tensors.ModeTensor(name, sign) for name, sign in
                 (("plus", True), ("minus", -1.0), ("plus", sp.Integer(1)), ("plus", -1),
                  ("minus", 1), ("unknown", 1)))
    calls.extend(lambda sequence=sequence: tensors.plus.covariant(sequence, 0, 0)
                 for sequence in ([0], "0", (), (0,)*5) if sequence != ())
    bad_contractions = ((), (("Y", "", "mn"),), (("Y", "", "mn"), ("Y", "", "mm")),
                        (("unknown", "", ""), ("Y", "", "mm"), ("Y", "", "nn")),
                        (("Y", "", "m"), ("Y", "", "m")), (("Y", "aaaaaa", "mn"), ("Y", "", "mn")),
                        (("Y", "", "00"), ("Y", "", "00")), (("Y", "", "mn", "extra"), ("Y", "", "mn")))
    calls.extend(lambda value=value: tensors.contraction(value) for value in bad_contractions)
    calls.extend(lambda value=value: serialize(value) for value in (1.0, sp.Float(1), sp.oo, sp.nan))
    rejected = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(calls):
        raise ValueError("An unsupported quadratic-pole input was accepted")
    return {"rejected_inputs": rejected, "nested_index_validation_precedes_cache_lookup": True,
            "literal_fourth_reference_withheld_not_coefficient_fitted": True,
            "compact_variation_boundaries_not_global_flux_deletion": True,
            "four_dimensional_residue_not_finite_dimensional_matching": True,
            "flat_potential_quadratic_term_not_counted_twice": True,
            "pole_sign_not_quantum_stability_or_cutoff_verdict": True}


@cache
def build_report():
    prior, identities = prior_checks(), residuals()
    exact, proofs = affine.certify_residuals(identities), proof_checks()
    if not all(value is True for value in proofs.values()):
        raise ValueError("A curved quadratic pole audit check failed")
    scalar_count = sum(value.rows*value.cols if isinstance(value, sp.MatrixBase) else 1 for value in identities.values())
    return {"schema": 1, "claim": "P8-S6.61.VECTOR_QUADRATIC_POLE", "date": "2026-09-08",
            "status": "CURVED_QUADRATIC_RETAINED_VECTOR_UV_RESIDUE; ORIGINAL_P8_OPEN",
            "prior_sha256": prior, "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in source_files()},
            "formulation": "FORMULATION.md", "written_proofs": ["notes/pole.md", "notes/lapse.md", "notes/literature.md"],
            "exact_residuals": exact, "named_exact_check_count": len(identities), "checked_scalar_entries": scalar_count,
            "proof_checks": proofs,
            "conventions": "Euclidean Gamma_div^(Y^2)=integral sqrt(g) Q/(32*pi^2*epsilon), epsilon=(4-d)/2. Opposite spatial Fourier legs, ordered covariant derivatives, and compact-support action comparisons. Lorentzian matching needs the frozen Wick/sign conversion.",
            "accepted_pole_bilinears": serialize({order: kernel.pole(order) for order in (0, 2, 4)}),
            "independent_scalar_conformal_check": serialize(checks.scalar_conformal()),
            "auxiliary_metric_definition": "A=1+Y/m^2, g_tilde(lambda)=sqrt(det(1+lambda*Y/m^2)) (1+lambda*Y/m^2)^(-1) g. Q4=-[lambda^2]sqrt(g_tilde)/sqrt(g) a4_scalar(g_tilde); a4_scalar=R^2/72-Ricci^2/180+Riemann^2/180, modulo the covariant total Laplacian only. No four-dimensional Gauss-Bonnet reduction defines this basis.",
            "independent_coordinate_geometry": serialize(geometry.heat_kernel()),
            "local_auxiliary_metric_quadratic_density": serialize(bimetric.local_quadratic_density()),
            "withheld_literal_reference": serialize(checks.withheld_control()),
            "actual_mass_profiles": serialize(dict(zip(("alpha", "beta"), response.mass_coefficients()))),
            "actual_clock_lapse_local_coefficients": serialize({order: response.local_form(order) for order in (0, 2, 4)}),
            "matching_boundary": "The full frozen flat potential already contains Q0, which is checked and not added again. The new derivative/curvature Y^2 data have zero value and first variation on the clock. They do not alter S6.60's selected first-variation cancellation. A dimensionally continued quadratic counterterm and its finite evanescent contribution remain to be specified and computed.",
            "controls": controls(),
            "verdict": "Independent geometric derivation and exact flat/curved checks supply the four-dimensional quadratic retained-vector mass-insertion UV residue and its actual clock-lapse local action. Original P8 remains open.",
            "not_established": ["Finite selected-state contact/retarded kernel, compatible varied-state covariance or integrated feedback norm",
                                "Full mixed metric/mass second variation, coupled quantum cones, interactions, cutoff or omitted-sector/higher-loop control",
                                "Common UV parent, finite-gravity remainder, V/G/B or original P8 closure"],
            "verification_boundary": "Exact symbolic/rational coordinate and frame tensor calculus, independent Feynman and conformal comparisons, written determinant and compact-support arguments; not proof-assistant formalization"}


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The curved quadratic vector pole report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.61.VECTOR_QUADRATIC_POLE replay passed; finite response and original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
