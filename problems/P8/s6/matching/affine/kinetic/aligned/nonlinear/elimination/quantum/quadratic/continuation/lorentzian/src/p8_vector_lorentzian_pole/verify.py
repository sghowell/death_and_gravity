"""Read-only physical-signature quadratic pole and evanescent counterterm."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as affine
from p8_vector_quadratic import tensors
from p8_vector_quadratic_matching import geometry
from p8_vector_quadratic_matching import verify as parent

from . import bimetric, checks, frame, response, wick

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"physical-signature-quadratic-vector-pole.json"
PARENT_SHA = "c85e3f36ce0d8f2d98fe8169af17d478a30e57aaef5567e261a980117ed8cc6d"
sha, serialize = affine.sha, affine.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_vector_lorentzian_pole/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen quadratic dimensional-continuation certificate changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {"S6_62_fully_rebuilt": PARENT_SHA, "original_Lorentzian_clock_and_selected_state_unchanged": True}


@cache
def residuals():
    out = {**checks.frame_checks(), **checks.signature_checks(), **checks.physical_geometry(), **checks.flat()}
    for order in (2, 4):
        out.update(response.checks(order))
        for jet in (0, 1):
            out["real_Lorentzian_pole_"+str(order)+"_"+str(jet)] = sp.expand(sp.im(wick.loop_pole(order, jet)))
    return out


@cache
def proof_checks():
    history = checks.history_control()
    out = {"naive_time_substitution_changes_the_chosen_history": history["nonzero_history_difference"] == -4*response.u**2,
           "original_and_changed_half_time_values_distinct": (history["original_half_time_scale_factor"], history["changed_half_time_scale_factor"]) == (sp.Rational(25, 16), sp.Rational(9, 16)),
           "physical_fourth_order_bounce_coefficient_positive": history["fourth_pole_second_time_square_bounce"] == sp.Rational(128, 6561)}
    fields = tensors.plus.temporal+tensors.plus.spatial+tensors.minus.temporal+tensors.minus.spatial
    for order, value in ((2, frame.second_loop_pole()), (4, bimetric.loop_pole())):
        polynomial = sp.Poly(value, *fields)
        out["quadratic_mass_degree_"+str(order)] = all(sum(powers) == 2 for powers, _ in polynomial.terms())
        out["polynomial_in_symbolic_dimension_"+str(order)] = all(coefficient.is_polynomial(geometry.dimension) for coefficient in polynomial.coeffs())
        out["no_float_in_physical_signature_calculation_"+str(order)] = not value.atoms(sp.Float)
        zero = {field: 0 for field in fields}
        out["new_operators_vanish_to_first_order_on_clock_"+str(order)] = (
            value.subs(zero) == 0 and all(sp.diff(value, field).subs(zero) == 0 for field in fields))
    return {name: bool(value) for name, value in out.items()}


def controls():
    bad_orders = (True, False, sp.true, sp.false, 2.0, sp.Float(2), sp.Integer(2), "2", -1, 0, 1, 6, sp.Rational(1, 2))
    bad_jets = (True, False, sp.true, sp.false, 0.0, sp.Float(0), sp.Integer(0), "0", -1, 2, sp.Rational(1, 2))
    calls = []
    for value in bad_orders:
        for function in (wick.loop_pole, response.checks, response.finite_counterterm):
            calls.append(lambda function=function, value=value: function(value))
        for function in (response.raw, response.compact, response.operator):
            calls.append(lambda function=function, value=value: function(value, 0))
    for value in bad_jets:
        for function in (wick.loop_pole, response.raw, response.compact, response.operator):
            calls.append(lambda function=function, value=value: function(2, value))
    for value in (True, False, sp.true, sp.false, 0.0, sp.Float(0), sp.Integer(0), "0", -1, 6):
        calls.extend((lambda value=value: frame.eta(value),
                      lambda value=value: frame.ricci(value, 0),
                      lambda value=value: frame.riemann(0, value, 0, 0),
                      lambda value=value: frame.covariant(tensors.plus, (value,), 0, 0),
                      lambda value=value: frame.covariant(tensors.plus, (), value, 0)))
    for value in (True, False, sp.true, sp.false, 1.0, sp.Float(1), "1", sp.oo, -sp.oo, sp.zoo, sp.nan):
        calls.extend((lambda value=value: wick.continuation(value),
                      lambda value=value: wick.graded_actual(value, 2)))
    calls.extend((lambda: frame.contraction(()), lambda: frame.contraction((("Y", "", "mn"),)),
                  lambda: frame.covariant(tensors.plus, [0], 0, 0),
                  lambda: frame.covariant(tensors.plus, "0", 0, 0),
                  lambda: frame.covariant(tensors.plus, (0, 0, 0), 0, 0),
                  lambda: frame.covariant("unknown", (), 0, 0)))
    calls.extend(lambda value=value: serialize(value) for value in (1.0, sp.Float(1), sp.oo, sp.nan))
    rejected = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(calls):
        raise ValueError("An unsupported physical-signature input was accepted")
    return {"rejected_inputs": rejected, "nested_index_and_native_jet_validation_before_cache": True,
            "local_covariant_continuation_not_global_state_Wick_rotation": True,
            "bare_action_counterterm_opposite_to_loop_pole": True,
            "full_flat_potential_checked_not_readded": True,
            "finite_evanescent_piece_not_complete_finite_response": True,
            "local_higher_derivative_poles_not_exact_new_EFT_modes": True}


@cache
def build_report():
    prior, identities = prior_checks(), residuals()
    exact, proofs = affine.certify_residuals(identities), proof_checks()
    if not all(value is True for value in proofs.values()):
        raise ValueError("A physical-signature pole check failed")
    scalar_count = sum(value.rows*value.cols if isinstance(value, sp.MatrixBase) else 1 for value in identities.values())
    return {"schema": 1, "claim": "P8-S6.63.LORENTZIAN_QUADRATIC_POLE", "date": "2026-09-08",
            "status": "PHYSICAL_SIGNATURE_LOCAL_QUADRATIC_POLE_AND_EVANESCENCE; ORIGINAL_P8_OPEN",
            "prior_sha256": prior, "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in source_files()},
            "formulation": "FORMULATION.md", "written_proofs": ["notes/signature.md", "notes/clock.md", "notes/scope.md"],
            "exact_residuals": exact, "named_exact_check_count": len(identities), "checked_scalar_entries": scalar_count,
            "proof_checks": proofs,
            "local_signature_convention": "Physical eta=(-,+,+,+); Y_00=-Yt,Y_ij=Ys delta_ij for the mixed eigenvalue amplitudes. H_E^(j)=(-i)^(j+1)H_L^(j), Y_E^(j)=(-i)^j Y_L^(j), k_E=k_L, Q_loop,L=-Wick(Q_loop,E). Independent frame and negative-lapse coordinate calculations agree in arbitrary symbolic D.",
            "physical_loop_pole_dimensional_jets": serialize({order: {jet: wick.loop_pole(order, jet) for jet in (0, 1)} for order in (2, 4)}),
            "physical_clock_compact_dimensional_jets": serialize({order: {jet: response.compact(order, jet) for jet in (0, 1)} for order in (2, 4)}),
            "physical_normalized_Euler_dimensional_jets": serialize({order: {jet: response.operator(order, jet) for jet in (0, 1)} for order in (2, 4)}),
            "physical_bare_counterterm_finite_evanescent_operators": serialize({order: response.finite_counterterm(order) for order in (2, 4)}),
            "unchanged_history_control": serialize(checks.history_control()),
            "counterterm_sign_and_measure": "The Lorentzian bare-action counterterm is -Q_L,D/(32*pi^2*epsilon), not the same-sign subtraction functional. With epsilon=(3-D)/2, its normalized lapse variation is -E_L,3/(32*pi^2*epsilon)+2*(partial_D E_L,D)_3/(32*pi^2). All a^D measure derivatives remain included before the dimensional limit.",
            "grading_identity": "Every actual compact coefficient and normalized Euler operator obeys X_L,n,j(q)=-(-1)^(n/2)X_E,n,j(-q), n=2,4 and j=0,1. This is checked against independent physical-signature projection with full time-dependent mass-profile chain rules.",
            "preserved_scope": "The original physical clock, S6.55 selected state, S6.60 fixed profiles and all lower-order matching are unchanged. The new Y^2 derivative operators vanish to first order on the clock. No global complex-time history or Euclidean state preparation is invoked.",
            "controls": controls(),
            "verdict": "The physical-signature local quadratic mass-insertion pole and specified evanescent bare-counterterm variation are independently checked. The finite causal kernel and original P8 remain open.",
            "not_established": ["Finite dimensionally regulated bare/contact/retarded response, compatible varied-state covariance or integrated feedback bounds",
                                "Finite matching/Wilson-coefficient data, coupled quantum stability/cones, interactions or cutoff",
                                "Common UV parent, finite-gravity remainder, V/G/B or original P8 closure"],
            "verification_boundary": "Exact tensor, coordinate-curvature, Feynman, grading and normalized-variation comparisons, with written local signature and history arguments; not proof-assistant formalization"}


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The physical-signature quadratic pole report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.63.LORENTZIAN_QUADRATIC_POLE replay passed; finite kernel and original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
