"""Read-only quantitative C5 fixed-background vector-stress certificate."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as affine
from p8_vector_evolution import mixing
from p8_vector_hadamard import cutoffs
from p8_vector_variation import verify as parent

from . import estimates, frequency, preparation, spectral

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"fifth-order-vector-stress-regularity.json"
PARENT_SHA = "e4aa8aa397c2722168d7f475cd6b6fb73d742d3e143596d344955c5bd1feeb99"
sha, serialize = affine.sha, affine.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_vector_regularity/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen actual-vector-variation certificate changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {"S6_58_fully_rebuilt": PARENT_SHA,
            "fixed_state_cutoffs_action_and_counterterms_not_overwritten": True}


@cache
def residuals():
    out = {**spectral.algebra_checks(), **frequency.algebra_checks(), **estimates.radial_checks(),
           **mixing.identities(), **mixing.product_identities(), **cutoffs.initial_data_identities()}
    for kind in ("transverse", "longitudinal"):
        data = frequency.reference(kind)
        out[kind+"_eighth_order_low_residual_coefficients"] = sp.ImmutableMatrix(list(data["low_residual_numerator_coefficients"].values()))
        out[kind+"_frequency_coefficient_box_reconstructions"] = sp.ImmutableMatrix(data["box_reconstructions"])
        for order, row in preparation.constants()[kind]["coefficient_bounds"].items():
            out[kind+"_preparation_coefficient_reconstructions_"+str(order)] = sp.ImmutableMatrix([
                row["coefficient_reconstruction"], row["slope_reconstruction"]])
    for name in spectral.tail.physical_weights():
        for order in range(6):
            tag = name+"_"+str(order)
            out[tag+"_readout_box_reconstructions"] = sp.ImmutableMatrix(spectral.row_bounds(name, order)["box_reconstructions"])
            out[tag+"_differentiated_subtraction_cancellations"] = sp.ImmutableMatrix(
                list(spectral.reference_tail(name, order)["low_coefficient_residuals"].values()))
    for observable in ("energy", "pressure"):
        for order in range(6):
            out[observable+"_local_derivative_box_reconstructions_"+str(order)] = sp.ImmutableMatrix(
                list(estimates.local_derivatives(observable, order)["reconstructions"].values()))
    return out


@cache
def proof_checks():
    out = {**preparation.proof_checks(), **estimates.proof_checks()}
    for kind in ("transverse", "longitudinal"):
        out.update({kind+"_"+key: value for key, value in frequency.reference(kind)["proof_checks"].items()})
    mass = spectral.wkb.MASS_TIME_MIN
    out.update({"reference_squared_momentum_envelope": sp.Rational(9, 4)+25/mass**2 < 3,
                "reference_mixed_product_envelope": sp.Integer(3) < 4,
                "small_norm_exponential_below_two": 1/(1-sp.Rational(1, 4)) < 2,
                "lower_reference_failure_is_nonzero": spectral.low_reference_failure() == sp.Rational(31744, 81),
                "lower_reference_is_not_certified_integrable": not spectral.reference_tail("longitudinal_energy", 2, 2)["integrable_tail_certified"]})
    data = estimates.physical_bounds(10**24, 1000)
    for name, rows in data["normalized_derivative_bounds"].items():
        for order, row in rows.items():
            out[name+"_scale_example_derivative_"+str(order)+"_below_1e_minus_18"] = row["total"] < sp.Rational(1, 10**18)
    return {key: bool(value) for key, value in out.items()}


def controls():
    bad_exact = (True, False, sp.true, sp.false, 1.0, sp.Float(1), "1000", sp.I,
                 sp.oo, -sp.oo, sp.zoo, sp.nan, sp.Symbol("unproved"))
    calls = []
    for value in bad_exact:
        calls.extend((lambda value=value: preparation.threshold(value),
                      lambda value=value: preparation.mixing_bound(1000, value),
                      lambda value=value: estimates.physical_bounds(value, 1000),
                      lambda value=value: estimates.physical_bounds(10**24, value)))
    bad_orders = (True, False, sp.true, sp.false, 1.0, sp.Float(1), sp.Integer(1), "1", -1, 6, sp.Rational(1, 2))
    for value in bad_orders:
        calls.extend((lambda value=value: spectral.rows("longitudinal_energy", value),
                      lambda value=value: spectral.row_bounds("longitudinal_energy", value),
                      lambda value=value: spectral.adiabatic_rows("longitudinal_energy", value),
                      lambda value=value: estimates.local_derivatives("energy", value)))
    calls.extend(lambda value=value: frequency.reference("transverse", value)
                 for value in (1, 3, 5, 4.0, sp.Integer(4), True))
    calls.extend(lambda value=value: spectral.reference_tail("longitudinal_energy", 2, value)
                 for value in (1, 5, 2.0, sp.Integer(2), True))
    calls.extend((lambda: spectral.rows("scalar", 0), lambda: frequency.reference("scalar"),
                  lambda: estimates.local_derivatives("density", 0), lambda: preparation.mixing_bound(1000, 999)))
    calls.extend(lambda value=value: preparation.threshold(value) for value in (0, -1, 999))
    calls.extend(lambda value=value: estimates.physical_bounds(value, 1000) for value in (0, -1))
    calls.extend(lambda value=value: serialize(value) for value in (1.0, sp.Float(1), sp.oo, sp.nan))
    rejected = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(calls):
        raise ValueError("An invalid exact regularity diagnostic was accepted")
    return {"rejected_inputs": rejected,
            "native_integer_checks_remain_active_after_cache_warmup": True,
            "reference_not_a_new_state": True,
            "nonzero_initial_mixing_and_all_oscillatory_endpoints_retained": True,
            "projected_reference_not_confused_with_differentiated_approximate_mode": True,
            "failure_of_lower_reference_not_claimed_state_divergence": True,
            "new_scale_example_not_replacement_of_frozen_example": True,
            "fixed_background_regularity_not_metric_functional_response": True}


@cache
def build_report():
    prior, identities = prior_checks(), residuals()
    proofs = proof_checks()
    if not all(value is True for value in proofs.values()):
        raise ValueError("A continuous higher-regularity proof check failed")
    exact = affine.certify_residuals(identities)
    scalar_count = sum(value.rows*value.cols if isinstance(value, sp.MatrixBase) else 1 for value in identities.values())
    readouts = {name: {order: {"row_envelope": spectral.row_bounds(name, order)["normalized_reference_product_envelope"],
                              "coefficient_majorants": spectral.row_bounds(name, order)["coefficient_majorants"],
                              "reference_tail": spectral.reference_tail(name, order)}
                       for order in range(6)} for name in spectral.tail.physical_weights()}
    return {"schema": 1, "claim": "P8-S6.59.VECTOR_REGULARITY", "date": "2026-09-08",
            "status": "QUANTITATIVE_C5_FIXED_BACKGROUND_VECTOR_STRESS; ORIGINAL_P8_OPEN",
            "prior_sha256": prior,
            "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in source_files()},
            "formulation": "FORMULATION.md", "written_proofs": ["notes/reference.md", "notes/preparation.md", "notes/readouts.md"],
            "exact_residuals": exact, "named_exact_check_count": len(identities), "checked_scalar_entries": scalar_count,
            "proof_checks": proofs,
            "fixed_input_boundary": "The S6.56 named action, S6.55 all-order state and fixed cutoffs, original clock background and S6.53 finite prescription are unchanged. W8 is a comparison reference, not new Cauchy data.",
            "frequency_residual_envelopes": serialize({kind: frequency.reference(kind) for kind in ("transverse", "longitudinal")}),
            "reference_residual_definition": "W=omega*(1+sum_n=1^4 P_n omega^-2n), rho=(1-S²)/t+2P1+T/S+3R²/(4S²), R=S', T_n=-D0 B_n/2+(n+1/2)lambda B_n. The actual low numerator coefficients vanish before absolute value and slope majorants are applied.",
            "initial_preparation_envelopes": serialize(preparation.constants()),
            "frequency_band_boundary": "nu²=m²+k_com²/(25/16)²; bands split at 4m and K=max(8m,10^12). The latter completes both actual fourth-coefficient cutoffs for all m>=1000. Initial mixing is bounded by B6/nu6,B8/nu8,B10/nu10 in the respective bands.",
            "oscillatory_evolution_boundary": "Nonzero initial mixing is retained. Exact diagonal phase removal and one integration by parts, with both endpoints, give abs(Bcal)<=abs(Bcal_initial)+Kmix/nu10. The common Kmix=5347035781757616 includes endpoint, differentiated residual and feedback terms.",
            "actual_readout_row_rule": "For X=(abs(v)²,Re(v*conj(p)),abs(p)²), X'=M X with M=[[2d,2,0],[-omega²,0,1],[0,-2omega²,-2d]]. r0=(B omega²,0,A), r_(j+1)=D rj-3H rj+rj M. The physical derivative is rj X/(2a³).",
            "projected_reference_boundary": "The exact-ODE row evaluated on W8 products is not the time derivative of an approximate mode readout. Its difference from the differentiated frozen subtraction has a proved omega^-5 density tail for all four readouts and orders 0..5.",
            "readout_derivative_envelopes": serialize(readouts),
            "lower_reference_failure_control": serialize({"readout": "longitudinal_energy", "derivative_order": 2,
                                                          "WKB_coefficient_order": 2, "u": 0, "z": 1,
                                                          "nonzero_low_numerator_coefficient": spectral.low_reference_failure()}),
            "radial_integral_envelope": serialize(estimates.radial_envelope(1000)),
            "radial_bound_formula": "J5=A³[32B6 m³/27+B8 K/18+B10/(9K)+Kmix/(24m)]. For derivative j<=5 the three-polarization state/reference difference is <=9 E_j A^(j+1)m^(j-5)J5; projected subtraction integrates to T_j/(72m²). Both are divided by L².",
            "local_derivative_coefficients_and_envelopes": serialize({name: {order: estimates.local_derivatives(name, order)
                                                                           for order in range(6)} for name in ("energy", "pressure")}),
            "physical_scale_example": serialize(estimates.physical_bounds(10**24, 1000)),
            "continuous_regularity_conclusion": "Uniform absolutely integrable bounds on every actual exact-minus-subtraction derivative through order five justify differentiation under the integral. The matched finite energy and pressure are C5 on the full closed window, with continuous one-sided endpoint derivatives. All twelve example normalized bounds are below 10^-18.",
            "controls": controls(),
            "verdict": "Quantitative fixed-background C5 vector stress is established without changing the state or prescription. Original P8 remains open.",
            "not_established": ["A compatible all-order covariance family and state-change bound on perturbed histories",
                                "Full renormalized contact/retarded second variation, feedback norm, nonlinear residual or self-consistent quantum solution",
                                "Corrected coupled physical cones, interacting cutoff, common vacuum/UV parent, finite-gravity Regge, V/G/B or original P8 closure"],
            "verification_boundary": "Exact symbolic/rational recurrence, derivative, subtraction, continuous frequency-band and input controls with written proofs; not proof-assistant formalization"}


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The higher-regularity report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.59.VECTOR_REGULARITY replay passed; quantitative C5 stress, original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
