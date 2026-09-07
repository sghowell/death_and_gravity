"""Read-only replay of the actual physical preparation-cost certificate."""
import argparse
import hashlib
import json
from fractions import Fraction
from functools import cache
from pathlib import Path

import sympy as sp
from p8_variable_prepared import verify as prior

from . import construction, core, gramian, independent, spectral

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"physical-preparation-cost.json"
PREPARED_SHA = "a8d89126980dd8ec458b3161bc8b397f74211736f02e64d4e9e29b95e43a0ade"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def serialize(value):
    if isinstance(value, (Fraction, sp.Rational)):
        return str(value)
    if isinstance(value, sp.Basic):
        if (value.has(sp.Float, sp.oo, sp.zoo, sp.nan) or value.is_number is not True
                or value.is_real is not True or value.is_finite is not True):
            raise TypeError("Only exact finite real report expressions are admitted")
        return str(value)
    if isinstance(value, dict):
        return {key: serialize(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [serialize(item) for item in value]
    if value is None or isinstance(value, (str, bool, int)):
        return value
    raise TypeError("Only exact finite report data are admitted")


def source_files():
    return (sorted(ROOT.glob("src/p8_preparation_cost/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(prior.REPORT) != PREPARED_SHA:
        raise ValueError("The frozen prepared-sector certificate changed")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return {"S6_23_prepared_and_recursive_response_variable_parent_contract": PREPARED_SHA}


def controls():
    construction.source_cost(1)
    spectral.calibration()
    calls = [lambda value=value: construction.source_cost(value)
             for value in (True, 1.0, sp.Float("1"), sp.oo, sp.nan, -1)]
    calls += [lambda value=value: spectral.lowpass_bounds(value)
              for value in (True, 0.1, -1, 101, sp.oo, "1")]
    calls += [lambda: construction.hermite_polynomial((1, 0, 0), 1),
              lambda: construction.hermite_polynomial((1.0, 0, 0, 0), 1),
              lambda: construction.hermite_polynomial((1, 0, 0, 0), 0),
              lambda: construction.hermite_polynomial((1, 0, 0, 0), 5),
              lambda: gramian.enclosure(order=True),
              lambda: gramian.enclosure(order=11),
              lambda: gramian.enclosure(steps=63),
              lambda: gramian.enclosure(precision=127),
              lambda: independent.Interval(0.0),
              lambda: independent.positive_gramian([[(0, 0)]*4 for _ in range(4)]),
              lambda: spectral.minimum_cost_envelope(100, 1, 1, 1, 1),
              lambda: spectral.minimum_cost_envelope(100, 1, 1, 2, 10, 2),
              lambda: spectral.finite_delta_error(0, 1, 1),
              lambda: spectral.finite_delta_error(sp.Rational(1, 10**8), 1, 1),
              lambda: spectral.finite_delta_error(sp.Rational(1, 10**9), 1, 1, target="raw_g"),
              lambda: spectral.source_units(0, 1),
              lambda: spectral.source_units(1, 0),
              lambda: spectral.sinc_remainder(1, sp.Rational(1, 2)),
              lambda: serialize(0.1), lambda: serialize(sp.Float("0.1")),
              lambda: serialize(sp.Symbol("x")), lambda: serialize(sp.I)]
    count = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError):
            count += 1
    if count != len(calls):
        raise ValueError("An invalid exact-domain or evidence input was admitted")
    return {"rejected_inputs": count, "zero_source_zero_target_admitted": True,
            "source_K_zero_not_imported_from_wider_ancestor_domain": True,
            "Omega_zero_lowpass_exactly_zero": True,
            "H0_2_minimum_at_boundary_budget_not_assumed": True,
            "computed_band_moments_and_optimizer_not_claimed": True,
            "continuous_K_not_inferred_from_a_point_Gramian": True,
            "parameterized_cost_function_does_not_certify_a_loading_map": True}


@cache
def build_report():
    pins = prior_checks()
    audits = [ROOT/"tests"/"test_cost_source_independent_audit.py",
              ROOT/"tests"/"test_cost_gramian_independent_audit.py"]
    if any(not path.is_file() for path in audits):
        raise ValueError("Both separately authored scientific audits are required before promotion")
    literal, loading, fourier = core.checks(), construction.checks(), spectral.checks()
    ancestry_bridge = construction.parent_checks()
    independent_loading = construction.independent_checks()
    norms = independent.checks()
    if core.calibration().keys() != norms.keys() or any(
            core.calibration()[key] != sp.Rational(value) for key, value in norms.items()):
        raise ValueError("The independent continuous energy reconstruction disagrees")
    validated = gramian.certificate()
    separate_matrix = independent.positive_gramian(validated["scaled_outer_gramian"])
    spectrum = spectral.calibration()
    for key, expected in (
        ("gramian_lower_input", core.calibration()["uniform_Frobenius_gramian_lower"]),
        ("gramian_upper_input", core.calibration()["uniform_Frobenius_gramian_upper"]),
        ("regular_light_norm_lower_input", construction.calibration()["regular_light_L2_lower"]),
        ("regular_even_norm_lower_input", construction.calibration()["regular_even_L2_lower"])):
        if spectrum[key] != expected:
            raise ValueError("A spectral calibration was not proved by its companion module")
    if (spectrum["Hermite_control"]["source_L2_upper"] != construction.calibration()["source_L2"]
            or spectrum["Hermite_control"]["source_second_derivative_L2_upper"]
            != construction.calibration()["source_second_derivative_L2"]):
        raise ValueError("The Fourier derivative-tail input changed")
    return {
        "schema": 1, "claim": "P8-S6.29.COST", "date": "2026-09-07",
        "status": "ACTUAL_FULL_SOURCE_COST_AND_SPECTRAL_INFIMUM_BOUNDS; ORIGINAL_S6_P8_OPEN",
        "prior_sha256": pins,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": ["notes/gramian.md", "notes/construction.md", "notes/optimization.md"],
        "source_audit": "notes/sources.md",
        "exact_residuals": {
            "literal_full_generator_and_energy_chart": dict.fromkeys(literal["residuals"], "0"),
            "literal_source_Hermite_and_boundary_jets": dict.fromkeys(loading["residuals"], "0"),
            "additional_parent_source_bridges": dict.fromkeys(ancestry_bridge, "0"),
            "Fourier_kernel_KKT_and_target_algebra": dict.fromkeys(fourier["residuals"], "0")},
        "strict_continuous_margins": serialize({
            "full_source_map_energy": literal["strict_margins"],
            "Hermite_and_true_light_source_cost": loading["strict_margins"],
            "spectral_and_finite_delta_cost": fourier["margins"]}),
        "validated_midpoint_Gramian": serialize(validated),
        "independent_rational_Gramian_LDL": serialize(separate_matrix),
        "literal_rational_to_Arb_coefficient_bridges": serialize(core.coefficient_bridges()),
        "continuous_K_energy_and_singular_value_bounds": serialize(core.calibration()),
        "independent_Fraction_energy_replay": serialize(norms),
        "Hermite_and_source_cost": serialize(construction.calibration()),
        "independent_Fraction_source_Hermite_replay": serialize(independent_loading),
        "spectral_envelopes_and_exact_remainder": serialize(spectrum),
        "controls": controls(),
        "unchanged_physical_problem": {
            "action": "S6.20 full two-tensor action, positive Einstein terms and actual g matter frame",
            "source": "sigma=tau^2 Pi/M^2; literal external TT source; physical derivative-state input(0,2,0,0)",
            "preparation_domain": "J=(-1/50,-1/100), K in[1,4], all four physical data zero on the left",
            "limiting_loading_map": "Full punctured delta0 equation into true Psi-normalized four-state endpoint frame",
            "source_class": "Real zero-extended H0^2(J); signed; compact support may meet J boundary; fixed independently of delta",
            "target": "Lv in true Frobenius coordinates, not raw Y or raw physical g",
            "frequency": "0<=Omega<=100 in u=T/tau; actual sigma Fourier energy, measure domega/(2pi)"},
        "proof_chain": {
            "midpoint": "64closedTaylorsteps,order18,256bitArb; continuous interval Taylor remainder and four strict LDL pivots",
            "continuous_momentum": "duration-resolved energy bound exp(11(t-s)/20), K derivative<1/5000, inputnorm<1/2000",
            "uniform_Gramian": "9/(4e12) I <= C C* <=16/25 I on full K interval",
            "construction": "degree11 physical f Hermite jets0..5, exact g/source recovery, C1 zero-extended H0^2 source",
            "light_cost": "actual symplectic source kernels give >19||v|| and >900|v_even| for nonzero targets/components",
            "spectrum": "sinc operator norm<=1/3, tail>=2I/3; exact L2 optimizer/dual, H0^2 infimum equality for strict budgets",
            "calibrated_budget": "unitv,S1e6: unconstrained spectral optimizer has norm^2<=2e12/3; infimum>722/3 and<=4e12/9",
            "finite_rank_interface": "nine uncomputed actual moments; sinc degree8 kernel error<=1/(3*11!)",
            "finite_delta": "fixed source:42e+delta(8600||v||+12600000S); actual prepared target uses8400; exampledelta1e-17,S1e6<1e-3"},
        "not_established": [
            "No numerically computed band moments, spectral optimizer or sharp minimum",
            "No automatic H0^2 attainment or feasibility at the exact minimum L2 budget",
            "No claim the actual controls are mostly in the calibrated temporal band",
            "No fixed peak/derivative budget for approximating the spectral infimum",
            "No raw physical derivative error substituted for the bounded normalized endpoint error",
            "No nonlinear stress/backreaction or quantum preparation/production bound",
            "No generic sourced local EFT, rolling gap/cutoff, vacuum positivity, C/D operator matching or original P8 closure"],
        "verification_boundary": "Written full operator, interval-ODE, continuous-parameter, Hermite, Hilbert-space and Fourier proofs; exact symbolic/Fraction replays and independent source/Gramian audits. Finite KKT tests are not substituted for the functional theorem; not FORMALIZED.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The physical preparation-cost certificate differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.29.COST: physical preparation-cost replay passed; original S6/P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
