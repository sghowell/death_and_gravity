"""Read-only exact replay of the free-Maxwell conformal-strip QSEI gate."""

import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8a import verify as geometric_prior
from p8a_see_qsei import sampling as prior_clock
from p8a_see_qsei import verify as clock_prior

from . import bounds, focusing, functional, independent, stress

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"maxwell-qsei.json"
PINS = {"A1": "f719c98068848003e296260d3ff6b03d32acad40662b9c9980f61c39175f209e",
        "A12": "f72698170ee7ca0468b42ceaa098e69e1a26054c6ffdbf6321f1968eaefd6584"}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def serialize(value):
    if isinstance(value, dict):
        return {key: serialize(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [serialize(item) for item in value]
    if isinstance(value, bool):
        return value
    return str(value)


@cache
def prior_checks():
    for name, module in (("A1", geometric_prior), ("A12", clock_prior)):
        if sha(module.REPORT) != PINS[name]:
            raise ValueError("a pinned geometric or general-clock input changed: "+name)
    # A.12 replays its source-hashed lineage, including A.1. This new field
    # gate imports no scalar numerical QSEI or scalar physical prescription.
    clock_prior.validate_report(json.loads(clock_prior.REPORT.read_text()), clock_prior.build_report())
    geometric_prior.validate_report(json.loads(geometric_prior.REPORT.read_text()), geometric_prior.build_report())
    return PINS.copy()


def exact_checks():
    groups = {"Maxwell_reference_and_FK_variation": stress.identities(),
              "exact_conformal_and_proper_functional": functional.identities(),
              "curvature_envelope_scaling": bounds.scaling_identities(),
              "conditional_SEE_and_A1_dictionary": focusing.identities(),
              "pinned_general_proper_clock": prior_clock.clock_identities()}
    for name, group in groups.items():
        if any(sp.simplify(value) != 0 for value in group.values()):
            raise ValueError("a free-Maxwell exact identity failed: "+name)
    return {name: dict.fromkeys(group, "0") for name, group in groups.items()}


def ibp_coefficients():
    t = sp.Symbol("t", real=True)
    h = sp.Function("H")(t)
    jets = tuple(sp.diff(h, t, j) for j in range(4))
    return {name: {",".join(map(str, powers)): str(value)
                   for powers, value in sp.Poly(expression, *jets).terms() if value}
            for name, expression in functional.ibp_coefficients(h, t).items()}


def source_coefficient_controls():
    result = {}
    for name, values in focusing.calibration()["signed_source_controls"].items():
        coarse = values["separate_source_penalty"]
        quantum = sp.simplify(values["Q0"]-coarse)
        raw = sp.simplify(values["raw_constant_before_nonnegative_coarsening"]-quantum)
        result[name] = {"Q2_pi_squared": sp.simplify(sp.pi**2*values["Q2"]),
                        "Q0_quantum_pi_squared": sp.simplify(sp.pi**2*quantum),
                        "coarsened_source": coarse, "raw_source": raw,
                        "upward_coarsening_gap": sp.simplify(coarse-raw)}
    return result


def checked_constants():
    old = json.loads(geometric_prior.REPORT.read_text())
    ratio = sp.Rational(old["parameters"]["partition_ratio"])
    if ratio != sp.Rational(3, 4) or old["parameters"]["rho0_times_tau_squared"] != "0":
        raise ValueError("the frozen A1 tail family does not match the exclusion control")
    fractions = independent.replay()
    primary = {"stress_jet_coefficients": stress.jet_coefficients(),
               "proper_IBP_coefficients": ibp_coefficients(),
               "radiation_control": functional.radiation_control(),
               "envelope_controls": bounds.calibration(),
               "source_controls": source_coefficient_controls(),
               "focusing_incompatibility": focusing.calibration()["Hmax_times_tau_at_most_one"]}
    if primary["focusing_incompatibility"]["A1_zeta_zero_gradient_lower"] != 3/ratio:
        raise ValueError("the old A1 gradient floor was not preserved")
    for name, value in primary.items():
        if serialize(value) != fractions[name]:
            raise ValueError("independent Fraction Maxwell replay disagrees: "+name)
    return {**serialize(primary), "independent_Fraction_replay": fractions}


def build_report():
    inputs = prior_checks()
    identities, constants = exact_checks(), checked_constants()
    paths = sorted(ROOT.glob("src/p8a_maxwell/*.py"))+sorted(ROOT.glob("tests/*.py"))
    paths += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    return {
        "schema": 1, "claim": "P8-A.16", "date": "2026-09-06",
        "status": "FREE_MAXWELL_CONFORMAL_ALL_HADAMARD_H2_QSEI_CERTIFIED; CONDITIONAL_FOCUSING_AND_P8_OPEN",
        "prior_sha256": inputs,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in paths},
        "formulation": "FORMULATION.md", "written_proof": "notes/proof.md", "source_audit": "notes/sources.md",
        "exact_residuals": identities, "derived_constants": constants,
        "theorem": {
            "field": "one free real physical Maxwell field, two polarizations, energy=(electric^2+magnetic^2)/2",
            "geometry": "smooth spatially flat conformal strip I_eta times R^3, a>0; comoving timelike geodesics",
            "target_states": "all physical Maxwell Hadamard states on that strip; no target homogeneity, quasifree, zero-mean or SEE assumption",
            "global_target_state_extension_assumed": False,
            "reference": "actual conformal transport of the restricted Minkowski Maxwell vacuum; zero flat reference stress",
            "renormalization_family": "hbar/(2880pi^2) times [62 H3_FK + beta_M I_FK], beta_M real constant and explicit",
            "finite_beta_inherited_from_scalar_gamma": False,
            "optional_named_specialization": "zero type D means beta_M=0; not silently fixed for the generic theorem",
            "clock": "dt=a deta, F=a^(-3/2) f(t(eta)); no approximation",
            "absolute_bound": "integral f^2 E_target >= -hbar/(8pi^2) integral |f''-2Hf'+(3H^2/4-3Hdot/2)f|^2 + integral f^2 E_conf",
            "samplers": "all real smooth compact samplers; real H2_0 extension on compact proper subintervals, both zero boundary traces",
            "optimal_coefficient_claimed": False,
            "radiation_control": "R=0 removes finite beta, not the anomaly; absolute IBP zeroth coefficient=4601/(1280t^4)",
            "counterterm_source_boundary": "printed Markowicz H1 not imported; independently conserved FK inverse-metric R^2 variation defines beta_M",
        },
        "conditional_geometric_dictionary": {
            "equation": "G_FK+Lambda g_FK=-kappa(T_Maxwell+T_other), kappa>0",
            "Ricci": "R_UU=Lambda-kappa(E_Maxwell+E_other)",
            "constants": "Q2=kappa B2; Q0=kappa B0+max(Lambda,0)+kappa max(-ell_other,0), with E_other>=ell_other",
            "separate_cosmological_Einstein_and_extra_source_choices_retained": True,
            "new_SEE_solution_asserted": False,
            "initial_pointwise_Ricci_requirement_removed": False,
            "all_A1_extended_normal_windows_verified": False,
            "new_K_threshold_or_global_causal_hypotheses_verified": False,
            "naive_small_quantum_constant_example_focuses": False,
            "negative_control": "Hmax*tau<=1 forces |K|tau<=3, below the frozen zeta=0 A1 gradient floor4 even at zero quantum ratios",
        },
        "verification_boundary": [
            "A1 and A12 report hashes checked; both replayed read-only with their source-hashed lineage",
            "Symbolic identities cross-checked by Fraction-only curvature-contraction and source-sign polynomial engine",
            "Independent coordinate audit checks every diagonal curvature stress, conservation, physical sign, radiation integration and strip measure",
            "Distributional positivity, local algebra identification, Sobolev density and causal geometry are written proofs, not proof-assistant formalizations",
        ],
        "not_established": [
            "interacting QED, charged or massive realistic fields, boundaries or nontrivial topological flux sectors",
            "an exact new Maxwell semiclassical Einstein solution or controlled EFT validity",
            "a concrete cosmological focusing example, maximal extension or incompleteness conclusion",
            "a QEIs-only singularity theorem without initial geometric/pointwise hypotheses",
            "arbitrary boosted/null-geodesic or other-hypersurface bounds",
            "completion of original P8a or P8",
        ],
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("the free-Maxwell certificate differs from exact read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8(a) A.16: free-Maxwell all-Hadamard conformal H2 QSEI replay passed; cosmological focusing and P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
