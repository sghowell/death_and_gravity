"""Read-only physical vector energy and finite curved mode comparison."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as affine
from p8_vector_curvature import verify as parent

from . import comparison, energy, wkb

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"physical-energy-mode-difference.json"
PARENT_SHA = "4f6bcb4fd4d93d4df5deef50c35c7a83f340a6280625001860a8cc1bf0f75669"
sha, serialize = affine.sha, affine.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_vector_state/*.py"))
            +sorted(ROOT.glob("tests/*.py"))+sorted(ROOT.glob("*.md"))
            +sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen local vector curvature certificate changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {"S6_49_fully_rebuilt": PARENT_SHA,
            "original_all64_action_physical_frame_and_clock_retained": True,
            "no_frozen_action_or_finite_counterterm_changed": True}


def residuals():
    out = {**energy.checks(), **wkb.checks(), **wkb.physical_chain_rule_checks(), **comparison.checks()}
    d = energy.canonical()
    out["zero_momentum_canonical_energy_chart_limit"] = sp.factor((d["rhoL"]-d["rhoT"]).subs(energy.q, 0))
    out["zero_momentum_canonical_pressure_chart_limit"] = sp.factor((d["pL"]-d["pT"]).subs(energy.q, 0))
    return out


def controls():
    bad = (True, False, sp.true, sp.false, 1.0, sp.Float(1), "1", sp.I,
           sp.oo, -sp.oo, sp.zoo, sp.nan, sp.Symbol("unproved"), 0, -1)
    calls = [lambda value=value: comparison.physical_bounds(value, 1000) for value in bad]
    calls += [lambda value=value: comparison.physical_bounds(10**12, value) for value in (*bad, 999)]
    calls += [lambda value=value: serialize(value) for value in (1.0, sp.Float(1), sp.oo, sp.nan)]
    calls += [lambda value=value: wkb.box_bound(value) for value in
              (1/(1+wkb.u**4), 1/(1-wkb.u), sp.Float(1), sp.Symbol("unknown"))]
    calls += [lambda: wkb.frequency("scalar")]
    rejected = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(calls):
        raise ValueError("An invalid scale, mode, coefficient envelope or report value was accepted")
    return {"rejected_inputs": rejected,
            "clock_mass_variation_cannot_be_omitted": energy.zero_order_pole()["omitting_clock_variation_pole_difference"] != 0,
            "second_order_reference_has_nonzero_order_inverse_frequency_squared_residual":
                all(wkb.frequency(kind)["P4"] != 0 for kind in ("transverse", "longitudinal")),
            "physical_lapse_energy_not_canonical_oscillator_Hamiltonian": True,
            "all_momenta_bounded_continuously_not_only_sampled": True,
            "finite_Gaussian_preparation_not_called_all_order_Hadamard": True,
            "evolution_difference_not_full_renormalized_energy_or_VGB": True}


@cache
def build_report():
    pins = prior_checks()
    identities = residuals()
    exact = affine.certify_residuals(identities)
    proofs = comparison.proof_checks()
    if not all(value is True for value in proofs.values()):
        raise ValueError("A continuous mode or all-momentum energy proof failed")
    canonical = energy.canonical()
    return {"schema": 1, "claim": "P8-S6.50.VECTOR_STATE", "date": "2026-09-08",
            "status": "PHYSICAL_VECTOR_ENERGY_AND_FINITE_CURVED_MODE_DIFFERENCE_CERTIFIED; ORIGINAL_P8_OPEN",
            "prior_sha256": pins,
            "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in source_files()},
            "formulation": "FORMULATION.md", "written_proofs": ["notes/energy.md", "notes/comparison.md"],
            "exact_residuals": exact, "named_exact_check_count": len(identities),
            "checked_scalar_entries": sum(value.rows*value.cols if isinstance(value, sp.MatrixBase) else 1 for value in identities.values()),
            "proof_checks": proofs,
            "literal_component": "Unchanged S6.42 retained canonical vector at fixed light fields, actual rolling physical FLRW metric, zero Gaussian mean and vanishing source first variation. The physical lapse variation retains nonzero a_N and b_N. Other determinants and affine contact/measure terms are not supplied.",
            "actual_clock_mass_direction": serialize(energy.actual_coefficients()),
            "physical_canonical_energy_and_pressure": serialize({name: canonical[name] for name in
                ("z", "omega2", "rateT", "rateL", "rhoT", "rhoL", "pT", "pL")}),
            "local_pole_control": serialize(energy.zero_order_pole()),
            "local_pole_boundary": "d=4-2epsilon_DR: the zero-point Gamma residue is -m0^4/(64pi² epsilon_DR). Including two transverse and one longitudinal mode gives 3+3a_N/2+9b_N/2=3+20/(9h)=C+C_N, not C. No finite or derivative counterterm matching is inferred.",
            "domain_and_state_preparation": "u in [-1/2,1/2], every real comoving momentum, m=m0*tau>=1000. Positive W=omega*(1+P2/omega²+P4/omega^4); exact modes match f=(2W)^(-1/2)exp(-i integral W du) and f' at u=-1/2. Wronskian i and positive Gaussian polarization covariance. No all-order Hadamard claim.",
            "exact_WKB_reference_coefficients": serialize({kind: {name: wkb.frequency(kind)[name] for name in
                ("U", "P2", "P4", "B2", "B4", "A", "B", "S", "reference_log_rate", "residual")}
                for kind in ("transverse", "longitudinal")}),
            "continuous_WKB_envelopes": serialize({kind: wkb.bounds(kind) for kind in ("transverse", "longitudinal")}),
            "reference_and_residual_bounds": "omega/2<W<3omega/2, abs(W'/W)<4, abs(f''+(omega²-U)f)<=C*abs(f)/omega^4, C=2000000. These follow from an exact rational residual and continuous polynomial-box bounds, not an unbounded asymptotic expansion.",
            "exact_evolution_comparison": "nu_k²=m²+k_com²/(25/16)². The variation-of-constants column-sum exponent J_k<=2C/nu_k^5. Both squared-mode and squared-physical-derivative differences are <=8C/nu_k^5 times their reference values. The physical reference energy per polarization is <=2omega/a_s³.",
            "all_momentum_integral_and_units": "integral d³k/(2pi)³ /nu_k^4=(25/16)³/(8pi*m). Summing all three polarizations gives abs(Delta rho)<=12C/m and abs(Delta p)<=108C/(5m), uniformly on I in tau=1 units. Physical densities are tau^-4 times these. Each difference is absolutely integrable; the individual unsubtracted integrals are not assigned finite values.",
            "physical_scale_example": serialize(comparison.physical_bounds(10**12, 1000)),
            "controls": controls(),
            "verdict": "The actual physical vector lapse/pressure forms and a finite all-momentum exact-mode/WKB evolution difference on the rolling bounce interval are established. Original P8 remains open.",
            "not_established": ["Full reference adiabatic subtraction, covariant finite mass/metric counterterm matching, all-order Hadamard admissibility or full renormalized quantum energy",
                                "All-time state tails, other field loops, higher-loop error, corrected bounce/constraints/cones, finite higher matching operators or interacting cutoff",
                                "Lorentz-invariant vacuum, finite-gravity Regge remainder, V/G/B admissibility or original P8 closure"],
            "verification_boundary": "Exact physical lapse/constraint variation, actual fixed-comoving WKB residual, continuous rational envelopes, Wronskian/CCR controls, variation-of-constants proof and analytic momentum integral; not proof-assistant formalization"}


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The physical vector energy comparison report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.50.VECTOR_STATE replay passed; finite mode difference, original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
