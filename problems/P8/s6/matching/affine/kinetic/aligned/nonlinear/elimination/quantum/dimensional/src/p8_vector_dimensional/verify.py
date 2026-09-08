"""Read-only dimensional vector poles and finite local matching diagnostics."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as affine
from p8_vector_subtraction import verify as parent

from . import adm, heat, local, mass, metric, replay

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"dimensional-local-vector-matching.json"
PARENT_SHA = "49cb4fc8230f472211f643f6b7ce15ce54804b6223eb404c79f38b07ed3282f4"
sha, serialize = affine.sha, affine.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_vector_dimensional/*.py"))
            +sorted(ROOT.glob("tests/*.py"))+sorted(ROOT.glob("*.md"))
            +sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen finite vector subtraction certificate changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {"S6_51_fully_rebuilt": PARENT_SHA,
            "original_action_physical_frame_and_state_preparation_retained": True,
            "no_frozen_finite_potential_or_lapse_jet_changed": True}


def residuals():
    return {**adm.variations(), **local.checks(), **metric.ordinary_matching()["identities"],
            **mass.checks(), **mass.frozen_flat_finite_match(), **replay.physical_checks(),
            **replay.gamma_checks(), **heat.checks()}


def controls():
    bad_orders = (True, False, sp.true, sp.false, 1.0, sp.Float(1), "1", -1, 3, sp.Integer(1), sp.nan)
    calls = [lambda value=value: local.radial_polynomial(sp.Integer(1), value) for value in bad_orders]
    bad_coefficients = (True, False, sp.true, sp.false, 1.0, sp.Float(1), "1", sp.oo, -sp.oo, sp.zoo, sp.nan, {"bad": 1})
    calls += [lambda value=value: local.radial_polynomial(value, 0) for value in bad_coefficients]
    calls += [lambda value=value: local.radial_polynomial(value, 0) for value in
              (1/(1+local.z), sp.sin(local.z), sp.exp(local.z))]
    calls += [lambda value=value: local.reference(value) for value in ("scalar", None, True)]
    calls += [lambda value=value: serialize(value) for value in (1.0, sp.Float(1), sp.oo, sp.nan)]
    rejected = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(calls):
        raise ValueError("An invalid order, coefficient, mode or report value was accepted")
    return {"rejected_inputs": rejected,
            "early_polarization_count_limit_loses_finite_term": "2+4b_N",
            "ordinary_second_order_radial_Ward_defect": serialize(metric.ordinary_matching()["ordinary_covariantly_subtracted_local_terms"][1]["radial_Ward_defect"]),
            "ordinary_fourth_order_radial_Ward_defect": serialize(metric.ordinary_matching()["ordinary_covariantly_subtracted_local_terms"][2]["radial_Ward_defect"]),
            "scalar_gradient_of_constrained_vector_not_dropped": True,
            "counterterm_variation_precedes_dimension_limit": True,
            "ordinary_finite_diagnostic_not_actual_clock_finite_matching": True}


@cache
def build_report():
    pins = prior_checks()
    identities = residuals()
    exact = affine.certify_residuals(identities)
    geometry = metric.variations()
    return {"schema": 1, "claim": "P8-S6.52.VECTOR_DIMENSIONAL", "date": "2026-09-08",
            "status": "ACTUAL_VECTOR_LOCAL_POLES_AND_ORDINARY_FINITE_MATCHING_DIAGNOSTIC_CERTIFIED; ORIGINAL_P8_OPEN",
            "prior_sha256": pins,
            "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in source_files()},
            "formulation": "FORMULATION.md", "written_proofs": ["notes/proof.md"],
            "exact_residuals": exact, "named_exact_check_count": len(identities),
            "checked_scalar_entries": len(identities),
            "literal_component": "Unchanged retained vector and physical clock of S6.42/50/51. First local lapse/pressure variations and their formal D=3-2epsilon_DR continuation, not a state on a noninteger-dimensional Hilbert space or a full affine determinant.",
            "dimensional_canonical_normalization": "gT²=a_s^(D-2), gL²=a_s^D*m²*q/(q+m²), rates H(D-2)/2 and H[(D-2)/2+z]. D-1 transverse polarizations and one longitudinal mode. The action variation, WKB coefficients and physical subtraction coefficients are independently replayed at D=3.",
            "radial_Gamma_prescription": "For an order-2n bracket sum c_j(D)z^j, G_n(D)=2sqrt(pi)*sum c_j(D)*(D/2)_j/Gamma(n+j-1/2). Relative to m^(4-2n)/(64pi²), P_n=(-1)^(2-n)G_n(3)/(2-n)!, E_n=-2*(-1)^(2-n)G_n'(3)/(2-n)!, radial finite=P_n[Harmonic(2-n)-ell]+E_n, ell=log(m²/mu²). This is the stated dimensional Gamma continuation, not an ordinary convergent integral assigned a finite value.",
            "actual_lapse_radial_poles_and_finite_parts": serialize(local.integrated()),
            "physical_pressure_radial_poles_and_finite_parts": serialize(local.integrated_pressure()),
            "radial_finite_part_boundary": "A component-wise radial finite part need not equal the covariantly matched stress component. Counterterm metric variation in D+1 dimensions can contribute an additional finite term. The zero-order actual lapse term agrees with the frozen S6.47 finite potential and N jet.",
            "dimension_dependent_local_heat_coefficients_and_finite_action": serialize(heat.coefficients()),
            "ordinary_counterterm_variation": serialize({name: geometry[name] for name in
                ("rho_pole_counterterm_in_D", "pressure_pole_counterterm_in_D", "energy_evanescent_counterterm",
                 "pressure_evanescent_counterterm", "independent_finite_heat_action_lapse_energy")}),
            "ordinary_covariantly_matched_finite_local_diagnostic": serialize(metric.ordinary_matching()["ordinary_covariantly_subtracted_local_terms"]),
            "ordinary_finite_matching_boundary": "The control sets a_N=b_N=0. Vary the four-dimensional covariant pole coefficients as invariants in D+1 dimensions before the limit. Subtract their epsilon coefficient from the radial finite components. The Ward identity is restored and an independent finite heat-action variation agrees. This is not the actual clock-dependent finite derivative matching.",
            "actual_clock_mass_covariant_pole_reconstruction": serialize(mass.reconstruction()),
            "actual_clock_mass_pole_boundary": "For alpha=a_N,beta=b_N, extra_P2=-(3H²+4H')alpha-(9H²+2H')beta; extra_P4=(beta-alpha)E/2+(alpha+beta)boxR/12, E=H'^2-6H²H'-2HH''. The covariant one-form/scalar-gradient coincidence reconstruction agrees with the independent radial result. Only local dimensional poles are used in the scalar EOM/trace argument; finite anomalies are not set to zero.",
            "finite_extension_boundary": "An evanescent continuation of mass counterterms can shift finite terms even when its four-dimensional pole is unchanged. Any actual clock-dependent finite extension must preserve S6.47's already frozen potential and lapse jets. That derivative matching remains open, not a user-choice blocker.",
            "local_boundary_terms": "Metric variations are compactly supported. boxR is a total divergence in this local Euler variation, not a noncompact infinite-time flux assumed to vanish; S6.49's boundary statement is unchanged.",
            "controls": controls(),
            "verdict": "Actual clock-sensitive local lapse poles, dimensionally normalized radial coefficients, preservation of the frozen flat finite lapse term, and the ordinary-Proca covariant finite matching diagnostic are established. Actual finite clock-curvature matching and original P8 remain open.",
            "not_established": ["Full finite clock-mass curvature counterterm matching or a matched finite actual-clock quantum energy/pressure bound",
                                "All-order Hadamard admissibility, all-time state control, other field/higher loops, corrected bounce/constraints/cones or interacting cutoff",
                                "Lorentz-invariant vacuum, finite-gravity Regge remainder, V/G/B admissibility or original P8 closure"],
            "verification_boundary": "Exact D-dimensional action variations, physical WKB replay, Gamma continuation, finite counterterm variations, independent local heat action and scalar-gradient pole reconstruction; written proof, not proof-assistant formalization"}


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The dimensional local vector matching report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.52.VECTOR_DIMENSIONAL replay passed; local matching controls, original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
