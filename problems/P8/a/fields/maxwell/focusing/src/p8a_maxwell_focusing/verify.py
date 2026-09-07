"""Read-only photon-QSEI history-based incompleteness certificate replay."""

import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8a_maxwell import verify as prior

from . import controls, envelope, focusing, history, independent

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"history-focusing.json"
PARENT_SHA = "f59d63524bf1d98998a928fbe5139eef45f3703ebd1e61cbdf470e738ce94b94"


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
    if sha(prior.REPORT) != PARENT_SHA:
        raise ValueError("the pinned actual Maxwell field/source inequality changed")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return {"A16": PARENT_SHA}


def exact_checks():
    groups = {"matched_sampler_history_and_index": history.identities(),
              "relative_jet_cost_and_clock_scaling": envelope.identities(),
              "physical_scale_and_conditional_margin": focusing.identities(),
              "smooth_complete_geometric_countercontrol": controls.identities()}
    for name, group in groups.items():
        if any(sp.simplify(value) != 0 for value in group.values()):
            raise ValueError("an exact history-focusing identity failed: "+name)
    return {name: dict.fromkeys(group, "0") for name, group in groups.items()}


def checked_constants():
    data = focusing.calibration()
    base, sourced = data["zero_source"], data["source_sigma_at_most_one"]
    cost = base["envelope"]
    primary = {"ratio": cost["ratio"], "contraction": base["contraction"],
               "delta": base["delta"], "caps": cost["past_caps"],
               "beta_absolute_cap": data["beta_absolute_cap"],
               "past_reference_loss": cost["past_reference_loss"],
               "future_reference_loss": cost["future_reference_loss"],
               "past_root": cost["past_root_after_ratio_factor"], "future_root": cost["future_root"],
               "past_cost": cost["past_cost"], "future_cost": cost["future_cost"], "total_cost": cost["total_cost"],
               "history_gain": base["history_gain"], "future_gradient": base["future_gradient"],
               "quantum_cost": base["quantum_cost"], "source_weight": cost["source_weight"],
               "zero_source_margin": base["strict_focusing_margin"],
               "sigma_one_margin": sourced["strict_focusing_margin"],
               "cost_coarsening_gap": data["cost_coarsening_gap"],
               "zero_source_margin_above_three_quarters": data["zero_source_margin_above_three_quarters"],
               "sourced_margin_above_one_third": data["sourced_margin_above_one_third"]}
    fractions = independent.replay(json.loads(prior.REPORT.read_text()))
    for name, value in (("constants", primary), ("cubic_moments", history.moments()),
                        ("radical_margins", envelope.radical_margins()),
                        ("geometric_countercontrol", controls.calibration())):
        if serialize(value) != fractions[name]:
            raise ValueError("independent Fraction history-focusing replay disagrees: "+name)
    return {"macroscopic_conditional_calibration": serialize(data),
            "complete_geometric_countercontrol": serialize(controls.calibration()),
            "cubic_moments": serialize(history.moments()),
            "strict_radical_margins": serialize(envelope.radical_margins()),
            "independent_Fraction_replay": fractions}


def build_report():
    inputs = prior_checks()
    identities, constants = exact_checks(), checked_constants()
    paths = sorted(ROOT.glob("src/p8a_maxwell_focusing/*.py"))+sorted(ROOT.glob("tests/*.py"))
    paths += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    return {
        "schema": 1, "claim": "P8-A.17", "date": "2026-09-06",
        "status": "FREE_MAXWELL_HISTORY_BASED_CONDITIONAL_INCOMPLETENESS_CERTIFIED; COSMOLOGICAL_APPLICATION_AND_P8_OPEN",
        "prior_sha256": inputs,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in paths},
        "formulation": "FORMULATION.md", "written_proof": "notes/proof.md", "source_audit": "notes/sources.md",
        "exact_residuals": identities, "derived_constants": constants,
        "theorem": {
            "spacetime": "global smooth spatially flat FLRW product I_t times R^3, a>0; not an unspecified local chart",
            "field": "one free real physical Maxwell field in an arbitrary Hadamard state",
            "prescription": "the explicit A16 real constant beta_M family; no scalar gamma inheritance or silent beta choice",
            "actual_equation": "G_FK+Lambda g_FK=-kappa(T_Maxwell+T_other), E_other>=ell on the required segments",
            "quantum_parameters": "delta=kappa hbar/(8pi^2 tau^2); sigma=tau^2[max(Lambda,0)+kappa max(-ell,0)]",
            "actual_past_requirement": "smooth [-r tau,0] exists; H<=-h/tau and |H^(j)|<=d_j/tau^(j+1), j0..3",
            "conditional_future_requirement": "only if a comoving normal reaches tau: |H^(j)(t)|<=c_j/(tau-t)^(j+1) for0<=t<tau",
            "future_caps_derived_from_initial_patch_or_QSEI": False,
            "initial_pointwise_SEC_or_Ricci_sign_assumed": False,
            "quantitative_contraction_history_assumed": True,
            "samplers": "matched cubic proper-time H2_0 sampler on [-r tau,tau], used only if that smooth segment exists",
            "exact_history_identity": "J[g]+K <=3(||g'||^2+||u'||^2)+Q_Max[f]-3 integral_past(u'-Hu)^2",
            "criterion": "3h+39h^2r/35-18/5-delta(C_P+C_F)-13(1+r)sigma/35>=0",
            "conclusion": "upper endpoint of I_t<=tau; no future timelike curve from t0 has length>tau; stipulated spacetime future timelike geodesically incomplete",
            "global_proof": "nonfocal comoving scale a/a0>0 gives J+K=3 integral(g'-Hg)^2>0 on any smooth segment reaching tau",
            "compact_Cauchy_theorem_imported_on_R3": False,
            "larger_extension_or_curvature_inextendibility_proved": False,
        },
        "named_calibration": {
            "r": "1/100", "h": "3/2", "past_and_future_caps": ["2", "4", "16", "96"],
            "beta_absolute_cap": "1", "delta_upper": "1/100000000", "sigma_upper": "1",
            "zero_source_margin_lower": "3/4", "sourced_margin_lower": "1/3",
            "geometric_not_actual_observed_or_radiation_calibration": True,
            "Planck_dictionary": "kappa=8piG implies delta=ell_P^2/(pi tau^2), c=1",
            "macroscopic_scale_compatibility": True,
            "observational_or_fundamental_EFT_validity_verified": False,
        },
        "countercontrol": {
            "definition": "positive normalized smooth convolution at widthtau/1000 of H=-(31/20tau)(1-p5), p5 quintic extended by0/1",
            "past_pointwise_SEC": "violated, R_UU>0",
            "geometric_history_and_relative_caps_hold": True,
            "future_timelike_and_null_complete": True,
            "actual_allowed_Maxwell_SEE_state_witness": False,
            "reference_vacuum_requires_source_penalty_above_calibration": True,
            "purpose": "geometric assumptions alone do not imply the conclusion; field/source inequality is necessary in this argument",
        },
        "verification_boundary": [
            "A16 report hash checked and full source-hashed lineage replayed read-only",
            "Exact sampler polynomials and all threshold constants independently reconstructed with Fraction from the pinned photon EED polynomial",
            "Independent covariant/source/endpoint and nonvacuity audit included unchanged",
            "Hadamard positivity and H2 density are inherited written field proofs; geometric quantifiers and completeness are proved in the new note",
            "No proof-assistant formalization or numerical spacetime evolution is substituted for the stated analytic proof",
        ],
        "not_established": [
            "a verified observed-universe or new exact Maxwell SEE solution satisfying a nonvacuous long-extension calibration",
            "future relative curvature bounds from initial data or the Maxwell QSEI alone",
            "interacting QED, arbitrary realistic matter, other hypersurfaces, boosted/null focusing or all renormalization coefficients",
            "curvature blow-up, C2 inextendibility, stability or fundamental EFT validity at the endpoint",
            "an incompleteness theorem without the explicit geometric contraction history and conditional future bounds",
            "completion of original P8a or P8",
        ],
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("the photon history-focusing certificate differs from exact read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8(a) A.17: photon-QSEI history-based conditional incompleteness replay passed; physical cosmological application and P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
