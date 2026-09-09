"""Read-only ordinary Proca C5 stress and declared fixed-profile replacement."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as original
from p8_constant_proca import verify as parent

from . import audit, estimates, profiles, readouts, subtraction

ROOT=Path(__file__).resolve().parents[2]
REPORT=ROOT/"certificates"/"fixed-ordinary-Proca-stress.json"
PARENT_SHA="73d13e9c26bc66ee8ed93453249e7944a0e741579a8a2942fe4f95631a7a358d"
sha,serialize=original.sha,original.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_proca_stress/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT)!=PARENT_SHA:
        raise ValueError("The frozen constant-mass physical Proca report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()),parent.build_report())
    return {"S6_81_fully_rebuilt":PARENT_SHA,
            "old_fixed_profile_nontransfer_boundary_preserved_before_explicit_replacement":True}


@cache
def build_report():
    prior=prior_checks()
    identities,gates=audit.residuals(),audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("An ordinary Proca stress/profile gate failed")
    exact=original.certify_residuals(identities)
    kinds=("transverse","longitudinal")
    return {"schema":1,"claim":"P8-S6.82.FIXED_ORDINARY_PROCA_STRESS","date":"2026-09-08",
            "status":"NEW_FIXED_ORDINARY_PROCA_C5_STRESS_AND_EXPLICIT_PROFILE_REPLACEMENT_CANCEL_RETAINED_GAUSSIAN_ONE_POINT; FULL_RESPONSE_MATCHING_AND_ORIGINAL_P8_OPEN",
            "prior_sha256":prior,"source_sha256":{str(p.relative_to(ROOT)):sha(p) for p in source_files()},
            "formulation":"FORMULATION.md",
            "written_proofs":["notes/readouts.md","notes/subtraction.md","notes/bounds.md",
                              "notes/profiles.md","notes/scope.md"],
            "exact_residuals":exact,"named_exact_check_count":len(identities),
            "checked_scalar_entries":sum(v.rows*v.cols if isinstance(v,sp.MatrixBase) else 1 for v in identities.values()),
            "proof_checks":gates,
            "actual_new_energy_rows_and_continuous_envelopes":serialize({
                kind:{j:{"row":readouts.rows(kind,j),"bound":readouts.row_bounds(kind,j)}
                      for j in range(6)} for kind in kinds}),
            "actual_differentiated_subtraction_and_integrable_reference_tails":serialize({
                kind:{j:{"adiabatic":subtraction.adiabatic_rows(kind,j),"tail":subtraction.reference_tail(kind,j)}
                      for j in range(6)} for kind in kinds}),
            "actual_lower_reference_failure_control":serialize(subtraction.reference_tail("longitudinal",4,2)),
            "actual_new_local_energy_derivative_envelopes":serialize({
                j:estimates.local_derivatives(j) for j in range(6)}),
            "continuous_new_stress_and_energy_change_bounds":serialize(estimates.physical_bounds(estimates.SCALE)),
            "actual_fixed_continuum_integrand_definition":serialize(profiles.integrands()),
            "explicit_new_fixed_profile_and_old_profile_difference":serialize(profiles.action()),
            "global_mixed_profile_and_local_contact_bounds":serialize(profiles.profile_bounds(estimates.SCALE)),
            "controls":audit.controls(),
            "verdict":"The new constant-mass ordinary Proca energy, differentiated subtraction and finite local matching are recomputed in the unchanged actual continuum state. Uniform C5 bounds for energy, unchanged pressure and the energy change are below 10^-770 at fixed m=1000 and L=10^400. An explicitly replaced, subsequently frozen scalar profile cancels the retained Gaussian one-point energy, pressure and conserved clock source. Both 21-entry global mixed profile-derivative triangles and the stated local contact bounds pass. This does not transfer the old quantum response or prove a full quantum solution.",
            "not_established":["The new Gaussian functional metric response, its spectrum or all-frequency quantum stability",
                "Mixed loops, omitted operators, heavy-threshold errors or an interacting Wilsonian cutoff",
                "Canonical higher-vertex bounds for this added profile correction from its coordinate C5 norm alone",
                "Full scalar-neighborhood cone, nonlinear spatial well-posedness or an actual global quantum bounce",
                "Finite Wilson matching, common-parent V/G/B, UV completion/exclusion or original P8 closure"],
            "verification_boundary":"Exact readout, ODE, differentiated subtraction, Ward and fixed-profile identities; continuous rational envelopes, full radial integrals and written proofs. Parent state/evolution/counterterm chain is natively rebuilt. Retained Gaussian one-point only, not functional response or proof-assistant formalization."}


def validate_report(expected,actual):
    if expected!=actual:
        raise ValueError("The ordinary Proca stress/profile report differs from read-only replay")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    actual=build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()),actual)
        print("P8 S6.82.FIXED_ORDINARY_PROCA_STRESS replay passed; full response, matching and original P8 OPEN")
    else:
        print(json.dumps(actual,indent=2,sort_keys=True))


if __name__=="__main__":
    main()
