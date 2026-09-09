"""Read-only retuned covariant margin and fresh classical-domain certificate."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as certificate
from p8_proca_spectator_repair import verify as parent

from . import audit, central, domain, flow, model, original

ROOT=Path(__file__).resolve().parents[2]
REPORT=ROOT/"certificates"/"retuned-covariant-margin-domains.json"
PARENT_SHA="13fed85f6b8c17191bb4bf2b9618da2c14fee818c7d233134fa738c1ec2cc4a5"
sha,serialize=certificate.sha,certificate.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_proca_retuned_margin/*.py"))+sorted(ROOT.glob("tests/*.py"))
        +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT)!=PARENT_SHA:
        raise ValueError("The frozen finite principal spectator-repair report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()),parent.build_report())
    return {"S6_95_fully_rebuilt":PARENT_SHA,
        "new_named_covariant_light_action_not_unchanged_spectator_addition":True,
        "new_actual_nearby_solution_has_its_own_constraint_and_conserved_charge":True,
        "no_old_quantum_state_profile_response_or_UV_matching_transferred":True}


@cache
def build_report():
    prior=prior_checks()
    identities,gates=audit.residuals(),audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A retuned-action, actual-solution or central-family gate failed")
    exact=certificate.certify_residuals(identities)
    signatures={name:{"numerator":model.old.polynomial_signature(value.numer),
        "denominator":model.old.polynomial_signature(value.denom)} for name,value in model.system().items()}
    boxes=domain.enclosures()
    return {"schema":1,"claim":"P8-S6.96.RETUNED_COVARIANT_MARGIN_DOMAINS","date":"2026-09-09",
        "status":"NEW_TOTAL_MARGIN_ONE_OVER_TWO_HUNDRED_PRESERVES_GLOBAL_CLASSICAL_CLOCK_AND_HAS_FRESH_ACTUAL_NEARBY_SUBLUMINAL_BOUNCE_PLUS_CENTRAL_POSITIVE_MATTER_CONTINUUM; FULL_RETARDED_INTERACTING_QUANTUM_V_G_B_AND_ORIGINAL_P8_OPEN",
        "prior_sha256":prior,"source_sha256":{str(p.relative_to(ROOT)):sha(p) for p in source_files()},
        "formulation":"FORMULATION.md",
        "written_proofs":["notes/action.md","notes/evolution.md","notes/enclosures.md","notes/central.md","notes/scope.md"],
        "exact_residuals":exact,"named_exact_check_count":len(identities),
        "checked_scalar_entries":sum(value.rows*value.cols if isinstance(value,sp.MatrixBase) else 1
                                     for value in identities.values()),
        "proof_checks":gates,
        "literal_new_covariant_action_and_fixed_phase_lapse_change":serialize(model.action()),
        "new_exact_rational_constraint_and_flow_signatures":signatures,
        "fresh_whole_real_and_complex_coordinate_enclosures":serialize({"real_boxes":boxes,
            "complex_modulus_annuli":{name:domain.complex_annulus(row) for name,row in boxes.items()},
            "actual_physical_Hubble_partial_quotient_boxes":domain.physical_hubble_partials()}),
        "complete_physical_proper_time_Hubble_acceleration":serialize(domain.acceleration()),
        "new_Picard_solution_conserved_charge_and_principal_health":serialize(
            {"Picard":flow.picard(),"initial_data":flow.initial(),"physical_health":flow.health()}),
        "preserved_global_clock_and_literal_deformation_budgets":serialize(original.data()),
        "central_positive_matter_continuum_and_weaker_margin_control":serialize(central.theorem()),
        "controls":audit.controls(),
        "verdict":"A separately named covariant action raises the existing total clock margin from 10^-6 to 1/200, adding only 4999/10^6 times (X-1)^2/h^2. Its first clock-trajectory variations vanish, and the original global classical bounce and two-chart subluminal principal result are explicitly bridged. A newly derived rational constraint/flow has a unique holomorphic solution on |u|<=10^-7 from N0=1+10^-6, z0=0 and R0=1, with its own nonzero conserved matter charge. Whole-box bounds give 0.986<c_clock^2<0.988, positive physical principal forms and 3.9997<dHphysical/dtau<4.001. Exact polynomial and Bernstein inequalities cover every positive-matter central datum in 0.9<=X<=1.1, giving 1/2<c_clock^2<1-1/4000 and bounce acceleration above three. A weaker 10^-4-margin fast-clock bounce is retained as a negative control. No old action, state, profile or matching claim is changed. Full retarded support, interacting/quantum control and UV completion are not established.",
        "not_established":["A complete retarded-support or microcausality theorem for the changed finite-frequency system",
            "A common quantitative time radius for the whole central family, whole off-center phase tube or global nearby continuation",
            "Nonlinear stability, canonical interaction/error bounds or an interacting temporal/spatial EFT validity band",
            "A new quantum state, Hadamard stress, profile cancellation, loops or self-consistent semiclassical backreaction",
            "Physical vacuum, amplitudes, dispersion and common-parent V/G/B matching, realistic-field P8(a) or original P8 closure"],
        "verification_boundary":"Native literal covariant scalar and volume/Hamiltonian derivatives, fixed-phase rational lapse constraint and fresh flow, exact conserved charge, full real/complex polynomial boxes, Cauchy/Picard inequalities, complete uncancelled physical-Hubble quotient derivatives, actual principal matrices through the inherited generic Euler-first derivation, independent central time-jet bridge, monotone central polynomials and positive Bernstein reconstructions. Written existence, chart-coverage and continuum proofs are pinned. Every parent is natively rebuilt, every source/proof/test byte pinned and every report field mutation rejected. No scientific-library patch, point-scan proof, reused old charge, inferred interacting matching or frozen action mutation. Not proof-assistant formalization."}


def validate_report(expected,actual):
    if expected!=actual:
        raise ValueError("The retuned covariant-margin report differs from read-only replay")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    actual=build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()),actual)
        print("P8 S6.96.RETUNED_COVARIANT_MARGIN_DOMAINS replay passed; full retarded, interacting quantum V/G/B and original P8 OPEN")
    else:
        print(json.dumps(actual,indent=2,sort_keys=True))


if __name__=="__main__":
    main()
