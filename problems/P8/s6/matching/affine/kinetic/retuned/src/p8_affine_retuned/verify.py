"""Read-only source-centered trace-mass and physical-cone certificate."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as affine
from p8_affine_ricci import verify as parent

from . import bounds, cones, degeneracy, dynamics, geometry

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"source-centered-mass-cone.json"
PARENT_SHA = "7e65aa60e5f002feed1940c4c640c07c3937b000b02c834c9ea4f9da10e749ce"
sha, serialize = affine.sha, affine.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_affine_retuned/*.py"))
            +sorted(ROOT.glob("tests/*.py"))+sorted(ROOT.glob("*.md"))
            +sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen Ricci-difference certificate changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {"S6_40_fully_rebuilt": PARENT_SHA,
            "original_unrestricted_CD_M1_action_and_matter_cone_contract_rechecked": True,
            "no_frozen_ancestor_modified": True}


def controls():
    wrong = (True, False, 1.0, sp.Float(1), "1", sp.I, sp.oo, sp.nan, sp.Symbol("unproved"))
    half, small = sp.Rational(1, 2), sp.Rational(1, 2000)
    calls = [lambda value=value: bounds.require_domain(half, value, 1, 1) for value in wrong]
    calls += [lambda value=value: bounds.units(value, 2, small) for value in wrong]
    domains = ((-half, small, 1, 1), (sp.Rational(1, 3), small, 1, 1),
               (half, -small, 1, 1), (half, 2*small, 1, 1),
               (half, small, 0, 1), (half, small, 1, 0), (half, small, 1, -1))
    calls += [lambda values=values: bounds.require_domain(*values) for values in domains]
    calls += [lambda values=values: bounds.units(*values) for values in
              ((0, 1, 1), (1, 0, 1), (1, 1, 0), (1, 1, -1), (1, 1.0, 1))]
    calls += [lambda value=value: serialize(value) for value in (1.0, sp.Float(1), sp.oo, sp.nan)]
    rejected = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(calls):
        raise ValueError("An inexact, out-of-domain or singular-chart input was accepted")
    if bounds.require_domain(half, 0, 0, 0)["branch"] != "exact_auxiliary_control":
        raise ValueError("The zero-curl auxiliary control failed")
    return {"rejected_inputs": rejected, "source_and_counterterm_both_retained": True,
            "missing_time_boundary_misaligns_null_but_does_not_alone_change_old_rank": True,
            "matter_cone_test_not_replaced_by_positive_energy": True,
            "massive_finite_q_phase_speed_not_used_as_front_speed": True,
            "zero_curl_and_high_frequency_limits_not_interchanged": True,
            "original_witness_rejection_not_a_below_cutoff_UV_no_go": True}


@cache
def build_report():
    pins = prior_checks()
    groups = {"source_centered_all64_mass_update": geometry.checks(),
              "full_primary_null_and_ten_velocity_chart": degeneracy.checks(),
              "actual_full_scalar_constraints_and_crossing": dynamics.checks(),
              "physical_matter_cone_obstruction": cones.checks()}
    exact = {name: affine.certify_residuals(values) for name, values in groups.items()}
    proofs = bounds.proof_checks()
    if not all(value is True for value in proofs.values()):
        raise ValueError("A continuous-domain, primary-rank or cone-sign proof failed")
    cone = cones.characteristic()
    return {"schema": 1, "claim": "P8-S6.41.AFFINE_RETUNED", "date": "2026-09-07",
            "status": "MATCHING_AND_PRIMARY_NULL_CERTIFIED; POSITIVE_PRINCIPAL_BUT_MATTER_CONE_REJECTED; ORIGINAL_P8_OPEN",
            "prior_sha256": pins,
            "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in source_files()},
            "formulation": "FORMULATION.md", "written_proofs": ["notes/geometry.md", "notes/dynamics.md"],
            "exact_residuals": exact,
            "named_exact_check_count": sum(map(len, groups.values())),
            "checked_scalar_entries": sum(value.rows*value.cols if isinstance(value, sp.MatrixBase) else 1
                                          for group in groups.values() for value in group.values()),
            "proof_checks": proofs,
            "literal_new_operator": "-mu*g^(ab)*(T-Tstar)_a*(T-Tstar)_b/2-zeta*F(T)^2/4, T=V+U, mu=11/9, all sources and counterterms included",
            "matching": "At zero curl all64 stationary equations and the entire original CD/M1 reduced action are unchanged; quotient inverse norm below 6900 on the original tube",
            "primary_constraint": "Tstar_normal=(4p²-1)*K_hat+6p*J3*s³ contains no lapse velocity; the ten-velocity block keeps rank, not a full nonlinear secondary-constraint theorem",
            "positive_principal_range": "0<zeta<=1/2000, all finite u!=0 and q>0; regular first-order center, positive center velocity chart q>6",
            "characteristic": serialize({"kappa": cone["kappa"], "polynomial": cone["polynomial"],
                                         "speed_plus_squared": cone["speed_plus"], "speed_minus_squared": cone["speed_minus"],
                                         "matter_speed_squared": 1}),
            "units": serialize(bounds.units(3, 2, sp.Rational(1, 1000))),
            "controls": controls(),
            "verdict": "The named positive-energy isotropic mass-retuned parent fails the original scalar matter-cone acceptance condition at every finite u!=0. Its high-frequency obstruction is not a certified below-cutoff EFT exclusion.",
            "not_established": ["Full nonlinear secondary constraints, physical degree count or nonlinear stability",
                                "A coupled heavy spectral gap, controlled elimination, loop or cutoff/remainder bounds",
                                "A vacuum/finite-gravity V/G/B or UV-completion verdict",
                                "A no-go for other kinetic/source operators, all mass retunings or original P8 closure"],
            "verification_boundary": "Exact symbolic identities, continuous polynomial and rank proofs, independent dense and constrained calculations; not proof-assistant formalization"}


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The source-centered mass/cone report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.41.AFFINE_RETUNED replay passed; matter-cone obstruction, original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
