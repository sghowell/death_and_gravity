"""Read-only source alignment, physical cones and canonical frequency replay."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as affine
from p8_affine_retuned import verify as parent
from p8_affine_retuned.bounds import units

from . import alignment, dynamics, modes, source_bounds

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"source-aligned-canonical-modes.json"
PARENT_SHA = "37c8c3035fd3e8789889e2532a8817680cf4992ffa0c1152a761da554b1f7982"
sha, serialize = affine.sha, affine.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_affine_aligned/*.py"))
            +sorted(ROOT.glob("tests/*.py"))+sorted(ROOT.glob("*.md"))
            +sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen source-centered mass/cone certificate changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {"S6_41_fully_rebuilt": PARENT_SHA,
            "all64_source_mass_and_original_CD_M1_matter_frame_rechecked": True,
            "no_frozen_ancestor_modified": True}


def controls():
    wrong = (True, False, 1.0, sp.Float(1), "1", sp.I, sp.oo, sp.nan, sp.Symbol("unproved"))
    small = sp.Rational(1, 2000)
    calls = [lambda value=value: modes.frequency_bounds(value, 0, 1) for value in wrong]
    calls += [lambda value=value: source_bounds.bound_values(value, 0, 0, 0) for value in wrong]
    calls += [lambda values=values: modes.frequency_bounds(*values) for values in
              ((0, 0, 1), (-small, 0, 1), (2*small, 0, 1), (small, 0, -1),
               (small, 0.0, 1), (small, 0, 1.0), (small, 0, sp.Symbol("k2", positive=True)))]
    calls += [lambda values=values: source_bounds.bound_values(*values) for values in
              ((-1, 0, 0, 0), (sp.Rational(11, 100), 0, 0, 0),
               (0, -1, 0, 0), (0, 0, -1, 0), (0, 0, 0, -1))]
    calls += [lambda value=value: serialize(value) for value in (1.0, sp.Float(1), sp.oo, sp.nan)]
    rejected = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(calls):
        raise ValueError("An inexact, unproved or out-of-domain input was accepted")
    if modes.frequency_bounds(small, 0, 0)["chart"] != "three_homogeneous_coordinate_vectors":
        raise ValueError("The homogeneous vector chart was lost")
    return {"rejected_inputs": rejected,
            "zero_curl_separate_auxiliary_rank_control": True,
            "zero_momentum_not_inverted_in_longitudinal_chart": True,
            "wrong_sign_and_fixed_physical_momentum_controls": True,
            "second_order_source_and_cubic_coupling_nonzero": True,
            "frequency_floor_not_stationary_gap_or_cutoff": True,
            "source_norm_not_retarded_or_action_remainder": True}


@cache
def build_report():
    pins = prior_checks()
    groups = {"complete_local_source_alignment": alignment.checks(),
              "full_quadratic_constraints_and_physical_cones": dynamics.checks(),
              "actual_canonical_vector_equations": modes.checks(),
              "nonlinear_ODE_source_and_curl_bounds": source_bounds.checks()}
    exact = {name: affine.certify_residuals(values) for name, values in groups.items()}
    proofs = {"modes_"+name: value for name, value in modes.proof_checks().items()}
    proofs.update({"source_"+name: value for name, value in source_bounds.proof_checks().items()})
    proofs["regular_first_order_center_without_Theta_inverse"] = dynamics.first_order()["no_Theta_inverse"]
    if not all(value is True for value in proofs.values()):
        raise ValueError("A canonical frequency, continuous source or crossing proof failed")
    conversion = units(3, 2, sp.Rational(1, 1000))
    return {"schema": 1, "claim": "P8-S6.42.AFFINE_ALIGNED", "date": "2026-09-07",
            "status": "SOURCE_ALIGNMENT_QUADRATIC_MATTER_CONES_AND_CANONICAL_FLOOR_CERTIFIED; NONLINEAR_SOURCE_BOUNDED; ORIGINAL_P8_OPEN",
            "prior_sha256": pins,
            "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in source_files()},
            "formulation": "FORMULATION.md", "written_proofs": ["notes/source.md", "notes/dynamics.md", "notes/modes.md"],
            "exact_residuals": exact,
            "named_exact_check_count": sum(map(len, groups.values())),
            "checked_scalar_entries": sum(value.rows*value.cols if isinstance(value, sp.MatrixBase) else 1
                                          for group in groups.values() for value in group.values()),
            "proof_checks": proofs,
            "literal_new_operator": "S6.41 centered mass unchanged; replace F(T) by F(W), W=T-B(phi,x)dphi, B=-15H(1+x)/(8h); zeta=zeta_physical/(M² tau²)",
            "nonlinear_source": "Sstar_normal=(4p²-1)*(K_hat-3H*s)+(3/2)*s*Q_lower; full background and first variation zero, second variation nonzero; no lapse velocity",
            "quadratic_action": "Exactly the unchanged physical CD/M1 action plus minimally coupled positive three-polarization Proca; full lapse/shift/temporal constraints retained",
            "matter_cones": "All three scalar and two transverse-vector principal characteristics exactly c²=1 in the original matter metric; original tensor block unchanged",
            "canonical_modes": serialize({"alpha": modes.canonical()["alpha"], "beta": modes.canonical()["beta"],
                                          "zeta_max": modes.ZETA_MAX, "uniform_squared_floor": "q+1/zeta-15>=q+1985"}),
            "source_bounds": "Original closed x tube: |Q|<(3/2)*(1+x)²/h², |Q_x|<3*|1+x|/h² away from zero distance, with non-strict equality at the clock; normal<=eta*|Delta_K|+(5/2)eta², electric<=eta*epsilon_K+(3|Delta_K|+10eta)*epsilon_s",
            "nonunit_example": serialize({"conversion": conversion,
                                          "frequency": modes.frequency_bounds(conversion["normalized_zeta"], sp.Rational(1, 2), sp.Rational(4375, 256)),
                                          "physical_frequency_squared_lower": 2998}),
            "controls": controls(),
            "verdict": "This named source-aligned candidate passes the original quadratic rolling matter-cone gate and has an actual canonical vector frequency floor. Nonlinear source bounds alone do not establish heavy elimination or UV admissibility.",
            "not_established": ["Full nonlinear secondary constraints, degree count or nonlinear stability",
                                "A controlled retarded heavy inverse, omitted-action or loop remainder, interacting cutoff or stationary mass gap",
                                "The adopted vacuum, finite-gravity and common-parent V/G/B UV conditions",
                                "Original P8 closure or a general metric-affine classification"],
            "verification_boundary": "Exact symbolic identities, continuous bounds and independent internal constrained/Euler calculations with written proofs; not proof-assistant formalization"}


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The aligned-source canonical-mode report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.42.AFFINE_ALIGNED replay passed; quadratic cones and canonical floor, original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
