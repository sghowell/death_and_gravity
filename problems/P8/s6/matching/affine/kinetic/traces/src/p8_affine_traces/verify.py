"""Read-only rank-complete constant two-trace curl certificate."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine_kinetic import verify as parent

from . import bridges, dynamics, geometry, rank_two

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"constant-two-trace-family.json"
PARENT_SHA = "32af5490520b8d9a40c4ad76e1f51560474f431f77b7c5c1772b414f6c0b438e"
sha, serialize = parent.sha, parent.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_affine_traces/*.py"))
            +sorted(ROOT.glob("tests/*.py"))+sorted(ROOT.glob("*.md"))
            +sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen specified-kinetic-screens report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {"S6_38_fully_rebuilt": PARENT_SHA,
            "S6_37_exact_action_and_original_CD_M1_rechecked": True,
            "original_target_and_adopted_contract_unchanged": True,
            "no_frozen_ancestor_modified": True}


def controls():
    wrong = (True, False, 1.0, sp.Float(1), "1", sp.I, sp.oo, -sp.oo, sp.nan, sp.Symbol("unproved"))
    calls = [lambda value=value: rank_two.classify(value, 1, 0) for value in wrong]
    calls += [lambda value=value: dynamics.require_domain(value, 0, 1, 1, 8) for value in wrong]
    calls += [lambda values=values: dynamics.require_domain(*values) for values in
              ((1, 0, 1, 0, 8), (1, 0, 1, 1, 0), (1, 0, 1, 1, -1),
               (1, 1, 1, 1, sp.Rational(2, 9)), (7+4*sp.sqrt(3), 1, 1, 0, 1))]
    calls += [lambda values=values: geometry.require_parameters(*values)
              for values in ((0, 0), (1.0, 0), (1, sp.oo))]
    calls += [lambda value=value: serialize(value) for value in (1.0, sp.Float(1), sp.oo, sp.nan)]
    rejected = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(calls):
        raise ValueError("A forbidden inexact, singular-chart or domain input was accepted")
    unchanged = rank_two.classify(0, 0, 0)["rank"] == 0
    null_branches = [dynamics.require_domain(root, 1, sign, sp.Rational(1, 10), 1)["branch"]
                     for root in (7-4*sp.sqrt(3), 7+4*sp.sqrt(3)) for sign in (-1, 1)]
    if not unchanged or any(value != "null_schur_ostrogradsky" for value in null_branches):
        raise ValueError("The unchanged or exceptional-rank controls failed")
    return {"rejected_inputs": rejected, "both_null_rays_and_both_coupling_signs": True,
            "zero_matrix_is_unchanged_auxiliary": True,
            "null_Schur_not_a_Proca_or_transverse_limit": True,
            "rolling_U_source_not_minus_V_source": True,
            "all_matrix_ranks_including_off_diagonal_forms": True,
            "no_frequency_cutoff_or_general_affine_no_go": True}


@cache
def build_report():
    pins = prior_checks()
    groups = {"literal_eight_trace_geometry_and_null_equations": geometry.checks(),
              "independent_complete_rank_one_dynamics": dynamics.checks(),
              "rank_two_longitudinal_inertia": rank_two.checks(),
              "full_connection_and_physical_constraint_interfaces": bridges.checks()}
    residuals = {name: parent.parent.certify_residuals(values) for name, values in groups.items()}
    proofs = {"rank_one_continuous_sign_and_Ostrogradsky": dynamics.proof_checks(),
              "full_quotient_and_original_background": bridges.proof_checks()}
    if not all(value is True for group in proofs.values() for value in group.values()):
        raise ValueError("A continuous domain or rank proof premise failed")
    return {
        "schema": 1, "claim": "P8-S6.39.AFFINE_TRACES", "date": "2026-09-07",
        "status": "ALL_NONZERO_CONSTANT_TWO_TRACE_CURL_FORMS_FAIL_LITERAL_PARENT_HEALTH; ORIGINAL_P8_OPEN",
        "prior_sha256": pins,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": ["notes/geometry.md", "notes/dynamics.md", "notes/rank-two.md", "notes/interfaces.md"],
        "exact_residuals": residuals,
        "named_exact_check_count": sum(map(len, groups.values())),
        "checked_scalar_entries": sum(value.rows*value.cols if isinstance(value, sp.MatrixBase) else 1
                                      for group in groups.values() for value in group.values()),
        "proof_checks": proofs,
        "family": "The unchanged unrestricted S6.37 action plus -Z_IJ F(CI).F(CJ)/4, C=(V,U), Z real constant symmetric 2x2",
        "geometry": serialize(geometry.calibration()),
        "rank_one": serialize(dynamics.calibration()),
        "rank_two": serialize(rank_two.calibration()),
        "interfaces_and_units": serialize(bridges.calibration()),
        "verdict": "Rank zero is the unchanged auxiliary action. Every nonzero rank-one or rank-two matrix has a physical negative kinetic direction or a nondegenerate unconstrained Ostrogradsky mode on the actual rolling solution.",
        "controls": controls(),
        "not_established": ["An unhealthy low-frequency EFT domain or quantitative frequency cutoff",
                            "A nonlinear open-tube inverse for the new kinetic family, especially its null-Schur rays",
                            "A no-go for field-dependent coefficients, other curvature/distortion contractions, mass tunings or additional fields",
                            "New vacuum, finite-gravity, loop, threshold or controlled matching V/G/B estimates",
                            "UV completion, original P8 closure or a new requirement on the completed scoped photon objective"],
        "verification_boundary": "Exact symbolic identities and written whole-domain rank, constraint and inertia proofs; not proof-assistant formalization",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The constant two-trace family report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.39.AFFINE_TRACES replay passed; constant two-trace family closed, original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
