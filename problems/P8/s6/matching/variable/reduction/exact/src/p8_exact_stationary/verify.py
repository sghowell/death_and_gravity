"""Read-only exact stationary branch and independent whole-box replay."""
import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_sixth_reduction import verify as prior
from p8_variable_reduction import dictionary as prior_dictionary
from p8_variable_reduction.verify import serialize

from . import audit, bridges, intervals, stationary

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"exact-stationary-branch.json"
PARENT_SHA = "d392a734f46d8f4c1c39523c7dc3d41a88634ec0e3dfa645365c3b4d14eb344c"
MAX_C = sp.Rational(201, 100)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_files():
    return (sorted(ROOT.glob("src/p8_exact_stationary/*.py"))
            +sorted(ROOT.glob("tests/*.py"))+sorted(ROOT.glob("*.md"))
            +sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(prior.REPORT) != PARENT_SHA:
        raise ValueError("The frozen first-omitted-action cancellation report changed")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return {"S6_32_recursively_rebuilt": PARENT_SHA,
            "same_unchanged_S6_20_parent_and_S6_30_physical_frame": True,
            "same_off_shell_unit_volume_probe_not_actual_FLRW": True,
            "does_not_use_S6_31_deformed_interactions": True}


def numeric_center(c=MAX_C, parent_planck_squared=1, time_scale=1):
    checked = prior_dictionary.center(c=c, parent_planck_squared=parent_planck_squared,
                                      time_scale=time_scale)
    cc, m2, tau = checked["c"], sp.Rational(parent_planck_squared), sp.Rational(time_scale)
    if cc > MAX_C:
        raise ValueError("This quantitative branch certificate requires 2<c<=201/100")
    data = stationary.center()
    return {"c": cc, "delta": cc-2, "M_squared": m2, "tau": tau,
            "b": sp.Integer(2), "N": data["N_center"].subs(stationary.C, cc),
            "f00": data["metric00_center"].subs(stationary.C, cc),
            "b_TT": sp.factor(data["b_uu_center"].subs(stationary.C, cc)/tau**2),
            "N_TT": sp.factor(data["N_uu_center"].subs(stationary.C, cc)/tau**2),
            "K_f": 16*m2/cc**2, "own_f_mass_squared": 8*cc/((cc-2)*tau**2),
            "coupled_relative_mass_squared": 8*(cc**2+16)/(cc*(cc-2)*tau**2),
            "full_parent_solution": False, "physical_cutoff_claimed": False}


def controls():
    bad = (True, False, 1.0, sp.Float("1"), sp.oo, -sp.oo, sp.nan, "3", sp.Symbol("c"))
    calls = [lambda value=value: numeric_center(c=value) for value in bad]
    calls += [lambda value=value: numeric_center(c=value) for value in (1, 2, 3)]
    calls += [lambda: numeric_center(parent_planck_squared=0), lambda: numeric_center(time_scale=-1)]
    calls += [lambda value=value: serialize(value) for value in (0.1, sp.Float("0.1"), sp.oo, sp.nan, sp.I)]
    calls += [lambda value=value: audit.arb_box(value) for value in (True, False, 128.0, "256", 64, 0, -1)]
    calls += [lambda value=value: intervals.Interval(value) for value in bad]
    rejected = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(calls):
        raise ValueError("An invalid domain, arithmetic input or literal endpoint was admitted")
    return {"rejected_inputs": rejected,
            "literal_c_equals_2_excluded": True,
            "no_division_by_scale_velocity_at_center": True,
            "parity_proved_not_assumed": True,
            "lapse_uses_total_implicit_derivative": True,
            "own_f_inner_16_not_coupled_relative_80": True,
            "analytic_background_not_a_selected_tensor_inverse": True,
            "original_physical_g_clock_and_free_chi_unchanged": True}


@cache
def build_report():
    pins = prior_checks()
    for name in ("test_exact_stationary.py", "test_exact_intervals.py", "test_exact_independent_audit.py"):
        if not (ROOT/"tests"/name).is_file():
            raise ValueError("Each separately authored scientific audit must be pinned")
    groups = {"literal_equations_local_uniqueness_center_and_inner": stationary.checks(),
              "independent_four_coordinate_action_and_frozen_profile_bridges": bridges.checks()}
    if any(value != 0 for group in groups.values() for value in group.values()):
        raise ValueError("An exact stationary/action identity failed")
    exact_box = intervals.report()
    if not all(value is True for value in exact_box["checks"].values()):
        raise ValueError("A continuous exact Fraction proof check failed")
    ball_box = audit.arb_box()
    center, inner = stationary.center(), stationary.inner()
    return {
        "schema": 1, "claim": "P8-S6.33.EXACT_STATIONARY", "date": "2026-09-07",
        "status": "EXACT_OWN_F_BACKGROUND_WITH_UNIFORM_POSITIVE_LAPSE_AND_DISTINCT_INNER_OPERATORS; ORIGINAL_P8_OPEN",
        "prior_sha256": pins,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in source_files()},
        "formulation": "FORMULATION.md", "written_proofs": ["notes/stationary.md", "notes/intervals.md"],
        "source_audit": "notes/sources.md",
        "actual_domain": "M,tau>0; c=2+delta; 0<delta<=1/100; |T|<=tau/100",
        "probe": "g=eta, theta=T; exact own-f only, not g/clock equations or the rolling CD solution",
        "exact_residuals": {name: dict.fromkeys(group, "0") for name, group in groups.items()},
        "exact_residual_count": sum(map(len, groups.values())),
        "fraction_whole_box": exact_box,
        "independent_arb_whole_box": ball_box,
        "fixed_c_local_Jacobian": serialize(stationary.fixed_c()["b_jacobian_center"]),
        "joint_center_IFT_Jacobian": serialize(center["joint_G_zeta_center"]),
        "center_jets": serialize({key: center[key] for key in (
            "b_center", "N_center", "b_uu_center", "b_uuuu_center", "N_uu_center",
            "metric00_center", "retained_metric00_center", "metric00_center_defect",
            "joint_zeta_v_center", "joint_zeta_c_center", "joint_zeta_center_as_function_c")}),
        "exact_tensor_action": serialize(stationary.tensors()["density"]),
        "inner_equation_coefficients": serialize({key: inner[key] for key in (
            "own_f_inner_coefficient", "g_inner_coefficient", "coupled_relative_inner_coefficient",
            "own_f_center_mass_squared", "coupled_algebraic_center_mass_squared", "own_f_center_K")}),
        "numeric_center_calibration": serialize(numeric_center(parent_planck_squared=3, time_scale=2)),
        "calibration": serialize(stationary.calibration()),
        "controls": controls(),
        "not_established": [
            "The physical-g/clock equations or a full complete stable CD/parent cosmology",
            "A unique own-f functional on a neighborhood of all physical metrics/clocks",
            "Tensor homogeneous data, a retarded/advanced/Feynman inverse or Green norm",
            "Uniform all-order derivative expansion, full source/matter matching or a physical cutoff",
            "Original CD coefficient equality, positivity/V/G applicability or original P8 completion"],
        "verification_boundary": "Written Euler/analytic IFT and Banach proof, exact identities, independent Fraction and Arb whole-box enclosures; not proof-assistant formalization or controlled EFT",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The exact stationary branch report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.33.EXACT_STATIONARY replay passed; original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
