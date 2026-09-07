"""Read-only complete CD/M1 auxiliary-affine dictionary certificate."""
import argparse
import hashlib
import json
from fractions import Fraction
from functools import cache
from pathlib import Path

import sympy as sp
from p8_m1 import verify as prior

from . import audit, bridges, connection, dictionary, lower

ROOT = Path(__file__).resolve().parents[2]
P8 = ROOT.parents[2]
REPORT = ROOT/"certificates"/"cd-auxiliary-lift.json"
PARENT_SHA = "6b0dcb44c0849050912f22546fc7c2ef4c55f148ebca32184b09b4abecf7ca7c"
PINS = {
    "FORMULATION.md": "73dd18f7a4b56a0205f9fbc6a4b09b213a3c280b3f0274414a255635b04e2b09",
    "s6/FORMULATION.md": "d00df35d5821b05baa61e897725a22ed1ac21d0536b74bdaece3ab0093215901",
    "certificates/witness-CD_matter.json": "caf8c8e688a7565b9d00f921c099a28da00f97522ed26ad182a8227eb80cd4dd",
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_files():
    return (sorted(ROOT.glob("src/p8_affine/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


def serialize(value):
    if isinstance(value, dict):
        return {str(key): serialize(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [serialize(item) for item in value]
    if isinstance(value, sp.MatrixBase):
        return serialize(value.tolist())
    if isinstance(value, bool) or value is None or isinstance(value, (str, int)):
        return value
    if value is sp.true or value is sp.false:
        return bool(value)
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, float):
        raise TypeError("Binary floats are not exact certificate entries")
    if isinstance(value, sp.Basic):
        if value.has(sp.Float, sp.oo, -sp.oo, sp.zoo, sp.nan):
            raise ValueError("Inexact or nonfinite symbolic certificate entry")
        return str(value)
    raise TypeError(f"Unsupported certificate entry: {type(value).__name__}")


@cache
def prior_checks():
    for relative, expected in PINS.items():
        if sha(P8/relative) != expected:
            raise ValueError(f"A frozen target or adopted contract changed: {relative}")
    if sha(prior.REPORT) != PARENT_SHA:
        raise ValueError("The frozen CD/M1 nonlinear/principal checkpoint changed")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return {**PINS, "S5_5_rebuilt": PARENT_SHA,
            "earlier_source_manifests_rechecked": True,
            "target_and_matter_not_reconstructed_or_reassigned": True,
            "bimetric_probe_not_a_prerequisite_of_this_new_candidate": True}


def certify_residuals(values):
    """Keep matrix identities distinct from the count of their scalar entries."""
    result = {}
    for name, value in values.items():
        if isinstance(value, sp.MatrixBase):
            if any(entry != 0 for entry in value):
                raise ValueError(f"A matrix identity failed: {name}")
            result[name] = {"shape": [value.rows, value.cols], "all_entries_exactly_zero": True}
        else:
            if value != 0:
                raise ValueError(f"An exact identity failed: {name}")
            result[name] = "0"
    return result


def controls():
    bad = (True, False, sp.true, sp.false, 1.0, sp.Float(1), "1", sp.I, sp.oo, -sp.oo, sp.zoo, sp.nan)
    calls = [lambda value=value: dictionary.tube_point(physical_X=value) for value in bad]
    calls += [lambda value=value: dictionary.tube_point(physical_X=value) for value in (0, -1, sp.Rational(1, 2), sp.Rational(6, 5))]
    calls += [lambda value=value: dictionary.tube_point(mass_squared=value) for value in (0, -1, True, 1.0)]
    calls += [lambda value=value: dictionary.tube_point(time_scale=value) for value in (0, -1, True, 1.0)]
    calls += [lambda value=value: audit.domain_at(0, value) for value in (True, 1.0, "-1", 0, -2)]
    calls += [lambda value=value: connection.require_domain(value, 1) for value in (True, 1.0, "1/2", 0, -1, sp.oo, sp.I, sp.sqrt(2)/4)]
    calls += [lambda value=value: connection.require_domain(sp.Rational(1, 2), value) for value in (True, 1.0, 0, -1, sp.nan)]
    calls += [lambda value=value: serialize(value) for value in (1.0, sp.Float(1), sp.oo, sp.nan)]
    rejected = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(calls):
        raise ValueError("An invalid input, algebraic locus or certificate value was accepted")
    if connection.exceptional()["quotient_rank"] != 57:
        raise ValueError("The three extra unsourced algebraic directions were lost")
    if lower.palatini_control()["printed_at_nonunit"] != -sp.Rational(135, 2):
        raise ValueError("The nonunit source-transcription control failed")
    return {"rejected_inputs": rejected, "printed_Q2_rejected_by_literal_action": True,
            "nonzero_Delta_not_mistaken_for_full_connection_rank": True,
            "compatible_exceptional_source_not_mistaken_for_uniqueness": True,
            "closed_endpoint_bounds_are_nonstrict": True,
            "fixed_coefficient_integral_not_a_retarded_time_inverse": True,
            "algebraic_inverse_not_a_mass_or_kinetic_health_certificate": True,
            "full_scalar_boundary_and_original_matter_retained": True,
            "no_unproved_kinetic_completion_imported": True}


@cache
def build_report():
    pins = prior_checks()
    groups = {
        "original_CD_principal_dictionary": dictionary.principal_identities(),
        "general_A1_zero_inverse": dictionary.generic_inverse_identities(),
        "lower_order_coefficient_matching": lower.identities(),
        "literal_normal_frame_boundary": lower.divergence_identities(),
        "full_connection_action_and_inverse": connection.checks(),
        "independent_signature_coefficients_and_domain": audit.identities(),
        "action_target_and_independent_premise_interfaces": bridges.identities(),
    }
    residuals = {name: certify_residuals(values) for name, values in groups.items()}
    checks = {"independent_domain_and_ODE": audit.checks(), "premise_discharge": bridges.checks()}
    if not all(value is True for group in checks.values() for value in group.values()):
        raise ValueError("A continuous coefficient or domain proof margin failed")
    return {
        "schema": 1, "claim": "P8-S6.37.AFFINE", "date": "2026-09-07",
        "status": "EXACT_SOURCE_PRESERVING_CLASSICAL_AUXILIARY_CD_M1_LIFT; KINETIC_UV_PARENT_OPEN; ORIGINAL_P8_OPEN",
        "prior_sha256": pins,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": ["notes/connection.md", "notes/dictionary.md", "notes/audit.md"],
        "domain": "all real u; closed 9/10<=X_repo<=11/10 inside open 4/5<X_repo<6/5; M^2,tau>0",
        "source_conventions": "g_source=-g_repo, x=-X_repo; Gamma^a_bc has derivative index last; unrestricted connection modulo four projective modes",
        "physical_source": "original free chi coupled to original physical g; no connection source or matter-metric reassignment",
        "exact_residuals": residuals,
        "named_exact_check_count": sum(map(len, groups.values())),
        "checked_scalar_entries": sum(value.rows*value.cols if isinstance(value, sp.MatrixBase) else 1
                                      for group in groups.values() for value in group.values()),
        "proof_checks": checks,
        "principal_lift": serialize(dictionary.lift()),
        "original_target_in_source_conventions": serialize(dictionary.target()),
        "fixed_lower_coefficient_ODE": serialize(lower.ode()),
        "parent_lower_coefficients": {
            "J3": "-(q+f_u)/(4*p*x)",
            "J2": "F_repo(u,-x)+x*q_u-3*x*(q+2*f_u)^2/(16*p^2)",
            "coefficient_integral_not_physical_history": True,
        },
        "literal_reduced_coefficients_with_generic_derivative_forcing": serialize(connection.coefficients()),
        "full_connection_rank_and_inverse": serialize(connection.calibration()),
        "independent_closed_tube_and_ODE_bounds": serialize(audit.calibration()),
        "all_60_inverse_row_bounds": serialize(connection.inverse_bound()["row_bounds"]),
        "physical_nonunit_calibration": serialize(dictionary.tube_point(physical_X=Fraction(2401, 2500), mass_squared=3, time_scale=2)),
        "pure_palatini_printed_formula_control": serialize(lower.palatini_control()),
        "controls": controls(),
        "solution_lifting": "Exact action equality modulo explicit covariant divergences, full 64 connection Euler equations and nondegenerate gauge quotient imply the same g/phi/chi equations and locally unique auxiliary lift on the open tube",
        "not_established": [
            "A propagating connection spectrum, heavy gap or healthy kinetic promotion",
            "A smooth vacuum extension or its controlled matching to this tube",
            "Loop, threshold, Regge/IR or dispersion remainders",
            "V, G, full B, positivity, UV completion or original P8 closure",
            "A new nonlinear-stability or affine-autoparallel completeness theorem"],
        "verification_boundary": "Exact finite-dimensional action algebra and separately written covariance, stationary-variation and smooth coefficient-ODE proofs; not proof-assistant formalization",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The auxiliary CD/M1 affine report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.37.AFFINE replay passed; auxiliary only, original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
