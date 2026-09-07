"""Read-only scalar microlocal principal/energy audit of the specified family."""

import argparse
import hashlib
import json
from fractions import Fraction
from functools import cache
from pathlib import Path

import sympy as sp
from p8_variable_constraints import verify as prior

from . import (
    bounds,
    domain,
    energy,
    independent,
    leading,
    physical,
    stueckelberg,
    symbol,
    weighted,
)
from . import rational as ra

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"variable-scalar-microlocal.json"
PRIOR_SHA = "9cb56291ce32499293a65ee1355e9b13d3abc4d429643f76c0dcf63e9b71a9b1"
CLAIM = "P8-S6.24.MICROLOCAL"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def serialize(value):
    if isinstance(value, ra.Rational):
        return "0" if not value else str(value.sympy())
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, sp.MatrixBase):
        return serialize(value.tolist())
    if isinstance(value, sp.Basic):
        if value.has(sp.Float, sp.oo, -sp.oo, sp.zoo, sp.nan):
            raise TypeError("Only exact finite report expressions are permitted")
        return str(value)
    if isinstance(value, dict):
        return {str(key): serialize(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [serialize(item) for item in value]
    if value is None or isinstance(value, (str, int, bool)):
        return value
    raise TypeError("Only exact JSON/report values are permitted")


def degrees(matrix):
    return [[None if value == -sp.oo else int(value) for value in row] for row in matrix]


@cache
def prior_checks():
    if sha(prior.REPORT) != PRIOR_SHA:
        raise ValueError("The frozen S6.22 actual scalar action/constraint certificate changed")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return {"S6_22_literal_scalar_action_and_regular_constraints_with_recursive_pins": PRIOR_SHA}


@cache
def independent_bridges():
    d = symbol.derive()
    records = []
    for fixture in independent.fixtures():
        at = tuple(fixture[key] for key in ("u", "c", "K"))
        count = 0
        for key in ("C", "inverse_C", "Mq", "Mv", "symplectic"):
            if [[x.evaluate(at) for x in row] for row in d[key]] != fixture[key]:
                raise ValueError("Independent Fraction/Taylor physical Cauchy replay failed")
            count += sum(map(len, fixture[key]))
        if d["det_C"].evaluate(at) != fixture["det_C"]:
            raise ValueError("Independent Fraction Cauchy determinant failed")
        records.append({"u": str(at[0]), "c": str(at[1]), "K": str(at[2]),
                        "exact_scalar_comparisons": count+1, "all_equal": True})
    l = leading.derive()
    coefficient_count = 0
    for record in domain.derive()["records"].values():
        terms = sp.Poly(record["numerator"], l["u"], l["c"]).as_dict()
        degree, values = independent.bernstein_from_time_lapse_terms(
            {power: Fraction(coefficient) for power, coefficient in terms.items()})
        if degree != record["degrees"] or values != record["coefficients"]:
            raise ValueError("Independent supplied-polynomial Bernstein conversion failed")
        coefficient_count += len(values)
    return {"literal_polarization_Fraction_Taylor_Cauchy_fixtures": records,
            "independent_supplied_polynomial_Bernstein_coefficients": coefficient_count,
            "independence_boundary": "The Fraction phase is independently polarized; the Bernstein conversion consumes supplied primary margin polynomials and is not an independent derivation of them",
            "total_exact_scalar_comparisons": sum(row["exact_scalar_comparisons"] for row in records)+coefficient_count}


def controls():
    c = leading.center()
    at = {c["c"]: 4}
    if c["gradient_without_final_boundary"][0, 0].subs(at) != 448 or c["gradient"][0, 0].subs(at) != 800:
        raise ValueError("The final pi*pi' boundary omission control failed")
    if c["Fpi_coefficient_at_center"] != 0 or c["Fpi_coefficient_derivative_at_center"].subs(at) != 352:
        raise ValueError("The center-first boundary omission was not detected")
    if sp.limit((c["c"]-2)*c["weight_second_log_jet"], c["c"], 2, dir="+") != -16:
        raise ValueError("The nonuniform joint-lapse connection guard failed")
    K = stueckelberg.action.K
    root = sp.sqrt(193)-13
    finite_det = 5*(K+60)*(K*K+26*K-24)/(19176*K**3)
    if sp.simplify(finite_det.subs(K, root)) != 0 or root.is_positive is not True:
        raise ValueError("The low-K Cauchy-chart pole guard failed")
    invalid = [lambda value=value: bounds.threshold(value) for value in (0, -1, 3, True, 0.01)]
    invalid += [lambda: physical.equation(0.0, 4), lambda: physical.equation(0, 2),
                lambda: ra.Rational(0.1), lambda: ra.Rational(True),
                lambda: independent.fixture(True, 4, 1),
                lambda: independent.fixture(0, 4, 0),
                lambda: serialize(sp.Float("0.1")), lambda: serialize(sp.oo)]
    rejected = 0
    for call in invalid:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(invalid):
        raise ValueError("An inexact or out-of-domain input bypassed validation")
    return {"c4_final_boundary_coefficient_at_center": "0", "c4_its_time_derivative": "352",
            "c4_gradient_with_omitted_boundary_NOT_correct": "448", "c4_correct_gradient": "800",
            "discarded_nonuniform_Psi_f_chart": "Frozen S6.22 documents apparent center5 and its time/momentum degree drop; it is not the new scalar cone",
            "correct_c4_center_helicity_speed_squared": "5/3",
            "finite_K_chart_zero": "sqrt(193)-13; positive and below1, not within the high-K gate",
            "joint_c_to_2_second_weight_jet_residue": "-16",
            "invalid_calls_rejected": rejected}


@cache
def build_report():
    previous = prior_checks()
    residuals = {"literal_Stueckelberg_action_and_boundaries": stueckelberg.checks(),
                 "leading_stationary_action_and_center_dictionary": leading.checks(),
                 "weighted_auxiliary_completion_and_inverse": weighted.checks()}
    if any(sp.factor(value) != 0 for group in residuals.values() for value in group.values()):
        raise ValueError("An exact scalar action or stationary-completion identity failed")
    generic, ec = symbol.checks(), energy.checks()
    positivity, cb = domain.checks(), bounds.derive()
    audit = ROOT/"tests"/"test_microlocal_independent_audit.py"
    if not audit.is_file():
        raise ValueError("The separately authored root audit must be present before freezing")
    sources = sorted(ROOT.glob("src/p8_variable_microlocal/*.py"))+sorted(ROOT.glob("tests/*.py"))
    sources += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    center = leading.center()
    return {
        "schema": 1, "claim": CLAIM, "date": "2026-09-07",
        "status": "EXACT_UNIFORM_LOCAL_SCALAR_PRINCIPAL_KINETIC_GRADIENT_AND_ENERGY; LIGHT_ONLY_MATCHING_AND_ORIGINAL_P8_OPEN",
        "prior_context_sha256": previous,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in sources},
        "formulation": "FORMULATION.md", "written_proof": "notes/review.md", "source_audit": "notes/sources.md",
        "exact_symbolic_residuals": {name: dict.fromkeys(group, "0") for name, group in residuals.items()},
        "exact_symbolic_residual_count": sum(map(len, residuals.values())),
        "generic_rational_symbol": {
            "parameter_field": "Exact Q(u,c,K); no interpolation or pointwise-in-time proof",
            "groups": {name: {"degree_bounds": degrees(group["degrees"]), "upper_degree": group["upper_degree"]}
                       for name, group in generic["groups"].items()},
            "zero_polynomial_degree_encoding": "null, not a nonfinite physical coefficient",
            "det_C_degree": generic["det_C_degree"],
            "leading_detC_times_det_kinetic_minus_one": "0",
        },
        "continuous_positive_principal_forms": serialize(positivity),
        "continuous_Bernstein_coefficient_count": sum(row["count"] for row in positivity.values()),
        "constructive_uniform_Cauchy_inverse": {
            "domain": cb["threshold_formula"],
            "leading_cubic_numerator": "c^4*d^18*(c*d^4-2)^2*(u^2-1)/344064",
            "leading_absolute_lower_factor_times_delta_squared": serialize(cb["leading_absolute_lower_factor_times_delta_squared"]),
            "lower_numerator_polynomial_coefficient_norms": serialize(cb["coefficient_norms"]),
            "rounded_bounds_n0_n1_n2": list(cb["rounded_coefficient_bounds"]),
            "relative_determinant_error_upper": serialize(cb["relative_determinant_error_upper"]),
            "not_cutoff": cb["interpretation"],
        },
        "exact_time_normalization_and_energy": {
            "normalized_Rq_degrees": degrees(ec["normalized_Rq_degrees"]),
            "normalized_Mv_degrees": degrees(ec["normalized_Mv_degrees"]),
            "order_k_energy_skew_identity_count": ec["energy_skew_identity_count"],
            "matrix_bound": "(1/4)I < U-normalized K_s,G_s <64I",
            "formula": ec["energy_rate_formula"],
            "estimate": "For each fixed positive delta_min, E(u)<=exp(C_delta*abs(u-u0))*E(u0), uniformly in all K above the stated threshold",
            "numeric_scope": "No optimized or expanded numerical C_delta is reported; boundedness follows from the explicit rational formulas on the compactified1/K box",
        },
        "center_physical_g_principal": serialize({key: center[key] for key in
            ("kinetic", "gradient", "helicity_speed_squared", "Fpi_coefficient_at_center",
             "Fpi_coefficient_derivative_at_center", "weight_second_log_jet")}),
        "independent_replay": independent_bridges(),
        "omission_and_scope_controls": controls(),
        "normalization_and_actual_action": {
            "physical_metric": "Same g, positive canonical sourced phi and free chi of frozen S6.20/S6.22",
            "curvature": "+---; R_B=-6(Hdot+2H²), Einstein action -M²R_B/2; no A-track sign import",
            "units": "u=T/tau,K=(tau*kcom)^2,overall action(Mtau)^2; physical squared speeds=eigenvalues(a² K_s^-1 G_s)",
            "regular_chart": "q=(pi=-e/K,xi=3(a³wz-PE)/(2aK),chi_B=z+wxi), no H division",
            "principal_weight": "Pi=U*pi,U=a²yP/2; actual time connections retained, not uniform atc2",
            "physical_smallness": "Energy uses fixed spatial derivative weights; no uniform unweighted raw-metric amplitude or source-response assertion",
        },
        "not_established": [
            "No finite-band light-only scalar matching or physical-source response; superluminal full principal cone is not that exclusion",
            "No cutoff, joint c->2 heavy-mode reduction, rolling spectral gap, instantaneous ground state or all-band positivity",
            "No global/nonlinear scalar health, loop or UV completion and no original DHOST/free-matter operator match",
            "No transfer to exploratory variable-lapse G1 or sibling tensor response without its own actual dictionary",
            "Completed scoped photon objective and frozen original linear classification unchanged; original P8 remains open",
        ],
        "verification_boundary": "Exact source-aware algebra, generic rational symbol identities, continuous Bernstein bounds, independent Fraction/Taylor Cauchy fixtures and root action audits; written uniform proof, not proof-assistant formalized",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("Scalar microlocal certificate differs from exact read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.24.MICROLOCAL: uniform local scalar principal/energy replay passed; light-only matching and original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
