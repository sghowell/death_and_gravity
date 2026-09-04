"""Read-only exact replay of S5.1, not a strong-coupling certificate."""

import argparse
import hashlib
import json
from fractions import Fraction
from functools import cache
from pathlib import Path

import sympy as sp
from p8.signs import even_polynomial, even_rational_positive, sturm_positive

from . import compact, lapse_series, nonlinear_d, scales

ROOT = Path(__file__).resolve().parents[2]
P8 = ROOT.parent
REPORT = ROOT/"certificates"/"d-nonlinear.json"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_prior_sources():
    """Check the checkpoint's pinned inputs, without rewriting its certificates."""
    prior = json.loads((P8/"certificates"/"classification.json").read_text())
    for relative, expected in prior["source_sha256"].items():
        if sha(P8/relative) != expected:
            raise ValueError(f"Linear checkpoint source changed: {relative}")
    return sha(P8/"certificates"/"classification.json")


def require_zero(residuals):
    bad = {key: value for key, value in residuals.items() if sp.cancel(value) != 0}
    if bad:
        raise ValueError(f"Nonzero exact S5 residuals: {bad}")
    return {key: "0" for key in residuals}


@cache
def build_report():
    prior_hash = check_prior_sources()
    m = nonlinear_d
    f, jets, series = m.functions(), m.lapse_jets(), lapse_series.derive()
    coefficient_data = {}
    factors = {}
    for group, orders in series.items():
        coefficient_data[group] = []
        for expr in orders:
            entries = {}
            for powers, value in lapse_series.coefficients(expr).items():
                den = sp.fraction(sp.cancel(value))[1]
                constant, factorization = sp.factor_list(den, m.u)
                if constant == 0:
                    raise ValueError("Zero denominator")
                encoded = []
                for base, exponent in factorization:
                    if sp.LC(sp.Poly(base, m.u)) < 0:
                        base = -base
                    label = str(base)
                    if label not in factors:
                        factors[label] = sturm_positive(even_polynomial(base, m.u))
                    encoded.append({"factor": label, "multiplicity": exponent})
                entries[powers] = {"coefficient": str(value), "at_bounce": str(value.subs(m.u, 0)),
                                   "nonzero_denominator_factors": encoded}
            coefficient_data[group].append(entries)
    J, lam = f["J"], f["Lambda"]
    sign_proofs = {
        "J_positive": even_rational_positive(J, m.u),
        "weighted_J_lower_half": even_rational_positive(m.d*J-sp.Rational(1, 2), m.u),
        "gamma_minus_Lambda_exceeds_half": even_rational_positive(-lam-sp.Rational(1, 2), m.u,
                                                                   upper=Fraction(1, 16)),
    }
    sc = scales.derive()
    # The two nonzero remainders below are positive, not exact-zero claims.
    zero_sc = {k: v for k, v in sc["residuals"].items() if k != "A3_variation_squared_bound"}
    require_zero({"J_upper_identity": 2/m.d-J-(4*m.u**6+12*m.u**4)/m.d**5,
                  "A3_variation_bound_identity": sc["residuals"]["A3_variation_squared_bound"]-36/m.d**2})
    sources = sorted(ROOT.glob("src/p8_s5/*.py")) + sorted(ROOT.glob("tests/*.py"))
    sources += sorted(ROOT.glob("*.md")) + sorted(ROOT.glob("notes/*.md"))
    return {
        "schema": 1, "claim": "P8-S5.1.D", "date": "2026-09-04",
        "status": "NONLINEAR_IDENTITIES_VERIFIED; PERTURBATIVE_CONTROL_OPEN",
        "scope": "D-only M0: exact nonlinear lapse-derivative removal, canonical Hamiltonian and stationary invariant expansion through fourth order",
        "prior_classification_sha256": prior_hash,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in sources},
        "exact_residuals": {"nonlinear_audit": require_zero(m.audit_checks()),
                            "generic_stationary_series": require_zero(lapse_series.generic_checks()),
                            "dimensionful_covariant_background": require_zero(scales.scaled_background_checks()),
                            "compact_nonlinear_normalization": require_zero(compact.checks()),
                            "scale_and_chart_identities": require_zero(zero_sc)},
        "lapse_jets": {key: list(map(str, values)) for key, values in jets.items()},
        "invariant_ordering": {"monomial_exponents": ["sigma", "rho", "shear2"],
                               "weights": [1, 1, 2], "lapse_orders": [1, 2, 3],
                               "Hamiltonian_orders": [0, 1, 2, 3, 4]},
        "stationary_coefficients": coefficient_data,
        "denominator_factor_proofs": factors,
        "all_time_sign_proofs": sign_proofs,
        "weighted_J_bounds": {"lower_strict": "1/2", "upper_weak": "2", "weight": "1+u^2",
                              "upper_gap": "(4*u^6+12*u^4)/(1+u^2)^4"},
        "lapse_hessian": {"expression": str(-2*J), "at_bounce": "-4",
                          "asymptotic_u_squared_times_J": str(sp.limit(m.u**2*J, m.u, sp.oo))},
        "scale_restoration": {"action_prefactor": str(sc["action_prefactor"]),
                              **sc["physical_functions"],
                              "E_ref": "1/(tau*sqrt(1+u^2))",
                              "E_curvature_over_E_ref": ["sqrt(2) inclusive", "sqrt(6) strict"],
                              "canonical_vertex_fixed_frequency_prefactor": "(M*tau)^(2-n)"},
        "finite_q_kinetic_check": {"domain": "abs(u)<=1/4, q=tau^2*k_physical^2>=64",
                                   "exact_kinetic": str(sp.factor(sc["gamma_kinetic_exact"])),
                                   "relative_error_upper": "8/(q-8) <= 1/7",
                                   "not_a_cutoff_or_complete_WKB_bound": True},
        "local_background_unit_compactification": {
            "x": "u/sqrt(1+u^2)", "ell": "tau*sqrt(1+u^2)",
            "P": str(compact.P), "P_bounds": ["1/4 strict", "1 inclusive"],
            "normalized_lapse_hessian": "-4*P, hence -4 <= h_NN < -1",
            "coordinate_warning": "Pointwise constant choice of local units, not a global time-dependent canonical transformation",
            "stationary_coefficient_bounds": compact.coefficient_report(),
            "principal_kinetic_normalizations": {
                "unitary": "abs(u)>=1/8: K_v/M^2=P/(2*x^2), 1/8<K_v/M^2<=65/2",
                "gamma": "abs(u)<=1/4, beta=b/ell: K_beta/M^2=2*P/Lambda^2, 1/2<K_beta/M^2<=8",
                "finite_q_gamma": "ell^2*k_physical^2>=68 implies the recorded kinetic relative error <=1/7",
                "gradient": "Same principal coefficients since the pinned D witness has K=G"},
            "physical_mode_reduction_still_required": True},
        "written_lemmas": ["notes/nonlinear-d.md", "notes/audit.md"],
        "not_established": ["fully momentum-reduced cubic or quartic vertices", "physical 2-to-2 amplitudes",
                            "uniform all-time weak-coupling hierarchy", "all-orders or loop control",
                            "nonlinear PDE stability", "M1 extension", "UV positivity or completion"],
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("S5.1 certificate differs from exact replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S5.1: nonlinear D identities and quartic invariant replay passed; control gate remains OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
