"""Read-only CD/M1 nonlinear preparation and principal-robustness replay."""

import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8.signs import even_rational_positive
from p8_s5.lapse_series import generic_checks

from . import nonlinear, regressions, robustness, series

ROOT = Path(__file__).resolve().parents[2]
P8 = ROOT.parents[1]
REPORT = ROOT / "certificates" / "cd-matter.json"
PRIOR = {
    "certificates/classification.json": "4e82b45d3d1daed0a5e4698f1dd39fa51e2eb09b6efcf3cf5caac3a12d45dec8",
    "certificates/witness-CD_matter.json": "caf8c8e688a7565b9d00f921c099a28da00f97522ed26ad182a8227eb80cd4dd",
    "s5/certificates/d-nonlinear.json": "5a1946d307a18ffd666c66208a689e23698cd18b2c7fa5ffb37100c29f314c59",
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prior_checks():
    for name, expected in PRIOR.items():
        path = P8 / name
        if sha(path) != expected:
            raise ValueError(f"Pinned prior certificate changed: {name}")
        if name.endswith(("classification.json", "d-nonlinear.json")):
            base = P8 if name.endswith("classification.json") else P8 / "s5"
            prior = json.loads(path.read_text())
            for relative, digest in prior["source_sha256"].items():
                if sha(base / relative) != digest:
                    raise ValueError(f"Pinned prior source changed: {base / relative}")
    return dict(PRIOR)


def require_zero(residuals):
    bad = {name: value for name, value in residuals.items() if sp.cancel(value) != 0}
    if bad:
        raise ValueError(f"Nonzero CD/M1 residuals: {bad}")
    return dict.fromkeys(residuals, "0")


@cache
def build_report():
    prior = prior_checks()
    f, jets = nonlinear.functions(), nonlinear.lapse_jets()
    sources = sorted(ROOT.glob("src/p8_m1/*.py"))+sorted(ROOT.glob("tests/*.py"))
    sources += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    coefficient_data = series.coefficient_report()
    compact = series.compact_jets()
    return {
        "schema": 1, "claim": "P8-S5.5.CD", "related_claims": ["P8-S5.5.cone"], "date": "2026-09-05",
        "status": "M1_NONLINEAR_PREPARATION_AND_PRINCIPAL_ROBUSTNESS; PHYSICAL_INTERACTION_CONTROL_OPEN",
        "prior_sha256": prior,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in sources},
        "exact_residuals": {
            "nonlinear_chart_and_seven_velocity_Legendre": require_zero(nonlinear.audit_checks()),
            "independent_literal_bounce": require_zero(regressions.bounce_checks()),
            "generic_stationary_lapse_series": require_zero(generic_checks()),
            "local_units_and_compact_coefficients": require_zero(series.compact_checks()),
            "principal_normal_form_and_covariant_deformation": require_zero(robustness.audit_checks()),
        },
        "nonlinear_chart": {
            "T": str(nonlinear.chart()["T"]), "omega": "-log(T)/4", "T_bounds": ["9/10", "11/10"],
            "domain": "u real, X=N^-2 in [9/10,11/10]",
            "primitive": "Psi(u,N)=integral from 1 to N of Psi_N(u,n) dn, a local coefficient primitive",
            "psi_N": str(nonlinear.chart()["psi_N"]),
            "A": "T^(1/4)", "P": "-T_u/(2*N*T^(3/4))-Psi",
            "U": "T^(-3/4)*(F+9*T_u^2/(16*N^2*T))-Psi_u/N",
            "matter_Hamiltonian": "N*T^(3/4)*(chi_dot+eta)^2/2+N*T^(-1/4)*z/2",
            "spatial_constraint_retained": "-2*hat_h_ij*hat_D_k*Pi^jk+pi_chi*partial_i chi=0",
            "physical_matter_metric": "original g; no redefinition of the physical causal cone",
        },
        "negative_control": regressions.boundary_negative_control(),
        "lapse_jets": {key: list(map(str, row)) for key, row in jets.items()},
        "invariant_ordering": {
            "variables": list(map(str, nonlinear.VARIABLES)), "weights": list(nonlinear.WEIGHTS),
            "pieces": "A+B*sigma+C*rho+L*eta+D*sigma^2+E*shear2+M*eta^2+Z*z",
            "jets_are_derivatives_not_Taylor_coefficients": True,
            "lapse_orders": [1, 2, 3], "Hamiltonian_orders": [0, 1, 2, 3, 4],
        },
        "stationary_invariant_coefficients": coefficient_data,
        "all_time_sign_proofs": {
            "weighted_J_lower": even_rational_positive(nonlinear.d*f["J"]-sp.Rational(1, 10), nonlinear.u),
            "weighted_J_upper": even_rational_positive(8-nonlinear.d*f["J"], nonlinear.u),
        },
        "local_units": {
            "ell": "tau*sqrt(1+u^2)", "x": "u/sqrt(1+u^2)", "y": "1/sqrt(1+u^2)",
            "domain": "-1<=x<=1, 0<=y<=1, x^2+y^2=1; endpoints bound the two tails",
            "scaled_invariants": "sigma_bar=ell*sigma, eta_bar=ell*eta, rho_bar=ell^2*rho, shear2_bar=ell^2*shear2, z_bar=ell^2*z (M=tau=1 convention)",
            "normalized_lapse_hessian": "-16 < hbar_NN=-2*d*J < -1/5",
            "bound_method": "exact polynomial coefficient l1 bounds and abs(1/hbar_NN)<5",
            "compact_jets": {key: list(map(str, row)) for key, row in compact.items()},
            "jet_absolute_bounds": {key: [str(series.polynomial_bound(value)) for value in row] for key, row in compact.items()},
            "coefficient_count": sum(len(row) for rows in coefficient_data.values() for row in rows),
            "action_prefactor": "(M*tau)^2 after phi=tau*u, chi=M*chi_bar, and coordinates x_phys=tau*x_bar",
            "scope": "local fixed units and invariant coefficients; no time-dependent canonical transformation or physical cutoff inference",
        },
        "principal_robustness": {
            "hypotheses": "K=[[A,B],[B,Q]] positive definite, G=[[C,D],[D,Q]], Q>0; regular two-derivative principal chart",
            "normalization": "k=A-B^2/Q, a=A-C, b=B-D, sigma=(a-2*B*b/Q)/k, r2=b^2/(k*Q)",
            "squared_speeds": "1-sigma/2 +/- sqrt(sigma^2/4+r2)",
            "strict_gradient_health_iff": "1-sigma-r2>0",
            "healthy_causal_iff": "b=0 and 0<=a<k",
            "mismatch_implies": "c_plus^2>1 for every b!=0, regardless of diagonal buffer a",
            "chart_invariant_r2": "chi_dot^2*E^2/(2*T^2*J)",
            "chart_invariant_sigma": "(T^2*J0-P_gamma-T*w*chi_dot*E)/(T^2*J)",
            "E": "Lambda-T*(1-3*delta)=-2*X*(2*F2X+X*A3) when A1=0",
            "pinned_bounce_mismatch_control": "K=[[6,1/20],[1/20,1/2]], G=diag(599/100,1/2), c^2=1+/-1/sqrt(1199)",
            "one_sided_covariant_family": {
                "delta_F": "epsilon*J_star(phi)*(X-1)^2/(4*(1+phi^2)^2)",
                "background": "exactly unchanged because delta_F and delta_FX vanish on X=1",
                "new_J": "(1+epsilon/d^2)*J_star", "speeds": ["1", "1/(1+epsilon/d^2)"],
                "epsilon_ge_zero": "healthy and causal principal modes; one is subluminal for epsilon>0",
                "minus_one_lt_epsilon_lt_zero": "healthy but superluminal at every finite time",
                "epsilon_eq_minus_one": "kinetic degeneracy at bounce",
                "canonical_tail_corrections": [str(value) for value in robustness.deformation()["tail_coefficients"]],
                "quantitative_tail_and_first_two_derivative_bounds": {
                    key: value for key, value in robustness.deformation_tail_bounds().items() if key != "residuals"},
                "preserves_old_linear_tube_and_tail_contract_for_fixed_nonnegative_epsilon": True,
                "not_a_computed_loop_correction_or_mismatch_margin": True,
            },
        },
        "written_proofs": ["notes/nonlinear.md", "notes/robustness.md"],
        "not_established": ["fully spatially constrained M1 physical cubic/quartic vertices",
                            "finite-k dispersion or interacting weak-coupling band", "M1 amplitudes or inclusive scattering",
                            "calculated radiative corrections or technical naturalness", "all-orders/loop control",
                            "nonlinear or BKL stability", "UV admissibility or completion", "full P8 completion"],
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("CD/M1 report differs from exact replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S5.5: CD/M1 nonlinear preparation and principal-cone robustness replay passed; physical interaction control OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
