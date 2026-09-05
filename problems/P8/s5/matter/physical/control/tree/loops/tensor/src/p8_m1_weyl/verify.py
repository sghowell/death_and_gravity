"""Read-only replay of a constant-Weyl² candidate's reduced tensor response."""

import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_m1_loops import verify as prior

from . import bounds, reduction

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"cd-weyl-tensor.json"
PRIOR_SHA = "da7e0b445245dd87369380c087dc8a6a45361b68e2b4711f2c6ffa989db42fb1"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prior_checks():
    if sha(prior.REPORT) != PRIOR_SHA:
        raise ValueError("Pinned S5.9.CD certificate changed")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return PRIOR_SHA


def residuals():
    groups = {"off_shell_operator_factorization": reduction.operator_checks(),
              "off_shell_and_branch_field_maps": reduction.field_map_checks(),
              "old_canonical_tensor_normalization": reduction.canonical_checks(),
              "original_fourth_order_residual": reduction.residual_checks(),
              "compact_coefficients_and_derivative_jets": reduction.compact_checks(),
              "finite_residual_coefficient_majorants": bounds.residual_majorant_checks()}
    for name, values in groups.items():
        if any(sp.expand(value) != 0 for value in values.values()):
            raise ValueError(f"A Weyl tensor response identity failed: {name}")
    return groups


def control_checks():
    values = reduction.controls()
    zero_controls = {"flat_background_reduced_correction", "speed_shift_zero_at_u_one"}
    for key, value in values.items():
        if (sp.simplify(value) == 0) != (key in zero_controls):
            raise ValueError(f"A Weyl tensor response control failed: {key}")
    return values


@cache
def build_report():
    previous = prior_checks()
    identities, controls = residuals(), control_checks()
    sources = sorted(ROOT.glob("src/p8_m1_weyl/*.py"))+sorted(ROOT.glob("tests/*.py"))
    sources += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    return {"schema": 1, "claim": "P8-S5.10.CD", "date": "2026-09-05",
            "status": "CONSTANT_WEYL_CANDIDATE_REDUCED_TENSOR_FINITE_WINDOW_CONTROL; FULL_CANDIDATE_AND_QUANTUM_CAUSALITY_OPEN",
            "prior_S5_9_CD_sha256": previous,
            "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in sources},
            "formulation": "FORMULATION.md", "written_proof": "notes/tensor-response.md",
            "primary_sources": {
                "perturbative_branch_and_order_reduction": "https://arxiv.org/pdf/gr-qc/9211002",
                "field_redefinitions_and_finite_order_equivalence": "https://arxiv.org/pdf/1709.09695",
                "reduced_equation_resummation_and_action_cautions": "https://arxiv.org/pdf/1710.01562"},
            "algebraic_residuals": {name: dict.fromkeys(values, "0") for name, values in identities.items()},
            "positive_and_omission_controls": {key: str(value) for key, value in controls.items()},
            "physical_tensor_equation": "gamma_ddot+[3H+8beta*(Hddot-2H*Hdot)]gamma_dot+[1-16beta*Hdot]*q*gamma=0; chosen representative, beta=cC/M^2",
            "off_shell_field_map": "gamma=y+2beta*E0[y]-8beta*H*y_dot+O(beta^2)",
            "branch_field_map": "gamma=y-8beta*H*y_dot+O(beta^2), only when E0[y]=O(beta)",
            "reduced_y_action": "M^2/8 integral a^3*((1+16beta*H^2)*y_dot^2-q*y^2), accurate through first order",
            "exact_original_operator_residual": "Efull[gamma_red]=4beta^2*(L0+2(H^2-Hdot))*L1*gamma_red, L1=8(Hddot-2H*Hdot)*dt-16Hdot*q",
            "compact_coefficients": {key: str(value) for key, value in reduction.compact_coefficients().items()},
            "all_time_coefficient_bounds": {"speed_shift": "64*abs(beta)/ell^2", "friction_shift": "320*abs(beta)/ell^3",
                                            "canonical_mass_shift": "abs(beta)*(64*q/ell^2+1920/ell^4)",
                                            "Hubble_derivative_jets": "abs(ell^(n+1)*H^(n))<=4*n!, verified recurrence through n=4 and Chebyshev identity"},
            "finite_window_energy_Duhamel_and_residual": bounds.report(),
            "operational_contract": {
                "candidate": "NEW S_CD/M1+cC*integral sqrt(-g)*C^2 with constant cC; old action/certificates not overwritten",
                "background": "exact old flat FLRW background retained because full first variation of C^2 vanishes there",
                "sector": "linear tensor sector only, both unit-norm TT polarizations; physical matter metric unchanged",
                "variable_and_time_units": "Y=(M/2)*a^(3/2)*gamma; sigma=(t-t0)/ell0; P=dY/dsigma; w0=ell0^2*W0; N_t=sqrt(abs(P)^2+w0*abs(Y)^2)",
                "matching": "an independent constant coefficient bound; no finite matching value inferred from isolated running",
                "initial_data": "same physical gamma,gamma_dot at centre t0, equivalently same old Y,P; not equal y data or a changed Gaussian state",
                "comparison": "exact solutions of the chosen reduced second-order representative versus the old tensor equation on one finite window",
                "sign_statement": "for nonzero cC the leading physical two-derivative speed shift changes sign as Hdot changes sign; not an ultraviolet front-velocity verdict",
                "residual_scope": "explicit O(beta^2) residual and first-order formal matching, not existence or closeness of an exact fourth-order branch",
            },
            "verification_boundary": [
                "S5.9 and its prior linear/free/tree lineage are pinned and replayed",
                "independent tests contract dynamic linearized curvature and vary the literal tensor action",
                "exact arithmetic checks operator identities, field maps, polynomial majorants and finite-window constants",
                "energy/Duhamel and perturbative interpretation are written proofs, not formalized PDE or UV theorems"],
            "not_established": [
                "scalar/matter constraints or complete candidate linear stability, nonlinear/tree or loop control",
                "a unique physical exact-fourth-order branch, control of arbitrary extra initial data, or a UV front velocity",
                "full quantum superluminality, a physical EFT ghost, a UV exclusion/completion or the full chi-loop correction",
                "bounds on finite matching, other operators, nonlocal/anomaly/state terms or unknown higher-order remainders",
                "global fixed-comoving evolution, a global vacuum, a corrected Gaussian-state comparison or backreaction solution control",
                "P8 completion or replacement of any old certificate"]}


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("Constant-Weyl tensor certificate differs from exact replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S5.10: constant-Weyl reduced tensor response passed; full candidate/quantum causality OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
