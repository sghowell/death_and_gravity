"""Read-only replay of the isolated CD/M1 matter-loop local-counterterm audit."""

import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_m1_tree import verify as prior

from . import background, heat_kernel

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"cd-matter-loop.json"
PRIOR_SHA = "839151b1ae103d4945b54ed703d9bb35e0aa4aa5ff3d173c20e7859767dd9787"
CLAIM = "P8-S5.9.CD"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prior_checks():
    if sha(prior.REPORT) != PRIOR_SHA:
        raise ValueError("Pinned S5.8.CD certificate changed")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return PRIOR_SHA


def residuals():
    groups = {"heat_kernel_decompositions": heat_kernel.decomposition_checks(),
              "real_scalar_pole_and_convention_signs": heat_kernel.convention_checks(),
              "illustrative_Einstein_field_redefinition_trade": heat_kernel.field_redefinition_checks(),
              "CD_curvature_and_covariant_variation": background.algebra_checks(),
              "tensor_principal_coefficient_only": background.tensor_principal_checks()}
    for group, values in groups.items():
        if any(sp.simplify(value) != 0 for value in values.values()):
            raise ValueError(f"A loop-audit algebraic residual failed: {group}")
    return groups


def control_checks():
    heat = heat_kernel.controls()
    geometry = background.controls()
    zeros = {"heat": ("conformal_R2_coefficient", "Ricci_flat_local_is_Euler"),
             "geometry": ("R_flat_radiation_R", "R_flat_radiation_I_R2_00",
                          "R_flat_radiation_I_R2_pressure", "full_A2_density_bounce")}
    for name, values in (("heat", heat), ("geometry", geometry)):
        for key, value in values.items():
            vanishes = sp.simplify(value) == 0
            if vanishes != (key in zeros[name]):
                raise ValueError(f"A loop-audit control failed: {name}/{key}")
    witness = heat_kernel.nonclosure_witness()
    if witness["acceleration_hessian"] != witness["a"] or witness["fourth_derivative_coefficient"] != witness["a"]:
        raise ValueError("Nonclosure witness lost its nonzero off-shell principal symbol")
    if witness["strict_class_acceleration_hessian"] != 0:
        raise ValueError("Strict-class restriction changed")
    return {"heat": {key: str(value) for key, value in heat.items()},
            "geometry": {key: str(value) for key, value in geometry.items()}}


@cache
def build_report():
    previous = prior_checks()
    algebra = residuals()
    controls = control_checks()
    sources = sorted(ROOT.glob("src/p8_m1_loops/*.py"))+sorted(ROOT.glob("tests/*.py"))
    sources += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    witness = heat_kernel.nonclosure_witness()
    return {
        "schema": 1, "claim": CLAIM, "date": "2026-09-05",
        "status": "ISOLATED_M1_LOCAL_ONE_LOOP_COUNTERTERM_AND_FROZEN_VARIABLE_NONCLOSURE; FULL_QUANTUM_CONTROL_OPEN",
        "prior_S5_8_CD_sha256": previous,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in sources},
        "formulation": "FORMULATION.md", "written_proof": "notes/one-loop.md",
        "primary_sources": {
            "heat_kernel": {"url": "https://arxiv.org/pdf/hep-th/0306138",
                            "equations": ["1.16", "1.18", "2.32", "3.6", "4.28"],
                            "sections": ["2.3 (local Lorentzian counterterms)", "7.1 (anomaly distinction)", "8.2 (nonlocal expansion)"]},
            "equations_of_motion_and_field_redefinitions": {"url": "https://www.numdam.org/article/AIHPA_1974__20_1_69_0.pdf",
                                                          "sections": ["6", "8"], "numerical_coefficients_imported": False},
            "perturbative_higher_derivative_interpretation": {"url": "https://arxiv.org/pdf/gr-qc/9211002", "sections": ["I", "II"]},
        },
        "algebraic_residuals": {group: dict.fromkeys(values, "0") for group, values in algebra.items()},
        "positive_and_omission_controls": controls,
        "real_minimal_scalar_local_A2": str(heat_kernel.scalar_a2()),
        "C2_E4_R2_BoxR_coefficients": {key: str(value) for key, value in heat_kernel.bulk_coefficients().items()},
        "convention_table": heat_kernel.convention_table(),
        "off_shell_fixed_variable_nonclosure": {
            "restriction": "phi=t, N=1, X_clock=1, arbitrary smooth a(t)>0",
            "bulk_density_mod_Euler_and_total_derivatives": str(witness["lagrangian"]),
            "acceleration_hessian": str(witness["acceleration_hessian"]),
            "fourth_order_Euler_Lagrange_coefficient": str(witness["fourth_derivative_coefficient"]),
            "strict_quadratic_DHOST_restriction_acceleration_hessian": "0",
            "equivalence_scope": "modulo total derivatives in frozen variables, not modulo derivative field redefinitions or CD equations of motion",
        },
        "exact_background_bounds_and_conditional_scaling": background.report(),
        "operational_contract": {
            "quantized_fields": "one real, massless, free, minimally coupled chi only; g and clock are external",
            "matter_frame": "unchanged original P8 physical metric, signature +---, R=-6*(Hdot+2H^2), Einstein action -M^2 R/2",
            "loop_order": "order hbar isolated chi determinant; no mixed, clock, graviton or other matter loops",
            "bulk_UV_scope": "local curvature-squared pole; no physical boundary/initial-surface counterterms or global determinant computed",
            "topological_quotient": "strictly four-dimensional bulk variation, compactly supported variations; dimensional Euler pole is retained before continuation",
            "scale_comparison": "mu and mu0 are constants; no substitution of a spacetime-dependent scale into the action",
            "comparison_denominator": "positive reference M^2/ell^2, not the actual vanishing Einstein density at the bounce",
            "conditional_smallness": "only specified constant logarithmic increment, or independently supplied finite cR bound; not a full quantum remainder",
            "S6_relationship": "higher EFT operators are permitted as separately named matching candidates; old finite DHOST truncation is not overwritten",
        },
        "verification_boundary": [
            "prior S5.8 is hash-pinned and replayed, not retroactively quantum-corrected",
            "the universal heat coefficient/local Lorentzian renormalization theorem is an external primary-source input",
            "exact algebra checks decomposition, pole normalization, off-shell symbol, all-time polynomial extrema and dimensions",
            "independent tests vary the lapse before gauge fixing and contract the linearized Weyl tensor before specializing",
            "field-redefinition and perturbative-EFT limitations are written arguments, not a complete quantum-equivalence or stability proof",
        ],
        "not_established": [
            "net beta functions or absence of cancellations in the fully quantized CD/M1 theory",
            "bounds on finite matching coefficients, complete renormalized stress tensor, anomaly or massless nonlocal/state-dependent terms",
            "a global Euclidean CD determinant, global vacuum, in/out S matrix, or infrared/boundary prescription",
            "a mass-gap or massive Schwinger-DeWitt expansion for the massless M1 field",
            "corrected coupled constraint reduction, kinetic positivity, propagation cone, physical ghost or quantum superluminality",
            "backreaction solution control, higher-derivative order-reduced evolution, higher-loop/all-order remainders or technical naturalness",
            "S6 UV matching/admissibility/completion, optimal physical cutoff, or full P8 completion",
        ],
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("CD/M1 isolated-loop certificate differs from exact replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S5.9: isolated M1 local-counterterm audit passed; full quantum control OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
