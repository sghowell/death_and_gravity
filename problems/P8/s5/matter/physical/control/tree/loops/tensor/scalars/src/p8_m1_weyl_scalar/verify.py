"""Read-only replay: same constant-Weyl candidate, coupled scalar response."""

import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_m1_weyl import verify as prior

from . import bounds, geometry, normalization, reduction

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"cd-weyl-scalars.json"
PRIOR_SHA = "05a66f42c280e099d73ef4629a0358cf84d1eba4af5eb692dc4ff83ecbd739d8"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prior_checks():
    if sha(prior.REPORT) != PRIOR_SHA:
        raise ValueError("Pinned S5.10.CD certificate changed")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return PRIOR_SHA


def residuals():
    groups = {"physical_scalar_Weyl_contraction": geometry.scalar_contraction_checks(),
              "physical_ADM_and_gauge_map": geometry.gauge_and_ADM_checks(),
              "regular_partial_Legendre_and_stationary_constraints": reduction.auxiliary_checks(),
              "full_and_branch_phase_reduction": reduction.field_map_checks(),
              "physical_lapse_shift_metric_reconstruction": reduction.reconstruction_checks(),
              "both_chart_rank_one_Hamiltonian": reduction.rank_one_checks(),
              "pinned_coupled_quadratic_and_clock_bridge": normalization.baseline_checks(),
              "full_canonical_volume_and_symmetric_boundary_map": normalization.generic_checks(),
              "exact_CD_lensing_rows": bounds.row_checks()}
    for name, values in groups.items():
        if any(sp.simplify(value) != 0 for value in values.values()):
            raise ValueError("A coupled scalar Weyl identity failed: "+name)
    return groups


def control_checks():
    groups = {"geometry": geometry.controls(), "phase_and_constraints": reduction.controls(),
              "full_canonical_normalization": normalization.controls(), "finite_band": bounds.controls()}
    zeros = {"conformal_scalar_metric_has_no_Weyl", "zero_matching_coefficient_has_zero_bound"}
    for name, values in groups.items():
        for key, value in values.items():
            if (sp.simplify(value) == 0) != (key in zeros):
                raise ValueError(f"A coupled scalar Weyl control failed: {name}/{key}")
    return groups


@cache
def build_report():
    previous = prior_checks()
    identities, controls = residuals(), control_checks()
    finite_band = bounds.report()
    sources = sorted(ROOT.glob("src/p8_m1_weyl_scalar/*.py"))+sorted(ROOT.glob("tests/*.py"))
    sources += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    return {"schema": 1, "claim": "P8-S5.11.CD", "date": "2026-09-05",
            "status": "CONSTANT_WEYL_CANDIDATE_FIRST_ORDER_COUPLED_SCALAR_PHASE_AND_FINITE_BAND_CONTROL; FULL_CANDIDATE_AND_QUANTUM_CAUSALITY_OPEN",
            "prior_S5_10_CD_sha256": previous,
            "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in sources},
            "formulation": "FORMULATION.md", "written_proof": "notes/scalar-response.md",
            "primary_sources": {
                "independent_scalar_Weyl_coefficient_and_metric_conventions": "https://arxiv.org/pdf/1703.08012",
                "perturbative_branch_and_order_reduction": "https://arxiv.org/pdf/gr-qc/9211002",
                "finite_order_field_redefinitions": "https://arxiv.org/pdf/1709.09695",
                "reduced_representative_resummation_cautions": "https://arxiv.org/pdf/1710.01562"},
            "exact_residuals": {name: dict.fromkeys(values, "0") for name, values in identities.items()},
            "positive_and_omission_controls": {name: {key: str(value) for key, value in values.items()}
                                                for name, values in controls.items()},
            "physical_scalar_Weyl_action": "(4cC/3)*integral dt*a^3*q^2*Lensing^2; one real unit-spatial-norm Fourier mode",
            "physical_lensing": "Lensing=(1-delta)*n-v_old+b_dot-H*b; zeta_phys=v_old+delta*n",
            "auxiliary_Hessian": "diag(-2J,-2q^2/3); determinant=4Jq^2/3; regular for J>0,q>0",
            "old_stationary_lensing_CD": "Lensing0=2delta*n0-2(v+H*b0), n0=-R/J,b0=-p/(2q)",
            "off_shell_lensing": "Lensing_off=Lensing0-Ep/(2q), Ep=p_dot+3H*p-2q(v+Lambda*n0)",
            "full_first_order_phase_map": "v_old=v+beta*(Ep-4q*Lensing0)/3, beta=cC/M^2, modulo O(beta^2) and stated boundary",
            "reduced_phase_Hamiltonian": "h0-(4beta/3)*q^2*Lensing0^2; both scalar pairs, full old gamma generator retained",
            "physical_reconstruction_first_order": {
                "v_old": "v-(4beta/3)*q*Lensing0",
                "lapse": "n0-(8beta/3)*delta*q^2*Lensing0/J",
                "shift": "b0+4beta*D_t Lensing0; old coupled phase flow and coefficient derivatives",
                "spatial_metric": "zeta_phys=v_old+delta*n_phys",
                "pre_mixed_auxiliary_metric_momentum": "p+3l*s",
                "physical_normalized_matter_density": "l+P_s+3l*v_old",
                "scope": "all identities through first order only; gamma canonical b0 is not corrected physical shift"},
            "normalized_ODE": "X'=([[-Omega,I],[-W,-Omega]]-(8epsilon0/3)*qfixed^2*Jcan*r*r^T)X",
            "full_inverse_canonical_map": "Q=a^-3/2*T^-1*Y, p=a^-3/2*T^T*(P-SY); rP=T*cp,rY=T^-T*cQ-S^T*rP",
            "finite_band_control": finite_band,
            "operational_contract": {
                "same_candidate": "NEW S_CD/M1+cC*integral sqrt(|g|)*C^2, constant coefficient, frozen physical matter metric",
                "background": "exact old flat FLRW background unchanged by the full first variation",
                "sector": "both coupled linear scalar channels, not a scalar-only clock truncation",
                "regularity": "auxiliary phase reduction at every finite time and q>0; positive velocity charts require the separate covering band",
                "units": "fixed centre sigma=(t-t0)/ell0, qfixed=ell0^2*kcom^2/a^2, epsilon0=cC/(M*ell0)^2; Nt^2=P^T P+Y^T W Y",
                "stationary_scope": "formal first-order stationary-action and derivative-field-map identity, not exact integration of new shift dynamics",
                "evolution_scope": "exact evolution of a specified reduced Hamiltonian, with an exact chosen first-order phase pullback only where explicitly stated",
                "data_distinctions": "equal reduced data, equal chosen pulled-old-phase data and equal physical metric data are not conflated",
                "matching": "finite |cC| is independently supplied; isolated running and nonlocal terms do not fix or bound it"},
            "verification_boundary": [
                "S5.10 and its full pinned scalar/free/tree/loop lineage are replayed",
                "direct four-index scalar curvature, independent gauge/field-map and canonical audits are tested and hashed",
                "symbolic equalities and exact rational positive Laurent majorants are replayed; no sampled interval proof",
                "stationary-action, finite-band energy/Dyson and Neumann arguments are written proofs, not formalized PDE or UV theorems"],
            "not_established": [
                "existence or closeness of an exact higher-derivative branch or control of arbitrary additional initial data",
                "bounds on unknown beta-squared action terms, complete matching, nonlocal/anomaly/state or mixed/graviton/clock loops",
                "full candidate nonlinear, interacting/tree, vector or extra-mode health",
                "quantum superluminality, a UV front velocity, a physical EFT ghost or a UV completion/exclusion",
                "equal physical metric-data evolution error merely from equal reduced or auxiliary phase data",
                "global fixed-comoving, unbounded-momentum, corrected vacuum or backreaction solution control",
                "completion of P8 or replacement of any prior certificate"]}


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("Coupled scalar Weyl certificate differs from exact replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S5.11: first-order coupled scalar Weyl phase and finite-band bounds passed; full candidate/quantum causality OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
