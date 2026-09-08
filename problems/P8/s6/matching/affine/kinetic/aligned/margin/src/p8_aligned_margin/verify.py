"""Read-only new-action clock-cone margin certificate; original P8 stays open."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as affine
from p8_vector_hadamard import transfer as vector_transfer
from p8_vector_hadamard import verify as parent

from . import action, bounds, dynamics

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"clock-cone-margin.json"
PARENT_SHA = "8ba1aa16cee4dbe890b233fac864a83de569402e06b26c9c18ae56ef763a9200"
sha, serialize = affine.sha, affine.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_aligned_margin/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen controlled Hadamard vector-state certificate changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {"S6_55_fully_rebuilt": PARENT_SHA,
            "frozen_action_certificates_not_overwritten": True,
            "vector_clock_operator_state_and_finite_matching_unchanged": True}


def residuals():
    return {**action.checks(), **dynamics.checks(), **bounds.checks()}


def controls():
    bad = (True, False, sp.true, sp.false, 1.0, sp.Float(1), "1", sp.I,
           sp.oo, -sp.oo, sp.zoo, sp.nan, sp.Symbol("unproved"), 0, -1)
    eps = sp.Rational(1, 10**6)
    calls = [lambda value=value: bounds.physical_bounds(value, 10**12, 1000)
             for value in (*bad, sp.Rational(1, 50))]
    calls += [lambda value=value: bounds.physical_bounds(eps, value, 1000) for value in bad]
    calls += [lambda value=value: bounds.physical_bounds(eps, 10**12, value) for value in (*bad, 999)]
    calls += [lambda value=value: serialize(value) for value in (1.0, sp.Float(1), sp.oo, sp.nan)]
    rejected = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(calls):
        raise ValueError("An invalid deformation, scale or report value was accepted")
    return {"rejected_inputs": rejected,
            "new_covariant_lower_scalar_not_reduced_matrix_edit": True,
            "both_scalar_charts_keep_actual_moving_boundary": True,
            "finite_q_velocity_chart_pole_not_phase_or_cutoff_pole": True,
            "negative_deformation_superluminal_control_not_allowed_model": True,
            "free_matter_luminal_direction_not_claimed_strictly_protected": True,
            "isolated_P_xx_not_full_quantum_cone_verdict": True,
            "nonzero_loop_tadpoles_not_silently_cancelled": True}


@cache
def build_report():
    pins, identities = prior_checks(), residuals()
    exact, proofs = affine.certify_residuals(identities), bounds.proof_checks()
    if not all(value is True for value in proofs.values()):
        raise ValueError("A continuous cone, chart or finite-jet bound failed")
    scalar_count = sum(value.rows*value.cols if isinstance(value, sp.MatrixBase) else 1 for value in identities.values())
    return {"schema": 1, "claim": "P8-S6.56.CLOCK_MARGIN", "date": "2026-09-08",
            "status": "NEW_BACKGROUND_PRESERVING_CLASSICAL_CLOCK_CONE_MARGIN_CANDIDATE; ORIGINAL_P8_OPEN",
            "prior_sha256": pins,
            "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in source_files()},
            "formulation": "FORMULATION.md",
            "written_proofs": ["notes/action.md", "notes/cones.md", "notes/finite-jet.md"],
            "exact_residuals": exact, "named_exact_check_count": len(identities),
            "checked_scalar_entries": scalar_count, "proof_checks": proofs,
            "literal_new_action": "Add Delta P=(M²/tau²)*epsilon_margin*(x+1)²/h² in the original physical matter metric, h=(1+u²)³,u=phi/tau. Require 0<epsilon_margin<=1/100. All other terms, vector alignment/mass/source functions, finite prescription and free-vector state remain unchanged. This is a new named candidate, not an overwrite of S6.42--S6.55.",
            "background_and_full_volume_variation": "Delta P and its first metric/clock variations vanish identically on x=-1. The exact classical complete bounce and free-matter history remain solutions. Differentiating the full density N*[(h-1+N^-2)/h]^(-3/4)*epsilon*(1-N^-2)²/h² gives the added lapse-square coefficient delta_J=4epsilon/h².",
            "literal_constraints_and_regular_Hamiltonian": "The unchanged temporal and shift equations leave only J->J+delta_J in the regular Hamiltonian C0+R²/(J+delta_J), with n=-R/(J+delta_J) and no inverse Theta. The vector quadratic block stays decoupled and unchanged.",
            "two_chart_principal_derivation": "The unitary chart has K=K_old+diag(delta_J/Theta²,0), G=K_old. A direct canonical Legendre calculation with -H*b*P_b and (a³q)'=H*a³q gives the gamma chart K=K_gamma_old+diag(delta_J/Lambda²,0), G=K_gamma_old. The beta' contribution is retained before the high-q limit, and an independent center time-jet calculation reproduces G(0)=[[6,1/20],[1/20,1/2]].",
            "physical_characteristics": "Both charts give det(G-c²K)/det(K)=(1-c²)[J/(J+delta_J)-c²]. J>0 and delta_J>0 imply positive K,G and positive semidefinite K-G. The clock mode is strictly subluminal, matter stays luminal; vector and tensor fronts remain one. Massive finite-q phase velocities are not used.",
            "chart_coverage_and_continuous_margins": "For |u|<=1/4, Lambda<=-1/5; for 1/4<=|u|<=1/2, abs(Theta)>=3/5. These charts cover the bounce interval. Theta is nonzero at every other finite time. Positive even polynomial coefficients prove J*h²>=1199/800 globally and J*h²<=3817516721573/107374182400<36 on |u|<=1/2. Thus c_clock²>=1199/(1199+3200epsilon); the compact margin and fractional kinetic increase have the recorded exact bounds.",
            "nonlinear_constraint_count_scope": "The full added Hamiltonian is -N*U*Delta P and introduces no auxiliary velocity/spatial derivative or momentum dependence. The rolling joint lapse/temporal Hessian becomes diag(-2[J+delta_J],-1). S6.45's smooth local implicit-function/Dirac argument persists with seven physical modes including free matter. No quantitative neighborhood radius, nonlinear energy or interacting cutoff follows.",
            "global_polynomial_bound_data": serialize(bounds.polynomial_bounds()),
            "center_data": serialize(dynamics.center()),
            "finite_potential_jet_diagnostic": serialize(bounds.potential_jet()),
            "finite_potential_interpretation": "At mu=m0 the fixed vector potential has B_xx=(B_NN+3B_N)/4=2392/(6561h²)>0, hence negative Lorentzian P_xx. Its ratio to the added positive 2epsilon/h² is <=R^4*1196/(576L²*6561epsilon). The isolated fixed-metric clock identity is K-G=4P_xx, not the full constrained quantum principal matrix. Nonzero loop background/first variations are not cancelled.",
            "physical_scale_example": serialize(bounds.physical_bounds(sp.Rational(1, 10**6), 10**12, 1000)),
            "unchanged_vector_first_variation_example": serialize(vector_transfer.physical_bounds(10**12, 1000)),
            "controls": controls(),
            "verdict": "This new named candidate preserves the exact classical history, local nonlinear degree count and vector state controls while giving an explicit positive clock-cone margin. It does not establish the quantum-corrected system or close P8.",
            "not_established": ["Quantum-corrected background and constraints, complete causal second functional variations and physical cones",
                                "Protection of the still-luminal matter/tensor/vector directions against arbitrary unsigned corrections",
                                "Interacting cutoff, other fields/loops, common vacuum/UV parent, finite-gravity Regge remainder or V/G/B",
                                "Nonlinear energy stability, global well-posedness or original P8 closure"],
            "verification_boundary": "Exact symbolic/rational controls and written action, nonlinear-constraint corollary, overlapping principal-chart and continuous-bound proofs; not proof-assistant formalization or a full quantum solution"}


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The new clock-cone margin report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.56.CLOCK_MARGIN replay passed; new classical margin candidate, original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
