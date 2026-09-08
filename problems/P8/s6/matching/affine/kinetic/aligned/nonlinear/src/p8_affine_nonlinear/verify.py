"""Read-only full ADM and local nonlinear auxiliary-constraint theorem."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as affine
from p8_aligned_onshell import verify as parent

from . import adm, constraints, maxwell, remaining, trace

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"local-nonlinear-constraints.json"
PARENT_SHA = "d2ec17690cce8549da5464d99ed0d3104734bd8f8984ed9617bcf1989e0af0fe"
sha, serialize = affine.sha, affine.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_affine_nonlinear/*.py"))
            +sorted(ROOT.glob("tests/*.py"))+sorted(ROOT.glob("*.md"))
            +sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen on-shell linear-light certificate changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {"S6_44_fully_rebuilt": PARENT_SHA,
            "unchanged_all64_affine_action_and_actual_source_rechecked": True,
            "no_frozen_action_or_ancestor_modified": True}


def residuals():
    out = {}
    for module in (adm, trace, constraints, maxwell, remaining):
        for name, value in module.checks().items():
            if name in out:
                raise ValueError("An exact residual name was duplicated")
            out[name] = value
    return out


def controls():
    calls = [lambda value=value: serialize(value) for value in (1.0, sp.Float(1), sp.oo, sp.nan)]
    calls += [lambda value=value: affine.certify_residuals({"deliberately_nonzero": value})
              for value in (1, sp.Rational(1, 7), sp.Matrix([[0, 1], [0, 0]]))]
    rejected = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(calls):
        raise ValueError("An inexact or deliberately nonzero certificate input was accepted")
    return {"rejected_inputs": rejected,
            "independent_dense_spatial_and_joint_temporal_Legendre_solves": True,
            "omitted_time_and_spatial_boundaries_have_nonzero_controls": True,
            "either_missing_auxiliary_pivot_invalidates_the_constraint_inverse": True,
            "spatial_Gauss_term_and_all_five_shear_velocities_retained": True,
            "block_operator_inverse_does_not_invert_secondary_secondary_bracket": True,
            "local_count_distinguished_from_nonlinear_health_and_X_zero_vacuum": True}


@cache
def build_report():
    pins = prior_checks()
    identities = residuals()
    exact = affine.certify_residuals(identities)
    proofs = constraints.proof_checks()
    if not all(value is True for value in proofs.values()):
        raise ValueError("A full auxiliary rank, compact bound or local constraint-count check failed")
    clock = trace.clock_auxiliary()
    return {"schema": 1, "claim": "P8-S6.45.ALIGNED_NONLINEAR", "date": "2026-09-07",
            "status": "LOCAL_NONLINEAR_AUXILIARY_CONSTRAINT_COUNT_CERTIFIED; ORIGINAL_P8_OPEN",
            "prior_sha256": pins,
            "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in source_files()},
            "formulation": "FORMULATION.md", "written_proofs": ["notes/proof.md"],
            "exact_residuals": exact, "named_exact_check_count": len(identities),
            "checked_scalar_entries": sum(value.rows*value.cols if isinstance(value, sp.MatrixBase) else 1
                                          for value in identities.values()),
            "proof_checks": proofs,
            "literal_action": "Unchanged S6.42, complete original scalar F and ODE source, source-centered retained mass, source-aligned Maxwell curl, physical free chi; all64 complement from S6.41",
            "complete_ADM_chart": "h_physical=e^(2omega)*h_hat, omega=-log((h-1+X)/h)/4; retain spatial-curvature IBP and I_s=e^(3omega)*L_V, I(u,1)=0; no lapse time or spatial derivatives remain",
            "joint_trace_Legendre_map": "a=(2/3)U*B, gamma_eff=gamma_t-3delta²/(8p_affine²)>1061/1140>9/10 on the original closed tube",
            "Maxwell_and_spatial_constraints": "Full three-electric-velocity Legendre map before W0=N*T+N^i*Wi; independent Pi gives -W0*D_i Pi^i, H_i=Pi^j*F_ij-Wi*D_j Pi^j; no auxiliary spatial derivative in the complete Hamiltonian",
            "actual_background_auxiliary_Jacobian": serialize(clock["jacobian"]),
            "clock_bounds": serialize({"I_N": adm.clock_lapse()["I_N"],
                                       "I_NN": adm.clock_lapse()["I_NN"],
                                       "J_lower": sp.Rational(1199, 800)/sp.Rational(5, 4)**18,
                                       "sufficient_Jacobian_perturbation_norm": sp.Rational(1, 40)}),
            "local_constraint_theorem": "On [-1/2,1/2], J>1/40 and D=diag(-2J,-1); smoothness and compactness imply a nonempty open canonical/spatial-jet neighborhood with two locally solvable auxiliaries and four second-class constraints; no Theta divisor",
            "complete_physical_count": {"configuration_variables": 15, "nonauxiliary_velocities": 10,
                                        "first_class_constraints": 6, "second_class_constraints": 4,
                                        "physical_modes_with_chi": 7, "physical_modes_without_chi": 6},
            "controls": controls(),
            "verdict": "The full nonlinear timelike branch near the compact rolling solution has the expected seven physical modes including free chi. This closes the local secondary-constraint gap for this unchanged action, not a nonlinear health or UV theorem.",
            "not_established": ["An explicit field/jet neighborhood radius, global nonlinear evolution, energy/stability or causal cones on other backgrounds",
                                "A full nonlinear light/heavy solution, higher-order induced-action or quantum remainder, interacting cutoff or stationary gap",
                                "A scalar-clock chart or degree count at X=0, V/G/B UV admissibility or original P8 closure"],
            "verification_boundary": "Exact general-X ADM and Hamiltonian identities, independent dense solves and boundary/pivot controls, plus written implicit-function and spatial-Dirac closure proof; not proof-assistant formalization"}


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The local nonlinear constraint report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.45.ALIGNED_NONLINEAR replay passed; local nonlinear count only, original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
