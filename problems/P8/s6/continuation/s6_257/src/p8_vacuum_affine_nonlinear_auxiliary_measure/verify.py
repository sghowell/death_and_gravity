"""Read-only whole retained nonlinear auxiliary and finite-measure report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_finite_window_growth import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates" / "polynomial-vacuum-affine-nonlinear-auxiliary-measure.json"
)
PARENT_SHA = "260d5c682b63eee5634810a421e969c01f1ef531558d78be2147020899ed50b8"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_nonlinear_auxiliary_measure/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError(
            "The frozen whole finite-window classical comparison parent changed"
        )
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    return {
        "S6_256_whole_current_finite_window_classical_growth_and_full_ancestry_rebuilt": PARENT_SHA,
        "all_current_functions_profiles_constants_and_heavy_source_retained": True,
        "classical_comparison_not_replacement_fixed_quantum_preparation": True,
        "second_class_cancellation_not_full_covariant_quantum_measure_or_P8": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError(
            "A full nonlinear auxiliary, canonical boundary or scope gate failed"
        )
    packets = audit.packets()
    canonical_keys = (
        "whole_current_nonlinear_trace_temporal_and_lapse_Hamiltonian",
        "whole_homogeneous_full_potential_lapse_Schur_bridge",
        "whole_actual_current_functions_profiles_and_regular_classical_datum",
        "whole_spatial_metric_vector_matter_and_Gauss_Legendre_blocks",
    )
    return {
        "schema": 1,
        "claim": "P8-S6.257.COMPLETE_RETAINED_NONLINEAR_CANONICAL_AUXILIARY_BRANCH_COTANGENT_BOUNDARY_AND_FINITE_SECOND_CLASS_MEASURE_WITH_UNCHANGED_PARENT",
        "date": "2026-09-14",
        "status": "EXACT_WHOLE_CURRENT_NONLINEAR_RETAINED_CANONICAL_AUXILIARY_AND_BOUNDARY_IDENTITIES_WITH_WRITTEN_LOCAL_BRANCH_AND_FINITE_REGULATOR_SECOND_CLASS_REDUCTION; NOT_FULL_GAUGE_AFFINE_COVARIANT_QUANTUM_MEASURE_FIXED_MEAN_ORIGINAL_V_G_B_OR_P8",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/canonical.md",
            "notes/cotangent.md",
            "notes/measure.md",
            "notes/centering.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_current_nonlinear_canonical_Hamiltonian_and_local_auxiliary_branch": serialize(
            {name: payload(packets[name]) for name in canonical_keys}
        ),
        "whole_cotangent_boundary_second_class_measure_and_source_centering": serialize(
            {
                name: payload(packets[name])
                for name in packets
                if name not in canonical_keys
            }
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The unchanged full retained unitary action has the displayed spatially inhomogeneous canonical Hamiltonian with all metric, vector, matter, curvature, gradient, fixed-profile, primitive, heavy-source and Gauss terms. Its lapse/normal-vector auxiliary Hessian is regular at the actual S256 classical datum, by the exact -N U/Gamma temporal pivot and -2J homogeneous lapse Schur bridge. The implicit function theorem supplies one local canonical finite-spatial-jet root branch, not an arbitrary inhomogeneous Cauchy theorem. The full extended point cotangent lift and primitive boundary momentum shift preserve canonical Liouville volume, with their complete time Hamiltonian terms and endpoint prescription. At a fixed finite regulator the entire nonzero-secondary-bracket second-class density cancels the auxiliary delta Jacobian on one regular branch. The remaining spatial gauge factors are retained and unevaluated. Nonlinear source pullback and the exact classical reference profile residual remain explicit; the old quantum-reference Gaussian is not replaced by the classical branch. No full affine/covariant/BRST quantum measure, physical subtraction, interacting state, quantum mean, controlled cutoff, UV/Regge or original P8 closure follows.",
        "not_established": [
            "The spatial gauge-fixing determinant, its regulator-preserved first-class algebra or a complete BRST measure",
            "A nonlinear quantum determinant prescription for the affine complement and projective gauge sector",
            "Continuum determinant factorization, ordering counterterms, physical covariant subtraction or a nonlinear Gaussian momentum integral",
            "An identity between the full unforced classical auxiliary reduction and the specified off-shell reference Gaussian",
            "The fixed interacting quantum state or mean, a uniform nonlinear inhomogeneous Cauchy neighborhood, global auxiliary-root uniqueness or macroscopic bounce",
            "Controlled Wilsonian matching and loop errors, physical UV scattering, gravitational IR/Regge remainder or original V/G/B/P8 closure",
            "Formalization of the written local implicit-function, finite-regulator delta-measure and continuum-scope arguments",
        ],
        "verification_boundary": "Exact whole-source and canonical identities support the written local branch and finite-regulator second-class theorem. Independent finite-cell, full-vector/matter Legendre and nonlinear-source diagnostics are not simulations of a quantum parent or its giant physical mass. Original SymPy is used for native/direct/ordinary/CLI replay; only captured full regression uses the audited exact-GCD adapter. Earlier frozen sources, profiles and the fixed reference preparation are unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The whole nonlinear auxiliary and measure report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.257 whole nonlinear canonical auxiliary and finite second-class measure replay passed; full quantum measure, fixed mean and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
