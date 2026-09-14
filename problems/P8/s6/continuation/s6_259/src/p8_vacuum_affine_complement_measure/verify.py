"""Read-only full affine complement and finite canonical-jet source report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_spatial_gauge import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-complement-measure.json"
PARENT_SHA = "23053dbcf675d2737666366ee1b73e5d1a2cdac9b14405fc5d60794b3f8306cc"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_complement_measure/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen whole nonlinear spatial-gauge parent changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    return {
        "S6_258_whole_spatial_gauge_ghost_reference_and_full_ancestry_rebuilt": PARENT_SHA,
        "same_current_parent_functions_profiles_sources_and_boundaries": True,
        "finite_complement_lift_not_a_replacement_fixed_quantum_preparation": True,
        "finite_constraint_density_not_completed_quantum_ordering_continuum_measure_or_P8": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError(
            "A whole affine complement, source, finite-jet or scope gate failed"
        )
    packets = audit.packets()
    algebra = tuple(
        name
        for name in packets
        if name
        in (
            "whole_56_affine_complement_Hessian_inverse_and_determinant",
            "whole_64_projective_configuration_and_gauge_Jacobians",
            "whole_56_exact_inertia_and_continued_Fresnel_phase",
            "whole_actual_64_connection_source_center_and_retained_trace",
        )
    )
    return {
        "schema": 1,
        "claim": "P8-S6.259.COMPLETE_CURRENT_56_AFFINE_COMPLEMENT_DETERMINANT_INERTIA_AND_FINITE_CANONICAL_JET_CONDITIONAL_SOURCE_MEASURE_BRIDGE_WITH_UNCHANGED_PARENT",
        "date": "2026-09-14",
        "status": "EXACT_WHOLE_CONNECTION_COMPLEMENT_PROJECTIVE_QUOTIENT_FINITE_FRESNEL_SOURCE_AND_JET_CONSTRAINT_DENSITY; NOT_FULL_QUANTUM_ORDERING_BRST_CONTINUUM_STATE_FIXED_MEAN_ORIGINAL_V_G_B_OR_P8",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/complement.md",
            "notes/projective.md",
            "notes/fresnel.md",
            "notes/jets.md",
            "notes/sources.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_64_connection_56_complement_source_determinant_and_inertia": serialize(
            {name: payload(packets[name]) for name in algebra}
        ),
        "whole_finite_Fresnel_jet_constraint_density_and_source_contacts": serialize(
            {name: payload(packets[name]) for name in packets if name not in algebra}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The unchanged current full 64-connection parent has the explicit regular projective quotient and 56-dimensional trace complement. All nine determinant blocks, the full inverse, actual source center and exact (26,30) inertia are retained on the strict timelike coefficient domain. The conditional finite Gaussian integral has a nonconstant flat-measure factor, positive cell-scale dependence and a continued Fresnel phase that is not the principal endpoint square root. Its normalized finite configuration insertion agrees with the absence of an additional complement density after the finite canonical-jet primary/secondary delta Jacobian cancels, while retaining all remaining derivative-coordinate constraints. Original connection sources are pulled back before reduction, leaving the full nonzero source-functional inverse contact and source-dependent remaining momentum equations. The actual nonlinear source center is averaged as a composite functional, not evaluated at the retained mean. The construction does not silently assume a regular unextended velocity Legendre transform or unchanged Weyl quantization under nonlinear coordinates. It is not a choice or proof of the complete nonlinear quantum time slicing, operator ordering, state, covariant/BRST continuum measure, physical subtraction, interacting fixed mean, controlled cutoff, nonlinear bounce, UV/Regge or original V/G/B/P8 closure.",
        "not_established": [
            "A complete interacting quantum time-slicing or operator-ordering prescription and equivalence of naively unchanged nonlinear Weyl symbols",
            "A full covariant/BRST continuum gauge measure or a quantum scalar-unitary gauge through a vanishing clock gradient",
            "Deletion of the remaining spatial/translation, lapse, metric or derivative-coordinate constraints or integration of the dynamical Proca field",
            "A vanishing original connection source contact or factorization of the nonlinear source-center expectation at the retained mean",
            "A physical continuum stress/counterterm prescription, interacting state, fixed quantum mean or bound on its correction",
            "A controlled Wilsonian cutoff/matching error, physical nonlinear bounce theorem, quantum gravitational limit, UV scattering or finite-gravity IR/Regge bound",
            "Formalization of the written finite-jet, canonical constraint-density, Fresnel limiting and continuum-scope arguments",
        ],
        "verification_boundary": "Exact full connection matrices, rational congruences and independent finite-parameter diagnostics support written finite conditional and constrained-measure statements. The first-order Einstein-Hilbert measures reviewed in the literature are not imported into this sourced unrestricted connection. Native/direct/ordinary/CLI use original SymPy; only captured full regression uses the audited exact-GCD adapter. Every earlier source, profile and fixed preparation remains unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete affine complement and finite-jet report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.259 whole affine complement and finite canonical-jet source measure replay passed; full quantum theory, fixed mean and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
