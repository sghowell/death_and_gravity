"""Read-only full current homogeneous nonreference Gaussian principal audit."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_physical_background_vertices import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT
    / "certificates"
    / "polynomial-vacuum-affine-coupled-principal-obstruction.json"
)
PARENT_SHA = "9b7e382f5ceb147457a8ee5d6eb5a5265834a2d7499ff5ec04ab5462195e927f"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_coupled_principal_obstruction/*.py"))
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
            "The frozen complete physical-background Gaussian parent changed"
        )
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    return {
        "S6_253_complete_physical_background_vertices_and_entire_ancestry_fully_rebuilt": PARENT_SHA,
        "all_frozen_scientific_bytes_and_same_QG2_H8A420_parent_unchanged": True,
        "same_full_product_preparations_both_profiles_and_finite_conditions": True,
        "no_finite_cutoff_or_on_shell_exclusion_substituted_for_original_P8": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A complete nonreference principal or scope gate failed")
    packets = audit.packets()
    return {
        "schema": 1,
        "claim": "P8-S6.254.COMPLETE_CURRENT_HOMOGENEOUS_NONREFERENCE_FOUR_MODE_GAUSSIAN_PRINCIPAL_SYMBOL_WITH_ENTIRE_HEAVY_SOURCE_ALL_EIGHT_MODES_AND_WRITTEN_TIME_DEPENDENT_UNBOUNDED_MOMENTUM_GROWTH_OBSTRUCTION",
        "date": "2026-09-13",
        "status": "EXACT_COMPLETE_HOMOGENEOUS_GAUSSIAN_FAST_SYMBOL_AND_WRITTEN_FIXED_NONREFERENCE_GROWTH_OBSTRUCTION; NOT_FINITE_CUTOFF_ON_SHELL_NONLINEAR_QUANTUM_UV_REGGE_ORIGINAL_V_G_B_OR_P8",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/reduction.md",
            "notes/principal.md",
            "notes/growth.md",
            "notes/neighborhood.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_current_parent_heavy_source_and_four_mode_constraint_reduction": serialize(
            {
                name: payload(packets[name])
                for name in (
                    "whole_current_background_and_heavy_source",
                    "whole_four_mode_joint_constraint_reduction",
                    "whole_time_dependent_physical_central_chart",
                )
            }
        ),
        "whole_eight_mode_principal_remainder_and_time_dependent_growth_normal_form": serialize(
            {
                name: payload(packets[name])
                for name in (
                    "whole_eight_mode_fast_symbol_and_exact_remainder",
                    "independent_whole_rank_one_determinant_diagnostic",
                    "whole_time_dependent_growth_normal_form",
                )
            }
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The unchanged full current parent, including its nonzero off-clock heavy source, gives a joint four-mode physical Gaussian block. Its exact time-dependent canonical transformation and whole-matrix rational remainder give a P^(3/2) fast symbol with characteristic lambda^4(lambda^4-r^2 L2^2 Yv/(8D^2J a^4)); retaining the four other physical modes gives lambda^12 times the same quartic. A written explicit smooth fast-cycle and slow-nilpotent normal form yields exponential-in-P^(3/2) growth on each fixed nonreference compact interval with strict margins. Arbitrarily small homogeneous off-shell lapse probes of the current reference realize that domain. Thus an all-momentum finite-derivative Gaussian stability neighborhood cannot be inferred from the healthy reference. The required momentum threshold is not evaluated against the EFT cutoff; no unforced nonlinear bounce, quantum corrected system, alternative parent, UV completion or original V/G/B/P8 is thereby excluded or completed.",
        "not_established": [
            "A field-redefinition-invariant off-shell Hessian obstruction beyond the fixed S253 affine-clock fluctuation chart",
            "The growing branch lies below a controlled physical Wilsonian cutoff",
            "A small unforced on-shell nonlinear solution realizes these prescribed off-shell backgrounds",
            "The full nonlinear or quantum corrected principal operator, interacting state, physical curved subtraction or omitted-loop norm",
            "An exclusion of the current finite-cutoff bounce, every alternative parent or any UV completion",
            "Quantum gravitational decoupling, physical UV scattering, finite-gravity IR/Regge or original V/G/B/P8 closure",
            "Formalization of the written time-dependent, continuum and neighborhood proofs",
        ],
        "verification_boundary": "Exact whole-action and whole-matrix identities plus independent source, Legendre, time-connection, spectrum and evolution diagnostics support the written scoped result. No stability verdict is drawn from an unclean frozen Hamiltonian, a chosen channel, or an unevaluated comparison with the cutoff. Native/direct/ordinary/CLI use original SymPy; only complete regression uses the audited exact-GCD adapter. Every prior scientific file remains unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete nonreference principal report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.254 complete nonreference Gaussian principal replay passed; finite-cutoff and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
