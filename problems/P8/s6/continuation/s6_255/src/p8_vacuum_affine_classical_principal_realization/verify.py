"""Read-only current local classical realization and on-shell principal report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_coupled_principal_obstruction import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT
    / "certificates"
    / "polynomial-vacuum-affine-classical-principal-realization.json"
)
PARENT_SHA = "34adb63ad201a137c2b510e7b7df9eb3c7d717f44a79089dcbcb211594c2f371"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_classical_principal_realization/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen entire four-mode principal parent changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    return {
        "S6_254_entire_current_principal_remainder_growth_and_full_ancestry_rebuilt": PARENT_SHA,
        "same_full_current_local_parent_and_every_frozen_scientific_byte": True,
        "classical_comparison_not_replacement_fixed_quantum_preparation": True,
        "no_finite_cutoff_quantum_mean_or_original_P8_closure_substituted": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError(
            "A current classical realization or on-shell scope gate failed"
        )
    packets = audit.packets()
    return {
        "schema": 1,
        "claim": "P8-S6.255.COMPLETE_CURRENT_HOMOGENEOUS_CLASSICAL_COMPARISON_CONSTRAINT_DATA_LOCAL_UNFORCED_REALIZATION_AND_ROLLING_HEAVY_ON_SHELL_PRINCIPAL_OBSTRUCTION_WITH_UNCHANGED_PARENT",
        "date": "2026-09-13",
        "status": "EXACT_CURRENT_CLASSICAL_CONSTRAINT_AND_ROLLING_BLOCK_WITH_WRITTEN_LOCAL_UNFORCED_REALIZATION_AND_ON_SHELL_UNBOUNDED_MOMENTUM_GROWTH; NOT_FIXED_QUANTUM_MEAN_CUTOFF_MACRO_BOUNCE_UV_REGGE_ORIGINAL_V_G_B_OR_P8",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/rolling.md",
            "notes/constraint.md",
            "notes/existence.md",
            "notes/principal.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_current_homogeneous_action_constraint_and_rolling_heavy_block": serialize(
            {
                name: payload(packets[name])
                for name in (
                    "whole_homogeneous_classical_action_and_regular_Euler_system",
                    "whole_rolling_heavy_Gaussian_block_and_Euler_contacts",
                )
            }
        ),
        "whole_current_local_classical_realization_and_on_shell_field_map_bridge": serialize(
            {
                name: payload(packets[name])
                for name in (
                    "whole_current_classical_constraint_data_and_local_realization",
                    "whole_on_shell_field_map_contact_bridge",
                )
            }
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The entire current local action has explicit nearby homogeneous classical constraint data with strict lapse and full reduced-velocity margins. The written implicit-lapse and full Euler argument gives local unforced comparison solutions with the rolling heavy source retained. Their complete four-mode Gaussian map keeps all shift, source and time contacts; the lower volume contacts vanish by the actual classical Euler equations, not by holding the heavy field artificially zero. The whole S254 eight-mode fast symbol and time-dependent growth proof survive on compact nonreference intervals of those solutions. On-shell field-map onepoint contacts vanish, so regular polynomially bounded physical phase changes do not remove that linearized unbounded-momentum obstruction. These classical data are not the fixed prepared quantum mean, the required momenta are not compared with a controlled cutoff, and no macroscopic physical bounce, UV/Regge or original V/G/B/P8 verdict follows.",
        "not_established": [
            "The fixed original quantum reference is a solution of the full mean equations or is replaced by this classical family",
            "A physical-scale bounce or a uniform-heavy-mass or macroscopic existence interval",
            "A growing branch below a controlled Wilsonian cutoff or absence of every finite-cutoff stability neighborhood",
            "The full quantum principal operator, physical curved subtraction, omitted-loop norm or interacting state",
            "Quantum gravitational decoupling, physical UV scattering, finite-gravity IR/Regge or original V/G/B/P8 closure",
            "Formalization of the written implicit-function, local Euler, Noether, continuum or on-shell functional proofs",
        ],
        "verification_boundary": "Exact complete-action and matrix identities, source-pinned full-current central bounds and independent finite diagnostics support a local classical comparison theorem. Diagnostic finite-parameter ODEs are not simulations of the giant physical mass or evaluated physical lifespan/cutoff bounds. Native/direct/ordinary/CLI use original SymPy; only captured full regression uses the audited exact-GCD adapter. Every earlier scientific file and fixed reference preparation is unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The current classical principal realization report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.255 current local classical realization replay passed; fixed quantum mean, cutoff and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
