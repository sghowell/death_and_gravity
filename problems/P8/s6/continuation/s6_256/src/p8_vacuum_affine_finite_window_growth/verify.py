"""Read-only evaluated classical finite-band growth and scope report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_classical_principal_realization import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-finite-window-growth.json"
PARENT_SHA = "b8fdddc40a30e37e4c77b509a9bb7442d139127b00f4b7cb9bb7f6f068165607"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_finite_window_growth/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen full-current classical realization parent changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    return {
        "S6_255_entire_classical_constraint_rolling_map_and_full_ancestry_rebuilt": PARENT_SHA,
        "all_current_functions_fixed_profiles_and_vacuum_constants_retained": True,
        "classical_comparison_not_replacement_fixed_quantum_preparation": True,
        "finite_subheavy_band_not_a_Wilsonian_or_original_P8_verdict": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError(
            "A whole finite-window coefficient, matrix or scope gate failed"
        )
    packets = audit.packets()
    background_keys = tuple(
        name
        for name in packets
        if name.startswith("whole_current")
        or name
        in (
            "whole_heavy_energy_and_evaluated_classical_continuation",
            "whole_initial_to_persistent_background_and_matrix_box_bridge",
        )
    )
    return {
        "schema": 1,
        "claim": "P8-S6.256.COMPLETE_CURRENT_CLASSICAL_EVALUATED_LIFETIME_AND_MASS_ADAPTED_ALL_MODE_FINITE_BAND_GROWTH_WITH_UNCHANGED_PARENT",
        "date": "2026-09-13",
        "status": "EXACT_WHOLE_CURRENT_COEFFICIENT_AND_FULL_MATRIX_ENCLOSURES_WITH_WRITTEN_EVALUATED_CLASSICAL_LIFETIME_AND_FINITE_BAND_GROWTH; NOT_WILSONIAN_CUTOFF_FIXED_QUANTUM_MEAN_NONLINEAR_BOUNCE_UV_REGGE_ORIGINAL_V_G_B_OR_P8",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/initial.md",
            "notes/lifetime.md",
            "notes/matrix.md",
            "notes/growth.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_current_function_bounds_and_evaluated_classical_lifetime": serialize(
            {name: payload(packets[name]) for name in background_keys}
        ),
        "whole_mass_adapted_physical_phase_matrix_and_finite_band_growth": serialize(
            {
                name: payload(packets[name])
                for name in packets
                if name not in background_keys
            }
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "For the unchanged current local classical action, the exact full-current constraint datum at u=0, a_hat=1, N=1+10^-6, Hhat=hbar=hbar_dot=0 and positive M1 root has an unforced solution on [0,10^-60]. Entire function and profile four-jet enclosures, the full constrained Euler system and a retained positive heavy energy prove the evaluated lifetime. For every comoving P in [10^64,2*10^64], the full eight-mode physical canonical propagator has norm at least 10^-300 exp(5*10^31). The finite-band full-matrix bound subtracts the positive massive oscillator only into a retained skew energy core, includes both TT and both transverse Proca modes, and keeps all shear and chart time connections. Both physical phase conversion norms are bounded. Physical wavenumbers and the principal growth timescale lie below the fixed heavy mass, without taking P or mass to infinity. No actual Wilsonian cutoff, omitted-operator or loop error budget, fixed quantum mean, nonlinear bounce, UV/Regge or original P8 closure follows.",
        "not_established": [
            "A controlled Wilsonian cutoff or matching bound on omitted operators, loops and nonlinear corrections",
            "Growth of the fixed interacting quantum reference mean or its replacement by the classical comparison",
            "A macroscopic physical bounce or a nonlinear instability or nonlinear Cauchy theorem",
            "All mode frequencies below the heavy mass: the retained heavy oscillator still has its mass frequency",
            "Quantum gravitational decoupling, physical UV scattering, finite-gravity IR/Regge or original V/G/B/P8 closure",
            "Formalization of the written coefficient analyticity, ODE continuation, energy cone and Fourier packet arguments",
        ],
        "verification_boundary": "Exact rational outward enclosures, full-action and matrix identities and independent finite-parameter diagnostics support the written finite-band classical theorem. Diagnostic ODEs are not numerical simulations at the giant physical mass. Original SymPy is used for native/direct/ordinary/CLI replay; only captured full regression uses the audited exact-GCD adapter. Earlier scientific files, fixed profiles and the reference preparation are unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The evaluated finite-window growth report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.256 evaluated full-mode finite-band classical growth replay passed; Wilsonian cutoff, fixed quantum mean and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
