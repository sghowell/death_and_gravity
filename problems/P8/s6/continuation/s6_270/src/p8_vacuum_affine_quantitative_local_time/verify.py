"""Read-only full finite tiny-time quantum domain and same-seed comparison."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_quantitative_phase_domain import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates/polynomial-vacuum-affine-quantitative-local-time.json"
PARENT_SHA = "f4bc0f6843458556163553c192d89003388362ab22e863a9cbfd7fd542065107"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_quantitative_local_time/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen finite bounce phase-domain report changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    return {
        "S6_269_entire_finite_bounce_phase_domain_and_ancestry_rebuilt": PARENT_SHA,
        "original_source_parameters_and_pure_nonzero_preparations_unchanged": True,
        "S261_refutation_and_S265_Gaussian_integrability_boundary_unchanged": True,
        "S269_Omega_star_is_fixed_balance_not_physical_dispersion": True,
        "finite_time_regulator_comparison_not_original_continuum_matching": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError(
            "A complete finite-time domain or quantum comparison gate failed"
        )
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.270.COMPLETE_EXPLICIT_FINITE_TIME_FULL_SOURCE_NONLINEAR_DOMAIN_AND_SAME_REFERENCE_QUANTUM_DYNAMICS_WITH_EVALUATED_CUTOFF_ORDERING_POSITIVITY_STATE_LEAKAGE_AND_TWO_REGULATOR_COMPARISON_WITH_UNCHANGED_PARENT",
        "date": "2026-09-14",
        "status": "EVALUATED_FINITE_TINY_TIME_FULL_SOURCE_DOMAIN_CALIBRATED_POSITIVE_VOLUME_STATE_LEAKAGE_AND_TWO_REGULATOR_COMPARISON; NOT_ORIGINAL_QUANTUM_MATCHING_CONTINUUM_LOOP_OR_V_G_B_P8",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/reference.md",
            "notes/domain.md",
            "notes/ordering.md",
            "notes/dynamics.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_full_time_source_and_reference_phase_domain": serialize(
            {name: payload(packets[name]) for name in names[:2]}
        ),
        "whole_complete_auxiliary_and_quantum_evolution_bounds": serialize(
            {name: payload(packets[name]) for name in names[2:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The same declared finite48-pair family has an explicit full-source nonlinear domain and exact regulated quantum evolution on |u|<=1e-2000, without replacing the original state, source, free reference or Hamiltonian by bounce jets. Every primitive, moving-chart, ghost, density, implicit-volume and cutoff-derivative contact remains. The first-Weyl-calibrated coherent volume is positive with operator distance fromI<2e-255 and heat SYMBOL error<1e-310. Same-seed evolution and coherent leakage are<1e-980. The two explicit radial cutoff regulators differ in state by<1e-1970 and readout mean by<1e-1230. No original interacting mean, arbitrary-cutoff theorem, continuum/physical matching, omitted-loop estimate or nonlinear global completeness follows. Original V/G/B/P8 remain OPEN.",
        "not_established": [
            "An original interacting volume mean or unregularized singular Hamiltonian",
            "Automatic Weyl OPERATOR norm matching from the heat SYMBOL estimate",
            "Bounds for arbitrary smooth cutoffs outside the two differentiated choices",
            "Uniform mode-count, torus-volume or physical regulator-removal limits",
            "A homogeneous quantum state or literal identity with the R3 state",
            "Original Wilsonian/finite-gravity matching, omitted loops or global completeness",
            "Exact phase-space support or a sharp joint law of noncommuting observables",
            "Formal verification of the written Banach and operator arguments",
        ],
        "verification_boundary": "Complete source bindings, all mixed C5 source bounds,187 full auxiliary derivative majorants, exact physical canonical generators, full invariant images and all radial/Cauchy/Leibniz contacts support written analytic and operator proofs. Complex phase arguments use the holomorphic continuation of the real adjoint, not a conjugate adjoint for complex coefficients; smooth fixed profiles use real-time derivatives. Diagnostics do not replace full PDE/CCR/regulator proofs. Native/direct/ordinary/CLI use original SymPy; only the captured full regression uses its audited exact-GCD adapter. Frozen ancestors and both historical errata remain unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete finite-time quantum comparison report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.270 full finite-time quantum domain and regulator comparison replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
