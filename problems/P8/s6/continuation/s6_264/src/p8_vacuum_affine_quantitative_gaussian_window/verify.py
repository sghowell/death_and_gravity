"""Read-only complete unchanged Gaussian preparation and limited physical band certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_finite_band_neighborhood import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates/polynomial-vacuum-affine-quantitative-gaussian-window.json"
)
PARENT_SHA = "e769415fc45aadf5c0ce9d157936146f5713a3f6cc75eb45a59fef71d3d02dc0"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_quantitative_gaussian_window/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen full-parent classical-window source changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    return {
        "S6_263_complete_classical_window_and_whole_ancestry_rebuilt": PARENT_SHA,
        "same_original_full_S251_preparation_and_physical_canonical_maps": True,
        "S261_false_physical_binding_is_not_renewed": True,
        "same_original_profiles_masses_H_Proca_states_and_determinants": True,
        "reference_covariance_not_interacting_quantum_mean_or_S263_classical_state": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError(
            "A complete preparation, finite-momentum, physical-row or covariance gate failed"
        )
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.264.COMPLETE_UNCHANGED_SAMPLED_SCALAR_TENSOR_GAUSSIAN_PREPARATION_FULL_FINITE_MOMENTUM_ENERGY_TRANSPORT_AND_BAND_LIMITED_PHYSICAL_LAPSE_VOLUME_COVARIANCE_WITH_UNCHANGED_PARENT",
        "date": "2026-09-14",
        "status": "COMPLETE_FIXED_SCALAR_TENSOR_PREPARATION_AND_PHYSICAL_BAND_COVARIANCE_WITH_FINITE_Weyl_VOLUME_JET; NOT_NONLINEAR_BRANCH_SUPPORT_INTERACTING_QUANTUM_MEAN_REGULATOR_CUTOFF_NONLINEAR_INHOMOGENEOUS_B_ORIGINAL_V_G_B_OR_P8",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/energy.md",
            "notes/transition.md",
            "notes/state.md",
            "notes/physical.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_fixed_preparation_energy_and_covariance_bounds": serialize(
            {name: payload(packets[name]) for name in names[:6]}
        ),
        "whole_full_physical_rows_band_covariances_and_jet_boundary": serialize(
            {name: payload(packets[name]) for name in names[6:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The SAME complete S251 scalar/tensor sampled preparation has a quantitative covariance bound at every P>=1e64 near the bounce, derived from full profile-dependent rational interval inequalities, complete finite-q kinetic and lower terms, both original canonical chart boundaries and the actual transported preparation Gramian. The scalar weighted covariance is at most1e21 P I; each tensor block is at most1e5 P I. On the original R3 band[1e64,2e64] and |u|<=1e-60, complete physical lapse variance is below1e-490, hat log-scale variance below1e-740, absolute symmetric cross covariance below1e-615 and both tensor Frobenius variance below1e-660. Each single linear lapse observable has a Gaussian spectral tail at threshold1e-230 bounded by2 exp(-1e29). The corrected finite second-order Weyl spatial-volume contact is below1e-489. All original preparation, profile, canonical, H/Proca and determinant data remain. These are fixed-reference linear covariance and finite-jet results, not nonlinear branch support, a joint noncommuting probability law, a spacetime supremum, an interacting quantum mean, a physical regulator/cutoff, nonlinear inhomogeneous stability or closure of original V/G/B/P8.",
        "not_established": [
            "The stationary interacting physical mean, complete quantum ordering/gauge regulator or vanishing state/endpoint/Ward defects",
            "A full nonlinear volume expectation from the second-order Weyl jet, or exact support of an unbounded Gaussian inside the nonlinear auxiliary branch",
            "Joint spectral probabilities for noncommuting physical observables, a uniform spacetime supremum or infinite-volume admissibility",
            "A Wilsonian cutoff or physical UV subtraction inferred from a mathematical band restriction",
            "Full interacting H/Proca/scalar/tensor onepoint response or omitted-loop bounds from the local scalar/tensor reference covariance",
            "Identification of the fixed reference with the S263 nearby classical family or a globally varied history from the original preparation slice",
            "Nonlinear inhomogeneous Cauchy stability, original UV/gravity-limit/IR-Regge or V/G/B/P8 closure",
            "Formal verification of the written energy, Gaussian preparation, interval and probability arguments",
        ],
        "verification_boundary": "Complete source-pinned algebra, exact interval inequalities and finite matrix/ODE/quadrature diagnostics support the specified unchanged-reference covariance result. The proof does not replace finite-P control by asymptotics. Every original time boundary and both physical kappa factors remain. Native/direct/ordinary/CLI use original SymPy; only captured full P8 regression uses the separately audited exact-GCD adapter. No frozen predecessor or original state is edited; S261 remains explicitly refuted where physically misidentified.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete quantitative Gaussian window report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.264 unchanged preparation and complete physical band covariance replay passed; nonlinear support, interacting mean, cutoff and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
