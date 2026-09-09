"""Read-only complete second-order physical matter density and pressure."""

import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as certificate
from p8_proca_finite_width import verify as parent

from . import audit, bridge, center, coefficients, controls, observable, spatial

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "actual-quadratic-physical-matter.json"
PARENT_SHA = "b952b4fafa6666729b42dc42a38745dda60f0bee08ec66068b3cc605b7f6fc15"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_proca_physical_matter/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen finite-width test-energy certificate changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_101_fully_rebuilt": PARENT_SHA,
        "same_named_retuned_constant_Proca_classical_action_and_physical_matter": True,
        "full_physical_observable_not_identified_with_the_previous_test_energy": True,
        "no_frozen_tadpole_or_quantum_state_reassigned": True,
    }


@cache
def build_report():
    prior = prior_checks()
    identities = audit.residuals()
    gates = {key: bool(value) for key, value in audit.gates().items()}
    if not all(value is True for value in gates.values()):
        raise ValueError("A complete physical matter reconstruction proof gate failed")
    exact = certificate.certify_residuals(identities)
    co = coefficients.data()
    ce = center.data()
    return {
        "schema": 1,
        "claim": "P8-S6.102.ACTUAL_QUADRATIC_PHYSICAL_MATTER",
        "date": "2026-09-09",
        "status": "COMPLETE_ACTUAL_PHYSICAL_CLOCK_NORMAL_MATTER_DENSITY_AND_PRESSURE_SECOND_VARIATION_WITH_COMPLETE_LAPSE_EXACT_CENTER_FOURTEEN_CHANNEL_HESSIANS_AND_CANONICAL_STATE_BRIDGE; RENORMALIZED_CONSERVED_SOURCE_FULL_QUANTUM_BACKREACTION_V_G_B_AND_ORIGINAL_P8_OPEN",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/model.md",
            "notes/observable.md",
            "notes/spatial.md",
            "notes/center.md",
            "notes/canonical.md",
            "notes/scope.md",
        ],
        "exact_residuals": exact,
        "named_exact_check_count": len(identities),
        "checked_scalar_entries": sum(
            v.rows * v.cols if isinstance(v, sp.MatrixBase) else 1
            for v in identities.values()
        ),
        "proof_checks": gates,
        "actual_retuned_lapse_jet_rows": serialize(
            {
                key: co[key]
                for key in (
                    "rows",
                    "old_background_row",
                    "added_background_row",
                    "new_background_row",
                    "Jnew",
                    "literal_retuned_Hamiltonian_increment",
                )
            }
        ),
        "actual_physical_normal_observable": serialize(
            {key: value for key, value in observable.data().items() if key != "checks"}
        ),
        "actual_center_phase_hessians": serialize(
            {
                "phase_order": spatial.CHANNELS,
                "normalization": "rho2=(1/2) Q^T M_rho Q at zero output momentum, coordinate average of pointwise density, not integrated physical-volume energy",
                "matrices": {
                    key: spatial.matrices()[key]
                    for key in ("rho", "pressure", "missing_if_lapse_truncated")
                },
                "independent_closed_invariant_formulas": {
                    key: ce[key] for key in ("rho", "pressure", "n1", "n2")
                },
            }
        ),
        "actual_canonical_state_bridge_and_indefinite_principal_control": serialize(
            {
                **{
                    key: value
                    for key, value in bridge.data().items()
                    if key != "checks"
                },
                **{
                    key: ce[key]
                    for key in (
                        "center_packet_map",
                        "actual_packet_density_Hessian",
                        "principal_packet_density_Hessian",
                        "negative_principal_configuration_test",
                    )
                },
            }
        ),
        "actual_nonzero_output_strip_controls": serialize(
            {
                point: spatial.nonzero_output_control(point)
                for point in (0, sp.Rational(-1, 100), sp.Rational(1, 100))
            }
        ),
        "controls": serialize(
            {
                **controls.summary(),
                "full_classical_energy_and_Ward_boundary": {
                    key: value
                    for key, value in controls.data().items()
                    if key != "checks"
                },
            }
        ),
        "verdict": "The actual fixed physical scalar matter density and isotropic pressure are reconstructed to second order in the declared constrained spatial phase. Their linear lapse dependence forces retention of the complete second-order lapse, unlike stationary Hamiltonian reduction at quadratic order. All 105 center phase pairs and 588 independent density/pressure/missing-lapse matrix comparisons agree. The exact existing scalar canonical map retains both matter boundary shifts. Explicit negative scalar principal, tensor-momentum and Proca-momentum second-variation controls prevent an automatic positive-square test-QEI transfer; they do not refute nonnegative full classical matter energy, establish instability, or yield an arbitrary-amplitude physical-energy no-go. The full renormalized source and original P8 remain open.",
        "not_established": [
            "A complete covariant renormalized physical stress, second-order dynamical mean response, Ward-completed source or self-consistent semiclassical solution",
            "A physical-matter QEI, negative full classical energy, arbitrary-amplitude quadratic validity, instability or an ultraviolet energy no-go",
            "All local bilinear kernels from the zero-output center Hessian, or quantum locality of the fixed spatial constraint gauge",
            "Full tensor/vector quantum state construction, interactions, omitted-operator and loop errors, heavy-cutoff bounds or actual vacuum/finite-gravity/common-parent V/G/B matching",
            "Original P8(b) classification or original P8 closure",
        ],
        "verification_boundary": "Native exact full scalar normal-frame expansion, fresh literal retuned constant-Proca lapse jets, complete quadratic spatial momentum/lapse reconstruction, independent closed fourteen-channel center Hessians, exact existing density-to-packet canonical bridge including both boundary shifts, nonzero-output controls at three exact times, physical normalization and second-order Ward bookkeeping. Written geometric, implicit-branch and scope arguments are pinned, not proof-assistant formalized. The Ward bookkeeping is not a completed quantum Ward identity. Every ancestor, source/proof/test byte and report field is independently replayed; scientific SymPy is unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The actual quadratic physical matter report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.102.ACTUAL_QUADRATIC_PHYSICAL_MATTER replay passed; full source, V/G/B and original P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
