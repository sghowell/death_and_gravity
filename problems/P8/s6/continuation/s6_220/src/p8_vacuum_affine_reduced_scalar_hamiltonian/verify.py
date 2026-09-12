"""Read-only regular QG1 scalar comparison and complete infrared Ward pullback."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_full_spatial_remainder import verify as parent
from p8_vacuum_affine_proca_gaussian import verify as gaussian

from . import audit, estimates, pullback, scalar

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates" / "polynomial-vacuum-affine-reduced-scalar-hamiltonian.json"
)
PARENT_SHA = "df188e09486c7666c0dcf9c811702d1be096849acf993c8efae18bfffc9532a6"
GAUSSIAN_SHA = "1175702d6980db090396863b807af0e767241daeb7d8dabf138ed11c67e1ae75"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_reduced_scalar_hamiltonian/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    for previous, expected in ((parent, PARENT_SHA), (gaussian, GAUSSIAN_SHA)):
        if sha(previous.REPORT) != expected:
            raise ValueError("A source-pinned reference metric input changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_219_complete_reference_Gaussian_metric_current_fully_rebuilt": PARENT_SHA,
        "S6_176_original_Gaussian_state_transitively_rebuilt": GAUSSIAN_SHA,
        "S6_182_fixed_QG1_coefficient_not_recomputed_stress": True,
        "no_frozen_scientific_proof_test_or_report_bytes_changed": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A scalar Hamiltonian or infrared pullback gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.220.REGULAR_RETUNED_SCALAR_COEFFICIENT_HAMILTONIAN_COMPLETE_WARD_INFRARED_PULLBACK_AND_COMPARISON_BOUNDS",
        "date": "2026-09-12",
        "status": "REGULAR_QG1_CLASSICAL_COEFFICIENT_PHASE_COMPARISON_AND_IR_SAFE_REFERENCE_SCALAR_RESPONSE; NOT_FULL_QUANTUM_CONSTRAINT_INVERSE_NONLINEAR_PARENT_BACKGROUND_STABILITY_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/scalar.md",
            "notes/retuning.md",
            "notes/ward.md",
            "notes/norms.md",
            "notes/phase.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "complete_scalar_action_fixed_retuning_and_regular_Hamiltonian": serialize(
            {
                "literal_ADM": payload(scalar.spatial_data()),
                "fixed_QG1": payload(scalar.retuning_data()),
                "Hamiltonian": payload(scalar.hamiltonian_data()),
            }
        ),
        "full_ordered_IR_pullback_clock_contact_and_comparison_bounds": serialize(
            {
                "Ward": payload(pullback.ward_data()),
                "clock": payload(pullback.clock_data()),
                "scalar_norms": payload(estimates.projection_data()),
                "phase": payload(estimates.phase_data()),
            }
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The complete literal nonzero-transfer scalar ADM action retains the original matter wave, spatial curvature, lapse mixing, shift and weighted boundary terms. The current fixed S182 QG1 profile contributes DeltaJ n^2+3Tcorr n v+9A v^2/2, not only a lapse-pivot change. Its coefficient-sector Legendre reduction gives b=pv/2 and n=Lcorr/(2Jnew), with Jnew>1/100 and no division by the vanishing Theta. Both weighted scalar phase pairs are retained. In the metric response the shift reconstruction itself contains inverse momentum. Full ordered source-density/ADM-chart and synchronous-detector identities cancel that pole, leaving only gradients paired with the shift or spatial transport. Both Fourier signs, endpoint preparations and the norm-one directional projector are checked. Removing the chart leaves a nonzero inverse-square control. The actual clock metric map adds its nonzero second-jet one-current contact. The complete S219 reference Gaussian metric response consequently extends to scalar inputs(n,zeta,b) with bound1e117 V03 U138, retaining derivative loss and no unweighted shift assumption. The current coefficient-sector phase map has detector bound400 and a finite explicitly defined source coefficient norm C13, giving1e120 C13 V05 U13,10; no numerical13-jet reference-stress estimate is invented. The exact classical coefficient-sector generator is bounded by1000(1+|P|^2)^2, yielding a unique compact-momentum retarded comparison propagator through the bounce. Its growing exponential is not a uniform continuum inverse or physical EFT cutoff. Independent constraint saddle solves, full order1024 profile derivatives, three-dimensional extrinsic-curvature means, density coordinate push and detector pull tests, and nonzero deleted-contact/noncoercivity controls supplement the written arguments. Full quantum constraints remain nonlocal and are not eliminated by this classical comparison. The genuinely coupled inverse, nonlinear sourced parent, quantum background/stability, heavy/physical cutoff and original V/G/B/P8 remain open.",
        "not_established": [
            "Elimination of the full nonlocal quantum lapse/shift constraints or the genuinely coupled scalar-mixed inverse",
            "A same-space contraction from the derivative-losing reference scalar response",
            "A numerical13-jet bound for fixed reference stress beyond the original five-jet estimate",
            "A literal global homogeneous solution from a measure-zero projector assignment",
            "A uniform continuum stable propagator or physical EFT cutoff from the compact-momentum comparison",
            "A finite-amplitude sourced-parent remainder, quantum background, other/heavy-loop matching or original V/G/B/P8 closure",
        ],
        "verification_boundary": "Exact literal ADM reduction, full fixed-profile lapse-volume Hessian, canonical constraints, full ordered Ward cancellation, both metric-chart contacts, projector identities and bound arithmetic are checked with independent numeric/variational fixtures. Sobolev extension, time-primitive, smooth-coefficient and matrix Volterra/Gronwall statements are written proofs, not FORMALIZED. Native/direct/ordinary/CLI use original SymPy; the independently audited exact-GCD adapter is confined to full regression. Frozen predecessors are unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The regular scalar Hamiltonian and infrared report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.220 regular QG1 scalar comparison and infrared pullback replay passed; full quantum inverse and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
