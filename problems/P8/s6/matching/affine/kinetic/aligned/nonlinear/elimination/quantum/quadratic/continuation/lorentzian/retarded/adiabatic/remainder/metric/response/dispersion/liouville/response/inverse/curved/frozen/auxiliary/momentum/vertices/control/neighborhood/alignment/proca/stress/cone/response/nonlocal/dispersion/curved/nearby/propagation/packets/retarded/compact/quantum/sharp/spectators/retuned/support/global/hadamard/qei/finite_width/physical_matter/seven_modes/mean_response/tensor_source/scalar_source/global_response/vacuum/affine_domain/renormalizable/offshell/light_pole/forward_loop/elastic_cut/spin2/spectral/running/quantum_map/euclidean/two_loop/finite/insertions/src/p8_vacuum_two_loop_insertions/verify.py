"""Read-only grouped full self-energy insertion bound with fixed subtractions."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_two_loop_finite import verify as parent

from . import audit, calibration, kernel, radial, routing, subtraction

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates" / "polynomial-vacuum-two-loop-self-energy-insertions.json"
)
PARENT_SHA = "d14b534b1f033896685b559e27e61777eb4ce45e75f3960a27102739c75d3e1f"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_two_loop_insertions/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen integrated finite-subsector report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_121_fully_rebuilt": PARENT_SHA,
        "same_actual_couplings_full_heavy_vertices_and_first_sheet": True,
        "same_once_fixed_inner_OS_and_inherited_outer_reference": True,
        "disjoint_finite_group_retained_without_double_counting": True,
        "no_complete_two_loop_pole_LSZ_or_P8_closure_inferred": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows = audit.residuals()
    gates = {k: bool(v) for k, v in audit.gates().items()}
    if not all(v is True for v in gates.values()):
        raise ValueError("A grouped self-energy insertion gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.122.POLYNOMIAL_VACUUM_TWO_LOOP_SELF_ENERGY_INSERTIONS",
        "date": "2026-09-09",
        "status": "COMPLETE_ACTUAL_GROUPED_SIXTY_FOUR_SELF_ENERGY_INSERTION_REFINEMENTS_FIXED_INNER_OS_LINEAR_INHERITED_OUTER_REFERENCE_COMPLEX_STRIP_DECAYING_REMAINDER_AND_INTEGRATED_COMPACT_FORWARD_BOUND; NOT_REMAINING_FORTY_REFINEMENTS_COMPLETE_TWO_LOOP_POLE_LSZ_HIGHER_LOOP_GRAVITY_V_G_B_OR_ORIGINAL_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/kernel.md",
            "notes/routing.md",
            "notes/subtraction.md",
            "notes/bounds.md",
            "notes/literature.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": len(rows),
        "proof_checks": gates,
        "grouped_topology_and_fixed_local_subtractions": serialize(
            payload(subtraction.data())
        ),
        "actual_complex_kernel_and_constant_remainder_split": serialize(
            payload(kernel.data())
        ),
        "all_three_channel_complex_momentum_routings": serialize(
            payload(routing.data())
        ),
        "convergent_radial_integral_and_actual_partial_bound": serialize(
            {
                "radial": payload(radial.data()),
                "actual": payload(calibration.data()),
                "zero_radial_enclosure": calibration.point(0),
                "unit_radial_enclosure": calibration.point(1),
            }
        ),
        "controls": serialize(audit.controls()),
        "verdict": "All 64 raw tadpole-insertion refinements are grouped with the same fixed inner on-shell subtraction. The actual relative self-energy multiplier splits into its fixed asymptotic constant and an absolutely integrable complex-strip remainder. Linear continuation of the inherited full regulated outer reference leaves twice the fixed multiplier times the same renormalized one-loop coefficient plus the explicitly bounded remainder. The absolute group error is below 10^-614 and below 2 times 10^-15 of the positive tree coefficient. Together with the disjoint 88 finite refinements, the partial sum remains below 3 times 10^-607 and below 10^-7 of tree. This is not the signed correction or the complete two-loop coefficient.",
        "not_established": [
            "The signed or exact value of this renormalized diagram family",
            "Integrated bounds and consistent local subtractions for the other 40 raw refinements",
            "Complete two-loop local counterterm accounting or light pole residue LSZ control",
            "A two-loop source-aware derivative-composite map or transfer to the truncated target",
            "All higher-loop errors all-energy contours finite-gravity IR/Regge Delta or a common controlled bounce parent",
            "Full V/G/B classification any excluded original P8 row or original P8 closure",
        ],
        "verification_boundary": "Native actual-kernel identities, explicit all-channel bilinear routing, grouped topology and local subtraction identities, anchored convergent radial integral and exact rational actual error comparisons. Continuous contour continuation, absolute domination, Cauchy estimates and renormalization grouping are source-pinned written proofs, not proof-assistant formalized, nonperturbative or peer reviewed. Every ancestor, own source and report field is rebuilt read-only with unmodified scientific SymPy.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The grouped insertion report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.122.POLYNOMIAL_VACUUM_TWO_LOOP_SELF_ENERGY_INSERTIONS replay passed; full V/G/B and original P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
