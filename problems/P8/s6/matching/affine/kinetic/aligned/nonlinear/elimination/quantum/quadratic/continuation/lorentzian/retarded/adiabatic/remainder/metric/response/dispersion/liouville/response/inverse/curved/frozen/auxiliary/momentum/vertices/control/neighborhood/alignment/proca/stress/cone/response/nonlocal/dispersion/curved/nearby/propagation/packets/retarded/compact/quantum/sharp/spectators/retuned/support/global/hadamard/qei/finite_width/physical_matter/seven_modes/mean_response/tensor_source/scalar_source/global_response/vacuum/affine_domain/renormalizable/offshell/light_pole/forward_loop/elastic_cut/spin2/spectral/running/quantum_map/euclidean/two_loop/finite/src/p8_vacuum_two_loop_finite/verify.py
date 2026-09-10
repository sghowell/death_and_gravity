"""Read-only integrated compact bound for the actual UV-finite two-loop subset."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_two_loop_denom import verify as parent

from . import audit, calibration, ordered, parametric, sectors, selection

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-two-loop-finite-subsector.json"
PARENT_SHA = "7bfc16a44e84a29b94ae4eeeef71b3529d62847f7e88dd13d2bb6ab3e53c5a6b"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_two_loop_finite/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(d):
    return {k: v for k, v in d.items() if k != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen complete two-loop denominator report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_120_fully_rebuilt": PARENT_SHA,
        "same_full_heavy_graphs_and_actual_parameters": True,
        "same_compact_denominator_and_crossing_derivative_bound": True,
        "same_external_assignment_Wick_weights": True,
        "no_subtraction_or_LSZ_completion_inferred": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows = audit.residuals()
    gates = {k: bool(v) for k, v in audit.gates().items()}
    if not all(v is True for v in gates.values()):
        raise ValueError("An integrated finite-subsector gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.121.POLYNOMIAL_VACUUM_TWO_LOOP_UV_FINITE_SUBSECTOR",
        "date": "2026-09-09",
        "status": "COMPLETE_ACTUAL_INTEGRATED_UV_FINITE_TWO_LOOP_FOUR_POINT_SUBSECTOR_COMPACT_FORWARD_BOUND_EXACT_HEPP_SECTORS_HEAVY_MASS_LOGARITHMS_AND_ABSOLUTE_ERROR_BELOW_ONE_E_MINUS_7_OF_TREE; NOT_SUBTRACTION_DEPENDENT_GRAPHS_COMPLETE_TWO_LOOP_LSZ_HIGHER_LOOP_HIGH_ENERGY_GRAVITY_V_G_B_OR_ORIGINAL_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/selection.md",
            "notes/sectors.md",
            "notes/parametric.md",
            "notes/bounds.md",
            "notes/literature.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": len(rows),
        "proof_checks": gates,
        "exact_finite_subsector_selection": serialize(payload(selection.data())),
        "all_anisotropic_box_sector_summaries": serialize(sectors.summaries()),
        "independent_elementary_and_Schwinger_integrals": serialize(
            {
                "ordered": payload(ordered.data()),
                "parametric": payload(parametric.data()),
            }
        ),
        "actual_coupling_weighted_integrated_bound": serialize(
            {
                "actual": payload(calibration.data()),
                "diagnostic_first_finite_box": calibration.polynomial(
                    "double_bubble", (0, 0, 2), 600
                ),
            }
        ),
        "controls": serialize(audit.controls()),
        "verdict": "All 88 individually UV-finite bare two-loop four-point refinements have absolutely integrable compact forward-disc majorants. Exact ordered sectors remove heavy mass powers and bound the remaining mass dependence by a quadratic logarithmic polynomial. With the fixed Gaussian loop normalization, Taylor factor, actual couplings and all external assignments, the sum has absolute forward-coefficient error below 3 times 10^-607 and below 10^-7 of the positive tree value. The signed correction is not determined. The 104 subtraction-dependent refinements, counterterm insertions and complete two-loop light pole/LSZ terms remain outside this bound.",
        "not_established": [
            "The signed or exact value of the UV-finite subset of two-loop amplitudes",
            "Integrated bounds for the 104 subtraction-dependent refinements or fixed finite counterterm insertions",
            "Complete two-loop light mass residue LSZ or forward coefficient",
            "All higher-loop errors all-energy complex contours or an ultraviolet completion",
            "Finite-gravity IR/Regge contour Delta_grav or a common controlled cosmological parent",
            "Exclusion or completion of an entire original P8 row full V/G/B or original P8 closure",
        ],
        "verification_boundary": "Native independent co-tree polynomials, every edge order and heavy-cutoff piece, exact mass-power cancellation, positive exponent and logarithmic bounds, independent anchored ordered integrals, Gaussian/Schwinger normalization and actual rational coupling sums. Continuous Tonelli domination, first-sheet analytic continuation and differentiation are explicit source-pinned written proofs, not proof-assistant formalized, nonperturbative or peer reviewed. Every ancestor, own source and report field is rebuilt read-only with unmodified scientific SymPy.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The integrated finite-subsector report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.121.POLYNOMIAL_VACUUM_TWO_LOOP_UV_FINITE_SUBSECTOR replay passed; full V/G/B and original P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
