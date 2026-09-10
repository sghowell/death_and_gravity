"""Read-only exceptional vacuum family and exact massive tree matching."""

import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as certificate
from p8_scalar_global_response import verify as parent

from . import (
    analytic,
    audit,
    calibration,
    decoupling,
    family,
    heavy,
    transition,
    uniform,
    vacuum,
)

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "exceptional-vacuum-tree-matching.json"
PARENT_SHA = "08e1bcdadb0df231e5ef5d854d0ee3f7da02ec996def12e7d1d8ec4e639acd8b"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_exceptional_vacuum/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key != "checks"}


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError(
            "The frozen global seven-mode relative-response report changed"
        )
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_107_fully_rebuilt": PARENT_SHA,
        "retuned_classical_coefficients_and_original_contract_rebuilt_through_ancestry": True,
        "new_named_classical_extension_does_not_edit_ancestors": True,
        "prior_absolute_counterterms_not_transferred_to_new_vacuum": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows = audit.residuals()
    gates = {key: bool(value) for key, value in audit.gates().items()}
    if not all(value is True for value in gates.values()):
        raise ValueError("An exceptional vacuum or matching proof gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.108.EXCEPTIONAL_VACUUM_TREE_MATCHING",
        "date": "2026-09-09",
        "status": "COMPLETE_DISTINCT_SMOOTH_EXACT_TUBE_AND_ANALYTIC_FINITE_JET_EXCEPTIONAL_CLASSICAL_VACUUM_FAMILIES_UNIFORM_FOUR_JET_BUDGET_FIXED_NONZERO_DHOST_CANONICAL_LIMIT_AND_EXACT_HEAVY_TREE_FOUR_POINT_MATCH_WITH_OMITTED_EXCHANGE_BOUND; NOT_FULL_V_G_B_COMMON_BOUNCE_PARENT_LOOP_OR_ORIGINAL_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/family.md",
            "notes/analytic.md",
            "notes/bounds.md",
            "notes/decoupling.md",
            "notes/amplitude.md",
            "notes/scope.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": sum(
            v.rows * v.cols if isinstance(v, sp.MatrixBase) else 1
            for v in rows.values()
        ),
        "proof_checks": gates,
        "smooth_classical_family": serialize(family.data()),
        "naive_transition_control": serialize(payload(transition.data())),
        "analytic_finite_jet_family": serialize(analytic.actual_family()),
        "actual_vacuum_jets": serialize(payload(analytic.local_jets())),
        "clock_tube_four_jet_bounds": serialize(
            {
                "curvature_and_DHOST": analytic.coefficient_error_bounds(),
                "full_scalar": analytic.lower_error_bounds(),
            }
        ),
        "uniform_co_scaled_family": serialize(payload(uniform.data())),
        "canonical_decoupling_limit": serialize(payload(decoupling.data())),
        "actual_labelled_scalar_contact": serialize(payload(vacuum.data())),
        "healthy_heavy_tree_matching": serialize(payload(heavy.data())),
        "exact_calibrations": serialize(
            {
                "first_even_order": calibration.order_point(),
                "next_even_order": calibration.order_point(1026),
                "complex_channel_remainder": calibration.heavy_remainder_bound(
                    heavy.CHANNEL_RADIUS
                ),
            }
        ),
        "controls": serialize(audit.controls()),
        "verdict": "A separately named smooth extension preserves the literal retuned classical action on an open clock tube and its necessary exceptional free-matter relation. A distinct real-analytic extension instead matches all clock X jets through n-1 with a uniform weighted C4 tube error below 10^-400 for every calibrated even order. Its explicit canonical M->infinity family retains fixed massive scalar, derivative and DHOST quartics. Independent labelled contacts match a healthy massive tree exchange through cubic channel order with an exact complex-polydisc omitted-term bound. Neither a full interacting dispersion proof nor a finite-gravity/common-bounce heavy parent is established.",
        "not_established": [
            "Scalar health on every off-clock transition, a uniform-in-n global tensor floor, a complete extended affine auxiliary dictionary or a healthy common bounce parent",
            "Equality of the analytic family with the exact open-tube action, all quantum vertices, old absolute tadpole profiles or reassignment of state/counterterms",
            "Full interacting vacuum cuts/absorptive density or V, finite-M gravity/Regge/IR contour error or G, complete common-parent B",
            "All-order loops, omitted operators beyond the exact tree-exchange remainder, nonlinear/quantum continuation, a UV completion or a universal no-go",
            "Original P8(b) classification or original P8 closure",
        ],
        "verification_boundary": "Native actual full scalar reconstruction, independently recomputed Ia completion and switch derivatives; exact vacuum jets, finite Taylor norm bounds, affine-in-order global envelopes, canonical co-scaling and Lorentzian labelled contractions; exact massive exchange and complex remainder identity. Continuous analytic, identity-theorem, limit and scope arguments are written/source-pinned, not proof-assistant formalized or peer reviewed. All ancestors, own source/proof/test bytes and report fields are read-only replayed with unmodified scientific SymPy.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The exceptional vacuum tree-matching report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.108.EXCEPTIONAL_VACUUM_TREE_MATCHING replay passed; full V/G/B, common parent and original P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
