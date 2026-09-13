"""Read-only complete conditional source and formal loop-filtration report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_heavy_parent_one_loop import verify as finite_input
from p8_vacuum_affine_heavy_scalar_parent import verify as source_input
from p8_vacuum_affine_two_mass_nonscalar_inverse import verify as parent

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-heavy-source-filtration.json"
PARENT_SHA = "619218777bc9d5e525a04d7276af6fa9fcf1f46826d7badcf35243859ba48b22"
SOURCE_SHA = "934b6c035fb5de0eb89ee1e1a03c2624c78a33985c0503bce883c41240b05905"
FINITE_SHA = "a09852126d79ca223af06555e10283328cb144a14945c127abee73b7ce8b9c05"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_heavy_source_filtration/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    for prior, expected in (
        (parent, PARENT_SHA),
        (source_input, SOURCE_SHA),
        (finite_input, FINITE_SHA),
    ):
        if sha(prior.REPORT) != expected:
            raise ValueError(
                "A frozen current conditional source or finite-extension input changed"
            )
        prior.validate_report(
            json.loads(prior.REPORT.read_text()), prior.build_report()
        )
    return {
        "S6_249_current_conditional_reference_and_complete_ancestry_fully_rebuilt": PARENT_SHA,
        "S6_238_whole_source_and_canonical_parent_fully_rebuilt": SOURCE_SHA,
        "S6_239_complete_stated_finite_extension_and_vacuum_loop_boundary_fully_rebuilt": FINITE_SHA,
        "all_frozen_original_and_separate_model_scientific_bytes_unchanged": True,
        "missing_curved_counterfunctional_not_inferred_from_flat_jet_prescription": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A complete source influence or formal filtration gate failed")
    packets = audit.packets()
    return {
        "schema": 1,
        "claim": "P8-S6.250.COMPLETE_CURRENT_CONDITIONAL_HEAVY_SOURCE_CTP_INFLUENCE_CANONICAL_CLOCK_JETS_PHYSICAL_FORCED_KG_BOUND_AND_FORMAL_COUNTERGRADED_SOURCE_GRAPH_LOWER_BOUNDS_WITH_MISSING_CURVED_RENORMALIZATION_RETAINED",
        "date": "2026-09-13",
        "status": "EXACT_CONDITIONAL_GAUSSIAN_SOURCE_INFLUENCE_AND_FORMAL_ANCESTRY_PRESERVING_LOOP_FILTRATION; NOT_COMPLETE_CURVED_COUNTERFUNCTIONAL_INTERACTING_STATE_NONZERO_RENORMALIZED_COEFFICIENT_OMITTED_LOOP_BOUND_NONLINEAR_UV_REGGE_ORIGINAL_V_G_B_OR_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/influence.md",
            "notes/graphs.md",
            "notes/counterterms.md",
            "notes/remaining.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_conditional_source_and_influence": serialize(
            {
                name: payload(data)
                for name, data in packets.items()
                if name != "full_countergraded_source_graph_filtration"
            }
        ),
        "formal_filtration_and_counterterm_boundary": serialize(
            payload(packets["full_countergraded_source_graph_filtration"])
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The complete physical source of the specified current conditional H action, with its fixed finite onepoint extension, has a canonical clock jet of order8 and coefficient t²exp(-A)(pi_dot-h00)^8/(64kappa³). Exact prepared Gaussian H integration retains both CTP density sources, their complete contour inverse and the full source-independent metric determinant. Its source influence starts degree16, causal coherent H response8, metric/light current15 and first response14, with explicit preparation and spatial-projection boundaries. The complete physical forced KG equation on the unchanged FLRW metric has an all-momentum energy bound with factor15T/sqrt(n). A formal connected countergraded graph identity gives E+2Ltotal>=16 for generic source-dependent nonheavy kernels, while a one-heavy mean first becomes topologically permitted at loop4. Complete source-ancestor countergraph contraction preserves weighted degree. This is conditional on background-jet and source-ancestry-preserving subtraction, not a construction of the missing curved counterfunctional or interacting initial state. The full Gaussian determinant and source-independent light/metric/M1 and mixed gravitational loops remain; graph bounds do not determine nonzero renormalized coefficients or omitted-loop sizes. No evaluated coupled nonlinear feedback, same-state bounce, quantum gravitational limit, physical UV/Regge or original V/G/B/P8 closure follows.",
        "not_established": [
            "A complete curved light/gravity/source counterfunctional, renormalized interacting clock amplitudes or interacting initial state",
            "Nonzero projected or renormalized graph coefficients, a convergent loop series or numerical omitted-loop bounds",
            "Removal of source-independent Gaussian, light/metric/M1, gauge or mixed gravitational quantum effects",
            "A bound for the complete nonlinear coupled feedback, same-state bounce, physical smallness or stability",
            "A physical gravitational quantum limit, exact UV S matrix, finite-gravity IR/Regge estimate or all-parent exclusion",
            "Original V/G/B/P8 closure or machine formalization of the continuum arguments",
        ],
        "verification_boundary": "Exact full-source, canonical-metric, complete CTP, density/operator variation, physical KG and countergraded graph identities support the written conditional Gaussian and formal filtration proofs. Independent unexpanded metric/source, ordered oscillator, Gaussian contraction, full graph-edge, counterterm ancestry and energy diagnostics test the formulas without constructing the missing interacting theory. Native/direct/ordinary/CLI use original SymPy; only full regression uses the audited exact-GCD adapter. All frozen scientific bytes remain unchanged. These arguments are not FORMALIZED.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete source influence and formal filtration report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.250 conditional source influence and formal filtration replay passed; interacting physics and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
