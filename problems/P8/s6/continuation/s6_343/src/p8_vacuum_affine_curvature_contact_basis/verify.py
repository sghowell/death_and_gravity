"""Read-only first curvature contact basis and selected-matter pole closure."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_factorized_bubble_radiation import verify as input4
from p8_vacuum_affine_heavy_parent_one_loop import verify as input6
from p8_vacuum_affine_local_tadpole_radiation import verify as input2
from p8_vacuum_affine_mixed_source_radiation import verify as input3
from p8_vacuum_affine_quadratic_radiation_cancellation import verify as input5
from p8_vacuum_affine_radiative_curvature_matching import verify as input1
from p8_vacuum_affine_triangle_box_radiation import verify as input0

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates/polynomial-vacuum-affine-curvature-contact-basis.json"
PARENT_SHA = "b73d47b8bb15941ebb82540312d1a681cea47295c0cef7f993aae7174aad1d04"
CURVATURE_SHA = "f8f9cd57c0240fe6008efd743d399c27dfc393a14d50501e8362289aae204db0"
LOCAL_SHA = "602dd0e80cc5fec993a20219f3423195439b0e3dba03e468794bc939b2627645"
MIXED_SHA = "faf709f6d72af441b9ff14d6ceb57fe573cf22d6b37a45d80ccf58205d5c9bc6"
BUBBLE_SHA = "174317dc92bc5b1c2956e6c5bebee4504f582325ae7e28b16f03d5de01f70206"
QUADRATIC_SHA = "79f298f3b8db27f6df41081e2f0ab80fd13eab8bc149699196c24eef4fda6cb7"
SOURCE_SHA = "a09852126d79ca223af06555e10283328cb144a14945c127abee73b7ce8b9c05"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_curvature_contact_basis/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    for module, digest in (
        (input0, PARENT_SHA),
        (input1, CURVATURE_SHA),
        (input2, LOCAL_SHA),
        (input3, MIXED_SHA),
        (input4, BUBBLE_SHA),
        (input5, QUADRATIC_SHA),
        (input6, SOURCE_SHA),
    ):
        if sha(module.REPORT) != digest:
            raise ValueError("A frozen curvature-basis or selected-pole input changed")
        module.validate_report(
            json.loads(module.REPORT.read_text()), module.build_report()
        )
    return {
        "S6_342_complete_selected_matter_radiation_frontier_rebuilt": PARENT_SHA,
        "S6_336_literal_curvature_contact_and_unassigned_matching_rebuilt": CURVATURE_SHA,
        "S6_337_complete_exact_D_local_TT_counterfunctional_rebuilt": LOCAL_SHA,
        "S6_338_all_mixed_source_radiation_and_pole_identities_rebuilt": MIXED_SHA,
        "S6_339_full_bubble_and_three_counterterm_radiation_rebuilt": BUBBLE_SHA,
        "S6_340_entire_selected_quadratic_external_cancellation_rebuilt": QUADRATIC_SHA,
        "S6_239_full_original_scalar_valence_inventory_rebuilt": SOURCE_SHA,
        "S279_remains_rejected_and_archived": True,
        "S275_S276_S277_and_scoped_P8a_unchanged": True,
        "all6_historical_physical_qualifications_unchanged": True,
        "all_original_primitive_and_matching_frontiers_retained": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("A curvature basis or selected-pole proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.343.FIRST_CURVATURE_CONTACT_BASIS_AND_SELECTED_MATTER_POLE_CLOSURE",
        "date": "2026-09-18",
        "status": "SCOPED_FIRST_CURVATURE_BASIS_AND_SELECTED_MATTER_TT_POLES; NOT_FULL_HARD_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/basis.md",
            "notes/poles.md",
            "notes/matching.md",
            "notes/calibration.md",
            "notes/scope.md",
            "notes/provenance.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_basis_and_selected_TT_poles": serialize(
            {name: payload(packets[name]) for name in names[:3]}
        ),
        "whole_strict_matching_and_original_calibrations": serialize(
            {name: payload(packets[name]) for name in names[3:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "At the first possible parity-even explicit-curvature derivative order, the local four-identical-scalar one-real-graviton contact has exactly one independent direction: the unassigned S336 chi*T. Full256-word Bose reduction and literal original massive-state contractions establish the stated basis. The selected minimal-matter radiative UV poles cancel with the fixed literal counterfunctionals, so no additional curvature-only physical-TT pole is needed in that representative. Neither fact assigns or bounds finite chi. The entire polynomial must match before a coefficient is extracted; higher derivatives, internal gravity, full curved matching, inclusive quantum probability and original V/G/B/P8 remain open.",
        "not_established": [
            "The finite value or a bound for the actual parent curvature coefficient chi",
            "A higher-derivative, parity-odd, off-shell or full curved-action basis",
            "Internal-graviton loops, a full gravity beta function or finite-gravity quantum decoupling",
            "A bound on an unknown curvature addition from the known minimal-matter remainder",
            "Full virtual or inclusive matching, exact LSZ, quantum unitarity, complex Regge, same-parent bounce, UV or original V/G/B/P8 closure",
        ],
        "verification_boundary": "Exact generic constrained a_i,H_ij Bose reduction of all256 ordered curvature words, symbolic rank1, all null-TT curvature traces and k-slots, plus independent literal256-component contractions on the unchanged original massive states. The complete exact-D local Wick and metric identities, mixed-source and bubble counterterms, UV-finite minimal triangle/box Gamma factors and quadratic cancellations are recomputed. Written index counting supplies the finite derivative-order completeness argument; tests do not replace that proof or determine unknown matching. Not kernel-formalized; all frozen ancestors unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete curvature-basis report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.343 first curvature basis and selected-matter poles replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
