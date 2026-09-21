"""Read-only complete selected scalar-loop curvature in one fixed jet convention."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_box_curvature_coefficient import verify as input4
from p8_vacuum_affine_core_curvature_coefficient import verify as input0
from p8_vacuum_affine_curvature_contact_basis import verify as input3
from p8_vacuum_affine_heavy_parent_one_loop import verify as input5
from p8_vacuum_affine_local_tadpole_radiation import verify as input1
from p8_vacuum_affine_mixed_source_radiation import verify as input2

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT
    / "certificates/polynomial-vacuum-affine-selected-matter-curvature-coefficient.json"
)
PARENT_SHA = "6b2abf4d56068636f2972200b4fc9a11411b51acd9e41612c81fd956ea7bf5e6"
LOCAL_SHA = "602dd0e80cc5fec993a20219f3423195439b0e3dba03e468794bc939b2627645"
MIXED_SHA = "faf709f6d72af441b9ff14d6ceb57fe573cf22d6b37a45d80ccf58205d5c9bc6"
BASIS_SHA = "dd49e3ecf06b03e905779e22625ad57147f3ab1d2b69f5d49896cea55338a89e"
JETS_SHA = "bcc2ef4ca3ccd84c523a0159e57a34289e99a35976a7f99201fd6384735c52bd"
SOURCE_SHA = "a09852126d79ca223af06555e10283328cb144a14945c127abee73b7ce8b9c05"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(
            ROOT.glob("src/p8_vacuum_affine_selected_matter_curvature_coefficient/*.py")
        )
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {k: v for k, v in data.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    for module, digest in (
        (input0, PARENT_SHA),
        (input1, LOCAL_SHA),
        (input2, MIXED_SHA),
        (input3, BASIS_SHA),
        (input4, JETS_SHA),
        (input5, SOURCE_SHA),
    ):
        if sha(module.REPORT) != digest:
            raise ValueError("A frozen selected scalar-loop matching input changed")
        module.validate_report(
            json.loads(module.REPORT.read_text()), module.build_report()
        )
    return {
        "S6_346_complete_old_polynomial_common_basis_coefficient_rebuilt": PARENT_SHA,
        "S6_337_full_exactD_local_tadpole_and_connection_prescription_rebuilt": LOCAL_SHA,
        "S6_338_all_three_source_classes_and_full_kernel_insertions_rebuilt": MIXED_SHA,
        "S6_343_first_real_TT_curvature_basis_and_selected_poles_rebuilt": BASIS_SHA,
        "S6_345_fixed_six_symmetrized_covariant_jets_rebuilt": JETS_SHA,
        "S6_239_full_scalar_source_inventory_and_fixed_finite_scheme_rebuilt": SOURCE_SHA,
        "S279_remains_rejected_and_archived": True,
        "S275_S276_S277_and_scoped_P8a_unchanged": True,
        "all6_historical_physical_qualifications_unchanged": True,
        "all_original_primitive_and_matching_frontiers_retained": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A selected scalar-loop coefficient proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.347.COMPLETE_SELECTED_SCALAR_LOOP_CURVATURE_IN_COMMON_JET_BASIS",
        "date": "2026-09-18",
        "status": "SCOPED_COMPLETE_SELECTED_SCALAR_LOOP_DEGREE6_CURVATURE_MATCHING; NOT_FULL_HARD_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/local.md",
            "notes/mixed.md",
            "notes/moment.md",
            "notes/aggregate.md",
            "notes/calibration.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_complete_inventory_local_and_mixed_matching": serialize(
            {name: payload(packets[name]) for name in names[:3]}
        ),
        "whole_exact_extra_and_total_known_coefficient": serialize(
            {name: payload(packets[name]) for name in names[3:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "Every selected formal scalar one-loop four-light class of the original classical limiting action is matched at homogeneous degree6 in the same six-word covariant jet lift. ExactD off-shell local tadpoles contribute zero curvature difference; all mixed-source classes give the computed chi_extra=g^2 E(n)/(16pi^2 kappa). A finite positive-weight bound proves0<chi_extra/A0<10^-604 at the original parameters. It is strictly smaller than the negative S346 core, so the complete known selected coefficient remains negative with magnitude/A0<10^-207. Independent extra parent chi, finite-gravity quantum decoupling, a full curved functional and original V/G/B/P8 closure remain open.",
        "not_established": [
            "The independent added parent-theory curvature coefficient chi or full physical matching",
            "A full curved-background effective or counterfunctional",
            "Internal-graviton loops or finite-gravity quantum decoupling",
            "A physical above-threshold amplitude approximation or local Taylor truncation-error bound",
            "Full inclusive matching, exact LSZ, quantum unitarity, complex Regge, same-parent bounce, UV or original V/G/B/P8 closure",
        ],
        "verification_boundary": "The exact graph-valence identity fixes the selected scalar inventory. All six local operators include every scalar label, metric seagull, connection and massive loop insertion in generic off-shell variables and exactD. Both mixed-source massive kernel insertions, EOM/delta terms and all three contraction classes are retained. Whole symbolic conversion identities, exact mass primitives, finite positive-integrand bounds, actual frozen kernels and original massive-vector checks determine the known coefficient. No on-shell EOM projection, floating fit or asymptotic-only bound substitutes for the calculation. All direct inputs and full source manifests are checked. Not kernel-formalized and not a finite-gravity or full nonlocal parent theorem.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete selected scalar-loop curvature report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.347 complete selected known scalar-loop common-basis curvature replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
