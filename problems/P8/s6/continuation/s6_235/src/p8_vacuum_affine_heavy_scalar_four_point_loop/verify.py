"""Read-only separate first-loop four-point representation and symmetric value matching."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_heavy_scalar_one_loop import verify as parent
from p8_vacuum_affine_heavy_scalar_tree_matching import verify as model_input

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates" / "polynomial-vacuum-affine-heavy-scalar-four-point-loop.json"
)
PARENT_SHA = "d6646d25d158bcfeb32bf07895e8fcc2ca08e05988919f803ce50aa5e45f58a5"
MODEL_SHA = "31716f2b701c0120c328dd8a2506f677eb30edeca956d531c57975e38d8d25cb"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_heavy_scalar_four_point_loop/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    for prior, expected in ((parent, PARENT_SHA), (model_input, MODEL_SHA)):
        if sha(prior.REPORT) != expected:
            raise ValueError(
                "A frozen model, first-loop prescription or matching frontier changed"
            )
        prior.validate_report(
            json.loads(prior.REPORT.read_text()), prior.build_report()
        )
    return {
        "S6_234_exact_heavy_integration_complete_light_hessian_and_explicit_first_loop_scheme_fully_rebuilt": PARENT_SHA,
        "S6_233_full_separate_classical_model_tree_and_first_elastic_cut_fully_rebuilt": MODEL_SHA,
        "original_S6_177_S6_182_target_and_S6_231_S6_232_physical_frontiers_transitively_rebuilt": True,
        "all_frozen_original_and_separate_model_scientific_bytes_unchanged": True,
        "V2S_T1_OS4_is_an_explicit_new_finite_condition_not_a_frozen_prescription_edit": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError(
            "A complete loop, UV, cut, symmetric sign or contact margin gate failed"
        )
    return {
        "schema": 1,
        "claim": "P8-S6.235.SEPARATE_V2S_T1_OS4_COMPLETE_RENORMALIZED_FIRST_FOUR_POINT_LOOP_UV_AND_CUT_MATCHING_WITH_SYMMETRIC_VALUE_CONTACT_AND_CLASSICAL_QUARTIC_MARGIN",
        "date": "2026-09-12",
        "status": "SEPARATE_COMPLETE_FIRST_LOOP_REPRESENTATION_AND_SYMMETRIC_VALUE_MATCHING; NOT_FULL_ANGULAR_HIGHER_COEFFICIENT_ALL_LOOP_ORIGINAL_AFFINE_PARENT_BOUNCE_UV_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/diagrams.md",
            "notes/counterterms-cut.md",
            "notes/zero-jet.md",
            "notes/symmetric.md",
            "notes/scope.md",
            "notes/literature.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "complete_first_four_point_integrals_UV_and_elastic_cut": serialize(
            payload(audit.amplitude.data())
        ),
        "independent_zero_jet_and_separate_symmetric_matching_contact": serialize(
            {
                "zero_jet": payload(audit.zero_jet.data()),
                "symmetric": payload(audit.symmetric.data()),
            }
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The full separate-model first four-light loop is represented by three bubble channels, all triangle contributions and six mass-ordered boxes with the original first-order external normalization. Complete contact, trilinear and heavy mass UV counterterms cancel its whole pole. Its forward imaginary part reproduces the full S233 rational first elastic coefficient; the light/heavy box routings are not interchanged. A separate full constant-background Hessian verifies the distinct off-shell zero-momentum jet. At the on-shell symmetric subthreshold point s=t=u=4/3, the unadjusted S234 loop value is negative and its magnitude exceeds the positive finely cancelled tree by more than10^190. This does not establish strong coupling or a UV no-go. The separately named V2S-T1-OS4 prescription sets deltaC_fin=-A1_base at that point. Complete positive parameter bounds prove0<deltaC_fin<10^-6 times24q, with q the entire classical heavy-square remaining quartic. The contact-only classical comparison remains coercive and positive; this is not a full quantum effective-potential or all-counterterm vacuum theorem. One finite first-loop value condition does not establish the remaining physical angular error, higher quantum coefficient matching, omitted-loop bounds, a continuum quantum UV construction or original affine/common-parent bounce and V/G/B/P8 closure.",
        "not_established": [
            "A uniform physical-angle real quantum remainder or controlled higher physical b20,b21,b40 matching",
            "An all-loop error, exact quantum S matrix, continuum UV measure, controlled resonance or high-energy arc",
            "A quantum effective-potential theorem or preservation of the old bare minimum by every finite counterterm",
            "Identification with the original affine/DHOST/Proca action, covariant common-parent bounce/state control or finite-gravity Regge completion",
            "A physical cutoff, full UV-parent exclusion or original V/G/B/P8 closure",
        ],
        "verification_boundary": "Exact routed algebra, independent full Hessians, parameter Jacobians, UV pole cancellation, complete first-cut normalization and rational margins accompany written subthreshold integral and matching proofs. High-precision quadratures are diagnostics, not continuum certificates. These arguments are not FORMALIZED. The base and new finite prescriptions, off-shell jet and on-shell point, and contact-only classical comparison and full quantum potential are explicitly separated. Native/direct/ordinary/CLI retain original SymPy; only full regression uses the audited exact-GCD adapter.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The separate four-point loop report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.235 separate first-loop symmetric matching replay passed; full angular quantum matching and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
