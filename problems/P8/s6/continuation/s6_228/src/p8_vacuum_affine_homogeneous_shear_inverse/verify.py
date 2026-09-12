"""Read-only actual homogeneous curved shear plus classical tensor inverse."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_covariant_current import verify as homogeneous_input
from p8_vacuum_affine_full_finite_local_reference import verify as parent
from p8_vacuum_affine_isolated_shear_resolvent import verify as shear_input

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates" / "polynomial-vacuum-affine-homogeneous-shear-inverse.json"
)
PARENT_SHA = "48c4549159564a4ece8cfd0eb274b24d5e76296f8bd4870b1ccfc2c1ef0d6347"
SHEAR_SHA = "c7c0eb49287a4456328fe86b20a8362b70f540309439ba98b74ad1f191e0c1e6"
HOMOGENEOUS_SHA = "16b91d124384e673e03edb32e2303e5572aeac2fe106593151bfb9b6e8550d8f"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_homogeneous_shear_inverse/*.py"))
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
        (shear_input, SHEAR_SHA),
        (homogeneous_input, HOMOGENEOUS_SHA),
    ):
        if sha(prior.REPORT) != expected:
            raise ValueError(
                "A source-pinned actual current, finite reference or shear inverse changed"
            )
        prior.validate_report(
            json.loads(prior.REPORT.read_text()), prior.build_report()
        )
    return {
        "S6_227_complete_original_finite_local_reference_and_boundary_fully_rebuilt": PARENT_SHA,
        "S6_223_original_finite_window_pole_plus_cut_shear_inverse_fully_rebuilt": SHEAR_SHA,
        "S6_194_actual_homogeneous_original_prescription_current_fully_rebuilt": HOMOGENEOUS_SHA,
        "S6_216_full_dimensional_pairs_S6_217_correct_phase_S6_189_local_tree_normalization_transitively_rebuilt": True,
        "original_state_prescription_and_all_frozen_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError(
            "An actual homogeneous matching, normal-form or inverse gate failed"
        )
    return {
        "schema": 1,
        "claim": "P8-S6.228.ACTUAL_HOMOGENEOUS_TRACEFREE_CURVED_NONLOCAL_QUANTUM_AND_CLASSICAL_TENSOR_PREPARED_CAUSAL_INVERSE",
        "date": "2026-09-12",
        "status": "ACTUAL_HOMOGENEOUS_FULL_NONLOCAL_SHEAR_AND_TREE_INVERSE; NOT_NONZERO_TRANSFER_FULL_S222_NUMERICAL_SMALLNESS_STABILITY_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/leading.md",
            "notes/remainder.md",
            "notes/local-tree.md",
            "notes/primitives.md",
            "notes/inverse.md",
            "notes/regularity.md",
            "notes/scope-validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "actual_homogeneous_leading_pair_full_dimensional_matching_and_original_local_tree": serialize(
            {
                key: payload(value)
                for key, value in audit.packets().items()
                if key != "actual_weak_log_normal_form_and_prepared_inverse"
            }
        ),
        "actual_complete_weak_log_remainder_and_ordered_prepared_inverse": serialize(
            payload(audit.volterra.data())
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The actual full constrained homogeneous shear pair gives C(d)=(2d^2+d-8)/(2d(d+2)), C3=13/30 and first dimension jet91/450 in every tracefree channel. The longitudinal1/30 is retained alongside transverse2/5. The corrected positive Kubo phase gives radial coefficient208/15 and off-diagonal coefficient52/5, matching four derivatives of13/(30tau) and the original F2 high logarithm. Proper-momentum rescaling matches the highest distribution throughout the original dimension neighborhood, including its finite fourth contact. The unchanged actual all-order state, complete original contacts and fixed finite action then give the written strong homogeneous normal form I4 T_Gamma=F2(Dt^2)+V_Gamma, with finite weak-log simultaneous-derivative majorants. Adding the actual tensor tree gives T_total=T_Gamma-16pi^2 kappa L and retains all local coefficients. For K=||K2||L1(0,1) and a valid complete actual remainder majorantC, the finite weight(4KC+1)^2 yields a strict bounded-causal Neumann contraction. The ordered inverse solves every smooth prepared homogeneous tracefree source; simultaneous time-shift estimates give smoothness and both unintegrated actual inverse identities. The original shear pole and complete preparation remain. The unweighted normalized-source bound is2exp((4KC+1)^2)K times||I4g||; the canonical source contributes16pi^2 kappa, andC itself retains kappa. No numericalC orK, small physical feedback, stability or half-lineL1 shear bound is asserted. This is the actual homogeneous quantum-plus-classical time problem per comoving volume, not a reference-only mass/state substitution, but it does not solve the nonzero-transfer scalar/clock/matter graph, nonlinear/finite-coupling parent problem, cutoff or original V/G/B/P8.",
        "not_established": [
            "A uniform nonzero-transfer spatial quantum inverse or the full S222 scalar/clock/matter graph",
            "Numerical values for the complete weak-log majorants or finite-window inverse-kernel L1 norm",
            "Small unweighted feedback, absence of growing modes, physical stability or a cutoff",
            "A nonlinear metric-neighborhood or finite-coupling common-parent remainder, original V/G/B or P8 closure",
        ],
        "verification_boundary": "Exact full-dimensional pair, Abel/lag matching, original local/tree normalization, ordered primitive identities and finite-weight constants are checked independently. The actual selected-state asymptotics and singular extension, weak-log normal form, causal Neumann convergence and time-shift regularity are written continuum proofs, not FORMALIZED. Finite matrix and polynomial tests are diagnostics, not continuum numerical certificates. Native/direct/ordinary/CLI use original SymPy; only full regression uses the audited exact-GCD adapter. Frozen predecessors remain unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The actual homogeneous shear inverse report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.228 actual homogeneous shear inverse replay passed; full spatial quantum inverse and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
