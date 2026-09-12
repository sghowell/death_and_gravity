"""Read-only original scalar partial unitarity and conditional matching disjunction."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_full_flat_tensor_response import verify as parent
from p8_vacuum_affine_quantum_retuning import verify as vacuum_input
from p8_vacuum_canonical_affine_decoupling import verify as tree_input

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-scalar-matching-scale.json"
PARENT_SHA = "9c2f9c6798e8b4a93af0343e8b0436abaea7428e2800d402020e6b25c1e20fc8"
TREE_SHA = "b85c821b59bc41489418c0c190a2cc573cd9852cb19254627c3e2150d1dd6446"
VACUUM_SHA = "f754f6d62795d5a1c0ade4c576395f1587fc35f9c9895bc8b7b15b2b2dec2e16"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_scalar_matching_scale/*.py"))
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
        (tree_input, TREE_SHA),
        (vacuum_input, VACUUM_SHA),
    ):
        if sha(prior.REPORT) != expected:
            raise ValueError(
                "A source-pinned original scalar tree, retained vacuum or full tensor comparison changed"
            )
        prior.validate_report(
            json.loads(prior.REPORT.read_text()), prior.build_report()
        )
    return {
        "S6_230_original_complete_flat_tensor_response_and_unchanged_frontier_fully_rebuilt": PARENT_SHA,
        "S6_177_actual_full_function_canonical_classical_decoupling_and_tree_amplitude_fully_rebuilt": TREE_SHA,
        "S6_182_original_retained_vacuum_and_unchanged_low_retuning_jets_fully_rebuilt": VACUUM_SHA,
        "original_state_prescription_and_all_frozen_scientific_bytes_unchanged": True,
        "no_old_GY14_quantum_amplitude_matching_imported": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError(
            "An original scalar normalization, unitarity or conditional matching gate failed"
        )
    return {
        "schema": 1,
        "claim": "P8-S6.231.ORIGINAL_MASSIVE_SCALAR_IDENTICAL_CHANNEL_UNITARITY_NECESSARY_CORRECTION_AND_CONDITIONAL_FULL_AMPLITUDE_DISPERSION_MATCHING_ERROR_TRADEOFF",
        "date": "2026-09-12",
        "status": "EXACT_ORIGINAL_TREE_NECESSARY_UNITARITY_CORRECTION_AND_CONDITIONAL_MATCHING_DISJUNCTION; NOT_PHYSICAL_CUTOFF_QUANTUM_DECOUPLING_FINITE_GRAVITY_UV_NO_GO_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/tree.md",
            "notes/normalization.md",
            "notes/unitarity.md",
            "notes/dispersion.md",
            "notes/scale-scope.md",
            "notes/literature.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "complete_original_tree_identical_partial_waves_and_first_elastic_coefficient": serialize(
            payload(audit.amplitude.data())
        ),
        "necessary_unitarity_correction_and_conditional_physical_matching_tradeoff": serialize(
            {
                "unitarity": payload(audit.unitarity.data()),
                "dispersion": payload(audit.dispersion.data()),
            }
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The same S182 retuned parent has the complete S177 fixed-canonical nongravitational massive scalar tree with mass1, lambda10^-600 and gamma1024*10^-800; no low retuning jet or higher independent vertex alters this connected four-point tree. The labelled invariant amplitude, exactly one final identical-particle factor and both normalized identical-channel factors give t0=beta(a+b/3)/(32pi), t2=beta*b/(240pi), and the full first elastic coefficient beta[(a+b/3)^2+4b^2/45]/(32pi). The exact S=1+2it circle, including positive inelastic channels, gives |Re t0|<=1/2. The original positive t0 is strictly increasing and dominates all its other tree waves; its unique nominal real-part ceiling occurs at COM energy between10^133 and2*10^133. At10^134 exact unitarity requires relative full t0 correction greater than1-10^-5. A real nonzero tree is not itself exactly unitary even below that ceiling. Separately, an actual scalar S matrix satisfying the explicitly unproved crossing, pole, positivity and arc conditions and a full angular L2 error eta on the annulus[E^2/2,E^2] must have physical b2/(4lambda)>9(1-eta)^2(E/10^125)^8, provided E^2/2>=16. Hence a relative b2 tolerance delta requires1+delta>9(1-eta)^2(E/10^125)^8. AtE10^125, eta<=1/2 and b2<=8lambda contradict the original elastic annulus lower bound b2>9lambda, despite an original tree amplitude below10^-46. At least one premise or tolerance must fail; this does not select which. The fixed S230 tensor-pole frequency lower10^398 is in the same anchored units but supplies no quantum decoupling control, physical cutoff or tensor-pole repair. All loops, heavy states and higher operators remain in the unproved full-amplitude matching error. No old finite-amplitude quantum matching, finite-gravity partial-wave bound, physical UV exclusion, nonlinear bounce control or original V/G/B/P8 closure is inferred.",
        "not_established": [
            "Existence and required analytic/arc properties of the actual nongravitational quantum S matrix, or a controlled full quantum gravitational decoupling limit",
            "A bound on physical b2 matching, the full angular amplitude error, real loop counterterms, higher loops, heavy thresholds or omitted operators",
            "A physical cutoff, validity of the conditional tensor equation at its complex poles, a tensor-pole modification or a UV-parent no-go",
            "Finite-gravity forward or Regge control, unrestricted curved S222 graph estimates, nonlinear common-parent bounce control, original V/G/B or P8 closure",
        ],
        "verification_boundary": "Original massive tree, Legendre projections, identical bubble normalization, exact circle, positive mass/angle polynomials, rational scale margins and dispersion coefficient are replayed independently. The continuum optical theorem, partial-wave projection, analytic dispersion implication and full-error norm argument are written proofs with explicit premises, not FORMALIZED. No computed exact zero proves the unestablished S-matrix or matching hypotheses. Native/direct/ordinary/CLI use original SymPy; only full regression uses the audited exact-GCD adapter. Frozen predecessors remain unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The original scalar conditional matching report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.231 original scalar unitarity and conditional matching replay passed; physical cutoff and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
