"""Read-only original-channel logarithmic bound and weighted reference inverse."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_curved_scalar_reference import verify as parent
from p8_vacuum_affine_isolated_shear_resolvent import verify as shear_input

from . import analytic, audit, curved, weighted

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates" / "polynomial-vacuum-affine-weighted-channel-resolvent.json"
)
PARENT_SHA = "5c4a3649bc435475d6bfd65901af47b63da2543faa9bc9bcd6364e00fac9a02d"
SHEAR_SHA = "c7c0eb49287a4456328fe86b20a8362b70f540309439ba98b74ad1f191e0c1e6"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_weighted_channel_resolvent/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    for prior, expected in ((parent, PARENT_SHA), (shear_input, SHEAR_SHA)):
        if sha(prior.REPORT) != expected:
            raise ValueError(
                "A source-pinned original channel or actual curvature input changed"
            )
        prior.validate_report(
            json.loads(prior.REPORT.read_text()), prior.build_report()
        )
    return {
        "S6_225_actual_curvature_coordinates_and_ordered_reference_fully_rebuilt": PARENT_SHA,
        "S6_223_original_shear_pole_plus_cut_fully_rebuilt": SHEAR_SHA,
        "S6_86_original_trace_first_sheet_gap_and_radial_factor_transitively_rebuilt": True,
        "S6_222_actual_full_quantum_force_graph_transitively_rebuilt_not_inverted_here": True,
        "original_state_prescription_and_all_frozen_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError(
            "An exterior logarithm, weighted inverse or reference-scope gate failed"
        )
    return {
        "schema": 1,
        "claim": "P8-S6.226.ORIGINAL_CHANNEL_FULL_COMPLEX_LOG_COERCIVITY_AND_ORDERED_WEIGHTED_BOUNDED_CURVATURE_REFERENCE_CORRECTIONS",
        "date": "2026-09-12",
        "status": "UNIFORM_WEIGHTED_BOUNDED_CHANNEL_REFERENCE_CORRECTIONS; NOT_ACTUAL_FULL_QUANTUM_INVERSE_MATCHING_UNWEIGHTED_SMALLNESS_STABILITY_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/analytic.md",
            "notes/weight.md",
            "notes/inverse.md",
            "notes/curved.md",
            "notes/control.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "original_full_complex_channel_coercivity": serialize(payload(analytic.data())),
        "weighted_channel_corrections_curvature_composition_and_units": serialize(
            {"channels": payload(weighted.data()), "curved": payload(curved.data())}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "Both original fixed scalar factors obey the common exterior bound Re A_i(p)>=(13/80)log(|p|/(4m^2))-4 on the complete first sheet outside |p|=4m^2. The proof retains the actual finite coefficients and uses complete cut-bank, inner-circle and uniform outer-circle bounds, with the separate global shear real-part lower bound-172/225. For lambda=sigma+i omega and everyq>=0, the exact joint-frequency identity gives |lambda^2+q|>=sigma^2. Consequently sigma_M=2m exp(32+8M) makes the original diagonal causal inverse uniformly bounded by5/(32+13M) on exponentially weighted L2_tH^r, for every realr and all time/spatial frequencies. Every measurable time-dependent channel matrix V(t) with norm at mostM admits the ordered causal inverse (I+Kdiag V)^-1 Kdiag for Fdiag+V, with norm at most5/(32+8M). Neither channel mixing nor coefficient multiplication is commuted through the memory kernel. The original shear pole and complete initial distributional boundary are retained. Composing with the actual S225 curvature-coordinate inverses gives the ordered inverse of Bc*(Fdiag+V)Bc, bounded without spatial derivative loss by125/[4(32+8M)(sigma_M-20)^2(sigma_M-108)^2] on the unchanged conformal slab. The physical density is restored on the right, retaining64pi^2 kappa a^4. Removing the exponential weight costs exp(sigma_M T). An explicit separate constant correction of norm below208/315 creates a simple positive-real reference pole atlambda=2m while the weighted theorem still applies, preventing an inference of stability. Independent radial/closed complex checks, two-bank exterior and near-wave-cone diagnostics, exact joint-frequency identities, noncommuting causal finite solves and a Neumann tail control supplement the written harmonic and functional proof. The actual mass/state/contact/local-remainder/tree/matter correction has not been proved to belong to this bounded middle class. Full S222 inversion, physical stability and original V/G/B/P8 remain open.",
        "not_established": [
            "Identification or compatible weighted bounds for the actual full curved mass/state/contact/local/classical/matter remainder",
            "A complete actual S222 auxiliary/clock graph inverse or finite-coupling/nonlinear common-parent control",
            "Unweighted small response, physical damping, absence of growing poles, stability or a physical cutoff",
            "Original V/G/B or P8 closure; no original pole, finite term, state or preparation has been changed",
        ],
        "verification_boundary": "Exact factor/branch coefficients, joint-frequency identity, positive rational margins, chosen weights, all Neumann/coordinate/density constants and the separate growing-pole control are checked. The harmonic minimum principle, Fourier-Laplace Plancherel, causal graph extension and bounded-operator inverse are written proofs, not FORMALIZED. Numerical complex probes and finite causal matrices are independent diagnostics, not replacements for those continuum arguments. Native/direct/ordinary/CLI retain original SymPy; only full regression uses the audited exact-GCD adapter. Frozen predecessors are unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The weighted channel reference report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.226 weighted channel reference replay passed; actual full quantum inverse and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
