"""Read-only actual prepared bounded-Fourier scalar/clock/matter inverse."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_full_spatial_remainder import verify as spatial_input
from p8_vacuum_affine_homogeneous_shear_inverse import verify as parent
from p8_vacuum_affine_quantum_forced_constraints import verify as constraints_input

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-band-coupled-inverse.json"
PARENT_SHA = "88da57811f211f16c36c38936bec26f6f616de515a86b93438ebde91e54da2a5"
CONSTRAINTS_SHA = "05e36391c6636fa167c5569d79f706f04029548d9517c0686fa29aba70b15a40"
SPATIAL_SHA = "df188e09486c7666c0dcf9c811702d1be096849acf993c8efae18bfffc9532a6"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_band_coupled_inverse/*.py"))
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
        (constraints_input, CONSTRAINTS_SHA),
        (spatial_input, SPATIAL_SHA),
    ):
        if sha(prior.REPORT) != expected:
            raise ValueError(
                "A source-pinned actual homogeneous inverse, full forced constraint or complete spatial current changed"
            )
        prior.validate_report(
            json.loads(prior.REPORT.read_text()), prior.build_report()
        )
    return {
        "S6_228_actual_homogeneous_strong_remainder_and_prepared_inverse_fully_rebuilt": PARENT_SHA,
        "S6_222_complete_forced_constraints_contacts_and_original_graph_fully_rebuilt": CONSTRAINTS_SHA,
        "S6_219_actual_complete_spatial_current_original_regulator_and_Ward_fully_rebuilt": SPATIAL_SHA,
        "S6_220_full_classical_S6_216_dimensional_pairs_S6_176_reference_parent_bridge_and_original_trace_shear_inverses_transitively_rebuilt": True,
        "original_state_prescription_and_all_frozen_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError(
            "An actual coupled principal, matching, normal-form or prepared inverse gate failed"
        )
    return {
        "schema": 1,
        "claim": "P8-S6.229.ACTUAL_FULL_SCALAR_CLOCK_MATTER_PREPARED_CAUSAL_INVERSE_AND_S222_FEEDBACK_ON_SMOOTH_BOUNDED_FOURIER_SUPPORT",
        "date": "2026-09-12",
        "status": "ACTUAL_SMOOTH_BOUNDED_EXTERNAL_MOMENTUM_COUPLED_INVERSE; NOT_UNRESTRICTED_GRAPH_NUMERICAL_SMALLNESS_STABILITY_NONLINEAR_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/coordinates.md",
            "notes/ward.md",
            "notes/matching.md",
            "notes/remainder.md",
            "notes/inverse.md",
            "notes/recovery.md",
            "notes/scope-validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "actual_full_classical_Ward_contacts_residual_adjoint_and_dimensional_scalar_matching": serialize(
            {
                key: payload(value)
                for key, value in audit.packets().items()
                if key != "ordered_mixed_inverse_and_actual_smooth_ball_graph"
            }
        ),
        "actual_mixed_normal_form_ordered_inverse_and_smooth_ball_S222_graph": serialize(
            payload(audit.inverse.data())
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The actual full sourced classical scalar action in prepared coordinates(n,zeta,b,sigma)=(eta',w+Heta,c'+qeta,r+ell eta) has mixed row orders(4,4,4,2) and principal diag(-6delta^2,0,0,-1), with the matter invariant essential. Both ordered Ward legs and the distinct nonlinear clock nn contact contain at most one derivative on either leg and no inverse external momentum after full cancellation. The complete dimensional scalar pair gives C(d),B(d) and the fixed-physical leading matrix[[4,-4/3],[-4/3,8/5]], including its fixed-physical first dimension jet. Proper full-dimensional matching retains the original fourth finite constants and the complete original factors. Uniform finite-transfer mode expansions on each fixed external momentum ball retain the orderzero residual phase, infinite internal loop, all-order prepared state, full contacts and classical retuning. They prove a strong mixed-row normal form with complete weak-log simultaneous derivative majorantsCj(Pmax). With the original finite-window trace/shear inverse kernels and b0=max(15625/6144,16Kmax/(9gamma),1), the finite weight(4b0C+1)^2 gives a strict bounded-causal contraction and ordered inverse(I+B0V)^-1 B0 diag(I4,I4,I4,I2). Prepared time-shift regularity retains derivatives of the variable eta pivot as well as the kernel. The weighted residual adjoint and original zero germ recover both lapse and shift equations. Full force, output density and rank-two direct auxiliary terms then give both identities for the ACTUAL S222 scalar/clock/matter feedback inverse on smooth prepared bounded-Fourier histories, including a general mathematical phase drive. The unweighted adapted-source bound is2 exp(weight)b0 times the row-primitive source norm; constants are finite but unevaluated and all kappa factors remain. This is an actual coupled reference linear result, not only a flat reference factorization, but it is not a uniform all-momentum bound, unrestricted or maximal X/Y graph inverse, numerical smallness, stability, nonlinear common-parent control, physical cutoff or original V/G/B/P8 closure.",
        "not_established": [
            "Uniform Pmax-to-infinity estimates or a bounded inverse on the unrestricted original completed X_r/Y_r graph",
            "Maximal graph characterization, graph-density, spatial Schwartz-output preservation or arbitrary unprepared initial data",
            "Numerical complete remainder and inverse norms, small physical feedback, absence of growing modes, stability or a physical cutoff",
            "Nonlinear and finite-coupling common-parent control, original V/G/B or P8 closure",
        ],
        "verification_boundary": "Exact complete classical principal, both Ward legs and clock contact, fixed-physical dimensional jets, residual adjoint, finite-weight constants and causal-order controls are replayed. The actual state expansion, original singular extension, uniform compact-ball strong normal form, causal Neumann convergence, smooth regularity and Schur graph recovery are written continuum proofs, not FORMALIZED. Independent finite matrices and polynomial tests are diagnostics, not numerical continuum certificates. Native/direct/ordinary/CLI use original SymPy; only full regression uses the audited exact-GCD adapter. Frozen predecessors remain unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The actual smooth Fourier-ball coupled inverse report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.229 actual smooth Fourier-ball coupled inverse replay passed; unrestricted graph, stability and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
