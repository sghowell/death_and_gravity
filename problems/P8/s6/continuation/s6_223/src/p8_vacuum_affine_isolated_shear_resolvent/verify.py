"""Read-only isolated shear pole-plus-cut and finite-window causal inverse."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_proca_rank_one_inverse import verify as trace_input
from p8_vacuum_affine_flat_tensor_cut import verify as cut_input
from p8_vacuum_affine_local_tensor_response import verify as finite
from p8_vacuum_affine_quantum_forced_constraints import verify as parent

from . import audit, kernel, normalization, spectral

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates" / "polynomial-vacuum-affine-isolated-shear-resolvent.json"
)
PARENT_SHA = "05e36391c6636fa167c5569d79f706f04029548d9517c0686fa29aba70b15a40"
FINITE_SHA = "567a8a84e981c637b000ef671c350b63c54f4fc357ca60eaf1ce1de9c1293b02"
CUT_SHA = "e7a37017276d0b5474d81a295ab3f2574f69b883ffee34b032a0dadf9566ceec"
TRACE_SHA = "0f04fd277040ad21525eb1d9102e6af45c4d695fac3f69d0d37710a7ab8a0db6"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_isolated_shear_resolvent/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {k: v for k, v in data.items() if k not in ("checks", "gates")}


@cache
def prior_checks():
    inputs = (
        (parent, PARENT_SHA),
        (finite, FINITE_SHA),
        (cut_input, CUT_SHA),
        (trace_input, TRACE_SHA),
    )
    for prior, expected in inputs:
        if sha(prior.REPORT) != expected:
            raise ValueError(
                "A source-pinned physical cut, finite coefficient or force-graph input changed"
            )
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    # Explicit read-only validation of the reused inputs as well as the linear chain.
    for prior, _ in inputs[1:]:
        prior.validate_report(
            json.loads(prior.REPORT.read_text()), prior.build_report()
        )
    return {
        "S6_222_force_dependent_constraints_and_ordered_graph_fully_rebuilt": PARENT_SHA,
        "S6_189_original_fixed_finite_tensor_Hessian_fully_rebuilt": FINITE_SHA,
        "S6_199_actual_normalized_full_Proca_cut_fully_rebuilt": CUT_SHA,
        "S6_86_old_trace_factor_fully_rebuilt_but_gap_not_transferred": TRACE_SHA,
        "original_state_preparation_finite_prescription_and_frozen_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError(
            "An isolated shear normalization, pole or causal inverse gate failed"
        )
    return {
        "schema": 1,
        "claim": "P8-S6.223.ISOLATED_FIXED_SHEAR_FACTOR_UNIQUE_SUBTHRESHOLD_ZERO_POLE_PLUS_CUT_AND_FINITE_WINDOW_CAUSAL_INVERSE",
        "date": "2026-09-12",
        "status": "ISOLATED_REFERENCE_POLE_PLUS_CUT_FINITE_WINDOW_INVERSE; NOT_HALF_LINE_L1_FULL_COUPLED_PHYSICAL_POLE_STABILITY_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/normalization.md",
            "notes/first-sheet.md",
            "notes/cut.md",
            "notes/dispersion.md",
            "notes/kernel.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "physical_finite_normalization_and_exact_first_sheet_spectral_data": serialize(
            {
                "normalization": payload(normalization.data()),
                "spectral": payload(spectral.data()),
            }
        ),
        "ordinary_causal_kernel_domain_primitive_and_mass_scaling": serialize(
            payload(kernel.data())
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The original finite local shear Hessian, independently combined with the full normalized Proca cut, gives F2(0)=-1/30 and A2=-F2 with a positive radial weight. Unlike the older trace factor, A2 has exactly one simple subthreshold zero p=-r m^2. Exact rational geometric-series enclosures prove0.5798<r<0.5799 and16<R/m^2<17 for the positive inverse weight R=1/A2'(-r m^2); the residue of1/F2 is negative. Strict half-plane sign, real monotonicity and positive cut imaginary part exclude all other first-sheet zeros. The inverse is the negative pole term plus the negative positive-density cut integral, with total static moment30 and next moment675/(14m^2). The continuum kernel alone belongs toL1 on the half-line, but the nonzero undamped pole means the total kernel does not. The total ordinary kernel isL1 on each finite window, has primitive in[-60,0], and gives both causal inverse identities on the stated zero-past graph domain including its initial boundary. For prepared C1 sources the bound is60T||f'||. The forward distribution is explicitly constructed from a bounded continuous sine kernel, avoiding a divergent local/cut separation. Independent curvature, rational enclosure, complex radial and high-precision pole-plus-cut tests supplement the written contour and oscillatory proofs. No pole is deleted and no state or finite prescription is changed. This is an isolated flat reference factor, not a full tree-plus-loop pole, physical instability, curved matrix inverse, finite-coupling remainder, cutoff or original V/G/B/P8 closure.",
        "not_established": [
            "A full nonzero-transfer curved matrix normal form or compatible-space coupled quantum inverse",
            "A physical pole, ghost, instability or spectral positivity of the full tree-plus-loop propagator",
            "Half-lineL1 of the full shear kernel, a uniform numerical finite-window operator norm, pole deletion or order reduction",
            "Finite-coupling/nonlinear parent remainder, quantum background/stability, heavy/physical cutoff or original V/G/B/P8 closure",
        ],
        "verification_boundary": "Exact physical normalization, geometric root/residue enclosure, cut constants and kernel arithmetic are certified. Analyticity, Cauchy contour deformation, oscillatory integrability and causal distribution graph inversion are written proofs, not FORMALIZED. High-precision quadrature is an independent diagnostic, not the sign or existence proof. Native/direct/ordinary/CLI use original SymPy; only full regression uses the separately audited exact-GCD adapter. Frozen scientific inputs remain unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The isolated shear pole-plus-cut report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.223 isolated shear finite-window inverse replay passed; full coupled inverse and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
