"""Read-only full-parent tiny classical bounce and finite-band energy certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_physical_volume_correction import verify as previous

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates" / "polynomial-vacuum-affine-finite-band-neighborhood.json"
)
PARENT_SHA = "71791d5566730c91d4e5b205887318d9f130d6b7f445f03e3010739bc78ff83a"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_finite_band_neighborhood/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    if sha(previous.REPORT) != PARENT_SHA:
        raise ValueError("The frozen corrected physical-parent binding changed")
    previous.validate_report(
        json.loads(previous.REPORT.read_text()), previous.build_report()
    )
    return {
        "S6_262_correct_actual_physical_binding_and_whole_ancestry_rebuilt": PARENT_SHA,
        "S261_false_physical_binding_is_not_renewed": True,
        "same_full_sources_profiles_preparation_and_canonical_boundaries": True,
        "distinct_classical_family_not_fixed_interacting_quantum_mean": True,
        "finite_band_not_Wilsonian_cutoff_or_original_P8_closure": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError(
            "A complete source, physical-turn or finite-band energy gate failed"
        )
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.263.FULL_CURRENT_SIGNED_TINY_CLASSICAL_PHYSICAL_BOUNCE_AND_COMPLETE_SIXTEEN_PHASE_FINITE_BAND_REFERENCE_ENERGY_WITH_UNCHANGED_PARENT",
        "date": "2026-09-14",
        "status": "FULL_CURRENT_TINY_HOMOGENEOUS_CLASSICAL_BOUNCE_AND_COMPLETE_FINITE_BAND_LINEAR_ENERGY; NOT_INTERACTING_QUANTUM_MEAN_REGULATOR_CUTOFF_NONLINEAR_INHOMOGENEOUS_B_ORIGINAL_V_G_B_OR_P8",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/trajectory.md",
            "notes/physical.md",
            "notes/reference.md",
            "notes/energy.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_complete_sources_classical_neighborhood_and_physical_turn": serialize(
            {name: payload(packets[name]) for name in names[:6]}
        ),
        "whole_full_reference_energy_and_finite_band_variational_bounds": serialize(
            {name: payload(packets[name]) for name in names[6:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The unchanged full parent has a signed tiny family of complete constrained homogeneous classical solutions on the symmetric clock slab[-10^-60,10^-60]. The complete initial lapse equation sets the positive M1 rate without a bare-root replacement, with abs(epsilon)<=10^-230. Fifth source jets, full constrained Euler derivatives and the retained heavy energy keep all trajectories in the regular branch. The corrected physical factor R**(-1/2) gives a unique local physical turn with proper-time Hubble derivative greater than3. The entire16-phase linear generator, including both coupled scalars, the heavy mode, all three Proca polarizations and both tensors, has root-energy growth below2 on P in[10^64,2*10^64] in the single unchanged reference energy. All reference time and original canonical endpoint maps remain. The larger-displacement S256 growth example is unchanged. These are exceptionally local classical and fixed-band linear results, not an interacting quantum mean, a Wilsonian cutoff, a nonlinear inhomogeneous Cauchy theorem or closure of original V/G/B/P8.",
        "not_established": [
            "A stationary interacting physical quantum mean, complete quantum ordering/gauge regulator or vanishing Ward/state/endpoint defects",
            "A Wilsonian cutoff or controlled all-sector UV matching inferred from the selected mathematical momentum band",
            "A nonlinear inhomogeneous Cauchy theorem, stable evolution on a macroscopic interval or an all-momentum same-space bound",
            "Refutation of the unchanged larger-lapse-displacement S256 finite-band growth result",
            "Physical curved subtraction, complete omitted-loop estimates or Gaussian support in the nonlinear auxiliary branch",
            "Original physical UV scattering, quantum gravity limit, finite-gravity IR/Regge or V/G/B/P8 closure",
            "Formal verification of the written source, ODE continuation, energy and physical-turn arguments",
        ],
        "verification_boundary": "Exact whole-source and canonical identities, deterministic rational interval enclosures and independent finite comparison fixtures support the stated local classical and finite-band linear result. Fixtures do not replace the full parent or its preparation. The S261 physical-binding refutation is retained. Native/direct/ordinary/CLI use original SymPy; only captured full P8 regression uses the audited exact-GCD adapter. No frozen predecessor is modified.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete finite-band neighborhood report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.263 tiny classical physical-bounce and complete finite-band linear-energy replay passed; quantum mean, cutoff and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
