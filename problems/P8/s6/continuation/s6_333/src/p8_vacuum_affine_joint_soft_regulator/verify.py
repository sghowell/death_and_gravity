"""Read-only joint infrared and dimensional additional-soft regulator report."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_borel_soft_conversion import verify as input1
from p8_vacuum_affine_physical_virtual_soft_pairing import verify as input4
from p8_vacuum_affine_radiative_angular_finite import verify as input3
from p8_vacuum_affine_radiative_soft_state_transfer import verify as input2
from p8_vacuum_affine_remaining_log_cutoff import verify as input0

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates/polynomial-vacuum-affine-joint-soft-regulator.json"
PARENT_SHA = "13e795b571664878ba3bc5d41e3b823a897a8f8ed94af1433c9514cb2391920c"
BOREL_SHA = "642539fc5235de3a4a79099a5f105785cd89535e14ec3a46636f5957cfe3bea8"
STATE_SHA = "32396aad07a5c1521a282fe6161c0a8e99ff6e312f55df4937de8f5501e30655"
ANGULAR_SHA = "5e18135d1a49c2846047f2a71c3de00f7c82b84a131d8a6f924890308f4ccd99"
PHASE_SHA = "40c04cadfa1a1c542c16eed48f6df12bcfe0605ad69b39dbb70fcefaefe6b3ed"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_joint_soft_regulator/*.py"))
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
        (input1, BOREL_SHA),
        (input2, STATE_SHA),
        (input3, ANGULAR_SHA),
        (input4, PHASE_SHA),
    ):
        if sha(module.REPORT) != digest:
            raise ValueError("A frozen joint soft regulator input changed")
        module.validate_report(
            json.loads(module.REPORT.read_text()), module.build_report()
        )
    return {
        "S6_332_full_remaining_log_and_common_space_density_split_rebuilt": PARENT_SHA,
        "S6_331_positive_Borel_regulated_increment_and_cloud_transfer_rebuilt": BOREL_SHA,
        "S6_325_original_Born_energy_moduli_and_regulator_identity_rebuilt": STATE_SHA,
        "S6_301_full_fixed_ball_formula_and_phase_remainder_rebuilt": ANGULAR_SHA,
        "S6_296_original_dimensional_phase_and_Taylor_bound_rebuilt": PHASE_SHA,
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
        raise ValueError("A joint soft regulator proof gate failed")
    packets = audit.packets()
    names = tuple(packets)
    return {
        "schema": 1,
        "claim": "P8-S6.333.JOINT_INFRARED_AND_DIMENSIONAL_SOFT_REGULATOR_BOUND",
        "date": "2026-09-18",
        "status": "SCOPED_QUANTITATIVE_JOINT_SOFT_MARK_LIMIT; NOT_HARD_INTERACTING_V_G_B_P8_CLOSURE",
        "prior_sha256": prior,
        "source_sha256": {
            str(path.relative_to(ROOT)): sha(path) for path in source_files()
        },
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/source.md",
            "notes/kernels.md",
            "notes/moments.md",
            "notes/regulator.md",
            "notes/uniform.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "whole_original_source_kernels_and_moments": serialize(
            {name: payload(packets[name]) for name in names[:3]}
        ),
        "whole_quantitative_joint_regulator_bounds": serialize(
            {name: payload(packets[name]) for name in names[3:]}
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The specified additional-soft insertion has an explicit joint common-space L1 cutoff/regulator error <=[132000*a*eta*L_eta^2+400000*epsilon*a*x*L_x^2]/kappa, physically <[106000*eta*L_eta^2+320000*epsilon*x*L_x^2]/kappa^2. Its arbitrary joint and both iterated limits agree at fixed x and original hard data. This does not establish TV convergence of configuration laws, outer hard dimensional matching, interacting detector dynamics, all-N hard summability or original P8 closure.",
        "not_established": [
            "TV convergence of finite and infinite emission-configuration laws",
            "Outer hard dimensional/evanescent matching or finite radiative hard remainder",
            "Full detector probability, state-dependent dynamics or all-N hard sum",
            "Quantum unitarity and absolute complex Regge control",
            "Common-parent bounce, UV or original V/G/B/P8 closure",
        ],
        "verification_boundary": "Exact original phase/current identities, differentiated radial majorants and conditional beta moments support written quantitative joint-limit proofs. Independent radial and logarithmic quadratures plus majorant-model calibrations do not replace those proofs or the original physical current. Only the additional soft factor is dimensionally continued. Not kernel-formalized; original SymPy outside FULL.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The joint soft regulator report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.333 quantitative joint soft-regulator replay passed; original P8 remains OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
