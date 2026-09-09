"""Read-only exact actual nearby classical scalar cones and ray separation."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as original
from p8_proca_nearby_bounce import verify as parent

from . import audit, bounds, geometry, principal, rational

ROOT=Path(__file__).resolve().parents[2]
REPORT=ROOT/"certificates"/"actual-nearby-classical-cones.json"
PARENT_SHA="e1063cf789b86b4110efadeb634438c8824bc4e8faeeb276cf06ec172f2cd832"
sha,serialize=original.sha,original.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_proca_nearby_cones/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT)!=PARENT_SHA:
        raise ValueError("The frozen nearby classical bounce solution report changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()),parent.build_report())
    return {"S6_88_fully_rebuilt":PARENT_SHA,
            "same_literal_sourcefree_classical_action_and_actual_nearby_solution":True,
            "no_parent_quantum_profile_state_or_response_transfer":True}


@cache
def build_report():
    prior=prior_checks()
    identities,gates=audit.residuals(),audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError("An actual nearby classical cone or separation gate failed")
    exact=original.certify_residuals(identities)
    return {"schema":1,"claim":"P8-S6.89.ACTUAL_NEARBY_CLASSICAL_CONES","date":"2026-09-09",
            "status":"ACTUAL_NEARBY_CLASSICAL_BOUNCE_HAS_POSITIVE_SCALAR_PRINCIPAL_FORMS_LUMINAL_MATTER_AND_UNIFORMLY_SUPERLUMINAL_CLOCK_ON_ABS_U_LE_10_MINUS_7; CLASSICAL_RAY_SEPARATION_POSITIVE; FINITE_BAND_QUANTUM_MATCHING_AND_ORIGINAL_P8_OPEN",
            "prior_sha256":prior,"source_sha256":{str(p.relative_to(ROOT)):sha(p) for p in source_files()},
            "formulation":"FORMULATION.md",
            "written_proofs":["notes/reduction.md","notes/evolution.md","notes/enclosures.md","notes/separation.md","notes/scope.md"],
            "exact_residuals":exact,"named_exact_check_count":len(identities),
            "checked_scalar_entries":sum(v.rows*v.cols if isinstance(v,sp.MatrixBase) else 1 for v in identities.values()),
            "proof_checks":gates,
            "generic_full_geometric_scalar_phase_pairs":serialize({
                left+"_"+right:geometry.pair(left,right) for i,left in enumerate(geometry.CHANNELS)
                for right in geometry.CHANNELS[i:]}),
            "Euler_first_scalar_principal_and_physical_characteristic_basis":
                serialize({"full_Hessian_limits":principal.data(),"formula":principal.formula(),"physical":principal.physical()}),
            "literal_primitive_free_rational_coefficients":serialize(rational.coefficients()),
            "exact_rational_system_polynomial_signatures":rational.signatures(),
            "continuous_physical_cone_whole_box_enclosures":serialize(bounds.enclosures()),
            "actual_nearby_solution_box_bridge":serialize(bounds.solution_bridge()),
            "classical_comoving_characteristic_separations":serialize({
                str(value):bounds.characteristic_separation(value) for value in (bounds.TIME_WINDOW,bounds.TIME_WINDOW/2)}),
            "controls":audit.controls(),
            "verdict":"The actual S6.88 nearby classical bounce with central lapse 1+10^-6 has positive scalar principal kinetic and gradient matrices throughout |u|<=10^-7. A new full arbitrary-homogeneous geometric reduction retains the moving trace boundary and Euler time terms. Its literal primitive-free rational evolution reproduces the independent central lapse and speed data. Exact whole-polynomial box bounds and the actual-solution bridge give 2.5*10^-6<c_clock-1<8*10^-6 in the physical matter metric, while free matter is exactly luminal. Outgoing principal characteristics emitted at the same event have comoving separation greater than 99/(2*10^14) across the full interval. This is a classical ray result, not a finite-frequency signal or EFT/UV verdict.",
            "not_established":["A globally complete nearby bounce or nonlinear inhomogeneous/all-wavelength stability",
                "A regular Legendre chart at every finite momentum, wavepacket mode-conversion errors or a resolved finite-band detector signal",
                "A quantum state, fixed profile, cancellation, response or inverse on these different nearby backgrounds",
                "Canonical interaction estimates, mixed loops, omitted operators, thresholds or an interacting Wilsonian cutoff on this domain",
                "Finite Wilson matching, healthy full-candidate Minkowski vacuum, common-parent V/G/B, UV completion/exclusion or original P8 closure"],
            "verification_boundary":"Exact full geometric phase-pair reduction, high-q Euler-first canonical limits, literal primitive and fixed-phase rationalization, central and parity regressions, entire rational-polynomial coefficient enclosures, and written solution-bridge and characteristic comparison proofs. Exact intermediate rational polynomials are fully recomputed and their canonical monomial coefficient signatures pinned. Every parent is natively rebuilt. No floating inequality, sampled curve, library patch, finite-frequency or quantum transfer. Not proof-assistant formalization."}


def validate_report(expected,actual):
    if expected!=actual:
        raise ValueError("The actual nearby classical cone report differs from read-only replay")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    actual=build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()),actual)
        print("P8 S6.89.ACTUAL_NEARBY_CLASSICAL_CONES replay passed; finite-band, quantum, matching and original P8 OPEN")
    else:
        print(json.dumps(actual,indent=2,sort_keys=True))


if __name__=="__main__":
    main()
