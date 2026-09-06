"""Read-only beta1 parent screen: flat vacuum match, exact CD obstruction."""

import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_m1_weyl_scalar import verify as scalar_prior
from p8_matching import verify as matching_prior

from . import background, matching, vacuum

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"beta1-parent-screen.json"
PARENT_SHA = "e9aac0256e1a57e684b123cc5517863f2f1463a423a44cee11b59986665199f8"
SCALAR_SHA = "b00d124cf5c7ab712cb7f3198766f1aa542065d801741273db37855a5824d27c"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@cache
def prior_checks():
    checks = {}
    for name, module, expected in (("S6.2", matching_prior, PARENT_SHA),
                                   ("S5.11.CD", scalar_prior, SCALAR_SHA)):
        if sha(module.REPORT) != expected:
            raise ValueError(f"Pinned {name} certificate changed")
        module.validate_report(json.loads(module.REPORT.read_text()), module.build_report())
        checks[name] = expected
    return checks


def residuals():
    groups = {"full_lapse_retaining_Einstein_and_potential_variation": background.variation_checks(),
              "Noether_and_undivided_Bianchi": background.bianchi_checks(),
              "exact_regular_bounce_obstruction": background.bounce_checks(),
              "forbidden_NEC_algebraic_point_not_a_parent_solution": background.NEC_violating_point_checks(),
              "full_relative_FP_mass_structure": vacuum.potential_checks(),
              "positive_mass_basis_and_source_couplings": vacuum.mass_basis_checks(),
              "source_preserving_Schur_and_true_massive_integration": matching.source_matching_checks(),
              "full_conserved_source_projectors": matching.projector_checks(),
              "physical_on_pole_residue_squares": matching.physical_residue_checks(),
              "all_ten_FP_equations_and_massless_source_equations": matching.fierz_pauli_source_checks(),
              "curvature_squared_projector_normalization": matching.flat_curvature_checks(),
              "direct_four_index_curvature_contractions": matching.direct_curvature_checks(),
              "flat_rational_kernel_remainder": matching.remainder_checks()}
    for name, values in groups.items():
        if any(sp.simplify(value) != 0 for value in values.values()):
            raise ValueError("A beta1 parent-screen identity failed: "+name)
    return groups


def positive_checks():
    spectrum, match = vacuum.spectrum(), matching.tt_schur()
    values = {key: spectrum[key] for key in ("M_squared", "relative_mode_kinetic_coefficient",
                                            "m_FP_squared", "massless_source_residue", "massive_source_residue")}
    values.update({"clock_mass_squared": vacuum.MPHI2,
                   "flat_quadratic_Weyl_coefficient": match["c_C"],
                   "positive_ratio_root_margin": sp.Rational(3, 4),
                   "positive_kernel_weight_margin": background.MG2/(background.MG2+background.MF2),
                   "strict_CD_derivative_mismatch": background.bounce_equations()["minimum_Hdot_error"]})
    if any(value.is_positive is not True for value in values.values()):
        raise ValueError("A stipulated-domain positive expression is not certified")
    return values


def control_checks():
    groups = {"background": background.controls(), "vacuum": vacuum.controls(), "matching": matching.controls()}
    zeros = {"zero_interaction_loses_branch_factor", "zero_spring_has_no_heavy_gap",
             "zero_curvature_R2_tree_coefficient_not_radiative_closure"}
    for name, values in groups.items():
        for key, value in values.items():
            if (sp.simplify(value) == 0) != (key in zeros):
                raise ValueError(f"A beta1 parent-screen control failed: {name}/{key}")
    return groups


def serialize(value):
    if isinstance(value, dict):
        return {key: serialize(item) for key, item in value.items()}
    if isinstance(value, sp.MatrixBase):
        return [[str(item) for item in row] for row in value.tolist()]
    return str(value)


@cache
def build_report():
    previous = prior_checks()
    identities, positive, controls = residuals(), positive_checks(), control_checks()
    sources = sorted(ROOT.glob("src/p8_bimetric/*.py"))+sorted(ROOT.glob("tests/*.py"))
    sources += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    return {"schema": 1, "claim": "P8-S6.3.beta1", "date": "2026-09-06",
            "status": "SPECIFIED_BETA1_PARENT_HAS_POSITIVE_FLAT_VACUUM_AND_QUADRATIC_TREE_WEYL_MATCH_BUT_NO_REGULAR_CD_BOUNCE; GENERAL_S6_MATCHING_OPEN",
            "prior_certificate_sha256": previous,
            "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in sources},
            "formulation": "FORMULATION.md", "written_proof": "notes/parent-screen.md",
            "primary_source_audit": "notes/sources.md",
            "primary_sources": {"full_bimetric_FLRW_Bianchi_and_equations": "https://arxiv.org/pdf/1111.1655",
                                "massive_mode_and_correct_elimination_cautions": "https://arxiv.org/pdf/1807.05011",
                                "source_and_wrong_field_elimination_cautions": "https://arxiv.org/pdf/1303.6940"},
            "exact_residuals": {name: dict.fromkeys(values, "0") for name, values in identities.items()},
            "positive_expressions_on_stipulated_domain": serialize(positive),
            "positive_and_omission_controls": serialize(controls),
            "model": {"signature": "+---; R_FLRW=-6*(Hdot+2H^2)",
                      "action": "-MG2*integral sqrt|g| Rg/2 -MF2*integral sqrt|f| Rf/2 -m4*integral sqrt|g| sum beta_n e_n(sqrt(g^-1*f)) +S_can[g]",
                      "parameters": "MG2>0, MF2>0, NU=m4*beta1>0; beta0=-3beta1,beta4=-beta1,beta2=beta3=0",
                      "matter": serialize(vacuum.canonical_matter()),
                      "physical_frame": "g is the prescribed matter metric; no derivative metric redefinition or added old DHOST action",
                      "FLRW_domain": "flat common spatial slicing; finite positive a,b,Ng,Nf and smooth positive square-root branch",
                      "matter_extension_of_obstruction": "any classical positive field-metric minimally coupled g-only matter, arbitrary smooth potential"},
            "flat_vacuum_spectrum": serialize(vacuum.spectrum()),
            "flat_source_preserving_match": serialize(matching.tt_schur()),
            "flat_conserved_source_response": "(P2-P0/2)/(M_squared*D) + MF2*P2/(MG2*M_squared*(D+m_FP_squared))",
            "tree_matching_scope": "only flat quadratic four-derivative pure-metric action modulo Euler/time boundaries: cC=MF2^2/(4NU), cR=0; physical sources retained",
            "flat_kernel_remainder": serialize(matching.remainder_record()),
            "undivided_dynamical_Bianchi": "Ng*b_dot-Nf*a_dot=0; beta1>0 forbids the algebraic branch; no division by H or a_dot",
            "exact_bounce_obstruction": serialize(background.bounce_equations()),
            "conclusion": "On this regular parent branch H_g=0 implies y=1 and -2*(MG2+MF2)*Hdot_g=rho+p>=0; CD requires Hdot_g=4/tau^2>0",
            "robust_window_mismatch": "On [-tau/2,tau/2], simultaneous endpoint |H_parent-H_CD|<8/(5tau) and uniform |Hdot_parent-Hdot_CD|<48/(25tau^2) are impossible on the stipulated regular branch in physical g cosmic time",
            "verification_boundary": [
                "Both previous certificates and their replayable source lineages are checked; no previous result is replaced",
                "Exact symbolic two-lapse, mass/source, curvature and rational remainder equalities are replayed",
                "Independent covariant/lapse/FP source tests are hashed and run separately by pytest",
                "The positive-branch bounce contradiction and absolute-value remainder estimate have written proofs, not formalized global PDE proofs"],
            "not_established": [
                "A bimetric match to the rolling CD/DHOST coefficients, couplings, modes or matter state",
                "A general bimetric, higher-curvature, NEC-violating or UV-completion exclusion",
                "Nonlinear/cosmological curvature-squared matching from the flat quadratic result",
                "Loop closure or bounds on finite R2/Weyl/nonlocal coefficients; tree cR=0 does not erase S5.9",
                "A physical cutoff, interaction weak coupling, all-mode parent health on a rolling state, or positivity at finite gravity",
                "An exact parent realization of the S5.10/S5.11 candidate or its finite-band control from the vacuum pole expansion",
                "P8(a), S6, or full P8 completion"]}


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("Beta1 parent-screen certificate differs from exact replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.3.beta1: positive flat spectrum/source match and exact regular CD-bounce obstruction replay passed; general matching OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
