"""Read-only exact chosen auxiliary-parent applicability certificate replay."""

import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_composite_vacuum import verify as prior

from . import extension, independent, matter, model, vacuum

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"auxiliary-parent.json"
CONTRACT = ROOT.parents[1]/"FORMULATION.md"
PRIOR_SHA = "b06bfa4c45aa4e7ad5b0bd4a2c5c928d9616534b69e4c765a2d752966c78462e"
CONTRACT_SHA = "d00df35d5821b05baa61e897725a22ed1ac21d0536b74bdaece3ab0093215901"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@cache
def prior_checks():
    if sha(prior.REPORT) != PRIOR_SHA:
        raise ValueError("Pinned S6.12 context report changed")
    if sha(CONTRACT) != CONTRACT_SHA:
        raise ValueError("The adopted S6 conditional vacuum/matching contract changed")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return {"S6_12_source_preserving_branch": PRIOR_SHA, "adopted_S6_contract": CONTRACT_SHA}


def controls():
    v, m, e = vacuum.controls(), matter.controls(), extension.controls()
    expected_v = (-2, -2, 0, -2, -2, 0, -1, sp.Rational(-1, 4), sp.Rational(2, 3), sp.Rational(9, 32))
    if tuple(v.values()) != expected_v:
        raise ValueError("A no-flat, disconnected-link, coefficient or exact-cancellation control failed")
    expected_m = (0, sp.Rational(5, 9), sp.Rational(-1, 6), sp.Rational(-1, 2), -1/model.B)
    if tuple(m.values()) != expected_m:
        raise ValueError("A physical-metric or sourced Lorentz-constraint control failed")
    if (e["positive_vacuum_mass_squared_G_F_q_1"], e["omitting_beta4g_restores_nonzero_g_Euler_00"],
        e["negative_link_changes_mass_sign"]) != (2, -2, -2):
        raise ValueError("The separately named beta4 extension control failed")
    return {"chosen_parent_no_flat_and_cancellation": {key: str(value) for key, value in v.items()},
            "actual_source_metric_and_lorentz_constraint": {key: str(value) for key, value in m.items()},
            "separately_named_beta4_extension": {key: str(value) for key, value in e.items()}}


@cache
def build_report():
    previous = prior_checks()
    residuals = {"full_e_f_u_Euler_and_undivided_source_map": model.checks(),
                 "canonical_actual_matter_vierbein_dictionary": matter.checks(),
                 "exact_auxiliary_vacuum_and_restricted_coefficients": vacuum.checks(),
                 "separate_beta4_flat_and_FP_countercontrol": extension.checks()}
    if any(sp.simplify(value) != 0 for group in residuals.values() for value in group.values()):
        raise ValueError("An auxiliary-parent, physical source or beta4 extension identity failed")
    rational = independent.checks()
    omissions = controls()
    ext = extension.derive()
    relative_kinetic = ext["G"]*ext["F"]/(ext["G"]+ext["F"])
    if ext["relative_mass_squared"].is_positive is not True or relative_kinetic.is_positive is not True:
        raise ValueError("The extended flat quadratic mass/kinetic positivity failed")
    sources = sorted(ROOT.glob("src/p8_trimetric/*.py"))+sorted(ROOT.glob("tests/*.py"))
    sources += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    return {
        "schema": 1,
        "claim": "P8-S6.13.AUXILIARY",
        "date": "2026-09-06",
        "status": "CHOSEN_AUXILIARY_PARENT_FAILS_REGULAR_CONSTANT_FLAT_VACUUM_PREREQUISITE; NAMED_BETA4_EXTENSION_SURVIVES_THIS_SCREEN; P8_OPEN",
        "prior_context_sha256": previous,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in sources},
        "formulation": "FORMULATION.md", "written_proof": "notes/proof.md", "source_audit": "notes/sources.md",
        "exact_residuals": {name: dict.fromkeys(group, "0") for name, group in residuals.items()},
        "independent_Fraction_replay": rational,
        "checked_negative_and_positive_controls": omissions,
        "specified_parent": {
            "source_model": "arXiv:1804.04671 section 2.1 choice m_h²=0 and beta_n^g=beta_n^f=0 for n>=2; not the paper's general trimetric action",
            "vierbein_dictionary": "g=e^T eta e,f=v^T eta v,h=u^T eta u; eta=diag(1,-1,-1,-1); real invertible positive-orientation frames",
            "root_domain": "Parent real square-root/symmetrization constraints are additionally required; ambient noncommuting variation fixtures are not assumed to solve them",
            "physical_matter_frame": "Actual h=u^T eta u only, with no direct e/v matter action; h is auxiliary but remains the full prescribed physical source metric",
            "gravitational_action": "-(G/2)sqrt|g|R[g]-(F/2)sqrt|f|R[f]-2 det(u)[B+tr(u^-1 Q)], G,F>0,Q=p_g e+p_f v",
            "chosen_action_omissions": "No h EH term and no separate g/f cosmological terms; the latter would be source beta4g,beta4f, excluded by n>=2=0",
            "matter": "epsilon S_m[h,psi]; canonical dictionary L_m=det(u)[Y/2-V], Y=partial psi h^-1 partial psi",
            "parameters": "B,p_g,p_f constant real; no-flat result needs either link nonzero, not B!=0 or link positivity",
            "units": "G,F have mass dimension 2; B,p_g,p_f,beta4g,beta4f,V dimension 4; frames and epsilon dimension 0; psi dimension 1",
            "epsilon_dictionary": "Source matter-strength bookkeeping parameter, not the old composite asymmetry epsilon=1/Y",
        },
        "full_variation_and_source": {
            "gradient_convention": "delta L=tr(E_e^T delta e+E_v^T delta v+E_u^T delta u), J_u=delta L_m/delta u, with compact variations or EH boundary terms",
            "e_equation": "E_e=G det(e) eta e Einstein[g]^{up,up}-2p_g det(u)u^-T",
            "f_equation": "E_v=F det(v) eta v Einstein[f]^{up,up}-2p_f det(u)u^-T",
            "u_equation": "E_u=-2det(u)[s u^-T-u^-T Q^T u^-T]+epsilon J_u; s=B+tr(u^-1 Q)",
            "matter_equation": "epsilon delta S_m/delta psi=0 retained separately",
            "source_trace_reversal": "T=u J_u^T u/(2det u), t=tr(u^-1 T),calT=(t/3)u-T; C=Q+(B/3)u-epsilon calT",
            "undivided_equivalence": "E_u=2det(u)u^-T[C^T-tr(u^-1 C)u^T]u^-T, and E_u=0 iff C=0 without dividing B",
            "source_convention": "T is the explicitly normalized inverse-vierbein source matrix, not an assumed identification with another physical-stress convention",
        },
        "exact_chosen_parent_prerequisite_failure": {
            "domain": "Any constant nonsingular vierbein triple with at least one nonzero link; hence in particular the stipulated regular common Lorentz-invariant Minkowski vacuum",
            "identity": "Einstein[g]=Einstein[f]=0; u^T E_e=-2p_g det(u)I and u^T E_v=-2p_f det(u)I",
            "invertible_nonzero_residual": "det(E_e)=16p_g^4 det(u)^3 for a nonzero g link, and similarly f; either equation contradicts a full flat solution",
            "independence": "No auxiliary equation, B division, Q inverse, matter truncation or FP formula is used; any source acting only on u leaves the contradiction unchanged",
            "constant_potential": "-epsilon V det(u) merely shifts B_eff=B+epsilon V/2 in the exact full action and cannot repair E_e/E_v",
            "both_zero_links": "Outside the premise; B_eff=0 can leave free Einstein metrics and undetermined u, not an invertible auxiliary-parent construction",
            "excluded_boundary": "Singular u is not a regular physical metric or allowed inverse chart",
            "interpretation": "A failure of this chosen parent's additional S6 flat-vacuum prerequisite, not a refutation of the source's ghost-free degree-of-freedom construction or exclusion of curved vacua",
        },
        "exact_elimination_and_action_nonidentity": {
            "extra_domain_for_elimination": "B!=0 and det(Q)!=0; formulas using p_f/p_g additionally require p_g!=0, unlike the undivided theorem",
            "vacuum_solution": "u0=-3Q/B; V_eff=-54det(Q)/B³; sqrt|h0|=81det(Q)/B⁴",
            "geometric_coefficients": "beta_n=-54p_g^(4-n)p_f^n/B³=C r^n, C=-54p_g^4/B³,r=p_f/p_g",
            "weights": "Actual u0 weights a=-3p_g/B,b=-3p_f/B have b/a=r; source equation 3.3 uses their opposite common sign in the metric square",
            "not_old_beta2_parent": "beta1 beta3-beta2²=0; old bare beta2 value is -1 and its S6.12 vacuum-shifted value is -1/4; endpoint-only beta0/beta4 shifts do not change this minor",
            "full_constant_V_resummation": "Where B_eff!=0, u=-3Q/B_eff and L_eff=54det(Q)/B_eff³; relative to fixed h0 volume, L_eff=(2B/3)[B/B_eff]³ sqrt|h0|",
            "truncated_cancellation_control": "First-order coefficient 2B/3-epsilon V vanishes at epsilon V/B=2/3, but the full coefficient is 9B/32 and u=3u0/4",
            "FP_caution": "A cancelled truncated geometric potential would remove its whole relative interaction, but this is not an actual flat parent vacuum; do not compute a physical Minkowski mass around a nonsolution",
        },
        "actual_matter_frame_and_constraint": {
            "canonical_gradient": "J_u=det(u)[(Y/2-V)u^-T-eta u w w^T], w=h^-1 k,k=partial psi",
            "canonical_source": "t=Y/2-2V; calT=u[-(Y/12+V/6)I+w k^T/2]",
            "local_formal_branch": "For fixed e,v,k, regular u0 and B!=0, the epsilon=0 u-Jacobian of C is (B/3)I; finite-dimensional analytic IFT supplies a local branch for canonical analytic matter, without a computed finite radius",
            "first_order": "u1=(3/B)calT(u0); h1=u0^T eta u1+u1^T eta u0=[3kk^T-(Y0/2+V)h0]/B",
            "actual_vs_leading": "Physical h=h0+epsilon h1+O(epsilon²); frozen h0 is not the exact matter frame; constant V agrees with h=[B/B_eff]²h0",
            "source_lorentz_constraint": "With both links nonzero, p_f(e^T eta v-v^T eta e)=epsilon[e^T eta calT(u0)-calT(u0)^T eta e]+O(epsilon²), with analogous second constraint",
            "omission_control": "At e=diag(1,2,3,4),v=I,B=-6,p_g=p_f=1,V=0,k=(1,1,0,0),Y0=5/9, the antisymmetric epsilon coefficient 01=-1/6 and h1_01=-1/2",
            "fixture_scope": "A local off-shell source/constraint dictionary witness, not a full cosmological solution",
            "metric_vierbein_boundary": "The sourced effective vierbein theory cannot generically be replaced by the old metric-only composite action; no finite-source, temporal-cutoff or original free-M1 error estimate follows from these formal coefficients",
        },
        "separately_named_beta4_countercontrol": {
            "new_action": "Add -2beta4g det(e)-2beta4f det(v), genuinely changing the specified parent without changing the auxiliary equation",
            "exact_flat_data": "p_g=p_f=q>0,B=-6q,beta4g=beta4f=-q,V=0,e=v=u=I; with constant V use B=-6q-epsilon V/2",
            "flat_checks": "All three full Euler matrices and total vacuum density vanish; omitting either compensating beta4 restores the respective nonzero equation",
            "exact_source_free_solution": "u0=(e+v)/2; L_eff,int=-q det(e+v)/4+2q(det e+det v)",
            "full_relative_invariant": "e=r(I+Delta),v=r(I-Delta),r=(e+v)/2; L_int=4q det(r)[e2(Delta)+e4(Delta)] with quadratic FP term -2q det(r)[trDelta²-(trDelta)²]",
            "physical_clock": "h=g=f=eta at the vacuum, no extra constant time rescaling",
            "TT_kernel": "-[G D h_g²+F D h_f²+q(h_g-h_f)²]/8, D=partial_T²+kphys²",
            "mass_and_residue": "Generalized mass eigenvalues 0,q(1/G+1/F); relative kinetic coefficient GF/(G+F)>0 for G,F,q>0",
            "equal_Einstein_dictionary": "gamma=(h_g+h_f)/2,Delta_TT=(h_g-h_f)/4 have coefficients 2G,8G and relative mass 2q/G",
            "unequal_Einstein_guard": "Physical average source coordinate is not the massless eigenmode unless G=F; no source-decoupling verdict imported",
            "scope": "A real extended-parent flat quadratic countercontrol, not a rolling CD/M1 match, all-mode stability or UV/cutoff result",
        },
        "primary_source_boundary": {
            "source": "https://arxiv.org/pdf/1804.04671",
            "model": "Sections 2.1-2.3, equations 2.1-2.17",
            "elimination": "Sections 3-4, equations 3.1-3.8 and 4.11-4.16",
            "metric_vierbein_caution": "Section 5.1 explicitly retains first-order symmetrization differences",
            "scope": "No generic novelty claim, no import of full UV or a numerical cutoff from the title's ghost-free completion statement",
        },
        "verification_boundary": "Exact symbolic replay, independent Fraction full-matrix jets and coefficientwise polynomial identities, root covariant audits and written proof; not proof-assistant formalized",
        "not_established": [
            "No-vacuum exclusion across added beta4 terms, other links/potentials, derivative matter couplings or other trimetric/completion actions",
            "A contradiction of the source's degree-of-freedom construction, or a requirement that every consistent cosmology possess a Minkowski vacuum",
            "A flat-spectrum verdict from a truncated cancellation or a mass formula evaluated on a nonsolution",
            "Original beta2 composite-action identity, frozen physical matter-metric identity or old CD/M1 matching",
            "A full rolling solution, quantitative finite-source remainder, strong-coupling/cutoff hierarchy, nonlinear or quantum health of the extension",
            "A required physical-time trajectory between vacuum and bounce",
            "Closure of S6, either original P8 track, or P8",
        ],
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("Auxiliary-parent applicability certificate differs from exact replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.13.AUXILIARY: chosen-parent vacuum prerequisite and named beta4 countercontrol replay passed; general matching and P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
