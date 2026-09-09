"""Read-only actual scalar double-null QEI FLRW incompleteness certificate."""

import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine.verify import serialize
from p8a_double_null import verify as prior

from . import affine, audit, control, costs, independent, thermal

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "double-null-qei-flrw-incompleteness.json"
PARENT_SHA = "7f37845da53396f5bbc171c43572a07eb10c493c135247b58fe255784c50b15b"
CLAIM = "P8-A.23.DOUBLE_NULL_QEI_FLRW_INCOMPLETENESS"
REPORT_KEYS = (
    "schema",
    "claim",
    "date",
    "status",
    "prior_sha256",
    "source_sha256",
    "formulation",
    "written_proofs",
    "exact_residuals",
    "named_exact_check_count",
    "checked_scalar_entries",
    "derived_constants",
    "proof_checks",
    "controls",
    "null_theorem",
    "physical_plane_to_geometry",
    "actual_null_index",
    "complete_geometry_control",
    "actual_thermal_past",
    "verification_boundary",
    "not_established",
    "verdict",
)


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def source_files():
    return (
        sorted(ROOT.glob("src/p8a_null_focusing/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


@cache
def prior_checks():
    if sha(prior.REPORT) != PARENT_SHA:
        raise ValueError(
            "The pinned physical double-null QEI and scalar ancestry changed"
        )
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return {
        "A22_and_all_scalar_ancestry_fully_rebuilt": PARENT_SHA,
        "actual_physical_conformal_QEI_and_reference_used": True,
        "null_line_QEI_or_state_homogeneity_assumed": False,
    }


def residuals():
    return {
        **affine.data()["checks"],
        **costs.data()["checks"],
        **control.data()["checks"],
        **thermal.data()["checks"],
    }


def exact_checks():
    rows = residuals()
    for name, v in rows.items():
        if sp.simplify(v) != 0:
            raise ValueError("An actual null-focusing identity failed: " + name)
    return dict.fromkeys(rows, "0")


def checked_constants():
    data = {
        "costs": {k: v for k, v in costs.data().items() if k != "checks"},
        "complete_geometry_control": {
            k: v for k, v in control.data().items() if k not in ("checks", "A")
        },
    }
    if serialize(data) != serialize(independent.replay()):
        raise ValueError(
            "The independent Fraction null-focusing reconstruction differs"
        )
    return {**serialize(data), "independent_Fraction_reconstruction_agrees": True}


def build_report():
    parent = prior_checks()
    rows = exact_checks()
    proofs = {k: bool(v) for k, v in audit.proof_checks().items()}
    if not all(proofs.values()):
        raise ValueError("An actual null-focusing proof gate failed")
    result = {
        "schema": 1,
        "claim": CLAIM,
        "date": "2026-09-09",
        "status": "COMPLETE_CONDITIONAL_GLOBAL_FLRW_NULL_AFFINE_INCOMPLETENESS_FROM_ACTUAL_NONMINIMAL_SCALAR_DOUBLE_NULL_QEI_WITH_MACROSCOPIC_BUDGETS; GENERAL_SPACETIME_AND_ORIGINAL_P8_OPEN",
        "prior_sha256": parent,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/plane.md",
            "notes/index.md",
            "notes/constants.md",
            "notes/control.md",
            "notes/thermal.md",
            "notes/scope.md",
        ],
        "exact_residuals": rows,
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": len(rows),
        "derived_constants": checked_constants(),
        "proof_checks": proofs,
        "controls": {
            "rejected_inputs": audit.rejected_inputs(),
            "complete_geometry_and_conformal_vacuum_not_claimed_small_source_SEE": True,
            "actual_thermal_past_not_claimed_to_prove_future_caps": True,
        },
        "null_theorem": {
            "global_geometry": "M=I_t times R^3, smooth positive flat FLRW; a(0)=1, H(0)<=-2/tau, tau>0",
            "affine_clock": "lambda(t)=integral_0^t a(s)ds, K=a^-1 partial_t+a^-2 partial_z; x=lambda/tau and A(x)=a(t(lambda))",
            "actual_past": "x in [-1/1000,0],9/10<=A<=2, absolute first-through-fourth affine x derivatives <=(3,16,1000,2000000)",
            "conditional_future": "IF affine endpoint exceeds 2tau, x in [0,2] has 1/2<=A<=2 and absolute first-through-fourth derivatives <=(50,10000,2000000,400000000)",
            "actual_SEE": "G_FK+Lambda*g_FK=-kappa*(Tscalar+Tother), kappa>0; Lambda arbitrary fixed and drops from null contraction",
            "actual_field": "one real massless conformal scalar xi=1/6 in any Hadamard target state; fixed A20 scalar scheme beta_S",
            "finite_plane_state_cap": "relative w_physical<=Phi_*^2>=0 on the fixed coordinate plane sample |z|<=tau, lambda in [-tau/1000,2tau]; future part conditional on extension",
            "other_source": "Tother_KK>=ell on the same finite plane sample; no other-source premise is needed for pure scalar",
            "budgets": "delta=kappa*hbar/(8*pi^2*tau^2),zeta=kappa*Phi_*^2/3,sigma=tau^2*kappa*max(-ell,0); delta*(1+abs(beta_S))<=1e-12,zeta<=1e-6,sigma<=1/10",
            "conclusion": "future null affine endpoint <=2tau with the stated normalization; null geodesic incompleteness of the stipulated global FLRW spacetime",
            "pointwise_NEC_SEC_or_homogeneous_quantum_state_is_a_premise": False,
            "future_caps_inferred_from_actual_past_or_effective_Newton_constant": False,
        },
        "physical_plane_to_geometry": {
            "actual_sampler": "f(t,z)=v(lambda(t))*j_tau(z); j is the A22 normalized clamped profile, spatial first norm squared 3/tau^2",
            "affine_sampler": "past cubic step on [-tau/1000,0], joined C1 to h=1-3(lambda/(2tau))^2+2(lambda/(2tau))^3 on [0,2tau]",
            "actual_measure": "dvol_plane=a*dt*dz=dlambda*dz; integral j_tau^2 dz=1",
            "actual_affine_operators": "P=a*v''-2a'*v'+2(a'^2/a-a'')*v; U=a*v'-2a'*v; V=v'-2a'/a*v, primes here d/dlambda",
            "full_quantum_cost": "hbar/(8*pi^2)*[(5/9)integral P^2+(1/3)*(3/tau^2)integral a^-4*U^2]",
            "full_state_cost": "Phi_*^2/3*[integral V^2+(3/tau^2)integral a^-4*v^2]",
            "actual_reference": "hbar/(2880*pi^2)*[-4a'^2*a''/a+beta_S*(48a'^2*a''/a+84a''^2+72a'*a'''+12a*a'''')]",
            "source_reduction": "actual SEE makes TOTAL null stress geometric and spatially homogeneous; lower bound on Tother gives the needed scalar comparison even for an inhomogeneous target",
            "spatial_width_removed_or_null_line_distribution_pulled_back": False,
        },
        "actual_null_index": {
            "outgoing_sphere": "initial radius tau; coordinate radius tau+eta, eta'=a^-2; positive screen scale b=a*(tau+eta)",
            "initial_expansion": "theta0=2*(H0+1/tau)<=-2/tau",
            "exact_Jacobi_equation": "b''/b=a''/a",
            "index": "I=integral_0^(2tau)[2h'^2+2(a''/a)h^2]dlambda",
            "positive_screen_identity": "I+theta0=2*integral b^2*((h/b)')^2>=0 for h(0)=1,h(2tau)=0",
            "sourced_dimensionless_upper": "tau*I<=6/5+104/7875+delta*(C0+abs(beta_S)*Cbeta)+zeta*Cstate+sigma*Csource",
            "strict_worst_case_margin": "2-tau*I>=2724477176468412151/4784062500000000000>1/2",
            "contradiction": "I<-theta0 contradicts positive outgoing screen scale if affine extension exceeds 2tau",
            "maximal_inextendibility_curvature_blowup_or_fundamental_EFT_endpoint_claimed": False,
        },
        "complete_geometry_control": {
            "actual_scale": "A(x)=1-exp(-100x)*Q(x), Q=2x+204x^2+(31240/3)x^3+(1064160/3)x^4",
            "domain": "x>=-1/1000 with a smooth extension slightly further to the past; positive and future complete",
            "future_uniform_lower": "247371/312500>1/2, A<=1, A tends to 1",
            "all_past_and_future_caps_verified_by_explicit_polynomial_bounds": True,
            "anchor_jets": "matches (1-6x)^(1/3) through fourth order, A0=1,A'_0=-2",
            "quantum_state_control": "actual conformal vacuum has relative Wick square zero",
            "allowed_small_source_SEE_witness": False,
        },
        "actual_thermal_past": serialize(
            {k: v for k, v in thermal.data().items() if k != "checks"}
        ),
        "verification_boundary": "Native actual affine transport and scalar reference, exact outgoing Jacobi/index identity and boundary, full finite spatial-width sampler integrals, exact cost and gate arithmetic, independent Fraction reconstruction, continuous complete-geometry bounds and actual scalar SEE past/affine conversion. The QEI-to-source argument, Sobolev limits and global geometric endpoint proof are written and pinned, not proof-assistant formalized or independently peer reviewed.",
        "not_established": [
            "General non-FLRW Penrose theorem, arbitrary nonminimal coupling, mass, interactions or optimal constants",
            "A single-null QEI, state-independent nonminimal bound, momentum cutoff or state homogeneity requirement",
            "Future geometry or state caps from past data, or a new QEI-only thermal endpoint mechanism",
            "Curvature blowup, maximal inextendibility, fundamental EFT control or original P8(b) V/G/B UV classification",
            "Original P8(a) or P8 as a whole closure",
        ],
        "verdict": "The actual nonminimal conformal scalar double-null QEI implies null affine incompleteness in the stated global FLRW class, with explicit macroscopic quantum, field and other-source budgets and a strictly positive index margin. Finite spatial smearing and the physical reference are retained; actual SEE supplies the geometry reduction without homogeneous-state assumptions. A complete geometry satisfies all geometric caps, while an actual thermal SEE state realizes the short past but not the conditional future premises. Original P8 remains open.",
    }
    assert tuple(result) == REPORT_KEYS
    return result


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete double-null FLRW theorem certificate differs from replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), result)
        print(
            "P8 A.23 double-null QEI FLRW null incompleteness replay passed; original P8 OPEN"
        )
    else:
        print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
