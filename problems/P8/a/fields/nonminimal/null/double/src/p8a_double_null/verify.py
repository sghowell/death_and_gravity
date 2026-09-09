"""Read-only physical double-null scalar inequality certificate."""

import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine.verify import serialize
from p8a_null_cap import verify as prior

from . import audit, boost, curved, flat, independent

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "physical-double-null-qei.json"
PARENT_SHA = "ada5315159f41bf96c455987f4779dffae2fa6bfef3fca0921e60761fccc4e88"
CLAIM = "P8-A.22.PHYSICAL_DOUBLE_NULL_QEI"
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
    "physical_flat_inequality",
    "normalized_product_bound",
    "curved_conformal_transport",
    "verification_boundary",
    "not_established",
    "verdict",
)


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def source_files():
    return (
        sorted(ROOT.glob("src/p8a_double_null/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


@cache
def prior_checks():
    if sha(prior.REPORT) != PARENT_SHA:
        raise ValueError(
            "The pinned null obstruction and physical scalar ancestry changed"
        )
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return {
        "A21_and_A20_ancestry_fully_rebuilt": PARENT_SHA,
        "physical_nonminimal_stress_and_scalar_reference_reused": True,
        "single_null_cap_bound_inferred_from_parent": False,
    }


def residuals():
    return {
        **flat.data()["checks"],
        **boost.data()["checks"],
        **curved.data()["checks"],
    }


def exact_checks():
    rows = residuals()
    count = 0
    for name, value in rows.items():
        entries = list(value) if isinstance(value, sp.MatrixBase) else [value]
        for entry in entries:
            if sp.simplify(entry) != 0:
                raise ValueError("A physical double-null identity failed: " + name)
            count += 1
    return dict.fromkeys(rows, "0"), count


def checked_constants():
    primary = serialize(audit.constants())
    if primary != serialize(independent.replay()):
        raise ValueError("The independent Fraction product reconstruction differs")
    return {**primary, "independent_Fraction_reconstruction_agrees": True}


def build_report():
    inputs = prior_checks()
    rows, count = exact_checks()
    checks = {k: bool(v) for k, v in audit.proof_checks().items()}
    if not all(checks.values()):
        raise ValueError("A double-null analytic proof gate failed")
    result = {
        "schema": 1,
        "claim": CLAIM,
        "date": "2026-09-09",
        "status": "COMPLETE_EXPLICIT_NONOPTIMAL_PHYSICAL_DOUBLE_NULL_DIFFERENCE_QEI_WITH_ONE_SIDED_WICK_CAP_AND_CONFORMAL_FLRW_TRANSPORT; NULL_FOCUSING_AND_ORIGINAL_P8_OPEN",
        "prior_sha256": inputs,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/plane.md",
            "notes/boost.md",
            "notes/curved.md",
            "notes/scope.md",
            "notes/sources.md",
        ],
        "exact_residuals": rows,
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": count,
        "derived_constants": checked_constants(),
        "proof_checks": checks,
        "controls": {
            "rejected_inputs": audit.rejected_inputs(),
            "quantum_stress_not_effective_gravitational_stress": True,
            "no_transversely_smeared_field_square_substitution": True,
        },
        "physical_flat_inequality": {
            "field": "one real massless scalar in four-dimensional Minkowski space; 0<=xi<=1/2",
            "plane": "(t,0,0,z), physical null D=partial_t+partial_z, plane volume dt*dz",
            "states": "every Hadamard target, relative to the Minkowski vacuum",
            "samplers": "real C_c^infinity(R^2), or the stated compact H0^2 product limit",
            "physical_operator": "T_DD=(D Phi)^2-xi*D^2(Phi^2)",
            "integrated_square_identity": "(1-2xi)*integral f^2*(D Phi)^2+2xi*integral [D(f Phi)]^2-2xi*integral (Df)^2*Phi^2",
            "vacuum_plane_spectrum": "hbar/(8*pi^2)*integral_{omega>=|kz|}domega*dkz*exp(-iomega*Delta_t+ikz*Delta_z)",
            "inequality": "integral f^2*:T_DD: >= -hbar/(8*pi^2)*[(2/3)(1-xi)||f_tt||^2+(4xi/3)<f_tt,f_tz>+2xi||f_tz||^2]-2xi*integral w*(Df)^2",
            "one_sided_state_cap": "w<=Phi_*^2>=0 on the entire two-dimensional sampling region gives -2xi*Phi_*^2*||Df||^2",
            "derivative_positive_decomposition": "2xi*(f_tz+f_tt/3)^2+(2/9)(3-4xi)*f_tt^2",
            "absolute_Wick_cap_or_momentum_cutoff_inferred": False,
        },
        "normalized_product_bound": {
            "coordinates": "x_plus=t+z,x_minus=t-z; dt*dz=dx_plus*dx_minus/2",
            "sampler": "f=sqrt(2)*h_delta_plus(x_plus)*h_delta_minus(x_minus), h_delta=delta^-1/2*h(x/delta)",
            "actual_profile": "h(s)=sqrt(315/256)*(1-s^2)^2 on [-1,1], zero elsewhere; norm one, first norm squared 3, second 63/2",
            "boost": "U=b*partial_plus+b^-1*partial_minus,E=b*partial_plus-b^-1*partial_minus,Dprime=b*D; b^2=r*delta_plus/delta_minus",
            "conformal_boost_family": "Gamma(r)=35r+24/r+21/r^3, multiplying hbar/(8*pi^2*delta_plus^3*delta_minus)",
            "unique_family_minimizer_squared": "(12+9sqrt(29))/35; not the optimal physical QEI",
            "selected_rational_ratio": "r=4/3, Gamma=14117/192<74",
            "explicit_conformal_bound": "integral f^2*:T_DD: >= -14117*hbar/(1536*pi^2*delta_plus^3*delta_minus)-4*Phi_*^2/delta_plus^2",
            "single_null_limit": "fixed delta_plus and delta_minus->0 diverges; no finite null-line bound follows",
            "optimal_sampler_or_optimal_physical_QEI_claimed": False,
        },
        "curved_conformal_transport": {
            "domain": "smooth positive spatially flat FLRW a(t), conformal massless scalar xi=1/6 only, compact timelike-plane sample",
            "volume_and_null_vector": "dvol_plane=a*dt*dz=a^2*deta*dz; K=a^-2*(partial_eta+partial_z), affine with this normalization",
            "relative_physical_stress": "Delta T_KK=a^-6*Delta Tflat_DD, flat sampler F=a^-2*f",
            "operators": "L=f_tt-3H*f_t+(2H^2-2Hdot)*f; M=f_tz-2H*f_z; G=f_t+a^-1*f_z-2H*f",
            "quantum_cost": "hbar/(8*pi^2)*integral dvol_plane*[(5/9)a^-2*L^2+(2/9)a^-3*L*M+(1/3)a^-4*M^2]",
            "state_cost": "(1/3)*integral dvol_plane*a^-2*w_physical*G^2, bounded above by replacing w_physical with Phi_*^2",
            "reference_added_to_full_stress_lower_bound": "integral dvol_plane*f^2*Tconf_KK",
            "actual_reference": "Tconf_KK=hbar/(2880*pi^2*a^2)*[-4H^2*Hdot+beta_S*(12Hthird+72Hdot^2+36H*Hddot)]",
            "FK_null_Ricci": "R_KK=2Hdot/a^2; cosmological constant drops from null SEE",
            "single_ray_focusing_or_SEE_incompleteness_deduced": False,
        },
        "verification_boundary": "Native physical sampler identities, the full mass-shell and half-Fourier normalization, exact derivative matrix and product Gram, independent Fraction profile/cost reconstruction, affine conformal weights and actual scalar-reference contraction. Positive-type restriction, compact Sobolev density and curved transport are written source-pinned arguments, not proof-assistant formalized or independently peer reviewed.",
        "not_established": [
            "An optimal QEI, optimal sampler, finite single-null bound or state-independent nonminimal bound",
            "General curved nonconformal fields, interacting fields or a UV momentum cutoff",
            "A null singularity theorem or self-consistent SEE solution from this two-dimensional average",
            "Original P8(a) or full V/G/B UV classification P8(b) closure",
        ],
        "verdict": "Smearing the actual physical null stress over both null directions gives an explicit finite nonoptimal lower bound with a one-sided Wick upper cap in the stated coupling range. A normalized compact conformal example retains both widths and has coefficient 14117/(1536*pi^2), and the physical affine-null FLRW transport retains the full scalar reference. The single-null limit diverges and a two-dimensional average alone is not a null focusing theorem; original P8 remains open.",
    }
    assert tuple(result) == REPORT_KEYS
    return result


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete physical double-null certificate differs from replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), result)
        print(
            "P8 A.22 physical double-null QEI replay passed; null focusing and original P8 OPEN"
        )
    else:
        print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
