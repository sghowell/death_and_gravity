"""Read-only exact nonminimal conformal-scalar QSEI/cosmological replay."""

import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine.verify import serialize
from p8a_maxwell_thermal import verify as prior

from . import audit, cosmology, independent, thermal

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "nonminimal-conformal-qsei-cosmology.json"
PARENT_SHA = "d77f1a4ef293c53c70bb7fd540b473c557c36018dfbea558ea1af9afc147375c"
CLAIM = "P8-A.20.NONMINIMAL_CONFORMAL_QSEI_COSMOLOGY"
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
    "field_qsei",
    "conformal_source",
    "cosmological_theorem",
    "coherent_obstruction",
    "actual_thermal_history",
    "verification_boundary",
    "not_established",
    "verdict",
)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_files():
    return (
        sorted(ROOT.glob("src/p8a_nonminimal/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


@cache
def prior_checks():
    if sha(prior.REPORT) != PARENT_SHA or thermal.PARENT_SHA != PARENT_SHA:
        raise ValueError(
            "The pinned generic ODE and geometric-history ancestry changed"
        )
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return {
        "A19_fully_rebuilt_with_A18_A17_A16_ancestry": PARENT_SHA,
        "only_geometric_sampler_tube_and_generic_ODE_lemmas_reused": True,
        "photon_field_QSEI_anomaly_state_or_coupling_substituted": False,
    }


def exact_checks():
    groups = audit.exact_groups()
    output = {}
    count = 0
    for group_name, group in groups.items():
        output[group_name] = {}
        for name, value in group.items():
            values = list(value) if isinstance(value, sp.MatrixBase) else [value]
            if any(sp.simplify(entry) != 0 for entry in values):
                raise ValueError("A nonminimal scalar exact identity failed: " + name)
            count += len(values)
            output[group_name][name] = (
                {"all_entries_exactly_zero": True, "shape": list(value.shape)}
                if isinstance(value, sp.MatrixBase)
                else "0"
            )
    return output, count


def checked_constants():
    primary = {
        "costs": cosmology.costs(),
        "gate": cosmology.gate(
            cosmology.DELTA_MAX, cosmology.ZETA_MAX, cosmology.SIGMA_MAX
        ),
        "thermal": thermal.calibration(),
    }
    alternate = independent.replay()
    if serialize(primary) != serialize(alternate):
        raise ValueError(
            "The independent Fraction scalar cost/state reconstruction differs"
        )
    return {**serialize(primary), "independent_Fraction_reconstruction_agrees": True}


def build_report():
    inputs = prior_checks()
    residuals, count = exact_checks()
    proofs = {name: bool(value) for name, value in audit.proof_checks().items()}
    if not all(proofs.values()):
        raise ValueError("A nonminimal scalar proof/domain gate failed")
    result = {
        "schema": 1,
        "claim": CLAIM,
        "date": "2026-09-09",
        "status": "COMPLETE_NAMED_NONMINIMAL_CONFORMAL_SCALAR_STATE_DEPENDENT_QSEI_COSMOLOGICAL_GATES_COHERENT_ABSOLUTE_BOUND_OBSTRUCTION_AND_ACTUAL_THERMAL_SEE_PAST; GENERAL_NONMINIMAL_NULL_SHARPNESS_AND_ORIGINAL_P8_OPEN",
        "prior_sha256": inputs,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/field.md",
            "notes/anomaly.md",
            "notes/cosmology.md",
            "notes/thermal.md",
            "notes/scope.md",
            "notes/sources.md",
        ],
        "exact_residuals": residuals,
        "named_exact_check_count": sum(map(len, residuals.values())),
        "checked_scalar_entries": count,
        "derived_constants": checked_constants(),
        "proof_checks": proofs,
        "controls": {
            "rejected_inputs": audit.rejected_inputs(),
            "future_complete_geometry_plus_scalar_conformal_vacuum_is_not_a_small_budget_SEE_witness": True,
            "future_state_amplitude_not_inferred_from_past": True,
        },
        "field_qsei": {
            "field": "one real free massless scalar with xi=1/6 on smooth positive global flat FLRW",
            "targets": "all Hadamard states; no homogeneity, quasifree or global Minkowski extension requirement",
            "flat_positive_square_range": "0<=xi<=1/4; cosmological conformal transport only at xi=1/6",
            "flat_quantum_coefficient": "hbar*(3-4*xi)/(48*pi^2)",
            "actual_physical_bound": "integral E_omega*f^2 >= -7*hbar/(144*pi^2)*||L_H f||^2 -(1/3)*integral w_omega*|G_H f|^2 +integral E_conf*f^2",
            "proper_clock_L": "f''-2*H*f'+(3*H^2/4-3*H'/2)*f",
            "proper_clock_G": "f'-3*H*f/2",
            "w_omega": "physical relative Wick square against the conformal vacuum on the actual strip",
            "one_sided_state_cap": "w_omega<=Phi_*^2, Phi_*^2>=0; not an absolute positive Wick square",
            "extra_a_squared_in_flat_Wick_penalty_measure_retained": True,
            "compact_H2_zero_value_and_derivative_trace_extension": True,
            "optimal_or_null_bound_claimed": False,
        },
        "conformal_source": {
            "convention": "g_FK=(+---), R=6*(H'+2*H^2), G_00=-3*H^2",
            "actual_reference": "hbar/(2880*pi^2)*(H3+beta_S*I), tensors defined in notes/anomaly.md",
            "new_finite_beta_S_not_identified_with_photon_or_older_scalar_parameters": True,
            "physical_zero_type_D_scalar_density": "hbar*H^4/(960*pi^2)",
            "full_trace_conservation_and_EED_checked": True,
            "additional_gravitational_R_squared_source_silently_discarded": False,
            "reference_trace_anomaly_omitted": False,
        },
        "cosmological_theorem": {
            "SEE": "G_FK+Lambda*g_FK=-kappa*(T_scalar+T_other), E_other>=ell",
            "actual_past": "x=s/tau in [-1/100,0], with smooth endpoint neighborhoods",
            "actual_C3_tube": "tau*H(tau*x)+p/(p/2-x), p in [1/2,2/3], derivative suprema <=(1/100,1/2,3,16)",
            "derived_past_contraction": "H<=-19/(10*tau)",
            "distinct_past_and_future_caps": "d=(21/10,9,70,800), c=(4,128,16384,1048576)",
            "future_geometry_hypothesis": "IF a normal reaches tau, |H^(j)(s)|<=c_j/(tau-s)^(j+1) on its whole segment",
            "state_hypothesis": "one-sided physical Wick cap on actual past, and IF reaching tau on its whole future segment",
            "gates": "delta*(1+abs(beta_S))<=1e-8, zeta=kappa*Phi_*^2/3<=1/5000, sigma<=5",
            "delta": "kappa*hbar/(8*pi^2*tau^2)",
            "sigma": "tau^2*(max(Lambda,0)+kappa*max(-ell,0))",
            "strict_worst_margin": "63325013/350000000 >9/50",
            "conclusion": "stipulated global proper-clock endpoint <=tau; timelike geodesic incompleteness and all timelike curve lengths from zero <=tau",
            "initial_pointwise_SEC_is_an_independent_premise": False,
            "this_specific_radiation_tube_exhibits_initial_SEC_violation": False,
            "future_caps_derived_from_past_or_effective_Newton_positivity": False,
            "inextendibility_or_endpoint_EFT_control_in_every_larger_spacetime": False,
        },
        "coherent_obstruction": {
            "local_exact_field": "A*(3+(t/tau)^2+|X|^2/(3*tau^2))",
            "actual_global_field": "compact smooth Cauchy cutoff equals polynomial data on radius 2*tau, vanishes outside radius 3*tau; exact wave evolution",
            "sampled_EED": "A^2/tau^2*((4/3)*(t/tau)^2-2)<=-2*A^2/(3*tau^2) for |t|<=tau",
            "state": "actual positive finite-norm coherent displacement; vacuum plus smooth field outer product",
            "each_state_has_finite_energy": True,
            "family_has_uniform_energy_or_field_strength_bound": False,
            "conclusion": "every fixed nonzero real compact average in (-tau,tau) is unbounded below over this Hadamard family",
            "actual_self_consistent_SEE_counterexample": False,
        },
        "actual_thermal_history": {
            "new_occupation_state": "one scalar oscillator with n(k)=(exp(beta_T*k)-1)^-1",
            "energy": "Q=hbar*pi^2/(30*beta_T^4)",
            "physical_relative_Wick_square": "hbar/(12*a^2*beta_T^2)",
            "named_prescription": "beta_S=Lambda=T_other=0, no independent nonzero gravitational R^2 term",
            "own_SEE_coupling": "lambda=delta/360, not 31*delta/180",
            "anchor": "a0=1,y0=2,Q=12*(1-4*lambda)/(kappa*tau^2)",
            "exact_clock": "x=Phi(2)-Phi(y), Phi=1/(2*y)+sqrt(lambda)*atanh(sqrt(lambda)*y)/2",
            "exact_scale": "a^4=4*(1-4*lambda)/(y^2*(1-lambda*y^2))",
            "past_actual_C3_and_state_gates_verified": True,
            "actual_initial_zeta_squared": "(20/9)*delta*(1-4*lambda)<1/5000^2",
            "endpoint_x_strict_interval": ["1/16", "1/4"],
            "endpoint_zeta_squared": "50; the future small field cap fails near the endpoint",
            "future_state_cap_verified_on_every_shorter_segment": False,
            "new_QEI_only_endpoint_argument": False,
        },
        "verification_boundary": "Native full-component stress, exact flat spectral moments, both proper-clock measures, scalar reference conservation/trace, generic continuous-ODE identities and independent Fraction costs. Source-pinned analytic Hadamard, positivity, finite-propagation, Sobolev and global square proofs are written, not proof-assistant formalized or independently peer reviewed. Parent field coefficients are not reused.",
        "not_established": [
            "State-independent all-Hadamard averaged lower bound without a state qualification; the explicit coherent family rules it out here",
            "Optimal or null bound, generic nonconformal coupling, mass, interactions, arbitrary spacetime or observed-universe theorem",
            "Future geometric or field-strength bounds from the actual short thermal past or effective Newton positivity",
            "New QEI-only thermal endpoint argument, fundamental endpoint EFT control, or inextendibility in every larger spacetime",
            "Original P8(a), P8(b) UV classification, V/G/B matching, full interacting coupled-system stress, or original P8 closure",
        ],
        "verdict": "A separately normalized conformally coupled scalar has an exact state-dependent QSEI and a cosmological-scale global-FLRW incompleteness calibration with positive margin above 9/50 under explicit quantum, physical field-strength and source gates. A finite-energy coherent family excludes deleting all state qualifications. A new one-scalar thermal SEE branch realizes the required short past and small field gate, not the conditional future cap. General nonminimal/null/sharp and original P8 obligations remain open.",
    }
    assert tuple(result) == REPORT_KEYS
    return result


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The complete nonminimal scalar certificate differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 A.20 nonminimal scalar QSEI/cosmological replay passed; general nonminimal/null/sharp and original P8 OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
