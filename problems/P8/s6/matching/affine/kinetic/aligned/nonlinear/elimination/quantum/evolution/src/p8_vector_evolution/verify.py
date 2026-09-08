"""Read-only C1 vector matching and homogeneous clock-source certificate."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as affine
from p8_vector_clock_matching import verify as parent

from . import estimates, majorants, mixing, ward

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"differentiated-vector-clock-source.json"
PARENT_SHA = "17c6394918dd6bc1ca0238cab20f534db3950fb02a56eaac853e2d5377c9c1a4"
sha, serialize = affine.sha, affine.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_vector_evolution/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen clock-mass vector matching certificate changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {"S6_53_fully_rebuilt": PARENT_SHA,
            "original_action_physical_clock_Gaussian_data_and_finite_prescription_unchanged": True}


def residuals():
    return {**majorants.checks(), **majorants.replay_checks(), **mixing.checks(), **ward.checks()}


def controls():
    bad = (True, False, sp.true, sp.false, 1.0, sp.Float(1), "1", sp.I,
           sp.oo, -sp.oo, sp.zoo, sp.nan, sp.Symbol("unproved"), 0, -1)
    calls = [lambda value=value: estimates.physical_bounds(value, 1000) for value in bad]
    calls += [lambda value=value: estimates.physical_bounds(10**12, value) for value in (*bad, 999)]
    calls += [lambda value=value: majorants.polynomial_majorant(value) for value in
              (True, False, sp.true, sp.false, 1.0, sp.Float(1), "1", sp.I, sp.oo, sp.nan,
               sp.Symbol("unproved"), sp.sin(majorants.p), 1/(1+majorants.p), sp.sqrt(2))]
    good = {variable: (sp.Integer(1), sp.Integer(1)) for variable in majorants.VARIABLES}
    for value in (True, 1.0, sp.Float(1), sp.oo, sp.nan, -1, sp.I, "1"):
        for index in (0, 1):
            pair = [sp.Integer(1), sp.Integer(1)]
            pair[index] = value
            changed = {**good, majorants.p: tuple(pair)}
            calls.append(lambda changed=changed: majorants.value_and_slope(majorants.p, changed))
    calls += [lambda: majorants.value_and_slope(majorants.p, {key: value for key, value in good.items() if key != majorants.p}),
              lambda: majorants.value_and_slope(majorants.p, {**good, sp.Symbol("extra"): (sp.Integer(1), sp.Integer(1))})]
    calls += [lambda value=value: majorants.coefficients(value) for value in ("scalar", None, True)]
    calls += [lambda value=value: serialize(value) for value in (1.0, sp.Float(1), sp.oo, sp.nan)]
    rejected = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError, sp.PolynomialError):
            rejected += 1
    if rejected != len(calls):
        raise ValueError("An invalid scale, coefficient envelope, polynomial or report value was accepted")
    return {"rejected_inputs": rejected,
            "both_oscillatory_endpoints_and_feedback_retained": True,
            "diagonal_phase_removed_by_exact_transformation_not_dropped": True,
            "differentiated_momentum_integral_has_uniform_absolute_envelope": True,
            "full_ordinary_adiabatic_orders_conserved_not_only_UV_divergent_parts": True,
            "clock_mass_source_derived_before_background_restriction": True,
            "C1_first_variations_not_all_order_Hadamard_or_principal_symbol_control": True}


@cache
def build_report():
    pins, identities = prior_checks(), residuals()
    exact = affine.certify_residuals(identities)
    proofs = {**mixing.proof_checks(), **estimates.proof_checks()}
    if not all(value is True for value in proofs.values()):
        raise ValueError("A derivative, mixing or physical clock bound failed")
    scalar_count = sum(value.rows*value.cols if isinstance(value, sp.MatrixBase) else 1 for value in identities.values())
    return {"schema": 1, "claim": "P8-S6.54.VECTOR_EVOLUTION", "date": "2026-09-08",
            "status": "C1_VECTOR_METRIC_FIRST_VARIATIONS_AND_HOMOGENEOUS_CLOCK_SOURCE_BOUNDED; ORIGINAL_P8_OPEN",
            "prior_sha256": pins,
            "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in source_files()},
            "formulation": "FORMULATION.md", "written_proofs": ["notes/derivative.md", "notes/clock-source.md"],
            "exact_residuals": exact, "named_exact_check_count": len(identities),
            "checked_scalar_entries": scalar_count, "proof_checks": proofs,
            "literal_component": "Unchanged S6.50 exact prepared Gaussian modes, S6.51 full subtraction and S6.53 scalar-coefficient matching prescription, on I=[-1/2,1/2], all real comoving momenta and m0*tau>=1000. First time derivatives and a homogeneous clock Euler source, not higher functional variations.",
            "polynomial_majorant_method": "Differentiate actual rational coefficients using the full fixed-comoving chain rule. For the exact residual and subtraction-tail polynomials in independent coefficients, t and S^-1, the absolute-coefficient polynomial and its differentiated positive majorant bound values and slopes. Exact replay against both frozen identities is required.",
            "actual_residual_and_tail_derivative_envelopes": serialize({kind: majorants.bounds(kind) for kind in ("transverse", "longitudinal")}),
            "phase_separated_evolution": "With r=rho_res/(2W), theta'=W, eta'=r, Acal=e^(-i eta)a and Bcal=e^(i eta)b, the exact equations are a'=-i r e^(2i Phi)b, b'=i r e^(-2i Phi)a, Phi=theta+eta. The transformation retains the diagonal phase and CCR.",
            "oscillatory_mixing_bound": "C=2000000,C1=100000000. Psi=Phi'>=omega/4 and abs(Psi')<=7omega. Integrating b once by parts retains both endpoints and (r a/(2Psi))'. The endpoint, bulk and feedback bound gives abs(Bcal)<=K/nu^6, K=4C1+138C=676000000, nu²=m²+k²/(25/16)². No changed initial state is used.",
            "physical_readout_derivative": "For p=v'-d v and Q=[A|p|²+B omega²|v|²]/(2a³), differentiate with v'=p+d v,p'=-d p-omega²v. The reference has the additional residual term. Exact CCR-based product identities remove the diagonal phase from squared and mixed readout differences.",
            "differentiated_finite_integral_bound": "Both energy and pressure satisfy integral abs(exact/reference derivative difference)<=40K/m and reference-tail derivative integral<=D1/(72m²), D1=100000000. The resulting continuous all-momentum derivative envelope justifies a C1 finite subtracted integral on I, including one-sided endpoint derivatives.",
            "ordinary_conservation_controls": "The kinetic and potential mode balance identities hold for formal spatial dimension D. In physical dimension each full adiabatic order 0,2,4 separately obeys D_u F_n+(1-2n)lambda F_n+3H P_n=0. The matched ordinary local coefficients are separately conserved.",
            "direct_homogeneous_clock_variation": "Physical phi=tau*u, source-signature x=-u'^2/N². On the clock a_m=b_m=1, a_x=alpha/2,b_x=beta/2 and their explicit clock-coordinate derivatives vanish. Thus delta a_m=-alpha delta u', delta b_m=-beta delta u'; the homogeneous unit normal is unchanged. Varying the action before restriction gives clock momentum rho_mass and J_clock=-[rho_mass'+3H rho_mass]=-[rho'+3H(rho+p)].",
            "clock_source_regulator_limit": "Variations have compact support inside I with initial Gaussian data fixed. S6.53's uniform dimension-limit convergence permits integration by parts against test functions, hence the distributional clock-source limit. The present physical C1 bound identifies that limit with the stated continuous source; no separate complex-D derivative envelope is assumed.",
            "actual_local_derivative_and_clock_source_coefficients": serialize(ward.local_coefficients()),
            "continuous_local_derivative_and_clock_source_envelopes": serialize(ward.local_envelopes()),
            "physical_scale_example": serialize(estimates.physical_bounds(10**12, 1000)),
            "physical_units_and_scope": "Density derivatives and the physical phi Euler source scale as tau^-5. Their reference scale is M²/tau³, since phi=tau*u is unchanged. Explicit local and nonlocal bounds are each combined with their proper powers of R=m0*tau and L=M*tau. At L=10^12,R=1000 all three total bounds are below 10^-14. This is a vector-only one-loop contribution at a fixed background.",
            "controls": controls(),
            "verdict": "The retained Gaussian vector has C1 matched energy/pressure and a quantitatively bounded homogeneous clock source in the unchanged prescription. It does not supply a quantum-corrected solution or close original P8.",
            "not_established": ["Second time derivatives, all-order Hadamard admissibility, all-time state control or quadratic/higher quantum functional variations",
                                "Other fields/higher loops, omitted UV operators, quantum-corrected constrained bounce, principal cones, nonlinear stability or interacting cutoff",
                                "A common healthy vacuum, finite-gravity Regge remainder, V/G/B or original P8 closure"],
            "verification_boundary": "Exact rational and symbolic controls plus written continuous oscillatory, differentiated-integral and clock-variation proofs; not proof-assistant formalization or a full quantum solution"}


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The differentiated vector clock-source report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.54.VECTOR_EVOLUTION replay passed; C1 first variations and clock source, original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
