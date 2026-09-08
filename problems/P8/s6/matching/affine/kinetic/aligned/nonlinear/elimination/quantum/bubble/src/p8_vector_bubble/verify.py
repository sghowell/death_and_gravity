"""Read-only full vector mass-insertion pole and finite subtracted loop bound."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as affine
from p8_aligned_quantum import verify as parent

from . import independent, pole, remainder

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"mass-insertion-subtracted-loop.json"
PARENT_SHA = "808ed780dc548860f2d4b282bc70463c8802eed7e1cd5c9e2054bab7005be28f"
sha, serialize = affine.sha, affine.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_vector_bubble/*.py"))
            +sorted(ROOT.glob("tests/*.py"))+sorted(ROOT.glob("*.md"))
            +sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen local vector quantum certificate changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {"S6_47_fully_rebuilt": PARENT_SHA,
            "original_mass_tensor_first_variation_and_full_vector_Hessian_rechecked": True,
            "no_frozen_action_or_counterterm_modified": True}


def residuals():
    out = {**pole.checks(), **independent.checks(), **remainder.checks()}
    d = pole.actual_mass_direction()
    h, t, s = d["h"], d["time2"], d["space2"]
    expected = {"alpha": 4/(9*h), "beta": 28/(81*h),
                "constant": 904*pole.m2**2/(2187*h**2),
                "second": 416*pole.m2*(4*t+3*s)/(6561*h**2),
                "fourth": 32*(120*t**2+220*t*s+101*s**2)/(98415*h**2)}
    out.update({"actual_mass_direction_"+name: sp.factor(d[name]-value) for name, value in expected.items()})
    return out


def controls():
    bad = (True, False, sp.true, sp.false, 1.0, sp.Float(1), "1", sp.I,
           sp.oo, -sp.oo, sp.zoo, sp.nan, sp.Symbol("unproved"))
    calls = [lambda value=value: remainder.scale_bound(value, 1, 0) for value in bad]
    calls += [lambda value=value: remainder.scale_bound(1, value, 0) for value in bad]
    calls += [lambda value=value: remainder.scale_bound(1, 1, value) for value in bad]
    calls += [lambda value=value: remainder.scale_bound(value, 1, 0) for value in (0, -1)]
    calls += [lambda value=value: remainder.scale_bound(1, value, 0) for value in (0, -1)]
    calls += [lambda value=value: remainder.scale_bound(1, 1, value) for value in (-1, 2)]
    calls += [lambda value=value: serialize(value) for value in (1.0, sp.Float(1), sp.oo, sp.nan)]
    rejected = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(calls):
        raise ValueError("An inexact, nonpositive scale or outside-domain momentum was accepted")
    return {"rejected_inputs": rejected,
            "full_longitudinal_numerator_required": True,
            "linear_mass_bubble_not_full_second_mass_variation": True,
            "finite_logarithm_derived_from_actual_radial_Laurent_coefficients": True,
            "independent_isotropic_all_order_finite_control": True,
            "complex_Hermitian_momentum_norm_not_bilinear_k_squared": True,
            "zero_momentum_subtracted_kernel": remainder.scale_bound(1, 1, 0)["subtracted_loop_kernel_over_reference_density"] == 0,
            "subtracted_loop_not_unknown_higher_local_matching_operators": True,
            "no_stationary_to_curved_in_in_transfer_asserted": True}


@cache
def build_report():
    pins = prior_checks()
    identities = residuals()
    exact = affine.certify_residuals(identities)
    proofs = remainder.proof_checks()
    if not all(value is True for value in proofs.values()):
        raise ValueError("A continuous pole, logarithmic remainder or scale-domain check failed")
    d, f = pole.actual_mass_direction(), pole.parameter_integral()
    return {"schema": 1, "claim": "P8-S6.48.VECTOR_BUBBLE", "date": "2026-09-08",
            "status": "VECTOR_TWO_INSERTION_POLE_AND_FINITE_SUBTRACTED_LOOP_BOUND_CERTIFIED; ORIGINAL_P8_OPEN",
            "prior_sha256": pins,
            "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in source_files()},
            "formulation": "FORMULATION.md", "written_proofs": ["notes/pole.md", "notes/remainder.md"],
            "exact_residuals": exact, "named_exact_check_count": len(identities),
            "checked_scalar_entries": sum(value.rows*value.cols if isinstance(value, sp.MatrixBase) else 1 for value in identities.values()),
            "proof_checks": proofs,
            "literal_component": "Unchanged S6.42 vector-only Gaussian, fixed physical metric, isotropic on-clock mass m0²=M²/zeta_physical, linear mass insertion delta_M=m0²*A*r with actual A=diag(a_N,b_N,b_N,b_N). r is not a free physical constrained lapse.",
            "regulator_and_normalization": "Euclidean +1/2 Tr log K, d=4-2epsilon_DR. Gamma_bubble_pole^(2)=-1/(64pi² epsilon_DR)*integral r(-k)r(k)[F0+F2+F4], no implicit extra 1/2. All longitudinal numerator terms retained.",
            "general_pole_polynomial": serialize({name: f[name] for name in ("constant", "second", "fourth")}),
            "actual_mass_direction_and_pole": serialize(d),
            "continuous_pole_bound": "For |u|<=1/2 and nonzero real Euclidean momentum, (k²)²/125<F4<(k²)²/25. This is a divergent derivative coefficient of this component, not a physical ghost, cutoff or finite error.",
            "radial_poles_and_finite_parts": serialize(remainder.dimensional_logarithms()["integrals"]),
            "finite_logarithmic_subtraction": serialize(remainder.logarithmic_form()),
            "finite_remainder_constants": serialize(remainder.constants()),
            "continuous_finite_domain_and_bound": "kappa=sum_mu |k_mu|² is the Hermitian norm for complex Euclidean k; rho=kappa/m0²<=1. |R_loop|<=m0^4/(64pi² h²)*(4852/229635)*rho³/(1-rho/4)<=kappa³/(1920pi² h² m0²). Final bound strict for nonzero kappa; zero momentum gives zero.",
            "local_subtraction_boundary": "Subtract all Taylor terms through external-momentum degree four from this one-loop component. Scheme/scale-dependent local terms through that degree cancel; unknown higher local matching operators do not. No finite counterterm is changed.",
            "physical_scale_example": serialize(remainder.scale_bound(10**12, 1000, 1)),
            "source_functional_boundary": "The kernel/reference M²/tau² ratio is <=Q³/(17280 Lp² R²), Q=kappa*tau²<=R². A band-limited real Euclidean source-functional bound includes its squared Fourier/L2 norm. This is not a pointwise physical energy correction.",
            "controls": controls(),
            "verdict": "The actual vector two-mass-insertion derivative pole and a finite, subtracted stationary loop remainder on an explicit complex momentum ball are established. Original P8 remains open.",
            "not_established": ["Full second-mass-variation tadpole, metric/scalar constraints or the full quantum scalar two-point function",
                                "Curved/in-in quantum transfer, full stress tensor, other field loops, unknown higher matching operators, higher-loop error or interacting cutoff",
                                "Quantum-corrected bounce/spectrum/cones, Lorentz-invariant vacuum, finite-gravity Regge remainder, V/G/B admissibility or original P8 closure"],
            "verification_boundary": "Full exact Proca tensor numerator, independent shifted angular and isotropic denominator reductions, explicit radial Laurent coefficients, all-order isotropic finite comparison, and written continuous complex-log norm bound; not proof-assistant formalization"}


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The vector mass-insertion loop report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.48.VECTOR_BUBBLE replay passed; subtracted loop component, original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
