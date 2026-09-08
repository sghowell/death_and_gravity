"""Read-only vector-only local one-loop coefficient and reference-scale bound."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as affine
from p8_aligned_elimination import verify as parent

from . import bounds, kernel, potential

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"vector-local-one-loop.json"
PARENT_SHA = "4a7985892a70803382be53e0beb67b140a3ddcba0c81fd32423c29d267cd9980"
sha, serialize = affine.sha, affine.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_aligned_quantum/*.py"))
            +sorted(ROOT.glob("tests/*.py"))+sorted(ROOT.glob("*.md"))
            +sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen ordered elimination/source-remainder certificate changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {"S6_46_fully_rebuilt": PARENT_SHA,
            "actual_mass_tensor_all64_action_and_source_independent_vector_Hessian_rechecked": True,
            "no_frozen_action_or_finite_counterterm_modified": True}


def residuals():
    out = {**kernel.checks(), **potential.checks(), **potential.independent_sign_controls()}
    jets = potential.clock_jets()
    h = jets["h"]
    expected = {"pole_weight": (3, 20/(9*h), -4*(3645*h-236)/(2187*h**2)),
                "finite_weight": (-sp.Rational(5, 2), -22/(27*h), 2*(8019*h+4784)/(6561*h**2))}
    for name, values in expected.items():
        for key, value in zip(("value", "N_first", "N_second"), values, strict=True):
            out["actual_clock_"+name+"_"+key] = sp.factor(jets[name][key]-value)
    return out


def controls():
    bad = (True, False, 1.0, sp.Float(1), "1", sp.I, sp.oo, sp.nan, sp.Symbol("unproved"))
    calls = [lambda value=value: bounds.scale_bound(value, 1) for value in bad]
    calls += [lambda value=value: bounds.scale_bound(1, value) for value in bad]
    calls += [lambda value=value: bounds.scale_bound(value, 1) for value in (0, -1)]
    calls += [lambda value=value: bounds.scale_bound(1, value) for value in (0, -1)]
    calls += [lambda value=value: serialize(value) for value in (1.0, sp.Float(1), sp.oo, sp.nan)]
    rejected = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(calls):
        raise ValueError("An inexact or nonpositive physical-scale input was accepted")
    return {"rejected_inputs": rejected,
            "temporal_and_longitudinal_determinant_not_dropped": True,
            "dimensionally_continued_transverse_multiplicity_required_for_finite_part": True,
            "independent_cutoff_log_and_three_dimensional_sign_controls": True,
            "reference_pole_convention_not_silently_identified": True,
            "actual_clock_jets_not_a_full_quantum_constraint_Jacobian": True,
            "finite_MSbar_potential_not_a_UV_or_full_quantum_verdict": True}


@cache
def build_report():
    pins = prior_checks()
    identities = residuals()
    exact = affine.certify_residuals(identities)
    proofs = bounds.proof_checks()
    if not all(value is True for value in proofs.values()):
        raise ValueError("A local mass, continuous finite-potential or physical-scale check failed")
    data, jets = potential.coefficients(), potential.clock_jets()
    return {"schema": 1, "claim": "P8-S6.47.ALIGNED_QUANTUM", "date": "2026-09-07",
            "status": "VECTOR_ONLY_LOCAL_ONE_LOOP_COEFFICIENT_AND_SCALE_BOUND_CERTIFIED; ORIGINAL_P8_OPEN",
            "prior_sha256": pins,
            "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in source_files()},
            "formulation": "FORMULATION.md", "written_proofs": ["notes/proof.md"],
            "exact_residuals": exact, "named_exact_check_count": len(identities),
            "checked_scalar_entries": len(identities), "proof_checks": proofs,
            "literal_action": "Unchanged S6.42; vector-only Gaussian at fixed light fields, canonical Wc=sqrt(zeta_physical)*W, m0²=M²/zeta_physical, a=1/gamma_t and b=1/gamma_s",
            "complete_constant_coefficient_kernel": serialize(kernel.determinant()["K"]),
            "complete_determinant": serialize(kernel.determinant()["expected"]),
            "regulator_and_scheme": "Euclidean +1/2 Tr log K; d=4-2epsilon_DR, d-2 transverse plus one longitudinal mode; scaleless contact determinants vanish in dimensional regularization; modified minimal subtraction with explicit mu scale",
            "pole_and_finite_weights": serialize({"C": data["pole_weight"], "B": data["finite_weight"]}),
            "potential_normalization": "V_pole=-m0^4*C/(64pi² epsilon_DR), V_MSbar=m0^4*B/(64pi²), Lscale=log(m0²/mu²); mu*dV_MSbar/dmu=-m0^4*C/(32pi²) with background coefficients fixed",
            "actual_clock_coefficient_jets_at_mu_m0": serialize({name: jets[name] for name in ("pole_weight", "finite_weight")}),
            "continuous_domain_and_bound": "Original closed X tube: 8/9<a<19/18 and 35/36<b<19/18; C<4, and at mu=m0 |B|<=3971/1296<4; this is the zero-curvature/zero-mass-derivative local term, not a global constant-gradient solution",
            "physical_scale_example": serialize(bounds.scale_bound(10**12, 1000)),
            "reference_scale_boundary": "The ratio uses M²/tau² as a reference density, not a nonzero observable energy or an interacting cutoff. This local term is below 10^-14 for the example, not the complete loop correction.",
            "source_and_sign_boundary": "Source alignment does not remove the actual clock-dependent determinant. Pole sign follows the explicit Gamma integral plus independent cutoff-log and d=3 controls; no sign-normalized match to the reference's opposite printed epsilon definition is asserted.",
            "controls": controls(),
            "verdict": "The retained vector's first local one-loop coefficient, nontrivial clock jets and a scheme-explicit continuous reference-scale bound are established. The full curved-background and causal quantum matching problem remains open.",
            "not_established": ["Curvature, mass-derivative, nonlocal/in-in or higher-loop remainders, full quantum stress tensor or a renormalized bounce",
                                "Quantum-corrected constrained spectrum/cones, full affine-complement/light/gravity measure and loops, finite matching counterterms or interacting cutoff",
                                "A Lorentz-invariant vacuum, finite-gravity Regge remainder, V/G/B UV admissibility or original P8 closure"],
            "verification_boundary": "Exact full determinant, independent physical/dense modes, Gamma and cutoff controls, direct original-mass clock jets, and written continuous logarithm/scale proof in a stated one-loop scheme; not proof-assistant formalization"}


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The local vector one-loop report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.47.ALIGNED_QUANTUM replay passed; vector-only local term, original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
