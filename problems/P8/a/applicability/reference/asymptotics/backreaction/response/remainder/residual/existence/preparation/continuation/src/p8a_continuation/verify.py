"""Read-only replay of the frozen causal block, not a SEE continuation."""

import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8a_preparation import verify as prior

from . import bounds, decomposition, independent, resolvent

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"causal-resolvent.json"
PIN = "8dba8215c5a28e2a6379c09e4d452838449065c396e10888aa733b186799456f"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def serialize(value):
    if isinstance(value, dict):
        return {key: serialize(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [serialize(item) for item in value]
    if isinstance(value, bool):
        return value
    return str(value)


@cache
def prior_checks():
    if sha(prior.REPORT) != PIN:
        raise ValueError("the pinned A11 map/state certificate changed")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return {"A11": PIN}


def exact_checks():
    groups = {"full_equation_and_map_decomposition": decomposition.identities(),
              "causal_multiplier_poles_cut_and_units": resolvent.identities()}
    for name, group in groups.items():
        if any(sp.simplify(value) != 0 for value in group.values()):
            raise ValueError("a frozen-block exact identity failed: "+name)
    return {name: dict.fromkeys(group, "0") for name, group in groups.items()}


def checked_constants():
    primary = serialize(bounds.calibration())
    fractions = independent.replay(json.loads(prior.REPORT.read_text()))
    if primary != fractions:
        raise ValueError("the independent Fraction frozen-resolvent replay disagrees")
    return {"rational_gate": primary, "independent_Fraction_replay": fractions,
            "named_symbolic_parameters": serialize(resolvent.named_parameters())}


def build_report():
    inputs = prior_checks()
    identities, constants = exact_checks(), checked_constants()
    sources = sorted(ROOT.glob("src/p8a_continuation/*.py"))+sorted(ROOT.glob("tests/*.py"))
    sources += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    return {
        "schema": 1, "claim": "P8-A.13", "date": "2026-09-06",
        "status": "EXACT_FROZEN_CAUSAL_BLOCK_CERTIFIED; ACTUAL_ROLLING_CONTINUATION_AND_P8_OPEN",
        "prior_sha256": inputs,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in sources},
        "formulation": "FORMULATION.md", "written_proof": "notes/resolvent.md",
        "full_map_dictionary": "notes/decomposition.md", "source_audit": "notes/sources.md",
        "next_falsifiable_gate": "notes/next-step.md",
        "exact_residuals": identities, "derived_constants": constants,
        "theorem": {
            "nature": "exact frozen constant-coefficient linear resolvent and exact preconditioner decomposition of the full A11 map",
            "full_rolling_Frechet_derivative": False,
            "model_changed": False,
            "scheme": "same massless minimal field, transported actual state, ordinary radiation, lambda=2sqrt(2) A eta_star^2, gamma=0, epsilon=1",
            "causal_block": "[log(s)+beta-c/s^2]/2, beta=gamma_E-19/30-log(af/2), c=af^2/(30delta)",
            "exact_inverse": "2/[log(s)+beta-c/s^2], principal logarithm, Laplace line right of positive pole",
            "pole_count": "exactly three on the principal s sheet: W0,W1,W-1 of 2c exp(2beta)",
            "pole_locations": "one positive simple pole, one strictly left-half-plane conjugate pair",
            "residues": "2p_j/(1+W_j)",
            "cut_density": "2/[(log(r)+beta-c/r^2)^2+pi^2]",
            "full_kernel_nonnegative": True,
            "nonnegativity_proof": "convergent positive Volterra-Neumann series in an exponential weight beyond the growing pole",
            "C0_norm": "integral_0^L Kf; same supremum on smooth flat input",
            "reference_calibration": "af=5/2, delta=10^-14; af is a preconditioner choice, not a new physical solution",
            "positive_pole_bracket": "10^6<p0<2*10^6",
            "norm_lower": "[exp(10^6 L)-1]/17-8/3",
            "duration_norm_examples": "norm>1000 at L=10^-5; norm>10^28 at L=10^-4, only for the frozen comparator",
            "weighted_norm": "at sigma=2*10^6, norm<=240/769; unweighting costs exp(sigma L)",
            "bounded_forcing_condition": "the causal response is bounded iff the complete forcing moment integral_0^infinity exp(-p0 t) g(t) dt vanishes",
            "condition_imposed_on_actual_state": False,
            "whole_remaining_map_retained": ["baseline residual", "rolling Einstein coefficient divided by delta", "anomaly history", "auxiliary a^2 h q difference", "actual nonlinear mode-response derivative and local coefficient difference"],
        },
        "negative_controls": serialize(resolvent.controls()),
        "verification_boundary": [
            "A11 certificate and complete earlier A lineage replayed without writes",
            "Complete trace and fixed-point difference algebra checked with independent covariant audit",
            "Pole count and contour/inverse/positivity arguments are written analytic proofs, not numerical root sampling",
            "All rational brackets, growth examples and weighted bounds independently reconstructed with Fraction",
            "No rounded special-function evaluation, a priori rolling stability assumption or new finite counterterm used",
        ],
        "not_established": [
            "macroscopic continuation or an a posteriori error bound for the actual A11 SEE solution",
            "nonlinear runaway, instability, generic no-go or a necessary maximum step length",
            "the actual remainder forcing has a nonzero or zero growing-pole projection",
            "the initial energy constraint alone selects a no-growing solution",
            "physical validity at the frozen comparator pole timescale",
            "a no-runaway final boundary prescription, order-reduced model or quantum state reset",
            "cosmological incompleteness or completion of original P8a or P8",
        ],
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("the frozen causal-resolvent certificate differs from exact replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8(a) A.13: full-map decomposition and frozen causal resolvent passed; actual continuation and P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
