"""Read-only replay of the scoped coupled M1 hard-channel tree estimate."""

import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_m1_control import verify as prior

from . import free, majorant

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"cd-tree-control.json"
PRIOR_SHA = "a95282c853fdcce4e29e1e221a71e058d230b4bba466bbdb627f093e4ac3f66b"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prior_checks():
    if sha(prior.REPORT) != PRIOR_SHA:
        raise ValueError("Pinned S5.7.CD certificate changed")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return PRIOR_SHA


@cache
def build_report():
    previous = prior_checks()
    residuals = free.checks()
    if any(sp.simplify(value) != 0 for value in residuals.values()):
        raise ValueError("A coupled initial-state/canonical propagation residual failed")
    negative = {}
    for key, value in free.negative_controls().items():
        entries = tuple(value) if isinstance(value, (list, sp.MatrixBase)) else (value,)
        if not any(sp.expand(entry) != 0 for entry in entries):
            raise ValueError(f"A coupled free-mode omission control vanished: {key}")
        negative[key] = list(map(str, entries))
    scale_residuals = majorant.scale_checks()
    if any(sp.cancel(value) != 0 for value in scale_residuals.values()):
        raise ValueError("A homogeneous action/canonical scale identity failed")
    matter_controls = majorant.negative_controls()
    if any(sp.expand(value) == 0 for value in matter_controls.values()):
        raise ValueError("A matter-source omission control vanished")
    result = majorant.build()
    sources = sorted(ROOT.glob("src/p8_m1_tree/*.py"))+sorted(ROOT.glob("tests/*.py"))
    sources += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    return {
        "schema": 1, "claim": "P8-S5.8.CD", "date": "2026-09-05",
        "status": "SCOPED_NONEXCEPTIONAL_FINITE_TIME_M1_CUBIC_QUARTIC_TREE_CONTROL; LOOPS_AND_UV_OPEN",
        "prior_S5_7_CD_sha256": previous,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in sources},
        "formulation": "FORMULATION.md", "written_proof": "notes/tree-bound.md",
        "coupled_free_initialization_and_CCR_residuals": dict.fromkeys(residuals, "0"),
        "coupled_free_omission_controls": negative,
        "homogeneous_scale_restoration_residuals": dict.fromkeys(scale_residuals, "0"),
        "matter_source_and_mixed_boundary_omission_controls": {key: str(value) for key, value in matter_controls.items()},
        "independent_free_energy_exponential_enclosure": free.exponential_enclosure(),
        "exact_positive_series_majorants": result,
        "operational_contract": {
            "background": "fixed CD/M1, unchanged physical matter metric and covariant tube",
            "window": "abs(t-t0)<=ell0/100 at any finite t0; ell0=tau*sqrt(1+u0^2); a0=1",
            "external_band_in_fixed_ell0_units": ["10^11", "10^12"],
            "all_nonempty_proper_subset_momenta": "norm>=10^11; total signed spatial momentum zero",
            "free_scalar_columns": "R=W(left)^1/2; U0=(2R)^-1/2; P0=-i*(R/2)^1/2; Ydot0=P0-Omega(left)*U0",
            "free_evolution": "exact coupled canonical Hamiltonian evolution, with two scalar and two TT initial-mode columns",
            "initial_state_meaning": "local Gaussian complex structure; not generally the ground state of E-P^T*Omega*Y",
            "canonical_phase_map": "Q=a^-3/2*T^-1*Y; p=a^-3/2*T^T*(P-S*Y); gamma swap and both mixed matter shifts retained",
            "observable": "connected cubic 1-to-2/2-to-1 and hard quartic-order 2-to-2 transition blocks",
            "quartic_order": "one H4 plus two time-ordered H3 insertions, including both scalar and both tensor internal columns",
            "quantum_ordering": "tree terms only; contractions within one vertex and other loop/order-dependent effects excluded",
            "norm": "Schur bounds in fixed-total-momentum fibers of continuum one/two-particle spaces on R^3",
            "target_block_norm": "1/1000", "scale": "one sufficient finite M*tau works at all center times",
            "parameters_chosen_for_this_scoped_gate_not_preregistered": True,
        },
        "verification_boundary": [
            "prior physical Hamiltonian and coupled free-control certificates are immutable inputs, replayed through S5.7",
            "positive-series induction bounds all retained physical kernel assignments without a complete expanded H4",
            "the matrix/free-ODE, finite-tree combinatoric and continuum Schur arguments are written proofs",
            "exact rational arithmetic certifies all operational scale and coefficient inequalities; no sampled-angle coverage claim",
        ],
        "not_established": [
            "an optimized or observationally realistic duration, a necessary cutoff, or control at a separately specified M*tau",
            "inclusive, soft, forward or zero-mode transitions; vacuum production or a full Fock-space operator norm",
            "a global vacuum, an in/out S matrix, or chart-independent/global fixed-comoving mode control",
            "higher tree orders, loops, radiative protection or a perturbative remainder bound to all orders",
            "nonlinear PDE, BKL or backreaction stability; UV matching/admissibility/completion; full P8 completion",
        ],
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("CD/M1 tree-control certificate differs from exact replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S5.8: scoped CD/M1 hard-channel cubic/quartic tree replay passed; loops and UV OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
