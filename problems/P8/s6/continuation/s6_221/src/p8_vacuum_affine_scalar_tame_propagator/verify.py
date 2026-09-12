"""Read-only uniform two-chart classical scalar propagator certificate."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8 import gamma as historical
from p8_affine import verify as certificate
from p8_vacuum_affine_reduced_scalar_hamiltonian import verify as parent

from . import audit, charts, coefficients, energy, majorants

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "certificates" / "polynomial-vacuum-affine-scalar-tame-propagator.json"
PARENT_SHA = "cfd123d4faed1ce2efddb229931f6af79df203a7ea3c5af091d9d0e1857db4d8"
HISTORICAL_SHA = "5e411e649be38bcdea745e8602840fd1e1fe774608798a6948ceb88dbfe35bfe"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_scalar_tame_propagator/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    if (
        sha(parent.REPORT) != PARENT_SHA
        or sha(Path(historical.__file__)) != HISTORICAL_SHA
    ):
        raise ValueError("A source-pinned classical scalar input changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {
        "S6_220_regular_full_retuned_scalar_coefficient_system_fully_rebuilt": PARENT_SHA,
        "earlier_core_gamma_complementary_principal_chart_source": HISTORICAL_SHA,
        "S6_182_fixed_QG1_and_original_five_stress_jets_unchanged": True,
        "no_frozen_scientific_proof_test_or_report_bytes_changed": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(v is True for v in gates.values()):
        raise ValueError("A complete chart or tame energy gate failed")
    return {
        "schema": 1,
        "claim": "P8-S6.221.UNIFORM_TWO_CHART_QG1_CLASSICAL_SCALAR_PHASE_PROPAGATOR_WITH_TWELVE_SPATIAL_DERIVATIVE_LOSS",
        "date": "2026-09-12",
        "status": "UNIFORM_POLYNOMIAL_MOMENTUM_CLASSICAL_COEFFICIENT_PROPAGATOR; NOT_SAME_SPACE_QUANTUM_INVERSE_SMALLNESS_NONLINEAR_STABILITY_CUTOFF_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/coefficients.md",
            "notes/charts.md",
            "notes/majorants.md",
            "notes/energy.md",
            "notes/comparison.md",
            "notes/scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "actual_coefficients_complete_charts_and_uniform_majorants": serialize(
            {
                "coefficients": payload(coefficients.data()),
                "charts": payload(charts.data()),
                "majorants": payload(majorants.data()),
            }
        ),
        "energy_and_all_transfer_tame_phase_comparison": serialize(
            payload(energy.data())
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "The actual unit-slab fixed-QG1 two-scalar coefficient-sector principal pair has speeds squared1 andF/Jc, with F>1/100 and Jc-F>1/1000. Earlier core gamma algebra is explicitly reused, not represented as new. The complete finite-transfer retuned Hamiltonian admits an outer v chart where |Theta|>=3/10 and a complementary b=-pv/(2q) chart where |E|>1/4. The weighted canonical boundary, all A,Tcorr,Jc terms and q'=-2Hq are retained. At physical |P|>=100 the central normalized auxiliary denominator lies strictly between1/8 and2. Both complete kinetic and principal gradient matrices lie between1e-8 I and1e4 I. Exact polynomial enclosures bound all48 full coefficient/remainder entries by1e18 using only actual profile jets through2. The complete energy identity, phase conversions and at most two switches yield ||U(t,s;P)||<=exp(1e29)(1+|P|^2)^6, while smaller momenta use the original regular first-order system. Thus the classical coefficient-sector comparison extends uniformly from H^(r+12) initial data and L1 H^(r+12) additive phase forcing to C_t H^r for every real r. Independent finite-q saddle solves, original canonical-flow versus complete Euler comparisons, coefficient and energy fixtures and deleted-boundary/finite-q-degeneracy controls supplement the written continuum proof. The enormous constant is not a smallness bound, twelve derivatives are lost, and the nonlocal quantum lapse/shift constraints are not solved. Full quantum inversion, physical forcing maps, nonlinear sourced-parent remainder, background/stability, physical cutoff/heavy control and original V/G/B/P8 remain open.",
        "not_established": [
            "A same-space scalar propagator, kappa-smallness or a contraction for the derivative-losing Gaussian response",
            "Elimination of full nonlocal quantum constraints or a physical-source-to-phase forcing map",
            "A numerical high-time-jet stress bound or a literal global homogeneous constraint from P=0 assignment",
            "Finite nonlinear parent remainder, quantum background/stability, heavy/physical cutoff or original V/G/B/P8 closure",
        ],
        "verification_boundary": "Exact complete finite-transfer reductions, actual coefficient inequalities, positive normalized denominator,48 polynomial enclosures, energy identity and conversion arithmetic are checked. Independent dynamics and variational fixtures supplement them. Chart patching, Gronwall, Duhamel, Plancherel and dominated Sobolev continuity are written proofs, not FORMALIZED. Native/direct/ordinary/CLI retain original SymPy; only separately audited full regression uses the exact-GCD adapter. Frozen predecessors remain unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The two-chart tame propagator report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.221 two-chart tame classical propagator replay passed; full quantum inverse and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
