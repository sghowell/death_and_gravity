"""Read-only full homogeneous metric local-response certificate."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as affine
from p8_vector_mass_remainder import verify as parent

from . import bounds, canonical, consistency, counterterms, proofs, radial
from . import controls as failures

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"finite-local-metric-response.json"
PARENT_SHA = "81a0419a304c689847d2005438080e584b7dba73d42db29780127da8c937d387"
sha, serialize = affine.sha, affine.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_vector_metric_local/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen finite prepared mass-response certificate changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {"S6_66_fully_rebuilt": PARENT_SHA, "original_clock_state_profiles_and_counterterm_prescription_retained": True}


def controls():
    calls = []
    for value in (True, False, sp.true, sp.false, 1.0, sp.Float(1), sp.Integer(1), "1", -1, 3):
        calls.extend((lambda value=value: radial.laurent("N", value),
                      lambda value=value: radial.matched("N", value),
                      lambda value=value: counterterms.density(value),
                      lambda value=value: bounds.physical_operator("energy", value)))
    for value in (True, False, sp.true, 0, 1, "n", "Zeta", "", None):
        calls.extend((lambda value=value: radial.laurent(value, 0),
                      lambda value=value: counterterms.operator(value, 0),
                      lambda value=value: consistency.coefficients("N", 0, value),
                      lambda value=value: canonical.data(value),
                      lambda value=value: bounds.physical_operator(value, 0)))
    for value in (True, False, sp.true, sp.Integer(0), 0.0, sp.Float(0), "0", -1, 2):
        calls.append(lambda value=value: counterterms.operator("N", 0, value))
    for value in (True, False, sp.true, 1.0, sp.Float(1), "1", sp.I, sp.oo, -sp.oo, sp.zoo, sp.nan, 0, -1):
        calls.extend((lambda value=value: bounds.operator_bound(value, 1000),
                      lambda value=value: bounds.operator_bound(10**24, value)))
    calls.extend(lambda value=value: bounds.operator_bound(10**24, value) for value in (1, 999))
    calls.extend(lambda value=value: serialize(value) for value in (1.0, sp.Float(1), sp.oo, sp.nan))
    rejected = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(calls):
        raise ValueError("An unsupported full metric local-response input was accepted")
    return {"rejected_inputs": rejected, "native_validation_before_cached_implementations": True,
            "second_mass_jets_and_mixed_metric_contacts_not_optional": True,
            "finite_potential_not_added_twice": True,
            "physical_stress_normalization_not_confused_with_action_Hessian": True,
            "finite_fourth_derivative_piece_not_a_particle_spectrum": True,
            "homogeneous_local_bound_not_a_nonlocal_or_no_loss_estimate": True}


@cache
def build_report():
    prior, identities = prior_checks(), proofs.residuals()
    certified, gates = affine.certify_residuals(identities), proofs.checks()
    if not all(value is True for value in gates.values()):
        raise ValueError("A full metric local-response audit gate failed")
    scalar_count = sum(value.rows*value.cols if isinstance(value, sp.MatrixBase) else 1 for value in identities.values())
    coefficients = {j: {output+source: consistency.coefficients(output, j, source)
                         for output in ("N", "Z") for source in ("N", "Z")} for j in (0, 1, 2)}
    return {"schema": 1, "claim": "P8-S6.67.FINITE_LOCAL_METRIC_RESPONSE", "date": "2026-09-08",
            "status": "FINITE_HOMOGENEOUS_LOCAL_METRIC_RESPONSE; ORIGINAL_P8_OPEN",
            "prior_sha256": prior, "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in source_files()},
            "formulation": "FORMULATION.md", "written_proofs": ["notes/vertices.md", "notes/matching.md", "notes/bound.md"],
            "exact_residuals": certified, "named_exact_check_count": len(identities), "checked_scalar_entries": scalar_count,
            "proof_checks": gates,
            "full_metric_vertices": serialize({sector: canonical.data(sector) for sector in ("T", "L")}),
            "finite_Euler_current_coefficient_rows": serialize(coefficients),
            "physical_stress_coefficient_majorants": serialize(bounds.envelopes()),
            "finite_local_scale_example": serialize(bounds.operator_bound(10**24, 1000)),
            "omitted_counterterm_failure_controls": serialize(failures.omitted_counterterm()),
            "second_mass_vertex_control": serialize(failures.second_mass_vertex()),
            "fourth_derivative_local_symbol_control": serialize(failures.top_derivative()),
            "dimensional_prescription": "Vary the full fixed C4 potential, ordinary fixed scalar-coefficient curvature terms, linear mass-deviation curvature terms and S6.63 quadratic derivative mass terms with the D-dimensional volume. Keep mass-weighted box R, D-1 polarizations and the first D jet. Combine radial finite and 2 partial_D of the common-normalization counterterm Hessian before D=3.",
            "source_and_output_domain": "Independent physical lapse N=1+n and logarithmic spatial scale a_hat=a exp(zeta), fixed scalar clock, spatial K=0, compactly supported time variations. Euler currents J_N=-a^D rho and J_Z=D N a^D p are varied before physical output normalization. Quantitative norms use I=[-1/2,1/2] and maximum source derivatives zero through four.",
            "controls": controls(),
            "verdict": "All six full homogeneous metric poles match the fixed subtraction. Its finite local two-current matrix is weighted self-adjoint, reproduces the full frozen potential and ordinary finite heat action, and yields a physical C4-to-C0 bound below 10^-38 at L=10^24,m=1000. The full metric state remainder and original P8 remain open.",
            "not_established": ["Finite nonlocal full metric response or arbitrary spatial momentum and initial state variation",
                                "No-loss feedback, quantum stability/cones, interactions or cutoff",
                                "Finite Wilson matching, a common UV parent, V/G/B or original P8 closure"],
            "verification_boundary": "Exact Hamiltonian, Riccati, radial, covariant-action and weighted-adjoint identities with rational continuous bounds and written arguments; not proof-assistant formalization"}


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The finite local metric-response report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.67.FINITE_LOCAL_METRIC_RESPONSE replay passed; full metric state response and original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
