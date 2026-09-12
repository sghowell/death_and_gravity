"""Read-only original full flat tree-plus-Proca tensor response."""

import argparse
import json
from functools import cache
from pathlib import Path

from p8_affine import verify as certificate
from p8_vacuum_affine_adiabatic_action import verify as subtraction_input
from p8_vacuum_affine_band_coupled_inverse import verify as parent
from p8_vacuum_affine_flat_tensor_cut import verify as cut_input
from p8_vacuum_affine_quantum_retuning import verify as vacuum_input

from . import audit

ROOT = Path(__file__).resolve().parents[2]
REPORT = (
    ROOT / "certificates" / "polynomial-vacuum-affine-full-flat-tensor-response.json"
)
PARENT_SHA = "076648eee2112f5b86302a072e11265350e155a92877f5b68b5007ab0cf09432"
SUBTRACTION_SHA = "98888c7a7fdcfe4fee142833a4220c62192b2834f736625192c8534201dca278"
VACUUM_SHA = "f754f6d62795d5a1c0ade4c576395f1587fc35f9c9895bc8b7b15b2b2dec2e16"
CUT_SHA = "e7a37017276d0b5474d81a295ab3f2574f69b883ffee34b032a0dadf9566ceec"
sha, serialize = certificate.sha, certificate.serialize


def source_files():
    return (
        sorted(ROOT.glob("src/p8_vacuum_affine_full_flat_tensor_response/*.py"))
        + sorted(ROOT.glob("tests/*.py"))
        + sorted(ROOT.glob("*.md"))
        + sorted(ROOT.glob("notes/*.md"))
    )


def payload(data):
    return {key: value for key, value in data.items() if key not in ("checks", "gates")}


@cache
def prior_checks():
    for prior, expected in (
        (parent, PARENT_SHA),
        (subtraction_input, SUBTRACTION_SHA),
        (vacuum_input, VACUUM_SHA),
        (cut_input, CUT_SHA),
    ):
        if sha(prior.REPORT) != expected:
            raise ValueError(
                "A source-pinned original vacuum, subtraction, full cut or coupled frontier changed"
            )
        prior.validate_report(
            json.loads(prior.REPORT.read_text()), prior.build_report()
        )
    return {
        "S6_229_actual_smooth_band_coupled_inverse_and_unchanged_primitive_frontier_fully_rebuilt": PARENT_SHA,
        "S6_193_original_full_noncommuting_adiabatic_current_and_finite_matching_fully_rebuilt": SUBTRACTION_SHA,
        "S6_182_original_retained_Minkowski_vacuum_retuning_and_unchanged_low_jets_fully_rebuilt": VACUUM_SHA,
        "S6_199_complete_original_massive_Proca_cut_normalization_fully_rebuilt": CUT_SHA,
        "S6_223_original_shear_factor_S6_226_exterior_bounds_S6_189_finite_action_and_S6_177_Einstein_normalization_transitively_rebuilt": True,
        "original_state_prescription_and_all_frozen_scientific_bytes_unchanged": True,
    }


@cache
def build_report():
    prior = prior_checks()
    rows, gates = audit.residuals(), audit.gates()
    if not all(value is True for value in gates.values()):
        raise ValueError(
            "An original full flat tensor normalization, pole, scale or causal graph gate failed"
        )
    return {
        "schema": 1,
        "claim": "P8-S6.230.ACTUAL_FULL_FLAT_VACUUM_TREE_PROCA_TT_COMPLEX_POLE_PAIR_LOW_ENERGY_BOUND_AND_ORIGINAL_ALL_MOMENTUM_CAUSAL_GRAPH_INVERSE",
        "date": "2026-09-12",
        "status": "ACTUAL_RETAINED_FLAT_MEAN_EQUATION_GROWING_POLES_AND_ALL_MOMENTUM_CAUSAL_INVERSE; NOT_CURVED_INSTABILITY_PHYSICAL_UV_NO_GO_CUTOFF_NONLINEAR_V_G_B_OR_ORIGINAL_P8",
        "prior_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files()},
        "formulation": "FORMULATION.md",
        "written_proofs": [
            "notes/bridge.md",
            "notes/count.md",
            "notes/bounds.md",
            "notes/dispersion.md",
            "notes/response.md",
            "notes/scale-scope.md",
            "notes/validation.md",
        ],
        "exact_residuals": certificate.certify_residuals(rows),
        "named_exact_check_count": len(rows),
        "checked_scalar_entries": audit.scalar_entry_count(),
        "proof_checks": gates,
        "actual_flat_vacuum_symbol_and_original_polynomial_bridge": serialize(
            payload(audit.normalization.data())
        ),
        "full_first_sheet_poles_scale_dispersion_and_causal_graph": serialize(
            {
                "analytic": payload(audit.analytic.data()),
                "response": payload(audit.response.data()),
            }
        ),
        "observable_and_scope": audit.observable(),
        "primitive_and_matching_frontier": serialize(
            {"primitive": audit.frontier(), "matching": audit.matching()}
        ),
        "controls": serialize(audit.controls()),
        "verdict": "At the separately specified original Minkowski vacuum, the complete canonical three-polarization covariance equation shows that the original derivative-marker subtraction equals the first three frequency Taylor terms, including the original second vertex. S193 matching therefore fixes the actual lower polynomial, rather than inferring it from the cut. With the original Einstein tree, Rvac1, unchanged low vacuum jets and finite m2R and Weyl terms, the physical force symbol is O_phys=p(C+pA2)/(16pi^2kappa), C=16pi^2kappa+5m^2/6 and p=lambda^2+P^2. For every C>2m^2/15 the original complete first sheet has exactly one simple nonreal zero in each p half-plane, both with negative real part, and no real zeros. The full real-boundary and infinity-arc winding is proved analytically. At actual m1000,kappa10^800, 10^796<|z|<kappa and Re sqrt(z)>10^393, whereas the physical low-disk propagator ratio differs from Einstein by less than197m^2/(30240kappa)<10^-796. The original full reciprocal retains the massless pole, complex pair and positive massive cut. Its two finite moment cancellations and exact residue bound250/(249C) give the original causal all-momentum force inverse norm below4T^2 exp(sqrt(kappa)T), without spatial derivative loss. The bounded paired forward distribution and Fourier cutoff limit give both identities on the explicitly specified flat TT causal graph. Smooth compact-time prepared real L2 TT forcing can excite growing modes; this is a retained-equation mathematical response obstruction, not an unprepared reset. The pole is between one hundredth and one times the sqrt(kappa) normalization, where the conditional approximation has not been certified. No pole prescription, state, regulator, finite counterterm or cutoff is changed. Physical UV exclusion, curved-bounce instability, unrestricted curved scalar/clock/matter inversion, nonlinear common-parent control and original V/G/B/P8 remain OPEN.",
        "not_established": [
            "A physical UV no-go, validity of the retained approximation at the complex poles, or a physical cutoff",
            "Absence of growing modes: the exact retained flat equation instead has prepared forced growth",
            "Curved CD-bounce instability or an unrestricted completed curved scalar/clock/matter S222 graph inverse",
            "Full scalar, graviton and mixed loops, nonlinear or finite-coupling common-parent control, original V/G/B or P8 closure",
        ],
        "verification_boundary": "Exact full covariance, subtraction, original finite-force signs, cut/moment algebra and coarse rational inequalities are replayed with independent diagnostics. The original adiabatic-to-vacuum bridge, full-sheet argument principle, global complex bounds, Cauchy representation, causal uniform inverse, graph limit and growing-packet construction are written continuum proofs, not FORMALIZED. Numerical quadrature, root and winding diagnostics are not continuum certificates. Native/direct/ordinary/CLI use original SymPy; only full regression uses the audited exact-GCD adapter. Frozen predecessors remain unchanged.",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError(
            "The original full flat tensor response report differs from read-only replay"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print(
            "P8 S6.230 original full flat TT response replay passed; retained growing poles found, physical UV and original P8 remain OPEN"
        )
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
