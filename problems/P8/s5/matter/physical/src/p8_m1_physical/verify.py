"""Read-only CD/M1 physical reduction replay; no normalization/control verdict."""

import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_m1 import series
from p8_m1.verify import prior_checks as invariant_prior_checks
from p8_physical.momentum import projector_checks
from p8_physical.verify import prior_checks as geometry_prior_checks

from . import background, lagrangian, momentum, quadratic, regressions, vertices

ROOT = Path(__file__).resolve().parents[2]
P8 = ROOT.parents[2]
REPORT = ROOT/"certificates"/"cd-interactions.json"
PRIOR = {
    "s5/matter/certificates/cd-matter.json": "6b0dcb44c0849050912f22546fc7c2ef4c55f148ebca32184b09b4abecf7ca7c",
    "s5/physical/certificates/interactions.json": "fec2b468d9144817026e923eeed1d7b573fdae068213b42a93cb2e4e77eb8c20",
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prior_checks():
    invariant_prior_checks()
    geometry_prior_checks()
    for relative, expected in PRIOR.items():
        path = P8/relative
        if sha(path) != expected:
            raise ValueError(f"Pinned input certificate changed: {relative}")
        prior = json.loads(path.read_text())
        base = path.parents[1]
        for source, digest in prior["source_sha256"].items():
            if sha(base/source) != digest:
                raise ValueError(f"Pinned input source changed: {base/source}")
    return dict(PRIOR)


def zeros(residuals):
    bad = {key: value for key, value in residuals.items() if sp.cancel(value) != 0}
    if bad:
        raise ValueError(f"Nonzero M1 physical residual: {bad}")
    return dict.fromkeys(residuals, "0")


def phase_examples():
    out = []
    kinds_list = (("s", "m", "m"), ("p", "m", "P"), ("m", "m", "t"),
                  ("s", "P", "P"), ("m", "P", "pi"), ("s", "m", "P"),
                  ("s", "s", "m", "m"), ("m", "m", "m", "m"),
                  ("s", "m", "t", "pi"), ("m", "m", "pi", "pi"),
                  ("p", "P", "p", "P"), ("m", "P", "t", "pi"))
    for chart, point in (("unitary", sp.Rational(1, 3)), ("gamma", sp.Integer(0))):
        for kinds in kinds_list:
            legs = regressions.fixture(kinds)
            result = vertices.hamiltonian_kernel(legs, point, chart)
            if sp.im(result["kernel"]) != 0 or not all(result["constraint_checks"].values()):
                raise ValueError("A physical fixture failed reality or spatial constraints")
            out.append({"legs": regressions.encode(legs), "chart": chart, "r": str(point),
                        "kernel": str(result["kernel"]), "constraints": result["constraint_checks"]})
    for point in (-1, 1):
        legs = regressions.fixture(("s", "m", "m"))
        result = vertices.hamiltonian_kernel(legs, point)
        out.append({"legs": regressions.encode(legs), "chart": "unitary", "r": str(point),
                    "compact_tail_limit": True, "kernel": str(result["kernel"]),
                    "constraints": result["constraint_checks"]})
    return out


def velocity_examples():
    out = []
    kinds_list = (("s_dot", "m_dot", "m_dot"), ("m", "t_dot", "t_dot"),
                  ("m_dot", "m_dot", "m_dot", "m_dot"), ("s", "m_dot", "s", "m_dot"),
                  ("s_dot", "m_dot", "t", "t_dot"))
    for chart, point in (("gamma", sp.Integer(0)), ("unitary", sp.Rational(1, 3))):
        for kinds in kinds_list:
            legs = regressions.fixture(kinds)
            result = lagrangian.kernel(legs, point, chart)
            corrections = [{"partition": [list(pair) for pair in item["partition"]],
                            **{key: str(value) for key, value in item.items()
                               if key not in ("partition", "tensor_polarizations")},
                            "tensor_polarizations": list(map(str, item["tensor_polarizations"]))}
                           for item in result["Legendre_corrections"]]
            out.append({"legs": regressions.encode(legs), "chart": chart, "r": str(point),
                        "kernel": str(result["kernel"]), "minus_H": str(result["minus_H"]),
                        "Legendre_corrections": corrections, "normalized": False})
    return out


@cache
def build_report():
    prior = prior_checks()
    inverse = projector_checks()["inverse"]
    exact = {
        "compact_background": zeros(background.compact_checks()),
        "York_inverse": zeros({f"{i}{j}": inverse[i, j] for i in range(3) for j in range(3)}),
        "York_symplectic_projection": zeros(momentum.symplectic_projector_checks()),
        "canonical_boundaries_and_local_q_drift": zeros(regressions.boundary_checks()),
        "stationary_lapse_bridge": zeros(vertices.stationary_bridge_checks()),
        "full_two_scalar_velocity_response": zeros(quadratic.identities()),
        "matrix_quartic_Legendre_completion": zeros(lagrangian.stationary_checks()),
        "symbolic_time_full_quadratic_bridge": zeros(regressions.quadratic_checks()),
    }
    phase, velocity = phase_examples(), velocity_examples()
    sources = sorted(ROOT.glob("src/p8_m1_physical/*.py"))+sorted(ROOT.glob("tests/*.py"))
    sources += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    return {
        "schema": 1, "claim": "P8-S5.6.CD", "date": "2026-09-05",
        "status": "M1_PHYSICAL_PHASE_AND_UNNORMALIZED_VELOCITY_REDUCTION; INTERACTION_CONTROL_OPEN",
        "prior_sha256": prior,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in sources},
        "exact_residuals": exact, "negative_controls": regressions.negative_controls(),
        "phase_kernel_examples": phase, "unnormalized_velocity_examples": velocity,
        "domain": {
            "background": "fixed CD/M1; same physical matter metric and covariant tube as P8-2.CD",
            "momentum": "arbitrary rational spatial momenta, total zero, every nonempty proper-subset sum nonzero",
            "tensors": "arbitrary nonzero real rational TT polarizations, both internal polarizations contracted",
            "phase_time": "r=u/(1+sqrt(1+u^2)) in [-1,1]; endpoints denote compact tail limits",
            "local_units": "ell=tau*sqrt(1+u^2), local a=1; q=ell^2*k_physical^2",
            "unitary_velocity": "Theta!=0; actual full scalar K positive definite",
            "gamma_velocity": "q*Lambda^2>J0=J+w^2/2 at every external and quartic internal momentum",
            "symbolic_velocity_calls": "rational expression only; positivity hypotheses not certified over its entire time domain",
            "zero_modes": "homogeneous tadpoles checked separately; homogeneous/backreaction dynamics not solved",
        },
        "construction": {
            "metric": "g_ij=(1+2*zeta)*delta_ij+gammaTT_ij",
            "mixed_boundary_inverse": "p_metric=p+3*l*s; pi_chi=l+P+3*l*zeta at local a=1",
            "mixed_generator": "F=3*a^3*l*zeta*s, partial_t F=0 because (a^3*l)'=0",
            "constraint": "partial_j Pi_bar^ij+Gamma^i_jk Pi_bar^jk-(pi_chi/2)*g^ij*partial_j s=0",
            "York_recursion_orders": [1, 2, 3],
            "physical_hamiltonian_orders": [3, 4],
            "stationary_lapse": "same pinned invariant quartic formula; only first two perturbative lapse-solution orders n1,n2 are needed for H4; coefficient derivatives through A4,L3,Q2 are retained",
            "canonical_generators_retained": ["metric expansion boundary", "-l*pi_chi", "-H*b*P_b in gamma chart"],
            "labelled_convention": "coefficient of product epsilon_i, no 1/n!; repeated fields have distinct leg labels",
            "Legendre_contact": "L4=-H4(P0)+H3_P^T*A^-1*H3_P/2; three unordered 2+2 partitions",
            "scalar_contact": "full two-by-two inverse Hessian, including off-diagonal terms",
            "tensor_contact": "sum over both TT basis tensors of (E:E)*left_pi*right_pi/4",
            "not_a_scattering_amplitude": True,
        },
        "finite_momentum_quadratic_distinction": {
            "unitary_H2": str(quadratic.symbolic("unitary")["density"]),
            "gamma_H2": str(quadratic.symbolic("gamma")["density"]),
            "gamma_K": [[str(value) for value in row] for row in quadratic.symbolic("gamma")["kinetic"].tolist()],
            "bounce_q8_K": [[str(value) for value in row] for row in quadratic.response(0, "gamma", 8)["kinetic"].tolist()],
            "bounce_velocity_pole": "q=6; phase Hamiltonian is regular there",
            "principal_limit": "q->infinity recovers [[J0/Lambda^2,-w/(2*Lambda)],[-w/(2*Lambda),1/2]]",
            "finite_q_matter_diagonal_is_not_generally_one_half": True,
            "no_application_of_principal_cone_lemma_to_finite_q_K": True,
        },
        "structural_all_time_phase_regularity": {
            "basis": "finite recursion in compact-polynomial jets, fixed nonzero transfer inverses and 1/J",
            "compact_variables": "x=2r/(1+r^2), y=(1-r^2)/(1+r^2); x^2+y^2=1",
            "J_bounds_from_pinned_S5_5": ["1/10", "8"],
            "compact_lapse_jet_absolute_bounds": {
                key: [str(series.polynomial_bound(value)) for value in row]
                for key, row in series.compact_jets().items()},
            "consequence": "each fixed nonexceptional-momentum phase coefficient extends continuously and is bounded on compact time",
            "no_inverse_Theta_or_Lambda_in_phase_construction": True,
            "not_a_uniform_momentum_or_velocity_bound": True,
            "no_expanded_symbolic_quartic_majorant_is_claimed": True,
        },
        "written_proofs": ["notes/reduction.md", "notes/audit.md"],
        "not_established": ["time-dependent coupled configuration/mode normalization",
                            "finite-k positive oscillator chart or free-mode propagation bound",
                            "nonexceptional interacting M1 tree-control window or cutoff hierarchy",
                            "cubic-exchange scattering amplitudes or inclusive/forward/soft limits",
                            "all-orders/loop control or radiative protection of the exceptional relation",
                            "global nonlinear Dirac/PDE, BKL, or backreaction stability",
                            "UV admissibility or completion", "full P8 completion"],
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("CD/M1 physical certificate differs from exact replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S5.6: CD/M1 physical phase and unnormalized velocity replay passed; interaction control OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
