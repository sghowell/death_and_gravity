"""Read-only S5.2 replay: reduced interactions, not a scattering/cutoff verdict."""

import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_s5 import compact
from p8_s5.verify import check_prior_sources

from . import lagrangian, momentum
from .vertices import Leg, hamiltonian_kernel, tensor_basis

ROOT = Path(__file__).resolve().parents[2]
S5 = ROOT.parent
REPORT = ROOT/"certificates"/"interactions.json"
PRIOR_SHA = "5a1946d307a18ffd666c66208a689e23698cd18b2c7fa5ffb37100c29f314c59"
TRI = ((10, 20, 0), (0, 10, 30), (-10, -30, -30))
QUAD = ((10, 20, 0), (0, 10, 30), (20, -10, 10), (-30, -20, -40))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prior_checks():
    report = S5/"certificates"/"d-nonlinear.json"
    if sha(report) != PRIOR_SHA:
        raise ValueError("Published S5.1 certificate changed")
    prior = json.loads(report.read_text())
    for relative, expected in prior["source_sha256"].items():
        if sha(S5/relative) != expected:
            raise ValueError(f"Published S5.1 input changed: {relative}")
    if check_prior_sources() != prior["prior_classification_sha256"]:
        raise ValueError("Published classification certificate changed")
    return PRIOR_SHA


def zero(value):
    if sp.cancel(value) != 0:
        raise ValueError(f"Nonzero exact residual: {value}")
    return "0"


def fixture(kinds):
    waves = TRI if len(kinds) == 3 else QUAD
    return tuple(Leg(k, kind, tensor_basis(k)[0] if kind in ("t", "t_dot", "pi") else None)
                 for k, kind in zip(waves, kinds))


def encode(legs):
    return [{"wave": list(map(str, leg.wave)), "kind": leg.kind,
             "polarization": [list(map(str, row)) for row in leg.polarization]
             if leg.polarization else None} for leg in legs]


def quadratic_checks():
    wave = (10, 20, 30)
    opposite = tuple(-k for k in wave)
    q = sum(k*k for k in wave)
    x, J = compact.x, 2*compact.P
    theta, lam = 2*x, 1-2*(1-x*x)**3
    expected = {
        "unitary": {("s", "s"): 2*(lam**2*q*q/J-q),
                    ("s", "p"): -theta*lam*q/J, ("p", "p"): theta**2/(2*J)},
        "gamma": {("s", "s"): 2*theta**2*q*q/J,
                  ("s", "p"): theta*lam*q/J-2*x, ("p", "p"): (q*lam**2-J)/(2*J*q)}}
    out = {}
    for chart, entries in expected.items():
        for kinds, target in entries.items():
            legs = (Leg(wave, kinds[0]), Leg(opposite, kinds[1]))
            out[f"{chart}_{'_'.join(kinds)}"] = zero(hamiltonian_kernel(legs, chart=chart)["kernel"]-target)
    E = tensor_basis(wave)[0]
    norm = sp.trace(sp.Matrix(E)**2)
    for kinds, target in {("t", "t"): q*norm/4, ("pi", "pi"): 4/norm,
                          ("t", "pi"): 0, ("s", "t"): 0}.items():
        legs = tuple(Leg(k, kind, E if kind in ("t", "pi") else None)
                     for k, kind in zip((wave, opposite), kinds))
        out[f"tensor_{'_'.join(kinds)}"] = zero(hamiltonian_kernel(legs)["kernel"]-target)
    for chart, point in (("unitary", sp.Rational(3, 5)), ("gamma", 0)):
        for kind in ("s_dot", "t_dot"):
            legs = tuple(Leg(k, kind, E if kind == "t_dot" else None) for k in (wave, opposite))
            out[f"normalized_{chart}_{kind}"] = zero(lagrangian.kernel(legs, point, chart)["kernel"]-1)
    for kind in ("s", "p"):
        out[f"tadpole_{kind}"] = zero(hamiltonian_kernel((Leg((0, 0, 0), kind),))["kernel"])
    return out


@cache
def build_report():
    prior = prior_checks()
    projectors = momentum.projector_checks()
    exact = {"York_inverse": [[zero(v) for v in row] for row in projectors["inverse"].tolist()],
             "quartic_Legendre": {key: zero(v) for key, v in lagrangian.stationary_checks().items()},
             "quadratic_bridge": quadratic_checks()}
    phase = []
    kinds_list = [tuple(["s"]*(n-t)+["t"]*t) for n in (3, 4) for t in range(n+1)]
    kinds_list += [("s", "p", "p"), ("s", "pi", "pi"), ("p", "t", "pi"),
                   ("s", "p", "t", "pi"), ("p", "p", "pi", "pi")]
    for kinds in kinds_list:
        legs = fixture(kinds)
        out = hamiltonian_kernel(legs)
        phase.append({"legs": encode(legs), "chart": "unitary", "time": "symbolic x in [-1,1]",
                      "kernel_and_bound": compact.bound_coefficient(out["rational_kernel"]),
                      "constraints": out["constraint_checks"]})
    velocity = []
    for chart, point in (("gamma", sp.Integer(0)), ("unitary", sp.Rational(3, 5))):
        for kinds in (("s_dot", "s_dot", "s_dot"), ("s", "t_dot", "t_dot"),
                      ("s", "s", "s_dot", "s_dot"), ("s", "s_dot", "t", "t_dot"),
                      ("t", "t", "t_dot", "t_dot")):
            legs = fixture(kinds)
            out = lagrangian.kernel(legs, point, chart)
            corrections = [{"partition": [list(pair) for pair in item["partition"]],
                            **{key: str(value) for key, value in item.items() if key != "partition"}}
                           for item in out["Legendre_corrections"]]
            velocity.append({"legs": encode(legs), "chart": chart, "x": str(point),
                             "rational_kernel": str(out["rational_kernel"]),
                             "external_normalization": str(out["normalization"]),
                             "kernel": str(out["kernel"]), "rational_minus_H": str(out["rational_minus_H"]),
                             "Legendre_corrections": corrections})
    sources = sorted(ROOT.glob("src/p8_physical/*.py"))+sorted(ROOT.glob("tests/*.py"))
    sources += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    return {
        "schema": 1, "claim": "P8-S5.2.D", "date": "2026-09-04",
        "status": "PHYSICAL_KERNEL_REDUCTION_VERIFIED; SCATTERING_AND_CONTROL_OPEN",
        "scope": "D-only M0, cubic/quartic labelled Fourier kernels, arbitrary nonexceptional rational spatial momenta and rational TT polarizations",
        "prior_s5_1_sha256": prior,
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in sources},
        "exact_residuals": exact, "York_determinant": str(projectors["determinant"]),
        "ordering": {"momentum_constraint_orders": [1, 2, 3], "interaction_orders": [3, 4],
                     "labelled_leg_convention": "coefficient of product epsilon_i; no 1/n! inserted",
                     "Legendre_correction": "three unordered 2+2 partitions; includes scalar and both tensor channels"},
        "phase_kernel_examples_with_all_time_bounds": phase,
        "canonically_normalized_velocity_examples": velocity,
        "domains": {
            "momentum": "total sum zero; every nonempty proper subset sum nonzero; homogeneous/backreaction modes excluded",
            "phase": "x in [-1,1], limiting endpoints included; no inverse Theta or Lambda",
            "unitary_velocity": "x^2>=1/65 is the covering exterior; K=J/Theta^2",
            "gamma_velocity": "x^2<=1/17 and each external/internal q>=68 is sufficient; K=q*J/(q*Lambda^2-J)",
            "time_derivative": "fixed local ell; D_w=(1-x^2)d_x-2*x*q*d_q-2*w*x, w=0 unitary, w=1 gamma",
            "normalization": "Q_c=sqrt(2*K)*Q and Qdot=(Qcdot-DlogZ*Qc)/Z; tensors Z=sqrt(E:E)/2",
        },
        "written_lemmas": ["notes/reduction.md", "notes/audit.md"],
        "not_established": ["zero-mode or forward/soft limits", "cubic-exchange scattering diagrams",
                            "on-shell 2-to-2 amplitudes", "uniform momentum/frequency cutoff hierarchy",
                            "all-orders or loop control", "nonlinear PDE stability", "M1", "UV positivity or completion"],
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("S5.2 certificate differs from exact replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S5.2: physical scalar/tensor quartic kernel replay passed; scattering/control remain OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
