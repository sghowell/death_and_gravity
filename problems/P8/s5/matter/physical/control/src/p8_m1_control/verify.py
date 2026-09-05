"""Read-only replay of coupled M1 normalization and finite-window free control."""

import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_m1_physical import verify as prior

from . import audit, bounds, independent, oscillator
from . import model as m

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"cd-free-control.json"
PRIOR_SHA = "b961a0e07825ba356bec71123d27f47b3c7362d60fc206ee8554f770dc73934f"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prior_checks():
    if sha(prior.REPORT) != PRIOR_SHA:
        raise ValueError("Pinned S5.6.CD certificate changed")
    certificate = json.loads(prior.REPORT.read_text())
    for relative, expected in certificate["source_sha256"].items():
        if sha(prior.ROOT/relative) != expected:
            raise ValueError(f"Pinned S5.6.CD source changed: {relative}")
    prior.prior_checks()
    return PRIOR_SHA


def coefficient_record(value, regular=True):
    if regular:
        m.no_high_frequency_pole(value)
    even, odd = m.parity(value)
    parts = {}
    for name, coefficient in (("even", even), ("odd_times_l", odd)):
        payload = str(coefficient).encode()
        parts[name] = {
            "rational_expression_sha256": hashlib.sha256(payload).hexdigest(),
            "numerator_terms": len(coefficient.numer.terms()),
            "denominator_terms": len(coefficient.denom.terms()),
            "z_zero_denominator_not_identically_zero": coefficient.denom.as_expr().subs(m.z, 0) != 0,
        }
    return {"exact_arithmetic": "QQ(x,z)[l]/(l^2-(1-x^2)^11/100)", **parts}


@cache
def build_report():
    previous = prior_checks()
    generic = {key: audit.zero(value) for key, value in independent.generic_checks().items()}
    negatives = {}
    for key, value in independent.generic_negative_controls().items():
        entries = tuple(value) if isinstance(value, sp.MatrixBase) else (value,)
        if not any(sp.expand(entry) != 0 for entry in entries):
            raise ValueError(f"Omission control vanished: {key}")
        negatives[key] = list(map(str, entries))
    charts = {}
    for chart in ("unitary", "gamma"):
        data = oscillator.derive(chart)
        coefficients = {key: coefficient_record(data[key]) for key in
            ("mass11", "mass22", "mass12_factor", "connection_factor",
             "covariant_mass11", "covariant_mass22", "covariant_mass12_factor")}
        boundary = {key: coefficient_record(data[key], regular=False) for key in
            ("momentum_boundary11", "momentum_boundary22", "momentum_boundary12_factor")}
        residuals = [coefficient_record(value) for value in data["unwhitened_mass_residual"]]
        charts[chart] = {"exact_residuals": audit.compact_checks(chart),
                         "regular_mass_and_connection_coefficients": coefficients,
                         "principal_subtracted_unwhitened_mass": residuals,
                         "original_phase_momentum_boundary_coefficients": boundary}
    tails = []
    for point in (-1, 1):
        actual = oscillator.physical_matrices("unitary", point, 1000)
        tails.append({"x": point, "mass_plus_30I": audit.matrix_zeros(actual["mass"]+30*sp.eye(2)),
                      "connection": audit.matrix_zeros(actual["connection"]),
                      "tensor_mass_plus_30": audit.zero(oscillator.tensor()["mass"].subs(m.x, point)+30)})
    sources = sorted(ROOT.glob("src/p8_m1_control/*.py"))+sorted(ROOT.glob("tests/*.py"))
    sources += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    return {
        "schema": 1, "claim": "P8-S5.7.CD", "date": "2026-09-05",
        "status": "COUPLED_CANONICAL_NORMALIZATION_AND_LOCAL_FREE_ENERGY_CONTROL; M1_INTERACTION_CONTROL_OPEN",
        "prior_S5_6_CD_sha256": previous,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in sources},
        "formulation": "FORMULATION.md", "written_proof": "notes/normalization.md",
        "generic_coupled_canonical_and_energy_residuals": generic,
        "generic_omission_controls": negatives,
        "compact_coefficient_records": charts,
        "independent_original_cosmic_time_fixtures": [audit.point_check(*point) for point in audit.POINTS],
        "symbolic_finite_q_bounce": audit.bounce_check(), "compact_tail_checks": tails,
        "covering_arithmetic_margins": bounds.covering_checks(),
        "uniform_exact_majorants_and_free_energy_contract": bounds.build_bounds(),
        "domain": {
            "background": "the fixed CD/M1 witness, same covariant tube and physical matter metric",
            "compact_time": "x=u/sqrt(1+u^2), u=t/tau; x=+-1 denotes compact tail limits only",
            "local_length": "ell=tau*sqrt(1+u^2)",
            "momentum": "fixed nonzero comoving momentum; q=ell^2*kcom^2/a^2 and z=1/q",
            "positive_kinetic_cover": {"unitary": "abs(x)>=1/9", "gamma": "abs(x)<=1/4", "both": "q>=1000"},
            "local_window": "abs(t-t0)<=ell0/100; q0>=2*10^20; gamma if abs(x0)<=9/50, unitary otherwise",
            "solutions": "any exact coupled scalar or individual tensor free solution; nonzero solutions for energy ratios",
            "not_a_global_fixed_comoving_high_frequency_band": True,
        },
        "canonical_dictionary": {
            "Hamiltonian_measure": "a^3; H2=p^T*A*p/2+p^T*B*Q+Q^T*C*Q/2",
            "normalization": "T^T*T=A^-1; Y=a^(3/2)*T*Q; pi=a^(3/2)*T^-T*p",
            "time_generator": "F=3*H*I/2+Tdot*T^-1; Bhat=T*B*T^-1+F",
            "boundary_and_connection": "S=sym(Bhat); Omega=-skew(Bhat); P=pi+S*Y",
            "potential": "W=T^-T*C*T^-1-S^2+[S,Omega]-Sdot=(q*I+M)/ell^2",
            "Lagrangian": "(Ydot+Omega*Y)^T*(Ydot+Omega*Y)/2-Y^T*W*Y/2",
            "original_phase_inverse": "Q=a^-3/2*T^-1*Y; p=a^-3/2*T^T*(P-S*Y)",
            "free_energy": "E=(P^dagger*P+Y^dagger*W*Y)/2; Edot=Y^dagger*(Wdot+[Omega,W])*Y/2",
            "mass_dimension_derivative": "D_n=(1-x^2)*partial_x+6*x*z*partial_z-11*x*l*partial_l-n*x",
            "finite_q_mixing_retained": True, "instantaneous_eigenvectors_used": False,
            "frequency_dependent_canonical_boundary_is_not_an_interaction_bound": True,
        },
        "verification_boundary": [
            "the pinned physical Hamiltonian and all-time J bounds are imported certified inputs; their source chain is checked",
            "exact compact identities and coefficient denominator/L1 bounds establish the all-time chart estimates",
            "independent original-time fixtures test the implementation but do not replace the written all-time derivation",
            "matrix energy and elementary interval/Gronwall lemmas are written proofs, not Lean formalizations",
        ],
        "not_established": [
            "normalized cubic/quartic M1 kernel majorants or finite-time interacting tree control",
            "an interacting cutoff, an optimized high-frequency threshold or a usable hierarchy at a specified M*tau",
            "a preferred vacuum, a WKB expansion or fixed-comoving high-frequency control over the entire bounce history",
            "all-wavelength stability, infrared/zero-mode or global chart-stitching energy estimates",
            "on-shell scattering amplitudes, soft/forward/inclusive limits, loops or radiative protection",
            "nonlinear PDE/BKL/backreaction stability, UV matching/completion or full P8 completion",
        ],
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("CD/M1 free-control certificate differs from exact replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S5.7: CD/M1 coupled normalization and local free-energy replay passed; interaction control OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
