"""Cross-route checks, kept separate from the compact oscillator derivation."""

from functools import cache

import sympy as sp
from p8_m1_physical import quadratic

from . import independent, oscillator
from . import model as m

POINTS = (("unitary", sp.Rational(-3, 4)), ("unitary", sp.Rational(3, 4)),
          ("gamma", sp.Rational(-9, 40)), ("gamma", sp.Integer(0)),
          ("gamma", sp.Rational(9, 40)))


def zero(value):
    if sp.simplify(value) != 0:
        raise ValueError(f"Nonzero independent oscillator residual: {value}")
    return "0"


def matrix_zeros(value):
    return [[zero(entry) for entry in row] for row in value.tolist()]


def encoded_matrix(value):
    return [[str(entry) for entry in row] for row in value.tolist()]


@cache
def point_check(chart, u_value):
    # Hold the actual comoving momentum fixed during all derivatives.  These
    # fixtures check identities, not positivity over a time interval.
    raw = independent.cosmic_time_normalization(chart, u_value, sp.Integer(10000))
    local = oscillator.physical_matrices(chart, raw["x"], raw["local_q"])
    ell = raw["ell"]
    data = oscillator.derive(chart)
    mapping = {m.x: raw["x"], m.z: 1/raw["local_q"],
               m.l: ell*raw["jets"][0][quadratic.l]}
    at = lambda value: sp.simplify(value.subs(mapping, simultaneous=True))
    root = sp.sqrt(at(data["factors"]["ratio"]))
    off = root*at(data["momentum_boundary12_factor"])
    boundary = sp.Matrix([[at(data["momentum_boundary11"]), off],
                          [off, at(data["momentum_boundary22"])]])
    residuals = {
        "full_finite_q_potential": matrix_zeros(ell**2*raw["potential"]-local["potential"]),
        "antisymmetric_connection": matrix_zeros(ell*raw["connection"]-local["connection"]),
        "canonical_momentum_boundary": matrix_zeros(ell*raw["momentum_boundary"]-boundary),
        "physical_kinetic_whitening": matrix_zeros(raw["kinetic_residual"]),
    }
    return {"chart": chart, "u": str(u_value), "x": str(raw["x"]),
            "ell_at_tau_one": str(ell), "kcom_squared": "10000",
            "q_local": str(raw["local_q"]), "exact_residuals": residuals,
            "local_potential": encoded_matrix(local["potential"]),
            "local_connection": encoded_matrix(local["connection"])}


@cache
def bounce_check():
    q = sp.Symbol("q_bounce", positive=True)
    raw, target = independent.bounce_from_cosmic_time(q), independent.bounce_expected(q)
    data = oscillator.derive("gamma")
    at = lambda value: sp.cancel(value.subs({m.x: 0, m.z: 1/q, m.l: sp.Rational(1, 10)}, simultaneous=True))
    actual = {key: at(data["factors"][key]) for key in ("d1", "d2", "chi", "ratio")}
    actual.update({"potential11": q+at(data["mass11"]),
                   "potential22": q+at(data["mass22"]),
                   "potential12_factor": at(data["mass12_factor"]),
                   "connection_factor": at(data["connection_factor"]),
                   "connection_dot_factor": at(m.derivative(data["connection_factor"], 1)
                       +data["factors"]["log_root_ratio_derivative"]*data["connection_factor"])})
    determinant = actual["potential11"]*actual["potential22"]-actual["ratio"]*actual["potential12_factor"]**2
    residuals = {f"compact_vs_cosmic_{key}": zero(value-raw[key]) for key, value in actual.items()}
    residuals.update({f"cosmic_vs_known_answer_{key}": zero(raw[key]-target[key]) for key in actual})
    residuals["determinant_known_answer"] = zero(determinant-target["determinant"])
    low = independent.bounce_expected(sp.Integer(8))
    if not (low["d1"] > 0 and low["d2"] > 0 and low["potential11"] < 0 and low["determinant"] < 0):
        raise ValueError("The low-q positive-kinetic/indefinite-potential control failed")
    return {"positive_velocity_domain": "q>6", "exact_residuals": residuals,
            "known_answer": {key: str(value) for key, value in target.items()},
            "negative_control_q8": {key: str(low[key]) for key in
                ("d1", "d2", "potential11", "determinant")},
            "negative_control_meaning": "positive kinetic does not imply positive W; NOT a dynamical-instability theorem"}


def compact_checks(chart):
    data, f = m.coefficients(chart), m.factors(chart)
    target = sp.Matrix([[f["d1"]+f["chi"]**2*f["d2"], f["chi"]*f["d2"]],
                        [f["chi"]*f["d2"], f["d2"]]])
    out = {}
    for name, residual in (("full_Cholesky", data["alpha"]-target),
                           ("full_velocity_cross", data["alpha"]*data["B"]+data["beta"])):
        for i in range(2):
            for j in range(2):
                value = m.canonical(residual[i, j])
                if value != 0:
                    raise ValueError(f"Nonzero compact residual {chart}/{name}/{i}{j}")
                out[f"{name}_{i}{j}"] = "0"
    return out
