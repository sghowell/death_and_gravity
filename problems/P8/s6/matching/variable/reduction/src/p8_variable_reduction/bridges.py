"""Continuous bridges to the literal frozen parent and CD operator input."""
import hashlib
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8 import rational_candidates
from p8_variable_beta import background

from . import dictionary, stationary, tensor

P8 = next(path for path in Path(__file__).resolve().parents if path.name == "P8")
CD = P8/"certificates"/"witness-CD_matter.json"
CD_SHA = "caf8c8e688a7565b9d00f921c099a28da00f97522ed26ad182a8227eb80cd4dd"


def target_input():
    if hashlib.sha256(CD.read_bytes()).hexdigest() != CD_SHA:
        raise ValueError("The frozen CD covariant witness changed")
    raw = json.loads(CD.read_text())
    for name, digest in raw["source_sha256"].items():
        if hashlib.sha256((P8/"src"/"p8"/name).read_bytes()).hexdigest() != digest:
            raise ValueError(f"The frozen CD source changed: {name}")
    u, x = sp.symbols("u X", real=True)
    local = {"t": u, "X": x}
    functions = {key: sp.sympify(value, locals=local)
                 for key, value in raw["covariant_functions"].items()}
    return {"u": u, "X": x, "functions": functions,
            "source_hash_count": len(raw["source_sha256"]), "sha256": CD_SHA}


@cache
def checks():
    parent, reduced = background.derive(), stationary.profile()
    u = parent["u"]
    mapping = {stationary.V: u**2, stationary.C: parent["c"]}
    pairs = {"d": "d", "y": "y", "h_u": "h_u", "b1": "b1",
             "b4": "b4", "kbar": "kbar"}
    residuals = {f"literal_parent_{left}": sp.cancel(reduced[left].subs(mapping)-parent[right])
                 for left, right in pairs.items()}
    residuals.update({
        "literal_parent_h_squared": sp.cancel(reduced["h_squared"].subs(mapping)-parent["h"]**2),
        "actual_stationary_root": sp.cancel(reduced["r_cubed"]*reduced["b4"]+reduced["b1"]),
        "actual_kappa_coefficient": sp.cancel(reduced["kappa_bar"]-reduced["r_cubed"]/(4*reduced["b1"])),
        "clock_even_center": sp.diff(parent["kbar"], u).subs(u, 0),
        "kappa_even_center": sp.diff(reduced["kappa_bar"].subs(mapping), u).subs(u, 0),
        "stationary_root_even_center": sp.diff(reduced["r_cubed"].subs(mapping), u).subs(u, 0),
    })
    target = target_input()
    v, x, functions = target["u"], target["X"], target["functions"]
    spec = rational_candidates.specification("CD_matter")
    replacements = {rational_candidates.j.t: v, rational_candidates.X: x}
    for key in ("F2", "A1", "A3", "K"):
        residuals[f"frozen_CD_{key}"] = sp.cancel(functions[key]-spec[key].subs(replacements))
    f2, a1, a3 = (functions[key].subs(v, 0) for key in ("F2", "A1", "A3"))
    xi = a1-2*sp.diff(f2, x)
    gt = -2*f2+2*x*a1
    residuals.update({"whole_tube_CD_center_F2": f2+x/2,
                      "whole_tube_CD_center_F2_X": sp.diff(f2, x)+sp.Rational(1, 2),
                      "whole_tube_CD_center_A1": a1,
                      "whole_tube_CD_center_A3": a3-1/x,
                      "whole_tube_CD_center_Xi": xi-1,
                      "whole_tube_CD_center_I": sp.cancel(x*xi/gt-1)})
    co = dictionary.coefficients()
    cz = co["symbols"]
    m2 = sp.Symbol("Mstar_squared", positive=True)
    at_center = {cz["r"]: sp.Integer(2), cz["q"]: sp.S.Zero,
                 sp.diff(cz["kappa"], cz["phi"]): sp.S.Zero,
                 cz["M_squared"]: m2/5}
    residuals["formal_center_tensor_value"] = co["F2"].subs(at_center)+m2/2
    residuals["formal_center_curvature_X_derivative"] = co["F2_X"].subs(at_center)
    residuals["formal_center_Horndeski_A1"] = co["A1"].subs(at_center)
    residuals["full_tensor_kappa_normalization"] = (
        tensor.calibration()["kappa_center"]
        -stationary.calibration()["center_kappa"])
    # The symbols M,tau,c deliberately have identical names and assumptions.
    residuals = {key: sp.factor(value) for key, value in residuals.items()}
    if any(value != 0 for value in residuals.values()):
        raise ValueError("An actual-parent or open-tube CD bridge failed")
    return residuals


def domain_margins():
    dmax = sp.Rational(101, 100)
    theta_upper = sp.Rational(7, 220)
    return {"y_above_19over10": 2/dmax**4-sp.Rational(19, 10),
            "one_minus_theta_upper": 1-theta_upper,
            "root_ratio_below_9over4": sp.Rational(9, 4)**3-8/(1-theta_upper),
            "center_clock_upper": sp.Rational(2399, 100),
            "center_clock_lower": sp.Rational(799, 100)}
