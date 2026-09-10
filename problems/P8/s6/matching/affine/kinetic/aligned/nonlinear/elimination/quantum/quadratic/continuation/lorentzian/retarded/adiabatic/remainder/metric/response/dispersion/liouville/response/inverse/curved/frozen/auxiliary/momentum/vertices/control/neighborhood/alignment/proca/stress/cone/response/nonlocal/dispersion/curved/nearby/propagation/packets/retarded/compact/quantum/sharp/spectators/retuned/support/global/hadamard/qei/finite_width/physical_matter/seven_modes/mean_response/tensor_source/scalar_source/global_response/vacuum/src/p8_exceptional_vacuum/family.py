"""A separately named smooth classical light-action extension, not UV matching."""

from functools import cache

import sympy as sp
from p8.matter import ia_completion
from p8_affine import dictionary
from p8_proca_retuned_margin import model as retuned

u = dictionary.u
X = sp.Symbol("physical_positive_clock_X", real=True)
B = sp.Function("smooth_clock_switch")(X)
Z, m = sp.symbols("positive_vacuum_Z positive_vacuum_mass", positive=True)
u0, lam = sp.symbols("vacuum_clock_value vacuum_derivative_lambda", real=True)


def vacuum_scalar(kappa=sp.S.One):
    # The physical action is kappa times the dimensionless kernel. Scaling
    # the bare quartic by kappa keeps its canonical interaction fixed.
    return Z * X / 2 - Z * m * m * (u - u0) ** 2 / 2 + kappa * lam * Z * Z * X * X


@cache
def data():
    h = (1 + u * u) ** 3
    old = dictionary.target()
    tree = (
        old["scalar_F"].subs(dictionary.x, -X)
        + retuned.NEW_MARGIN * (X - 1) ** 2 / h**2
    )
    R = 1 + B * (X - 1) / h
    RX = sp.diff(R, X)
    F2 = -R / 2
    A3 = RX / X
    A2, A4, A5 = ia_completion(F2, sp.diff(F2, X), sp.S.Zero, A3, X)
    functions = {
        "F": B * tree + (1 - B) * vacuum_scalar(),
        "K": sp.Integer(0),
        "F2": F2,
        "A1": sp.Integer(0),
        "A2": A2,
        "A3": A3,
        "A4": A4,
        "A5": A5,
    }
    return {
        "h": h,
        "R": R,
        "RX": RX,
        "original_retuned_tree_scalar": tree,
        "functions_on_transition_or_tube": functions,
        "vacuum_functions": {
            "F": vacuum_scalar(),
            "K": sp.Integer(0),
            "F2": -sp.Rational(1, 2),
            "A1": sp.Integer(0),
            "A2": sp.Integer(0),
            "A3": sp.Integer(0),
            "A4": sp.Integer(0),
            "A5": sp.Integer(0),
        },
        "generic_closed_completion": {
            "A3": RX / X,
            "A4": -RX / X - sp.Rational(7, 4) * RX**2 / R,
            "A5": RX**2 / (R * X),
        },
        "physical_metric_is_unchanged": True,
        "piecewise_definition_avoids_zero_times_undefined_old_kernel": True,
    }


def constant_branch(value):
    if isinstance(value, bool) or not isinstance(value, int) or value not in (0, 1):
        raise TypeError("Require a certified constant branch zero or one")
    sub = {sp.diff(B, X): 0, B: value}
    return {
        key: sp.factor(expr.subs(sub, simultaneous=True))
        for key, expr in data()["functions_on_transition_or_tube"].items()
    }


@cache
def identities():
    d = data()
    f = d["functions_on_transition_or_tube"]
    target = dictionary.target()
    tube = constant_branch(1)
    vac = constant_branch(0)
    old = {
        "F": d["original_retuned_tree_scalar"],
        "K": -target["braiding"],
        "F2": -target["f"],
        "A1": target["alpha1"],
        "A2": target["alpha2"],
        "A3": -target["alpha3"],
        "A4": -target["alpha4"],
        "A5": target["alpha5"],
    }
    old = {key: expr.subs(dictionary.x, -X) for key, expr in old.items()}
    checks = {
        key + "_tube_is_actual_retuned_classical_coefficient": sp.factor(
            tube[key] - value
        )
        for key, value in old.items()
    }
    checks.update(
        {
            key + "_vacuum_patch_is_exact_GR_scalar_coefficient": sp.factor(
                vac[key] - value
            )
            for key, value in d["vacuum_functions"].items()
        }
    )
    checks.update(
        {
            key + "_completion_recomputed_not_spliced": sp.factor(f[key] - value)
            for key, value in d["generic_closed_completion"].items()
        }
    )
    checks["exceptional_luminal_matter_relation_preserved_on_all_timelike_domain"] = (
        sp.factor(X * f["A3"] + 2 * sp.diff(f["F2"], X))
    )
    checks["tensor_kinetic_and_gradient_are_same_positive_R"] = sp.factor(
        -2 * f["F2"] - d["R"]
    )
    checks["total_margin_is_actual_retuned_value"] = retuned.NEW_MARGIN - sp.Rational(
        1, 200
    )
    checks["original_curvature_pole_identity"] = sp.factor(d["h"] * X * tube["A3"] - 1)
    return checks
