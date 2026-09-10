"""A distinct rank-regular analytic switch for the same reduced clock jets."""

from functools import cache

import sympy as sp
from p8_affine import connection, dictionary
from p8_exceptional_vacuum import analytic as previous
from p8_exceptional_vacuum import family as original

u, X = original.u, original.X
N = previous.ORDER


@cache
def data():
    T = X**N / (X**N + (1 - X) ** N)
    bump = N * X**2 * sp.exp(-N * X**2)
    B = T + (1 - T) * bump
    S = (1 - T) * (1 - bump)
    h = (1 + u * u) ** 3
    R = 1 + B * (X - 1) / h
    RX = sp.diff(R, X)
    tree = original.data()["original_retuned_tree_scalar"]
    F = tree + sp.exp(-(u**4)) * S * (previous.vacuum_lower_function() - tree)
    return {
        "selected_integer_order": N,
        "T": T,
        "bump": bump,
        "B": B,
        "S": S,
        "R": R,
        "RX": RX,
        "K": sp.Integer(0),
        "A1": sp.Integer(0),
        "A2": sp.Integer(0),
        "F": F,
        "F2": -R / 2,
        "A3": RX / X,
        "A4": -RX / X - sp.Rational(7, 4) * RX**2 / R,
        "A5": RX**2 / (R * X),
        "tensor_factor_lower": sp.Rational(1, 2) + sp.Rational(1, 8 * N),
        "additional_quotient_factor_lower": sp.Rational(1, 4 * N),
        "all_real_u_domain": "-1/(4n)<X<6/5, every even n>=1024",
        "no_old_kinetic_affine_parent_transfer": True,
    }


@cache
def old_obstruction():
    n = sp.Integer(N)
    midpoint = sp.Rational(1, 4) + sp.Rational(3, 4) * sp.Rational(15, 16) ** n
    half = sp.Rational(1, 2) + sp.Rational(1, 2) * sp.Rational(3, 4) ** n
    p = connection.P
    R = sp.Symbol("positive_reduced_R", positive=True)
    C = sp.Symbol("positive_inverse_normalization", positive=True)
    return {
        "old_R_at_zero": sp.Integer(1),
        "old_R_at_quarter": midpoint,
        "old_R_at_half": half,
        "old_R_at_one": sp.Integer(1),
        "full_quotient_determinant": connection.quotient()["determinant"],
        "additional_quotient_factor": 2 * R - 1,
        "vacuum_c_pole_numerator": 2 * C - 1,
        "vacuum_quartic_double_pole_numerator": 2 * C * C - sp.Rational(1, 2),
        "bounds": {
            "old_midpoint_below_one_half": midpoint < sp.Rational(1, 2),
            "old_half_point_above_one_half": half > sp.Rational(1, 2),
        },
        "checks": {
            "actual_quotient_additional_factor": (8 * p * p - 1).subs(p, sp.sqrt(R) / 2)
            - (2 * R - 1),
            "regular_c_forces_unique_positive_normalization": (2 * C - 1).subs(
                C, sp.Rational(1, 2)
            ),
            "same_normalization_removes_quartic_double_pole": (
                2 * C * C - sp.Rational(1, 2)
            ).subs(C, sp.Rational(1, 2)),
        },
    }


@cache
def principal_checks():
    z = sp.Symbol("negative_source_X", negative=True)
    r = sp.Function("positive_R")(z)
    p = sp.sqrt(r) / 2
    c = (sp.sqrt(r) - r) / z
    f4 = (1 - r) / (2 * z * z)
    d = dictionary.principal(p, sp.diff(p, z), c, sp.diff(c, z), f4, z)
    rx = sp.diff(r, z)
    expected = {
        "f": r / 2,
        "Delta": 1,
        "alpha1": 0,
        "alpha2": 0,
        "alpha3": -rx / z,
        "alpha4": rx / z + sp.Rational(7, 4) * rx**2 / r,
        "alpha5": -(rx**2) / (r * z),
    }
    return {
        key + "_literal_principal_dictionary": sp.factor(d[key] - value)
        for key, value in expected.items()
    }


@cache
def local_checks():
    d = data()
    h = (1 + u * u) ** 3
    checks = {
        "switch_decomposition": sp.factor(d["S"] + d["B"] - 1),
        "vacuum_R": d["R"].subs(X, 0) - 1,
        "vacuum_RX": d["RX"].subs(X, 0),
        "vacuum_A3_removable_limit": sp.factor(
            sp.diff(d["RX"], X).subs(X, 0) + 2 * N / h
        ),
        "vacuum_A5_removable_limit": sp.factor(sp.diff(d["RX"] ** 2, X).subs(X, 0)),
        "clock_R": d["R"].subs(X, 1) - 1,
    }
    for j in range(6):
        checks["rational_step_vacuum_jet_" + str(j)] = sp.diff(d["T"], X, j).subs(X, 0)
        checks["complement_clock_jet_" + str(j)] = sp.diff(d["S"], X, j).subs(X, 1)
    return checks
