"""Literal full quadratic half-trace with one light covariance insertion."""

from functools import cache

import sympy as sp
from p8_vacuum_light_pole import normalization


@cache
def data():
    old = normalization.data()
    h = sp.symbols("formal_covariance_insertion")
    r, s, t = sp.symbols(
        "free_light_covariance11 free_light_covariance12 free_light_covariance22",
        real=True,
    )
    p11, p12, p22 = sp.symbols(
        "inserted_covariance11 inserted_covariance12 inserted_covariance22", real=True
    )
    P = sp.Matrix([[p11, p12], [p12, p22]])
    phi = sp.Matrix(sp.symbols("phi1 phi2", real=True))
    a, b, c = sp.symbols("positive_K11 K12 positive_K22", real=True)
    K = sp.Matrix([[a, b], [b, c]])
    G, L = sp.symbols("cubic_G quartic_L", real=True)
    shifted = old["half_trace"].subs(
        {r: r + h * p11, s: s + h * p12, t: t + h * p22}, simultaneous=True
    )
    actual = sp.diff(shifted, h).subs(h, 0)
    local_quartic = L * sum(P[i, i] * phi[i] ** 2 for i in range(2)) / 4
    J = sp.Matrix([z * z for z in phi])
    stationary = -G * G * sum(P[i, i] * (K.inv() * J)[i] for i in range(2)) / 4
    mixed = -G * G * sp.trace(P * sp.diag(*phi) * K.inv() * sp.diag(*phi)) / 2
    M, Tinsert = sp.symbols("heavy_mass_squared whole_regulated_inserted_tadpole")
    h_one_point = G * Tinsert / 2
    h_counter = -h_one_point
    reduced_mass = (L - G * G / M) * Tinsert / 2
    counter_mass = -G * h_counter / M
    z, c0, c1 = sp.symbols("s whole_local_constant whole_local_slope")
    affine = c0 + c1 * z
    checks = {
        "differentiate_full_frozen_half_trace": sp.factor(
            actual - sp.trace(P * old["actual_light_quadratic_Hessian_insertion"]) / 2
        ),
        "all_local_and_mixed_covariance_insertions_retained": sp.factor(
            actual - local_quartic - stationary - mixed
        ),
        "exactly_one_light_covariance_factor": sp.diff(shifted, h, 2),
        "heavy_one_point_reference_anchor": h_one_point + h_counter,
        "stationary_heavy_mass_piece_cancelled": sp.factor(
            reduced_mass + counter_mass - L * Tinsert / 2
        ),
        "whole_outer_affine_reference_cancels": sp.expand(
            affine - affine.subs(z, 1) - (z - 1) * sp.diff(affine, z).subs(z, 1)
        ),
    }
    for i in range(2):
        for j in range(2):
            checks[f"mixed_external_derivatives_{i}_{j}"] = sp.factor(
                sp.diff(mixed, phi[i], phi[j]) + G * G * P[i, j] * K.inv()[i, j]
            )
    return {
        "inserted_light_covariance": P,
        "literal_full_quadratic_insertion": actual,
        "local_quartic_tadpole": local_quartic,
        "stationary_heavy_tadpole": stationary,
        "mixed_heavy_light_trace": mixed,
        "H_one_point_increment": h_one_point,
        "fixed_H_one_point_counterterm": h_counter,
        "remaining_local_scalar_mass": sp.factor(reduced_mass + counter_mass),
        "placement_count": 1,
        "scope": "The full scalar half-trace differentiated once in its light covariance: all its local and mixed terms, with exactly one inserted light line in the nonlocal H-Phi bubble. Not every enlarged-model two-loop quadratic graph.",
        "checks": checks,
    }
