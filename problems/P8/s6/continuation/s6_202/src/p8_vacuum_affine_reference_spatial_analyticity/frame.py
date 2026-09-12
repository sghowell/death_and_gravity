"""Analytic full polarization frame with genuine complex Euclidean bounds."""

from functools import cache

import sympy as s

from . import domain

X, Y, Z = s.symbols("nx ny nz")
N = s.Matrix([X, Y, Z])
BASE = s.Matrix([0, 0, 1])
UNIT = X * X + Y * Y + Z * Z - 1


def reduce_unit(value):
    num, den = s.fraction(s.cancel(value))
    _, rem = s.div(num, UNIT, X, Y, Z)
    return s.cancel(rem / den)


@cache
def rotation():
    K = N * BASE.T - BASE * N.T
    return s.ImmutableMatrix(s.eye(3) + K + K * K / (1 + Z))


@cache
def constants():
    b = 2 * domain.DELTA
    eta = 2 * b + b * b
    ndev = (b + eta) / s.Rational(99, 100)
    safe = s.Rational(1, 10**5)
    return {
        "relative_vector_perturbation": b,
        "relative_bilinear_length_squared_defect": eta,
        "direction_displacement_upper": ndev,
        "safe_direction_displacement": safe,
        "rotation_minus_identity_upper": 2 * safe + 4 * safe**2 / (2 - safe),
    }


@cache
def data():
    R = rotation()
    c = constants()
    checks = {
        "complete_complex_orthogonal_frame": (R.T * R - s.eye(3)).applyfunc(
            reduce_unit
        ),
        "frame_sends_base_direction_to_new_direction": (R * BASE - N).applyfunc(
            reduce_unit
        ),
        "orientation_preserved": reduce_unit(R.det() - 1),
        "both_transverse_directions": (N.T * R[:, :2]).applyfunc(reduce_unit),
        "high_band_spatial_radius_factor": c["relative_vector_perturbation"]
        - 2 * domain.DELTA,
    }
    return {
        "construction": "For each fixed real k choose a real orthonormal frame with n0=-k/|k|. Let rho=sqrt(l.l) on the branch rho(0)=|k| and n=l/rho. Set K=n n0^T-n0 n^T and R=I+K+K^2/(1+n.n0). The two transverse vectors are R e1,R e2 and the longitudinal one is n.",
        "analyticity": "On the joint p-ball, nu<=2|k| gives ||p||/|k|<=2delta. Thus |rho^2/|k|^2-1|<=2b+b^2<1e-5, rho stays on its positive branch, and ||n-n0||<1e-5. The denominator1+n.n0 is bounded away from zero. This is a local analytic frame for each fixed real k, not an assumed global frame over the momentum sphere.",
        "complex_norm": "Bilinear n.n=1 does not imply Euclidean norm1. Use the explicit displacement: ||n||<2, ||K||op<=2e-5 and ||R-I||op<1e-4. Hence every polarization has genuine complex Euclidean norm<2.",
        "polarization_sum": "R^T R=I and det R=1 on n.n=1. Summing the two degenerate transverse modes and the longitudinal mode removes the arbitrary initial real transverse-frame choice. All physical polarizations remain.",
        "constants": c,
        "checks": checks,
        "gates": {
            "bilinear_length_defect_safe": c["relative_bilinear_length_squared_defect"]
            < s.Rational(1, 10**5),
            "length_modulus_above_point_nine_nine": 1
            - c["relative_bilinear_length_squared_defect"]
            > s.Rational(99, 100) ** 2,
            "direction_displacement_safe": c["direction_displacement_upper"]
            < c["safe_direction_displacement"],
            "frame_denominator_safe": 2 - c["safe_direction_displacement"]
            > s.Rational(199, 100),
            "rotation_norm_below_two": c["rotation_minus_identity_upper"]
            < s.Rational(1, 10000),
        },
    }
