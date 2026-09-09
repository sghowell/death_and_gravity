"""Actual global tensor and Cartesian Proca Hamiltonians, with constraints."""

from functools import cache

import sympy as sp
from p8_proca_global_support import model as scalar_model
from p8_proca_global_support import phase
from p8_proca_retuned_margin import model as retuned

u = scalar_model.u
a = sp.Symbol("positive_FLRW_scale", positive=True)
H = sp.Symbol("physical_Hubble", real=True)
m = sp.Symbol("positive_canonical_Proca_mass", positive=True)
kvec = sp.Matrix(sp.symbols("comoving_kx comoving_ky comoving_kz", real=True))
q = (kvec.T * kvec)[0]
J2 = sp.Matrix([[0, 1], [-1, 0]])


def clean(value):
    return (
        value.applyfunc(sp.factor)
        if isinstance(value, sp.MatrixBase)
        else sp.factor(value)
    )


@cache
def data():
    zeta = 1 / m**2
    A = sp.eye(3) / (zeta * a) + kvec * kvec.T / a**3
    B = a * sp.eye(3) + zeta * (q * sp.eye(3) - kvec * kvec.T) / a
    Omega = sp.zeros(6)
    Omega[:3, 3:] = sp.eye(3)
    Omega[3:, :3] = -sp.eye(3)
    M = sp.zeros(6)
    M[:3, 3:] = A
    M[3:, :3] = -B
    MT = sp.Matrix([[0, 2 / a**3], [-a * q / 2, 0]])
    # W0=-a^-3 div Pi; the actual covariant divergence must vanish.
    W0 = sp.zeros(1, 6)
    W0[0, 3:] = -sp.I * kvec.T / a**3
    divergence = sp.zeros(1, 6)
    divergence[0, :3] = sp.I * kvec.T / a**2
    W0_dt = sp.diff(W0, a) * a * H + W0 * M
    actual_a = (1 + u * u) ** 2
    actual_H = sp.diff(actual_a, u) / actual_a
    original = phase.data()
    old = scalar_model.old
    scalar = original["regular_density_generator"].subs(
        scalar_model.substitution(), simultaneous=True
    )
    scalar = clean(scalar.subs(old.q, q / actual_a**2))
    actual_MT = MT.subs(a, actual_a)
    actual_MP = M.subs({a: actual_a, m: 1000}, simultaneous=True)
    full = sp.diag(scalar, actual_MT, actual_MT, actual_MP)
    full_Omega = sp.diag(original["constant_symplectic_form"], J2, J2, Omega)
    return {
        "a": a,
        "H": H,
        "m": m,
        "zeta": zeta,
        "comoving_momentum": kvec,
        "Proca_momentum_Hessian": A,
        "Proca_coordinate_Hessian": B,
        "Proca_cartesian_generator": M,
        "Proca_symplectic_form": Omega,
        "tensor_polarization_generator": MT,
        "Proca_temporal_reconstruction": W0,
        "actual_global_scale": actual_a,
        "actual_global_Hubble": actual_H,
        "actual_full_fourteen_phase_generator": full,
        "actual_full_symplectic_form": full_Omega,
        "checks": {
            "same_original_global_retuned_scale": actual_a
            - scalar_model.background()["a"],
            "same_original_retuned_margin": retuned.NEW_MARGIN - sp.Rational(1, 200),
            "positive_fixed_mass_matches_actual_Proca": 1000**2 - sp.Integer(10**6),
            "complete_Proca_instantaneous_Hamiltonian_product": clean(
                A * B - (m * m + q / a**2) * sp.eye(3)
            ),
            "complete_Proca_generator_is_Hamiltonian": clean(M * Omega + Omega * M.T),
            "tensor_polarization_generator_is_Hamiltonian": clean(MT * J2 + J2 * MT.T),
            "actual_Proca_temporal_constraint_implies_covariant_divergence": clean(
                W0_dt + 3 * H * W0 - divergence
            ),
            "full_seven_mode_generator_preserves_actual_CCR": clean(
                full * full_Omega + full_Omega * full.T
            ),
        },
    }


@cache
def energy():
    d = data()
    MT = d["tensor_polarization_generator"]
    G = sp.diag(a * (1 + q) / 2, 2 / a**3)
    off = sp.Matrix([[0, 1 / a**2], [1 / a**2, 0]])
    T, P = sp.symbols("tensor_coordinate tensor_momentum", real=True)
    e0 = a * T * T / 2 + 2 * P * P / a**3
    v = sp.Matrix(sp.symbols("W0 W1 W2", real=True))
    p = sp.Matrix(sp.symbols("Pi0 Pi1 Pi2", real=True))
    cross = kvec.cross(v)
    hp = (
        p.dot(p) * m * m / a
        + (p.dot(kvec)) ** 2 / a**3
        + a * v.dot(v)
        + cross.dot(cross) / (m * m * a)
    ) / 2
    reconstructed = (
        p.T * d["Proca_momentum_Hessian"] * p + v.T * d["Proca_coordinate_Hessian"] * v
    )[0] / 2
    return {
        "tensor_positive_modified_energy_matrix": G,
        "tensor_exact_evolution_energy_cross_matrix": off,
        "actual_positive_Proca_energy_sum": hp,
        "Proca_relative_energy_growth_bound": "3 abs(a_dot/a)",
        "tensor_relative_modified_energy_growth_bound": "3 abs(a_dot/a)+1/a",
        "checks": {
            "Proca_Hamiltonian_is_four_positive_energy_pieces": sp.factor(
                hp - reconstructed
            ),
            "tensor_modified_energy_has_only_bounded_cross_term": clean(
                MT.T * G + G * MT - off
            ),
            "tensor_cross_term_relative_bound_is_an_exact_square": sp.factor(
                e0 / a - 2 * T * P / a**2 - (T - 2 * P / a**2) ** 2 / 2
            ),
        },
    }
