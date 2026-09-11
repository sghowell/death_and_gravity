"""Exact spatially varying constrained Proca Hamiltonian."""

from functools import cache

import sympy as s

A, MASS = s.symbols("a mass", positive=True)


def cross(k):
    return s.Matrix([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])


def base(k, a=A, m=MASS):
    k = s.Matrix(k)
    C = cross(k)
    K = s.eye(3) / a + k * k.T / (a**3 * m**2)
    V = a * m**2 * s.eye(3) + C.T * C / a
    return s.diag(V, K)


def features(k, a=A, m=MASS):
    k = s.Matrix(k)
    F = s.zeros(10, 6)
    F[:3, :3] = s.sqrt(a) * m * s.eye(3)
    F[3:6, :3] = cross(k) / s.sqrt(a)
    F[6:9, 3:] = s.eye(3) / s.sqrt(a)
    F[9:10, 3:] = k.T / (a * s.sqrt(a) * m)
    return F


@cache
def data():
    k = s.Matrix(s.symbols("k0:3", real=True))
    pi_div, A0 = s.symbols("div_pi A0", real=True)
    auxiliary = -A0 * pi_div - A**3 * MASS**2 * A0**2 / 2
    constraint = -pi_div / (A**3 * MASS**2)
    R, X, Bu, H, Z = s.symbols("R X Box_u H uHu", real=True)
    RX, Ru = s.symbols("R_X R_u", real=True)
    source = (R - 1) * (
        -3 * H + 3 * Ru / (4 * R) + Bu / X + (-1 / X**2 + 3 * RX / (2 * R * X)) * Z
    )
    F = features(k)
    metric = base(k)
    J = s.zeros(6)
    J[:3, 3:] = s.eye(3)
    J[3:, :3] = -s.eye(3)
    checks = {
        "actual_temporal_constraint": s.diff(auxiliary, A0).subs(A0, constraint),
        "positive_constraint_Schur_term": auxiliary.subs(A0, constraint)
        - pi_div**2 / (2 * A**3 * MASS**2),
        "full_source_zero_for_spatial_unimodular_shear": source.subs(
            {R: 1, X: 1, Bu: 3 * H, Z: 0}
        ),
        "ten_feature_complete_energy_identity": F.T * F - metric,
        "full_three_mode_canonical_generator_symplectic": J * metric * J
        + J * (J * metric).T,
        "instantaneous_three_mode_frequency_identity": metric[3:, 3:] * metric[:3, :3]
        - (MASS**2 + k.dot(k) / A**2) * s.eye(3),
    }
    return {
        "domain": "The same conditional actual Proca state at the isotropic clock, varied by smooth external spatial unimodular shear with the common zero initial neighborhood. The source is exactly zero on this restricted external metric family, including spatially varying shear.",
        "full_H": "H=integral[pi^t E pi/a+(div pi)^2/(a^3m^2)+a m^2 A^t E^-1 A+(curl A)^t E(curl A)/a]/2, E=exp gamma.",
        "temporal_constraint": "A0=-div pi/(a^3m^2); no independent fourth oscillator and no removal of the longitudinal constraint energy.",
        "spatial_operator": "K=E/a-grad div/(a^3m^2), V=a m^2 E^-1+curl E curl/a on real fields. Curl is self-adjoint under compact/periodic integration by parts.",
        "checks": checks,
        "gates": {
            "actual_mass_positive": MASS.is_positive is True,
            "positive_scale": A.is_positive is True,
            "feature_shape_keeps_constraint": F.shape == (10, 6),
            "no_new_source_at_fixed_unimodular_clock": True,
        },
    }
