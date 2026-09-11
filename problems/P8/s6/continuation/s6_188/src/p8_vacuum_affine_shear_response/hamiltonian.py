"""Complete constrained Proca Hamiltonian on a unimodular shear history."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise.modes import KAPPA, MASS

a, m = s.symbols("a mass", positive=True)
k = s.Matrix(s.symbols("k0:3", real=True))
k2 = k.dot(k)
cross = s.Matrix([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
K0 = s.eye(3) / a + k * k.T / (a**3 * m**2)
V0 = a * m**2 * s.eye(3) + cross.T * cross / a
omega2 = m**2 + k2 / a**2


def symmetric(prefix, traceless=False):
    x = s.symbols(prefix + "0:6", real=True)
    return s.Matrix(
        [
            [x[0], x[1], x[2]],
            [x[1], x[3], x[4]],
            [x[2], x[4], -x[0] - x[3] if traceless else x[5]],
        ]
    )


@cache
def data():
    Q = symmetric("Q")
    qk = Q * k
    magnetic = (k.dot(qk)) * Q - qk * qk.T
    adj = Q.adjugate()
    q = s.Symbol("momentum", nonnegative=True)
    longitudinal = (K0 * V0).applyfunc(s.cancel)
    pi = s.Symbol("constraint_divergence", real=True)
    A0 = s.Symbol("A0", real=True)
    # Hamiltonian before solving A0: -A0 div(pi)-a^3 m^2 A0^2/2.
    aux = -A0 * pi - a**3 * m**2 * A0 * A0 / 2
    solved = aux.subs(A0, -pi / (a**3 * m**2))
    checks = {
        "temporal_constraint_equation": s.diff(aux, A0) + pi + a**3 * m**2 * A0,
        "temporal_constraint_positive_reduced_energy": solved
        - pi * pi / (2 * a**3 * m**2),
        "all_direction_magnetic_cross_identity": s.simplify(
            cross.T * cross - (k2 * s.eye(3) - k * k.T)
        ),
        "full_anisotropic_magnetic_adjugate_identity": s.simplify(
            cross.T * adj * cross - magnetic
        ),
        "base_three_polarization_frequency_identity": s.simplify(
            longitudinal - omega2 * s.eye(3)
        ),
        "base_transverse_normalizer": K0[0, 0].subs({k[0]: 0, k[1]: 0, k[2]: q})
        - 1 / a,
        "base_longitudinal_normalizer": K0[2, 2].subs({k[0]: 0, k[1]: 0, k[2]: q})
        - (m * m + q * q / a**2) / (a * m * m),
        "base_transverse_potential": V0[0, 0].subs({k[0]: 0, k[1]: 0, k[2]: q})
        - a * (m * m + q * q / a**2),
        "base_longitudinal_potential": V0[2, 2].subs({k[0]: 0, k[1]: 0, k[2]: q})
        - a * m * m,
        "actual_mass_and_anchor": MASS - 1000,
    }
    return {
        "chart": "g_ij=-a^2(exp gamma)_ij; gamma(t) real symmetric tracefree, compact in time after the common Cauchy neighborhood. Homogeneous external shear, not general spatial response.",
        "full_parent_restriction": "The S6.187 exact pure tensor restriction has X=R=1 and full S=0. The ordinary canonical Proca action is retained; no scalar source surrogate is used.",
        "canonical_variables": "Z=(A_i,pi^i); pi^i=a^3 g_spatial^(ij)(A_j'-partial_j A0), A0=-div(pi)/(a^3 m^2)",
        "complete_Hamiltonian": "H_gamma=1/2[pi^T K_gamma pi+A^T V_gamma A], K_gamma=a^-1 exp(gamma)+kk^T/(a^3 m^2), V_gamma=a m^2 exp(-gamma)+a^-1 C_k^T exp(gamma) C_k",
        "complete_base_K": K0,
        "complete_base_V": V0,
        "frequency_boundary": "K_gamma V_gamma=[m^2+a^-2 k^T exp(-gamma) k]I by the adjugate identity and det exp(gamma)=1. This is an instantaneous Hamiltonian identity, not a time-dependent normal-mode gap or stability theorem.",
        "positivity": "All terms are positive quadratic forms and m>0. If ||gamma||op<=delta, exp(-delta) M0<=M_gamma<=exp(delta) M0, where M=diag(V,K).",
        "checks": checks,
        "gates": {"same_positive_mass": MASS > 0, "same_fixed_kappa": KAPPA == 10**800},
    }
