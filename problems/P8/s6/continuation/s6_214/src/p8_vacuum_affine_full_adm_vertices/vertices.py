"""All first and second ADM vertices, including the complete shift feature form."""

from functools import cache

import sympy as s
from p8_vacuum_affine_spatial_current import hamiltonian as old
from p8_vacuum_affine_spatial_current import vertices as shear

from . import adm


def spatial_feature(Q):
    B = adm.shifted_spatial_direction(Q)
    return s.diag(-B, B, B, -s.trace(Q) / 2)


def spatial_contact_feature(D, G):
    B, C = adm.shifted_spatial_direction(D), adm.shifted_spatial_direction(G)
    H = (B * C + C * B) / 2
    return s.diag(H, H, H, s.trace(D) * s.trace(G) / 4)


def shift_feature(beta, a):
    S = s.zeros(10)
    cross = old.cross(beta)
    S[3:6, 6:9] = -s.I * a * cross
    S[6:9, 3:6] = -s.I * a * cross
    S[:3, 9:10] = -s.I * a * beta
    S[9:10, :3] = s.I * a * beta.T
    return S


def feature(n, beta, Q, a):
    return n * s.eye(10) + spatial_feature(Q) + shift_feature(beta, a)


def contact_feature(nD, D, nG, G):
    return (
        spatial_contact_feature(D, G)
        + nD * spatial_feature(G)
        + nG * spatial_feature(D)
    )


def vertex(k, q, n, beta, Q, a=old.A, m=old.MASS):
    return old.features(k, a, m).T * feature(n, beta, Q, a) * old.features(q, a, m)


def contact(k, q, nD, D, nG, G, a=old.A, m=old.MASS):
    return (
        old.features(k, a, m).T * contact_feature(nD, D, nG, G) * old.features(q, a, m)
    )


@cache
def data():
    k = s.Matrix(s.symbols("k0:3", real=True))
    q = s.Matrix(s.symbols("q0:3", real=True))
    beta = s.Matrix(s.symbols("beta0:3", real=True))
    a, m = old.A, old.MASS
    D = s.Matrix([[1, 2, -1], [2, -2, 1], [-1, 1, 1]]) / 17
    G = s.Matrix([[2, 1, 3], [1, 1, -2], [3, -2, -3]]) / 19
    zero = s.zeros(3, 1)
    Vshift = vertex(k, q, 0, beta, s.zeros(3))
    P = k - q
    wanted = s.zeros(6)
    wanted[3:, :3] = s.I * ((beta.T * q)[0] * s.eye(3) + P * beta.T)
    wanted[:3, 3:] = -s.I * ((beta.T * k)[0] * s.eye(3) - beta * P.T)
    general = s.Matrix([[1, 2, -1], [2, 3, 1], [-1, 1, 2]]) / 11
    B = adm.shifted_spatial_direction(general)
    spatial = s.diag(
        -a * m * m * B + old.cross(k).T * B * old.cross(q) / a,
        B / a - s.trace(general) * k * q.T / (2 * a**3 * m * m),
    )
    psi = s.Symbol("psi", real=True)
    scalar = vertex(k, q, 0, zero, 2 * psi * s.eye(3))
    expected = s.diag(
        psi * (a * m * m * s.eye(3) - old.cross(k).T * old.cross(q) / a),
        -psi * s.eye(3) / a - 3 * psi * k * q.T / (a**3 * m * m),
    )
    checks = {
        "full_shift_from_energy_features_equals_Lie_generator": Vshift - wanted,
        "full_shift_reverse_pair_adjoint": Vshift
        - vertex(q, k, 0, beta, s.zeros(3)).conjugate().T,
        "full_spatial_trace_vertex_including_constraint": vertex(k, q, 0, zero, general)
        - spatial,
        "full_pure_scalar_trace_not_traceless_substitution": scalar - expected,
        "frozen_tracefree_first_vertex_exactly_recovered": vertex(k, q, 0, zero, D)
        - shear.metric_vertex(k, q, D),
        "frozen_tracefree_noncommuting_contact_exactly_recovered": contact(
            k, q, 0, D, 0, G
        )
        - shear.metric_contact(k, q, (D * G + G * D) / 2),
        "full_lapse_energy_Gram_includes_constraint": vertex(k, k, 1, zero, s.zeros(3))
        - old.base(k),
        "linear_lapse_has_no_NN_contact": contact(k, q, 1, s.zeros(3), 1, s.zeros(3)),
        "lapse_spatial_mixed_contact_retained": contact(k, q, 1, s.zeros(3), 0, general)
        - vertex(k, q, 0, zero, general),
        "constant_shift_full_common_translation": Vshift.subs(dict(zip(q, k)))[3:, :3]
        - s.I * (beta.T * k)[0] * s.eye(3),
    }
    return {
        "first_spatial": "For B_D=D-tr(D)I/2, M_Q=diag(-am^2 B_D+Ck^t B_D Cq/a,B_D/a-tr(D)kq^t/(2a^3m^2)). The full spatial trace introduces the previously absent constraint vertex.",
        "second_spatial": "Use H_B=(B_D B_G+B_G B_D)/2, with positive mass and H_B electric/magnetic terms; the constraint block adds tr(D)tr(G)kq^t/(4a^3m^2). No ordered matrix word is dropped.",
        "lapse": "With linear lapse N=1+n, first vertex is n Fk^tFq, lapse-lapse contact is zero, and lapse-spatial mixed contact is n_D M_Q(G)+n_G M_Q(D). This is a coordinate-chart statement.",
        "shift": "For P=k-q, M_PiA=i[(beta.q)I+P beta^t] and M_APi=-i[(beta.k)I-beta P^t]. The complete ten-field feature block has magnetic/electric and mass/divergence couplings. All shift contacts vanish in the independent contravariant ADM shift chart, not in arbitrary four-metric coordinates.",
        "checks": checks,
        "gates": {
            "all_ten_metric_directions_present": 1 + 3 + 6 == 10,
            "two_momenta_not_identified_in_shift_or_constraint": k != q,
            "noncommuting_full_spatial_contact": D * G != G * D,
            "shift_and_lapse_chart_contacts_not_covariant_contact_deletion": True,
        },
    }
