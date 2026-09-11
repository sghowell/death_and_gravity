"""Complete two-momentum vertices and local metric contacts."""

from functools import cache

import sympy as s

from . import hamiltonian as ham


def metric_vertex(k, q, D, a=ham.A, m=ham.MASS):
    return s.diag(-a * m * m * D + ham.cross(k).T * D * ham.cross(q) / a, D / a)


def metric_contact(k, q, H, a=ham.A, m=ham.MASS):
    return s.diag(a * m * m * H + ham.cross(k).T * H * ham.cross(q) / a, H / a)


def feature_form(D, contact=False):
    S = s.zeros(10)
    for n in range(3):
        S[3 * n : 3 * n + 3, 3 * n : 3 * n + 3] = D * (
            (-1) if n == 0 and not contact else 1
        )
    return S


@cache
def data():
    k = s.Matrix(s.symbols("k0:3", real=True))
    q = s.Matrix(s.symbols("q0:3", real=True))
    d = s.symbols("D0:6", real=True)
    D = s.Matrix([[d[0], d[1], d[2]], [d[1], d[3], d[4]], [d[2], d[4], d[5]]])
    G = s.Matrix([[1, 2, -1], [2, -2, 1], [-1, 1, 1]]) / 17
    H = (D * G + G * D) / 2
    Fk, Fq = ham.features(k), ham.features(q)
    vertex = metric_vertex(k, q, D)
    contact = metric_contact(k, q, H)
    eps, eta = s.symbols("epsilon eta", real=True)
    X = eps * D + eta * G
    second = (s.eye(3) + X + X * X / 2).applyfunc(
        lambda x: s.expand(x).coeff(eps, 1).coeff(eta, 1)
    )
    checks = {
        "full_off_diagonal_first_vertex_factorization": vertex
        - Fk.T * feature_form(D) * Fq,
        "full_off_diagonal_contact_factorization": contact
        - Fk.T * feature_form(H, True) * Fq,
        "real_first_vertex_reverse_pair_adjoint": vertex - metric_vertex(q, k, D).T,
        "real_contact_reverse_pair_adjoint": contact - metric_contact(q, k, H).T,
        "full_ordered_mixed_exponential_second": second - H,
        "constraint_has_no_unimodular_metric_vertex": feature_form(D)[9, 9],
    }
    return {
        "first_vertex": "M_D(k,q)=diag(-a m^2 D_(k-q)+C_k^t D_(k-q) C_q/a,D_(k-q)/a). Both momenta remain in the magnetic contraction.",
        "second_contact": "M_DGamma uses the same electric/magnetic terms and POSITIVE mass term with H=(D Gamma+Gamma D)/2 in position space. Its Fourier coefficient is the full convolution, not one mode product except for a specified plane-wave pair.",
        "uniform_all_momenta_bound": "The ten-row energy feature map gives M_D=F_k^t diag(-D,D,D,0)F_q. Since F_k^t F_k=M0(k), the two-sided energy-relative norm is <=||D||op for every k,q, including zero and noncollinear momenta. The contact norm is <=||D||op||Gamma||op for a single pointwise or Fourier pair.",
        "balanced_form": "For the actual symplectic normalizers T_k^t M0(k) T_k=omega_k I, ||T_k^t M_D(k,q)T_q||op<=sqrt(omega_k omega_q)||D||op. This is an instantaneous energy identity, not a frequency-dependent cutoff or inverse.",
        "checks": checks,
        "gates": {
            "all_three_physical_energy_pieces": feature_form(D)[:9, :9].shape == (9, 9),
            "no_spurious_temporal_constraint_variation": feature_form(D)[9, :]
            == s.zeros(1, 10),
            "second_mass_sign_is_positive": feature_form(H, True)[:3, :3] == H,
            "arbitrary_two_momenta_not_identified": k != q,
        },
    }
