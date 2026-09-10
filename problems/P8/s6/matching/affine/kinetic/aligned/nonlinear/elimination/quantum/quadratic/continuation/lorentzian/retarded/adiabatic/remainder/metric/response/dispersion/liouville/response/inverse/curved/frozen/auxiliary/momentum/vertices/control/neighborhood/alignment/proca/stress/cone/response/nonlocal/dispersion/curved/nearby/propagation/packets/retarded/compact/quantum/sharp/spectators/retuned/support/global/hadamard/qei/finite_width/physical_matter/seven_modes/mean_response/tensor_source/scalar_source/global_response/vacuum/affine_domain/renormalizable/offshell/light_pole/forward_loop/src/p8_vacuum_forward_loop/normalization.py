"""Literal quartic Hessian trace and all finite-index vertex contractions."""

from functools import cache
from itertools import combinations_with_replacement, product

import sympy as sp
from p8_polynomial_vacuum import gaussian


@cache
def data():
    d = gaussian.data()
    phi = sp.Matrix(sp.symbols("phi1 phi2", real=True))
    G, L = sp.symbols("cubic_G quartic_L", real=True)
    a, b, c = sp.symbols("positive_K11 K12 positive_K22", real=True)
    kd, ke, kf = sp.symbols("light_K11 light_K12 light_K22", real=True)
    r, s, t = sp.symbols(
        "free_light_covariance11 free_light_covariance12 free_light_covariance22",
        real=True,
    )
    inverse = sp.Matrix([[a, b], [b, c]]).inv()
    covariance = sp.Matrix([[r, s], [s, t]])
    W = sp.hessian(d["effective_quartic_action"], tuple(phi)) - sp.Matrix(
        [[kd, ke], [ke, kf]]
    )
    pairs = list(product(range(2), repeat=2))
    delta = lambda i, j: sp.Integer(i == j)

    def expected_vertex(i, j, k, l):
        return L * delta(i, j) * delta(i, k) * delta(i, l) - G * G * (
            delta(i, j) * delta(k, l) * inverse[i, k]
            + delta(i, k) * delta(j, l) * inverse[i, j]
            + delta(i, l) * delta(j, k) * inverse[i, j]
        )

    actual = sp.Matrix(
        [
            [
                sp.diff(d["effective_quartic_action"], phi[i], phi[j], phi[k], phi[l])
                for k, l in pairs
            ]
            for i, j in pairs
        ]
    )
    expected = sp.Matrix(
        [[expected_vertex(i, j, k, l) for k, l in pairs] for i, j in pairs]
    )
    Q = {
        (k, l): sp.Matrix(
            [[expected_vertex(i, j, k, l) for j in range(2)] for i in range(2)]
        )
        for k, l in pairs
    }
    Gamma = -sp.trace(covariance * W * covariance * W) / 4
    rows = {
        "actual_all_sixteen_quartic_tensor_entries": (actual - expected).applyfunc(
            sp.factor
        )
    }
    for i, j, k, l in combinations_with_replacement(range(2), 4):
        target = (
            -sum(
                sp.trace(covariance * Q[a, b] * covariance * Q[c, d])
                for a, b, c, d in ((i, j, k, l), (i, k, j, l), (i, l, j, k))
            )
            / 2
        )
        rows["literal_four_external_derivatives_" + "".join(map(str, (i, j, k, l)))] = (
            sp.factor(sp.diff(Gamma, phi[i], phi[j], phi[k], phi[l]) - target)
        )
    return {
        "actual_quartic_vertex_tensor": actual,
        "full_quadratic_Hessian_insertion": W,
        "one_loop_quartic_trace": Gamma,
        "three_channel_symmetry_factor": sp.Rational(1, 2),
        "four_dimensional_radial_channel_prefactor": 1 / (32 * sp.pi**2),
        "momentum_vertex": "A4=-lambda4+g[(M-s)^(-1)+(M-t)^(-1)+(M-u)^(-1)], including every nonlocal heavy contraction",
        "continuum_scope": "The written arbitrary-index quartic derivative formula Fourier transforms into the full momentum vertex. The exact two-site calculation checks all 16 tensor entries and all five independent fourth derivatives of the full one-loop trace. External self-energy legs are removed by the independently fixed mass-one/unit-residue scheme.",
        "checks": rows,
    }
