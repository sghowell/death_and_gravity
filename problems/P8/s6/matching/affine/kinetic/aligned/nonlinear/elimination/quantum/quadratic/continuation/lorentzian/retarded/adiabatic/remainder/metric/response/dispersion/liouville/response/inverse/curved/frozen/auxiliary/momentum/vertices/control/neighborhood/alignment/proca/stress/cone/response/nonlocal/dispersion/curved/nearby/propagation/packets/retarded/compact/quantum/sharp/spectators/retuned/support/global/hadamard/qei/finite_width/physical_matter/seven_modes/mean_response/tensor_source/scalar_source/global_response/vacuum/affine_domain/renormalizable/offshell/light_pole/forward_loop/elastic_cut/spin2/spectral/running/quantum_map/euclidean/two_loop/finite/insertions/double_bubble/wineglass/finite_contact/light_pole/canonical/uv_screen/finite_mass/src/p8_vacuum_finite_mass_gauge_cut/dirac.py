"""Exact four-dimensional Clifford, norm and Ward numerator identities."""

from functools import cache

import sympy as sp


@cache
def data():
    I = sp.I
    pauli = [sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -I], [I, 0]]), sp.diag(1, -1)]
    zero = sp.zeros(2)
    eye = sp.eye(2)
    gamma = [sp.BlockMatrix([[zero, eye], [eye, zero]]).as_explicit()]
    gamma += [
        sp.BlockMatrix([[zero, -I * t], [I * t, zero]]).as_explicit() for t in pauli
    ]
    q = sp.symbols("real_loop_q0:4", real=True)
    k = sp.symbols("external_shift_k0:4", real=True)
    m = sp.symbols("positive_fermion_mass", positive=True)
    slash = lambda p: sum((p[i] * gamma[i] for i in range(4)), sp.zeros(4))
    qk = tuple(q[i] + k[i] for i in range(4))
    N = m * sp.eye(4) - I * slash(q)
    Nk = m * sp.eye(4) - I * slash(qk)
    h = m * m + sum(v * v for v in q)
    hk = m * m + sum(v * v for v in qk)
    checks = {}
    for i, g in enumerate(gamma):
        checks[f"Hermitian_gamma_{i}"] = g - g.conjugate().T
        for j, f in enumerate(gamma):
            checks[f"Clifford_{i}_{j}"] = g * f + f * g - 2 * int(i == j) * sp.eye(4)
    checks["propagator_inverse_numerator"] = (
        N * (m * sp.eye(4) + I * slash(q)) - h * sp.eye(4)
    ).applyfunc(sp.expand)
    checks["propagator_norm_numerator"] = (
        N.conjugate().T * N - h * sp.eye(4)
    ).applyfunc(sp.expand)
    checks["Ward_resolvent_numerator"] = (
        N * I * slash(k) * Nk - N * hk + Nk * h
    ).applyfunc(sp.expand)
    e = sp.symbols("real_polarization_e0:4", real=True)
    checks["real_polarization_norm_numerator"] = (
        slash(e).conjugate().T * slash(e) - sum(v * v for v in e) * sp.eye(4)
    ).applyfunc(sp.expand)
    return {
        "gamma_matrices": gamma,
        "loop_components": q,
        "shift_components": k,
        "mass": m,
        "unshifted_propagator_numerator": N,
        "unshifted_positive_denominator": h,
        "norm_squared_identity": "S0(q).dagger S0(q) = I/(q^2+mF^2), for real Euclidean q",
        "complex_shift_bound": "||i gamma.P|| <= sum_mu |P_mu|; ||S0 i gamma.P|| <= 12/mF on the written routing domain",
        "Ward_identity": "S(q) i gamma.k S(q+k) = S(q)-S(q+k); the polynomial numerator identity extends to complex k",
        "scope": "Four-dimensional norm estimates apply to the absolutely convergent Taylor tail, not to the unsubtracted degree-zero UV term. The full gauge-invariant regulated diagram sum is treated before that term is removed.",
        "checks": checks,
    }
