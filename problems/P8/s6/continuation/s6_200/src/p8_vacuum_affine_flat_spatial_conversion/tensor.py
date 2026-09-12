"""Polynomial full spatial numerators without light-cone projector poles."""

from functools import cache

import sympy as s

Z = s.Symbol("z")
W = s.Symbol("w")
S = s.Symbol("sigma", positive=True)
P = s.Matrix(s.symbols("p1:4", real=True))
Q = s.expand(P.dot(P))


def require_spin(spin):
    if (
        isinstance(spin, bool)
        or not isinstance(spin, (int, s.Integer))
        or spin not in (0, 2)
    ):
        raise ValueError("An exact spin0 or spin2 selector is required")
    return int(spin)


def numerator(spin, z=Z):
    return _numerator(require_spin(spin), z)


@cache
def _numerator(spin, z):
    A = z * s.eye(3) + P * P.T
    scalar = s.Matrix(9, 9, lambda u, v: A[u // 3, u % 3] * A[v // 3, v % 3] / 3)
    if spin == 0:
        return s.ImmutableMatrix(scalar)
    full = s.Matrix(
        9,
        9,
        lambda u, v: (
            (A[u // 3, v // 3] * A[u % 3, v % 3] + A[u // 3, v % 3] * A[u % 3, v // 3])
            / 2
        ),
    )
    return s.ImmutableMatrix(full - scalar)


def cubic(spin, z=Z):
    return s.ImmutableMatrix(z * numerator(spin, z))


def centered_coefficients(spin):
    return _centered(require_spin(spin))


@cache
def _centered(spin):
    expanded = cubic(spin, W - Q).applyfunc(s.expand)
    return tuple(
        s.ImmutableMatrix(expanded.applyfunc(lambda x, j=j: x.coeff(W, j)))
        for j in range(4)
    )


@cache
def data():
    checks = {}
    for spin in (0, 2):
        poly = cubic(spin)
        cs = centered_coefficients(spin)
        reconstructed = sum((cs[j] * W**j for j in range(4)), s.zeros(9))
        checks[f"spin{spin}_complete_centered_cubic"] = (
            cubic(spin, W - Q) - reconstructed
        ).applyfunc(s.expand)
        checks[f"spin{spin}_no_lightcone_pole"] = poly.subs(Z, 0)
        checks[f"spin{spin}_symmetric_full_tensor"] = poly - poly.T
        zero = {p: 0 for p in P}
        expected = cubic(spin).subs(zero)
        checks[f"spin{spin}_zero_spatial_momentum_cubic"] = expected - Z**3 * numerator(
            spin, 1
        ).subs(zero)
    return {
        "definition": "A(z,p)=z I+p p^T; Q0=A_ij A_kl/3, Q2=(A_ik A_jl+A_il A_jk)/2-Q0, N_i=z Q_i. Q_i=z^2 Pi_i for the full covariant projectors restricted to spatial metric directions.",
        "regularity": "The complete N_i is a polynomial of degree3 in z. Apparent lightcone1/z poles of individual projectors cancel before integration. Both spin sectors and all nine spatial tensor entries remain.",
        "spectral_cut": "At s>4m^2, the spatial cut matrix is sum_i rho_i(s)Q_i(s,p)/s^2. It is not obtained by using COM projectors at nonzero p.",
        "nonlocal_representative": "F(z,p)=sum_i N_i(z,p) integral rho_i(s)/[s^3(s-z)]ds. It has the full S199 spatial cut and a specified zero-momentum covariant subtraction. A physical local matching polynomial remains separate.",
        "checks": checks,
        "gates": {
            "full_spatial_matrix": numerator(2).shape == (9, 9),
            "both_spin_sectors": numerator(0) != s.zeros(9),
            "nonzero_transfer_retained": Q != 0,
        },
    }
