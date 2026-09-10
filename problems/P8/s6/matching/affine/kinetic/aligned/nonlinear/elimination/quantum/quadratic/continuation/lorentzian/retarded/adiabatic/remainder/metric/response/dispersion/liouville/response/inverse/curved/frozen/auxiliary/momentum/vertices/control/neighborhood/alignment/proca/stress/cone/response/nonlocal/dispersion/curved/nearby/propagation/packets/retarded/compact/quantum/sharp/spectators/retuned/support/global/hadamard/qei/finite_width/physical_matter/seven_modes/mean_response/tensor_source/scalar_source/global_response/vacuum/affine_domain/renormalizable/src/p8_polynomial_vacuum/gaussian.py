"""Exact finite-regulator Gaussian integration and full Hessian Schur identity."""

from functools import cache
from itertools import permutations

import sympy as sp


@cache
def data():
    a, b, c = sp.symbols("positive_K11 K12 positive_K22", real=True)
    d, e, f = sp.symbols("light_K11 light_K12 light_K22", real=True)
    G, L = sp.symbols("cubic_G quartic_L", real=True)
    phi = sp.Matrix(sp.symbols("phi1 phi2", real=True))
    heavy = sp.Matrix(sp.symbols("H1 H2", real=True))
    K = sp.Matrix([[a, b], [b, c]])
    Km = sp.Matrix([[d, e], [e, f]])
    source = sp.Matrix([v * v for v in phi])
    Hstar = -G * K.inv() * source / 2
    UV = (
        (phi.T * Km * phi)[0] / 2
        + (heavy.T * K * heavy)[0] / 2
        + G * (heavy.T * source)[0] / 2
        + L * sum(v**4 for v in phi) / 24
    )
    effective = (
        (phi.T * Km * phi)[0] / 2
        + L * sum(v**4 for v in phi) / 24
        - G * G * (source.T * K.inv() * source)[0] / 8
    )
    shifted = heavy - Hstar
    completed = (shifted.T * K * shifted)[0] / 2 + effective
    fullH = sp.hessian(UV, tuple(phi) + tuple(heavy))
    A = fullH[:2, :2]
    B = fullH[:2, 2:]
    KK = fullH[2:, 2:]
    sub = dict(zip(heavy, Hstar))
    schur = (A - B * K.inv() * B.T).subs(sub, simultaneous=True)
    target = sp.hessian(effective, tuple(phi))
    # Literal 24-permutation determinant avoids generic symbolic pivot
    # choices without altering any scientific algebra implementation.
    literal_det = sum(
        (-1) ** sum(p[i] > p[j] for i in range(4) for j in range(i + 1, 4))
        * sp.prod(fullH[i, p[i]] for i in range(4))
        for p in permutations(range(4))
    )
    target_det = target[0, 0] * target[1, 1] - target[0, 1] * target[1, 0]
    return {
        "stationary_heavy": Hstar,
        "effective_quartic_action": effective,
        "positive_kernel_condition": "K>=mu^2 I; Km>=m^2 I at a fixed Euclidean regulator",
        "effective_action_lower": "S_eff>=Phi.Km.Phi/2 + (lambda4-3G^2/mu^2) sum(Phi_i^4)/24",
        "determinant_scope": "The Gaussian determinant of K is Phi-independent. Mixed heavy/light loops remain in the nonlocal light functional and are not omitted.",
        "checks": {
            "full_Gaussian_square_completion": sp.factor(UV - completed),
            "actual_stationary_heavy_equations": sp.Matrix(
                [sp.diff(UV, H) for H in heavy]
            )
            .subs(sub, simultaneous=True)
            .applyfunc(sp.factor),
            "actual_reduced_stationary_action": sp.factor(
                UV.subs(sub, simultaneous=True) - effective
            ),
            "heavy_kernel_is_field_independent": KK - K,
            "actual_offdiagonal_Hessian": B - G * sp.diag(*phi),
            "full_reduced_Hessian_is_Schur_complement": (schur - target).applyfunc(
                sp.factor
            ),
            "literal_full_block_determinant": sp.factor(
                literal_det.subs(sub, simultaneous=True) - (a * c - b * b) * target_det
            ),
        },
    }
