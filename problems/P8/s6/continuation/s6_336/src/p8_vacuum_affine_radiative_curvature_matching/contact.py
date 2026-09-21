"""Literal linear curvature vertex and its complete Bose-symmetric contact."""

from functools import cache
from itertools import combinations, permutations, product

import sympy as s
from p8_vacuum_affine_dimensional_gravity_radiation import tensor

ETA = s.diag(1, -1, -1, -1)
PAIRS = tuple(combinations(range(4), 2))


def dot(p, q):
    return (p.T * ETA * q)[0]


def exact_coefficient(value):
    if isinstance(value, (bool, float, str)) or not isinstance(value, (int, s.Basic)):
        raise TypeError("Require an explicit exact real curvature coefficient")
    value = s.sympify(value)
    if (
        value.has(s.Float)
        or value.is_real is not True
        or value.is_finite is not True
        or value.free_symbols
    ):
        raise ValueError("Require an explicit finite exact real curvature coefficient")
    return value


def _core(ps, k):
    total = s.zeros(4)
    for i, j in PAIRS:
        v = dot(k, ps[i]) * ps[j] - dot(k, ps[j]) * ps[i]
        total += v * v.T
    return total.applyfunc(s.factor)


def core(momenta, k):
    ps, k = tensor.kinematics(momenta, k)
    return _core(ps, k)


def transverse_core(momenta, k):
    ps, k = tensor.kinematics(momenta, k)
    basis = tensor.frame(k)
    full = _core(ps, k)
    return s.Matrix(
        2, 2, lambda i, j: (basis[i].T * ETA * full * ETA * basis[j])[0]
    ).applyfunc(s.factor)


def linear_curvature(k, h):
    """Return sqrt(kappa)*R^(1) for delta g=2h/sqrt(kappa); internal symbolic API."""
    kl = ETA * k
    return {
        (a, b, c, d): kl[b] * kl[d] * h[a, c]
        + kl[a] * kl[c] * h[b, d]
        - kl[b] * kl[c] * h[a, d]
        - kl[a] * kl[d] * h[b, c]
        for a, b, c, d in product(range(4), repeat=4)
    }


def curvature_pair(curvature, p, q):
    return sum(
        curvature[a, b, c, d] * p[a] * q[b] * p[c] * q[d]
        for a, b, c, d in product(range(4), repeat=4)
    )


def literal_vertex(ps, k, polarization):
    curvature = linear_curvature(k, polarization)
    return s.factor(
        s.Add(
            *(
                curvature_pair(curvature, ps[a], ps[b])
                for a, b, c, d in permutations(range(4))
            )
        )
        / 4
    )


@cache
def data():
    k = s.Matrix(s.symbols("k0:4"))
    p = s.Matrix(s.symbols("p0:4"))
    q = s.Matrix(s.symbols("q0:4"))
    z = s.Matrix(s.symbols("z0:4"))
    h = s.zeros(4)
    for (i, j), value in zip(
        ((row, col) for row in range(4) for col in range(row, 4)), s.symbols("h0:10")
    ):
        h[i, j] = h[j, i] = value
    curvature = linear_curvature(k, h)
    v = dot(k, p) * q - dot(k, q) * p
    swapped = dot(k, q) * p - dot(k, p) * q
    gauge = (ETA * k) * (ETA * z).T + (ETA * z) * (ETA * k).T
    pure = linear_curvature(k, gauge)
    scale = s.Symbol("soft_scale")
    rows = {
        "unrestricted_literal_curvature_dyad_identity": s.expand(
            curvature_pair(curvature, p, q) - (v.T * h * v)[0]
        ),
        "unrestricted_pair_Ward_identity": s.expand(dot(k, v)),
        "all256_literal_pure_gauge_curvature_components": s.Matrix(
            16, 16, lambda i, j: s.expand(pure[i // 4, i % 4, j // 4, j % 4])
        ),
        "pair_exchange_keeps_dyad": (v * v.T - swapped * swapped.T).applyfunc(s.expand),
        "exact_quadratic_soft_scaling": s.expand(
            curvature_pair(linear_curvature(scale * k, h), p, q)
            - scale**2 * curvature_pair(curvature, p, q)
        ),
        "complete_identical_scalar_assignment_count": s.Integer(
            len(tuple(permutations(range(4)))) - 4 * len(PAIRS)
        ),
    }
    return {
        "checks": rows,
        "gates": {
            "unrestricted_algebra_precedes_physical_projection": True,
            "full24_Bose_assignments_not_one_unsymmetrized_operator": True,
            "all_six_pairs_and_pure_gauge_components_retained": True,
            "coefficient_API_requires_explicit_exact_real_argument": True,
            "no_universal_basis_completeness_or_UV_claim": True,
        },
        "whole_comparison_operator": "Delta L=hbar*chi*sqrt(-g)*Phi^2*C_abcd*(nabla^a nabla^c Phi)*(nabla^b nabla^d Phi)/4. C is Weyl; on a null transverse-traceless plane wave its linear response equals Riemann. The coefficient chi is an unspecified real mass^-6 matching coordinate, not a selected new source.",
        "whole_canonical_vertex": "In g=eta+2h/sqrt(kappa), v_ij=(k.p_i)p_j-(k.p_j)p_i and T=sum_(i<j)v_ij v_ij^T. Delta M5^(1)=chi*T:epsilon/sqrt(kappa). Four identical-field assignments per unordered pair cancel the action's factor1/4. All metric/inverse/connection factors beyond linear curvature first enter at two gravitons.",
        "whole_Ward_Bose_and_soft_statement": "k.v_ij=0 without mass-shell assumptions. T is conserved, polynomial and Bose symmetric. It is homogeneous degree2 in emitted k at fixed hard momenta and uniformly O(omega^2) on the original recoil domain. The contact changes none of the leading through sub-subleading soft coefficients.",
    }
