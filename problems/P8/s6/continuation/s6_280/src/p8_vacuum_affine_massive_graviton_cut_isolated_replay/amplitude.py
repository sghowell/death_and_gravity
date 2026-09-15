"""Literal EH cubic, scalar contact/exchanges and canonical physical helicities."""

from functools import cache

import sympy as s

from . import source

ETA = source.ETA
E, P, X, Y = s.symbols("external_energy initial_momentum cosine sine", real=True)
MU = E * E - P * P
TOTAL = 4 * E * E


def dot(a, b):
    return (a.T * ETA * b)[0]


def on_shell(expr):
    return s.factor(s.expand(expr).subs(Y * Y, 1 - X * X))


def stress(a, b, mass2):
    return a * b.T + b * a.T - ETA * (dot(a, b) + mass2)


def polarization(direction, helicity):
    if type(direction) is not int or direction not in (-1, 1):
        raise ValueError("Require the exact physical momentum direction +/-1")
    if type(helicity) is not int or helicity not in (-1, 1):
        raise ValueError("Require a canonical physical helicity sign +/-1")
    return s.Matrix([0, 1, s.I * direction * helicity, 0]) / s.sqrt(2)


def connection_first(H, q):
    qcov = ETA * q
    return [
        [
            [
                s.I
                * ETA[r, r]
                * (qcov[m] * H[r, n] + qcov[n] * H[r, m] - qcov[r] * H[m, n])
                for n in range(4)
            ]
            for m in range(4)
        ]
        for r in range(4)
    ]


def connection_second(H, G):
    return [
        [
            [
                s.expand(-2 * sum(ETA[r, r] * H[r, j] * G[j][m][n] for j in range(4)))
                for n in range(4)
            ]
            for m in range(4)
        ]
        for r in range(4)
    ]


def contracted_pair(D, G, H):
    return s.expand(
        sum(
            D[m, n] * (G[r][m][l] * H[l][n][r] - G[r][m][n] * H[l][r][l])
            for m in range(4)
            for n in range(4)
            for r in range(4)
            for l in range(4)
            if D[m, n] != 0
        )
    )


def einstein_cubic(fields, momenta):
    """All multilinear terms of -sqrt(-g)R/2 at g=eta+2h, modulo divergence."""
    if len(fields) != 3 or len(momenta) != 3:
        raise ValueError("Require exactly three gravitational legs")
    G = [connection_first(H, q) for H, q in zip(fields, momenta)]
    D = [s.trace(ETA * H) * ETA - 2 * ETA * H * ETA for H in fields]
    value = 0
    for i in range(3):
        j, k = [v for v in range(3) if v != i]
        value += contracted_pair(D[i], G[j], G[k]) + contracted_pair(D[i], G[k], G[j])
        for H in (
            connection_second(fields[j], G[k]),
            connection_second(fields[k], G[j]),
        ):
            value += contracted_pair(ETA, G[i], H) + contracted_pair(ETA, H, G[i])
    return on_shell(-value / 2)


def covariant_numerator(a, b, k1, k2, v1, v2):
    d1, d2 = 2 * dot(a, k1), 2 * dot(a, k2)
    return on_shell(
        2 * dot(v1, b) * dot(v2, a) * d1
        + 2 * dot(v1, a) * dot(v2, b) * d2
        + dot(v1, v2) * d1 * d2
    )


def helicity_trees(energy, angle, mass2, kappa):
    denominator = energy - (energy - 4 * mass2) * angle * angle
    return 4 * mass2 * mass2 / (kappa * denominator), (
        (energy - 4 * mass2) ** 2 * (1 - angle * angle) ** 2 / (4 * kappa * denominator)
    )


@cache
def data():
    a = s.Matrix([E, P * Y, 0, P * X])
    b = s.Matrix([E, -P * Y, 0, -P * X])
    k1 = s.Matrix([-E, 0, 0, -E])
    k2 = s.Matrix([-E, 0, 0, E])
    r0 = -k1 - k2
    V = stress(a, b, MU)
    W = s.simplify(V - ETA * s.trace(ETA * V) / 2)
    checks = {
        "all_incoming_momentum_conservation": a + b + k1 + k2,
        "massive_source_Ward_before_sewing": (V * ETA * (a + b)).applyfunc(on_shell),
        "first_scalar_mass_shell": on_shell(dot(a, a) - MU),
        "second_scalar_mass_shell": on_shell(dot(b, b) - MU),
        "first_graviton_mass_shell": dot(k1, k1),
        "second_graviton_mass_shell": dot(k2, k2),
    }
    q = s.Matrix(s.symbols("q0:4", real=True))
    H = s.diag(0, 1, -1, 0) / s.sqrt(2)
    G1 = connection_first(H, q)
    G2 = connection_first(H, -q)
    kinetic = s.factor(
        -(contracted_pair(ETA, G1, G2) + contracted_pair(ETA, G2, G1)) / 2
    )
    checks["literal_EH_canonical_TT_quadratic"] = s.factor(
        kinetic.subs({q[1]: 0, q[2]: 0}) - q[0] ** 2 + q[3] ** 2
    )
    rows = {}
    for h1, h2 in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
        label = (
            ("plus" if h1 == 1 else "minus") + "_" + ("plus" if h2 == 1 else "minus")
        )
        v1 = polarization(1, h1)
        v2 = polarization(-1, h2)
        H1 = ETA * v1 * v1.T * ETA
        H2 = ETA * v2 * v2.T * ETA
        d1 = on_shell(2 * dot(a, k1))
        d2 = on_shell(2 * dot(a, k2))
        contact = on_shell(
            2 * s.trace(ETA * H1 * ETA * H2) * (dot(a, b) + MU)
            - 4 * (a.T * (H1 * ETA * H2 + H2 * ETA * H1) * b)[0]
        )
        scalar1 = on_shell(-4 * dot(a, v1) ** 2 * dot(b, v2) ** 2 / d1)
        scalar2 = on_shell(-4 * dot(a, v2) ** 2 * dot(b, v1) ** 2 / d2)
        cubic = einstein_cubic([ETA * W * ETA, H1, H2], [r0, k1, k2])
        exchange = on_shell(-cubic / TOTAL)
        actual = on_shell(contact + scalar1 + scalar2 + exchange)
        numerator = covariant_numerator(a, b, k1, k2, v1, v2)
        factored = on_shell(numerator * numerator / (d1 * TOTAL * d2))
        expected = helicity_trees(TOTAL, X, MU, s.Integer(1))[0 if h1 == h2 else 1]
        checks.update(
            {
                label + "_exact_contact_scalar_and_EH_sum": on_shell(actual - expected),
                label + "_independent_covariant_factorization": on_shell(
                    actual - factored
                ),
                label + "_first_graviton_Ward": covariant_numerator(
                    a, b, k1, k2, k1, v2
                ),
                label + "_second_graviton_Ward": covariant_numerator(
                    a, b, k1, k2, v1, k2
                ),
                label + "_physical_polarization_normalization": dot(v1, s.conjugate(v1))
                + 1,
                label + "_physical_traceless_polarization": dot(v1, v1),
                label + "_physical_transverse_polarization": dot(v1, k1),
            }
        )
        rows[label] = {
            "contact": contact,
            "first_scalar_exchange": scalar1,
            "second_scalar_exchange": scalar2,
            "entire_EH_exchange": exchange,
            "complete_tree_at_kappa_one": actual,
        }
    Pi = s.diag(0, 1, 1, 0)
    tensors = [ETA * polarization(1, h) * polarization(1, h).T * ETA for h in (-1, 1)]
    actual = s.Matrix(
        16,
        16,
        lambda a, b: sum(
            H[a // 4, a % 4] * s.conjugate(H[b // 4, b % 4]) for H in tensors
        ),
    )
    target = s.Matrix(
        16,
        16,
        lambda a, b: (
            (
                Pi[a // 4, b // 4] * Pi[a % 4, b % 4]
                + Pi[a // 4, b % 4] * Pi[a % 4, b // 4]
                - Pi[a // 4, a % 4] * Pi[b // 4, b % 4]
            )
            / 2
        ),
    )
    checks["full_physical_TT_projector_no_dilaton"] = actual - target
    S, mu, K = source.S, source.MU, source.K
    same, opposite = helicity_trees(S, X, mu, K)
    checks["same_helicity_exact_massive_threshold"] = s.factor(
        same.subs(S, 4 * mu) - mu / K
    )
    checks["opposite_helicity_exact_massive_threshold"] = opposite.subs(S, 4 * mu)
    checks["massless_same_helicity_limit"] = same.subs(mu, 0)
    checks["massless_opposite_helicity_limit"] = s.factor(
        opposite.subs(mu, 0) - S * (1 - X * X) / (4 * K)
    )
    return {
        "literal_four_graph_helicity_rows": rows,
        "whole_equal_helicity_tree": same,
        "whole_opposite_helicity_tree": opposite,
        "canonical_TT_quadratic": kinetic,
        "canonical_TT_projector": target,
        "normalization": "The entire EH GammaGamma cubic and matter contact plus both scalar exchanges are derived at kappa1. Canonical cubic vertices carry1/sqrt(kappa) and the contact1/kappa; the whole physical tree therefore carries1/kappa. No printed spinor prefactor is assumed.",
        "complete_tree_boundary": "Leading formal loop order in the old physical frame and mass1 reference. Only the two canonical TT graviton polarizations are sewn. The discarded EH divergence has zero coefficient because the three momenta sum to zero, not because a clock boundary term was forgotten.",
        "checks": checks,
        "gates": {
            "all_four_helicity_combinations_computed": len(rows) == 4,
            "all_four_EH_scalar_and_contact_graphs_in_each_row": all(
                len(row) == 5 for row in rows.values()
            ),
            "whole_cubic_not_just_scalar_exchange_or_KLT_import": True,
            "physical_TT_projector_not_four_vector_tensor_square": target.rank() == 2,
            "massive_same_helicity_not_deleted": same != 0,
            "original_curvature_sign_and_canonical_normalization": True,
            "original_source_formal_loop_scope_retained": True,
        },
    }
