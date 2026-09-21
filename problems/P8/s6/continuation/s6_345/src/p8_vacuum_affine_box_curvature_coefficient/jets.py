"""Literal null-TT metric response of the fully stated scalar-jet basis."""

from functools import cache
from itertools import combinations_with_replacement, permutations, product

import sympy as s
from p8_vacuum_affine_local_tadpole_radiation import vertices as original

from . import basis


def word_contact(word, G, H, a):
    metric = -2 * sum(
        H[i, j] * s.prod(G[ii, jj] for pos, (ii, jj) in enumerate(word) if pos != edge)
        for edge, (i, j) in enumerate(word)
    )
    if word[0] == word[1] == word[2] and word[0][0] != word[0][1]:
        i, j = word[0]
        return (
            metric
            - (3 * G[i, j] + a[j]) * (2 * a[j] * H[i, j] - a[i] * H[j, j])
            - (3 * G[i, j] + a[i]) * (2 * a[i] * H[i, j] - a[j] * H[i, i])
        )
    connection = 0
    for i in range(4):
        incidence = [
            (edge, second if first == i else first)
            for edge, (first, second) in enumerate(word)
            if first == i or second == i
        ]
        degree = sum((first == i) + (second == i) for first, second in word)
        if degree == 2 and len(incidence) == 2:
            (first, j), (second, k) = incidence
            rest = s.prod(
                G[ii, jj]
                for pos, (ii, jj) in enumerate(word)
                if pos not in (first, second)
            )
            connection -= (a[j] * H[i, k] + a[k] * H[i, j] - a[i] * H[j, k]) * rest
        elif degree > 2:
            raise ValueError(
                "Only the stated rank-three triple edge or rank-two words are implemented"
            )
    return metric + connection


def bose_contact(word, G, H, a):
    return s.expand(
        sum(
            word_contact(tuple((p[i], p[j]) for i, j in word), G, H, a)
            for p in basis.PERMS
        )
    )


@cache
def data():
    checks = {}
    pc = s.symbols("p0:4")
    kc = s.symbols("k0:4")
    pairs = tuple(combinations_with_replacement(range(4), 2))
    ev = s.symbols("e0:10")
    E = s.zeros(4)
    for (i, j), value in zip(pairs, ev):
        E[i, j] = E[j, i] = value
    metric = (1, -1, -1, -1)

    def connection(a, b):
        return sum(
            metric[d] * (kc[a] * E[b, d] + kc[b] * E[a, d] - kc[d] * E[a, b]) * pc[d]
            for d in range(4)
        )

    Z = {(a, b): s.expand(connection(a, b)) for a, b in pairs}

    def zij(a, b):
        return Z[tuple(sorted((a, b)))]

    for indices in product(range(4), repeat=3):
        raw = 0
        symmetric = 0
        for a, b, c in permutations(indices):
            raw -= (
                kc[a] * zij(b, c)
                + pc[a] * zij(b, c)
                + pc[c] * zij(a, b)
                + pc[b] * zij(a, c)
            )
            symmetric -= (3 * pc[a] + kc[a]) * zij(b, c)
        checks["literal_symmetrized_third_jet_" + "".join(map(str, indices))] = (
            s.expand((raw - symmetric) / 6)
        )
    # Independent agreement with the existing literal Galileon quartic contact.
    G, H = original.symbolic_grams()
    aa = s.symbols("soft0:4")
    new = -bose_contact(((0, 1), (1, 2), (3, 3)), G, H, aa) + bose_contact(
        ((0, 1), (1, 2), (2, 3)), G, H, aa
    )
    old = original.quartic_contact("Gal", G, H, aa, 0, (0, 0, 0, 0))
    checks["whole_original_Galileon_metric_and_connections"] = s.expand(new - old)
    # The actual null-TT contracted scalar connection, not an off-shell deletion.
    k = s.Matrix([1, 0, 0, 1])
    eta = s.diag(1, -1, -1, -1)
    for pol, E in enumerate(
        (
            s.diag(0, 1, -1, 0),
            s.Matrix([[0, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 0]]),
        )
    ):
        checks[f"null_TT_scalar_connection_trace_{pol}"] = 2 * E * k - k * s.trace(
            eta * E
        )
    return {
        "checks": checks,
        "gates": {
            "all64_third_jet_components_derived_from_connection": True,
            "metric_and_connection_terms_both_retained": True,
            "all24_scalar_label_assignments_retained": True,
            "existing_original_Galileon_normalization_checked": True,
            "traced_connection_removed_only_after_real_TT": True,
        },
        "whole_third_jet_identity": "For J_r(u)=u^a...u^b nabla_(a...b)Phi, the linear normalized connection variation obeys delta J3=-(3p.u+k.u)*[2k.u*epsilon(p,u)-(k.p)*epsilon(u,u)], with the common Fourier factor stripped consistently with the flat vertex. The direct component expansion of partialGamma, Gamma*partialpartial and all three index positions is checked before symmetrization for every triple of spacetime indices.",
        "whole_rank_two_connection": "For a Hessian on field i whose two slots contract with momenta p_j,p_k, the connection term is-[a_j H_ik+a_k H_ij-a_i H_jk]. A self-trace instead has the contracted scalar connection, which is zero only for the specified physical TT graviton.",
        "whole_triple_edge_contact": "For Gij^3, delta=-6Gij^2 Hij -(3Gij+a_j)(2a_j Hij-a_i Hjj) -(3Gij+a_i)(2a_i Hij-a_j Hii). Omitting the k-dependent connection piece would change the six-derivative curvature matching.",
    }
