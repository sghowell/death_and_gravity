"""All Bose-symmetric explicit-curvature contacts through six derivatives."""

from functools import cache
from itertools import combinations, combinations_with_replacement, permutations, product

import sympy as s
from p8_vacuum_affine_radiative_curvature_matching import contact

a0 = s.symbols("a0:3")
a = (*a0, -sum(a0))
h0 = s.symbols("h00 h01 h02 h11 h12 h22")
h = s.zeros(4)
for (i, j), value in zip(combinations_with_replacement(range(3), 2), h0):
    h[i, j] = h[j, i] = value
for i in range(3):
    h[i, 3] = h[3, i] = -sum(h[i, j] for j in range(3))
h[3, 3] = sum(h[i, j] for i in range(3) for j in range(3))
variables = (*a0, *h0)
pairs = tuple(combinations(range(4), 2))
perms = tuple(permutations(range(4)))


@cache
def q(i, j, k, l):
    return s.expand(
        a[i] * a[k] * h[j, l]
        - a[i] * a[l] * h[j, k]
        - a[j] * a[k] * h[i, l]
        + a[j] * a[l] * h[i, k]
    )


canonical = s.expand(sum(q(i, j, i, j) for i, j in pairs))


def predicted(i, j, k, l):
    if i == j or k == l:
        return 0
    if {i, j} == {k, l}:
        return 4 if (i, j) == (k, l) else -4
    common = {i, j}.intersection((k, l))
    if not common:
        return 0
    shared = next(iter(common))
    return -2 * (1 if i == shared else -1) * (1 if k == shared else -1)


@cache
def data():
    checks, gates = {}, {}
    averages, multiplicities = [], {}
    for word in product(range(4), repeat=4):
        full = s.expand(sum(q(*(perm[x] for x in word)) for perm in perms))
        expected = predicted(*word)
        checks["Bose24_word_" + "".join(map(str, word))] = s.expand(
            full - expected * canonical
        )
        multiplicities[expected] = multiplicities.get(expected, 0) + 1
        averages.append(full)
    monomials = sorted(
        set().union(*(set(s.Poly(expr, *variables).monoms()) for expr in averages))
    )
    coefficient_matrix = s.Matrix(
        [
            [s.Poly(expr, *variables).coeff_monomial(m) for m in monomials]
            for expr in averages
        ]
    )
    checks["all256_Bose_averages_rank_exactly_one"] = coefficient_matrix.rank() - 1
    checks["each_epsilon_row_annihilates_total_momentum"] = s.Matrix(
        [sum(h[i, j] for j in range(4)) for i in range(4)]
    )
    checks["soft_dot_total_momentum_null"] = sum(a)
    for i, j, k in product(range(4), repeat=3):
        checks[f"edge_conservation_{i}{j}{k}"] = s.expand(
            sum(q(i, l, j, k) for l in range(4))
        )
    for i, j, k, l in product(range(4), repeat=4):
        checks[f"linear_Riemann_Bianchi_{i}{j}{k}{l}"] = s.expand(
            q(i, j, k, l) + q(i, k, l, j) + q(i, l, j, k)
        )
    gates["canonical_not_zero_before_mass_shell_specialization"] = canonical != 0
    gates["Bose_projection_not_identity_on_a_single_monomial"] = (
        s.expand(q(0, 1, 0, 1) - canonical / 6) != 0
    )
    gates["disjoint_pair_unsymmetrized_control_nonzero"] = q(0, 1, 2, 3) != 0
    gates["shared_pair_omission_changes_projection"] = predicted(0, 1, 0, 2) == -2

    # Generic null TT frame: Lorentz-covariant vanishing of every trace and k slot.
    omega, eplus, ecross = s.symbols("omega eplus ecross", real=True)
    wave = s.Matrix([omega, 0, 0, omega])
    eps = s.zeros(4)
    eps[1, 1], eps[2, 2], eps[1, 2], eps[2, 1] = eplus, -eplus, ecross, ecross
    curvature = contact.linear_curvature(wave, eps)
    for slot in range(4):
        other = tuple(i for i in range(4) if i != slot)
        for labels in product(range(4), repeat=3):
            values = dict(zip(other, labels))
            check = s.S.Zero
            for mu in range(4):
                indices = tuple(mu if j == slot else values[j] for j in range(4))
                check += wave[mu] * curvature[indices]
            checks["null_TT_all_k_slots_" + str(slot) + "".join(map(str, labels))] = (
                s.expand(check)
            )
    for u, v in combinations(range(4), 2):
        other = tuple(i for i in range(4) if i not in (u, v))
        for labels in product(range(4), repeat=2):
            values = dict(zip(other, labels))
            check = s.S.Zero
            for mu in range(4):
                indices = tuple(mu if j in (u, v) else values[j] for j in range(4))
                check += contact.ETA[mu, mu] * curvature[indices]
            checks[
                "null_TT_all_metric_traces_"
                + str(u)
                + str(v)
                + "".join(map(str, labels))
            ] = s.expand(check)

    gates["exact_symbolic_upper_bound_and_nonzero_canonical_polynomial"] = (
        coefficient_matrix.rank() == 1 and canonical != 0
    )
    return {
        "checks": checks,
        "gates": gates,
        "whole_generic_canonical_polynomial": canonical,
        "whole_Bose_coefficient_multiplicities": multiplicities,
        "whole_symbolic_basis_rank": coefficient_matrix.rank(),
        "whole_basis_derivation": "Linear Ricci and its derivatives vanish in real null TT; two curvatures start at two gravitons. A single Weyl tensor needs four scalar momenta because every trace and k-slot vanishes. At most six total derivatives leaves exactly those four slots, so every term is epsilon(v_ij,v_kl), followed by full S4 Bose symmetrization. All256 words reduce to the same T=sum_(i<j)epsilon(v_ij,v_ij). Conservation gives sum a=0 and sum_j H_ij=0; no massless scalar limit is used.",
        "whole_representative": "hbar*chi*Phi^2*C_abcd*(nabla^a nabla^c Phi)*(nabla^b nabla^d Phi)/4 gives chi*T/sqrt(kappa), with the unchanged physical metric. The basis has dimension1 at six derivatives and0 below within the stated explicit-curvature parity-even four-identical-scalar class.",
        "whole_scope": "Only the one-null-physical-graviton contact of explicit-curvature parity-even local operators through six derivatives. Not higher derivative orders, parity-odd terms, other fields, off-shell metric response, or equality of curved actions. No finite chi is assigned.",
    }
