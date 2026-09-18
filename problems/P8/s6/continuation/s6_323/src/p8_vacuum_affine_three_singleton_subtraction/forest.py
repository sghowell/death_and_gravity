"""Original3767 external-chain forest and exact eight-class cumulant source."""

from functools import cache
from itertools import product

import sympy as s
from p8_vacuum_affine_all_multiplicity_tree_source import trees as e
from p8_vacuum_affine_complete_two_graviton_tree import trees as lower

from . import source


class HardCore(e.TreeEngine):
    @lower.instance_cache
    def current(self, mask, kind):
        if mask.bit_count() > 1 and all(
            self.legs[i][0] == "h" for i in range(len(self.legs)) if mask >> i & 1
        ):
            return (lower.ZERO if kind == "h" else s.S.Zero), 0
        return super().current(mask, kind)


class CentralCore(HardCore):
    @lower.instance_cache
    def current(self, mask, kind):
        if kind == "phi" and mask.bit_count() > 1:
            return s.S.Zero, 0
        return super().current(mask, kind)


@cache
def data():
    ps, hs = e.configuration()
    hs = tuple((q, s.ImmutableMatrix(A / 2)) for q, A in hs)
    assert all(sum(abs(x) ** 2 for x in A) <= 1 for q, A in hs)
    pars = source.original_parameters()
    assert (
        pars["heavy"] == s.Rational(10) ** 200 / 512 + 2
        and pars["kappa"] == s.Integer(10) ** 800
    )
    assert (
        sum(ps, lower.VECTOR_ZERO) + sum((q for q, A in hs), lower.VECTOR_ZERO)
        == lower.VECTOR_ZERO
    )

    @cache
    def momentum(mask):
        return s.ImmutableMatrix(
            sum((hs[i][0] for i in range(3) if mask >> i & 1), lower.VECTOR_ZERO)
        )

    @cache
    def line(leg, mask):
        if not mask:
            return s.S.One, 1
        p = ps[leg]
        Q = momentum(mask)
        D = (p + Q).dot(lower.ETA * (p + Q)) - 1
        assert D != 0
        value = 0
        count = 0
        block = mask
        while block:
            rest = mask ^ block
            J, n = line(leg, rest)
            fields = tuple(hs[i][1] for i in range(3) if block >> i & 1)
            V = lower.scalar_vertex(
                s.ImmutableMatrix(p + momentum(rest)),
                s.ImmutableMatrix(-p - Q),
                fields,
                1,
            )
            value -= V * J / (D * pars["kappa"] ** s.Rational(len(fields), 2))
            count += n
            block = (block - 1) & mask
        return s.factor(value), count

    assert all(line(leg, 7)[1] == 13 for leg in range(4))
    reference = HardCore(
        [("phi", p, 1) for p in ps[1:]] + [("h", q, A) for q, A in hs], **pars
    )
    actual, actual_count = reference.amplitude()
    assert actual_count == 3767
    total = 0
    count = 0
    kernels = {}
    direct_by_central = {i: s.S.Zero for i in range(4)}
    central_counts = (7, 19, 67, 307)
    by_central = {i: 0 for i in range(4)}
    for index, assignment in enumerate(product(range(5), repeat=3), 1):
        masks = tuple(
            sum(1 << i for i, owner in enumerate(assignment) if owner == leg)
            for leg in range(5)
        )
        assert sum(masks) == 7 and all(
            not (masks[i] & masks[j]) for i in range(5) for j in range(i)
        )
        shifted = tuple(
            s.ImmutableMatrix(p + momentum(mask)) for p, mask in zip(ps, masks[:4])
        )
        central_rays = tuple(hs[i] for i in range(3) if masks[4] >> i & 1)
        assert (
            sum(shifted, lower.VECTOR_ZERO)
            + sum((q for q, A in central_rays), lower.VECTOR_ZERO)
            == lower.VECTOR_ZERO
        )
        kernel = CentralCore(
            [("phi", p, 1) for p in shifted[1:]]
            + [("h", q, A) for q, A in central_rays],
            **pars,
        )
        coefficient, central_count = kernel.amplitude()
        kernels[masks] = coefficient
        nrays = len(central_rays)
        assert central_count == central_counts[nrays]
        chains = tuple(line(leg, mask) for leg, mask in enumerate(masks[:4]))
        nc = central_count * s.prod(item[1] for item in chains)
        count += nc
        by_central[nrays] += nc
        total += coefficient * s.prod(item[0] for item in chains)
        direct_by_central[nrays] += coefficient * s.prod(item[0] for item in chains)
    assert count == 3767 and by_central == {0: 1288, 1: 1368, 2: 804, 3: 307}
    assert s.factor(total - actual) == 0
    # Separate EGF and chain count oracle; this does not replace amplitude equality.
    x = s.Symbol("x")
    chain = 1 / (2 - s.exp(x))
    central_egf = s.exp(x) + 6 * s.exp(2 * x) / (2 - s.exp(x))
    for N in range(4):
        assert s.diff(chain, x, N).subs(x, 0) == (1, 1, 3, 13)[N]
        assert s.diff(central_egf, x, N).subs(x, 0) == central_counts[N]
    assert s.diff(central_egf * chain**4, x, 3).subs(x, 0) == 3767

    weights = tuple(q[0] for q, A in hs)
    rays = tuple(q / w for (q, A), w in zip(hs, weights))
    fields = tuple(A for q, A in hs)
    dot = lambda p, q: (p.T * lower.ETA * q)[0]
    contract = lambda p, A, q: (p.T * A * q)[0]

    @cache
    def S(leg, i):
        return s.factor(weights[i] * s.sqrt(pars["kappa"]) * line(leg, 1 << i)[0])

    @cache
    def F(leg, mask):
        count = mask.bit_count()
        return s.factor(
            s.prod(weights[i] for i in range(3) if mask >> i & 1)
            * pars["kappa"] ** s.Rational(count, 2)
            * line(leg, mask)[0]
        )

    @cache
    def T(leg, i, j):
        return s.factor(F(leg, (1 << i) | (1 << j)) - S(leg, i) * S(leg, j))

    @cache
    def U(leg):
        return s.factor(
            F(leg, 7)
            - sum(F(leg, 7 ^ (1 << i)) * S(leg, i) for i in range(3))
            + 2 * s.prod(S(leg, i) for i in range(3))
        )

    def leading_at(p, i, j):
        A, B = fields[i], fields[j]
        ni, nj = rays[i], rays[j]
        di, dj = 2 * dot(p, ni), 2 * dot(p, nj)
        Li, Lj = 2 * contract(p, A, p), 2 * contract(p, B, p)
        mij, mji = 4 * contract(p, A, nj), 4 * contract(p, B, ni)
        z = 2 * dot(ni, nj)
        M = lower.ETA * lower.density_inverse_coefficient((A, B)) * lower.ETA
        base = lower.determinant_coefficient((A, B)) - contract(p, M, p)
        R0 = -Li * Lj * z + Li * mji * dj + Lj * mij * di + base * di * dj
        a, b = weights[i], weights[j]
        return s.factor(a * b * R0 / (di * dj * (a * di + b * dj)))

    @cache
    def L(leg, i, j):
        return leading_at(ps[leg], i, j)

    for leg in range(4):
        for i in range(3):
            for j in range(i + 1, 3):
                assert s.factor(leading_at(-ps[leg], i, j) + L(leg, i, j)) == 0

    def masks_for(owners):
        return tuple(
            sum(1 << i for i, owner in enumerate(owners) if owner == leg)
            for leg in range(5)
        )

    def H(owners):
        masks = masks_for(owners)
        return pars["kappa"] ** s.Rational(masks[4].bit_count(), 2) * kernels[masks]

    parts = {
        name: s.S.Zero
        for name in (
            "central3",
            "central2",
            "central1_pair",
            "central1_leading",
            "external0_connected3",
            "external0_pair_remainder",
            "external0_pair_leading",
            "external0_leading3",
        )
    }
    parts["central3"] = s.prod(weights) * H((4, 4, 4))
    for singleton in range(3):
        i, j = [k for k in range(3) if k != singleton]
        for leg in range(4):
            owners = [4, 4, 4]
            owners[singleton] = leg
            parts["central2"] += weights[i] * weights[j] * S(leg, singleton) * H(owners)
    for central in range(3):
        i, j = [k for k in range(3) if k != central]
        for left, right in product(range(4), repeat=2):
            owners = [0, 0, 0]
            owners[central] = 4
            owners[i] = left
            owners[j] = right
            parts["central1_leading"] += (
                weights[central] * S(left, i) * S(right, j) * H(owners)
            )
        for leg in range(4):
            owners = [leg, leg, leg]
            owners[central] = 4
            parts["central1_pair"] += weights[central] * T(leg, i, j) * H(owners)
    for owners in product(range(4), repeat=3):
        parts["external0_leading3"] += s.prod(S(owners[i], i) for i in range(3)) * H(
            owners
        )
    for singleton in range(3):
        i, j = [k for k in range(3) if k != singleton]
        for pairleg, singleleg in product(range(4), repeat=2):
            owners = [pairleg, pairleg, pairleg]
            owners[singleton] = singleleg
            kernel = H(owners) * S(singleleg, singleton)
            parts["external0_pair_leading"] += L(pairleg, i, j) * kernel
            parts["external0_pair_remainder"] += (
                T(pairleg, i, j) - L(pairleg, i, j)
            ) * kernel
    for leg in range(4):
        parts["external0_connected3"] += U(leg) * H((leg, leg, leg))
    normalization = s.prod(weights) * pars["kappa"] ** s.Rational(3, 2)
    for nr, names in (
        (3, ("central3",)),
        (2, ("central2",)),
        (1, ("central1_pair", "central1_leading")),
        (
            0,
            (
                "external0_connected3",
                "external0_pair_remainder",
                "external0_pair_leading",
                "external0_leading3",
            ),
        ),
    ):
        assert (
            s.factor(
                sum(parts[name] for name in names)
                - normalization * direct_by_central[nr]
            )
            == 0
        )
    assert s.factor(sum(parts.values()) - normalization * actual) == 0

    checks = {
        "independent_original_root_inventory": s.Integer(actual_count - 3767),
        "all125_forest_assignments": s.Integer(len(kernels) - 125),
        "independent_forest_inventory": s.Integer(count - 3767),
        "whole_forest_equals_original_root": s.factor(total - actual),
        "whole_eight_class_sum_equals_original_root": s.factor(
            sum(parts.values()) - normalization * actual
        ),
    }
    for nr, names in (
        (3, ("central3",)),
        (2, ("central2",)),
        (1, ("central1_pair", "central1_leading")),
        (
            0,
            (
                "external0_connected3",
                "external0_pair_remainder",
                "external0_pair_leading",
                "external0_leading3",
            ),
        ),
    ):
        checks[f"central_group{nr}_exact_cumulant_identity"] = s.factor(
            sum(parts[name] for name in names) - normalization * direct_by_central[nr]
        )
    for leg in range(4):
        for i in range(3):
            for j in range(i + 1, 3):
                checks[f"line{leg}_pair{i}{j}_leading_odd"] = s.factor(
                    leading_at(-ps[leg], i, j) + L(leg, i, j)
                )
    for N in range(4):
        checks[f"central_EGF_coefficient{N}"] = (
            s.diff(central_egf, x, N).subs(x, 0) - central_counts[N]
        )
        checks[f"scalar_chain_EGF_coefficient{N}"] = (
            s.diff(chain, x, N).subs(x, 0) - (1, 1, 3, 13)[N]
        )
    for i, p in enumerate(ps):
        checks[f"external_mass_shell{i}"] = s.factor((p.T * lower.ETA * p)[0] - 1)
    for i, (q, A) in enumerate(hs):
        checks[f"external_null_shell{i}"] = s.factor((q.T * lower.ETA * q)[0])
        checks[f"physical_TT_transversality{i}"] = A * q
        checks[f"physical_TT_trace{i}"] = s.trace(lower.ETA * A)
    return {
        "checks": checks,
        "gates": {
            "all_four_central_class_inventories": by_central
            == {0: 1288, 1: 1368, 2: 804, 3: 307},
            "all_eight_cumulant_contributions_nonzero": sum(
                s.factor(value) != 0 for value in parts.values()
            )
            == 8,
            "original_couplings_not_diagnostic_parameters": pars
            == source.original_parameters(),
            "unit_bound_physical_TT_fields": all(
                sum(abs(value) ** 2 for value in A) <= 1 for q, A in hs
            ),
            "central_kernel_has_no_internal_Phi_or_pure_soft_propagators": True,
            "source_identity_not_sampled_uniform_estimate": True,
        },
        "whole_central_graph_counts": central_counts,
        "whole_forest_graph_counts_by_central_leaves": by_central,
        "whole_soft_canonical_power_stripped_parts": parts,
        "whole_points": ps,
        "whole_radiation": hs,
        "whole_source_identity": "The sum over125 assignments to four scalar chains and a central source equals the independent original3767 tree. Splitting normalized line currents into singles, connected pairs, their leading/remainder split and connected triples gives eight exactly equal algebraic contributions. These include signed lower products, not a partition into positive physical observables.",
    }
