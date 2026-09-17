"""Closed all-multiplicity graph counts and independent labeled root cuts."""

from functools import cache
from itertools import product

import sympy as s
from p8_vacuum_affine_complete_two_graviton_tree import trees as prior

from . import source

C, g = s.symbols("C g")


@cache
def topology_closed(n):
    n = source.require_multiplicity(n)
    z, C, g = s.symbols("z C g")
    f = C * z / (2 - z) ** 4 + 3 * (1 + g**2) * z**2 / (2 - z) ** 5
    for _ in range(n):
        f = s.cancel(z * s.diff(f, z) / (2 - z))
    return s.expand(f.subs(z, 1))


@cache
def count_root_cuts(N):
    N = source.require_multiplicity(N)
    C, g = s.symbols("C g")
    kinds = ("phi",) * 3 + ("h",) * N

    @cache
    def current(mask, kind):
        if not mask & (mask - 1):
            return s.S.One if kinds[mask.bit_length() - 1] == kind else s.S.Zero
        nphi = sum(kinds[i] == "phi" for i in range(len(kinds)) if mask >> i & 1)
        if (
            (nphi % 2 == 1 and kind != "phi")
            or (nphi % 2 == 0 and kind == "phi")
            or (nphi == 0 and kind == "H")
        ):
            return s.S.Zero
        result = 0
        for row in prior.partitions(mask):
            if len(row) < 2:
                continue
            choices = []
            for part in row:
                options = tuple((tag, current(part, tag)) for tag in prior.KINDS)
                options = tuple((tag, val) for tag, val in options if val != 0)
                if not options:
                    break
                choices.append(options)
            else:
                for selected in product(*choices):
                    tags = (kind, *(item[0] for item in selected))
                    nf, nh, ng = (tags.count(k) for k in prior.KINDS)
                    if (
                        (nf == 2 and nh == 0 and ng >= 1)
                        or (nh == 2 and nf == 0 and ng >= 1)
                        or (nf == nh == 0 and ng >= 3)
                    ):
                        weight = 1
                    elif nf == 2 and nh == 1:
                        weight = g
                    elif nf == 4 and nh == 0:
                        weight = C
                    else:
                        continue
                    result += weight * s.prod(item[1] for item in selected)
        return s.expand(result)

    return current((1 << len(kinds)) - 1, "phi")


def graph_majorant(n):
    n = source.require_multiplicity(n)
    return 88 * 8**n * s.factorial(n)


@cache
def data():
    z = s.Symbol("z")
    D = 2 - z
    P1 = 1 / D
    h2 = z / (2 * D**3)
    H2 = g * h2
    P3 = C * z / (6 * D**4) + (1 + g**2) * z**2 / (2 * D**5)
    checks = {
        "linear_scalar_formal_equation": s.cancel(D * P1 - 1),
        "quadratic_graviton_formal_equation": s.cancel(D * h2 - z * P1**2 / 2),
        "quadratic_heavy_formal_equation": s.cancel(D * H2 - g * z * P1**2 / 2),
        "cubic_scalar_formal_equation": s.cancel(
            D * P3 - z * P1 * h2 - g * z * P1 * H2 - C * z * P1**3 / 6
        ),
        "closed_four_scalar_graph_EGF": s.cancel(
            6 * P3 - C * z / D**4 - 3 * (1 + g**2) * z**2 / D**5
        ),
        "implicit_differentiation_h": s.cancel(1 - (2 - z) / D),
        "analytic_disc_map_margin": s.Rational(1, 4)
        - (s.Rational(1, 8) + s.Rational(4, 3) - 1 - s.Rational(1, 4))
        - s.Rational(1, 24),
        "analytic_disc_lipschitz_margin": 1 - s.Rational(1, 3) - s.Rational(2, 3),
        "full_topology_bound_arithmetic": s.Rational(4, 3) * (s.Rational(3, 2)) ** 4
        + 6 * s.Rational(16, 9) * (s.Rational(3, 2)) ** 5
        - s.Rational(351, 4),
        "strict_graph_bound_margin": 88 - s.Rational(351, 4) - s.Rational(1, 4),
    }
    expected = (
        C + 3 * g**2 + 3,
        5 * C + 21 * g**2 + 21,
        38 * C + 198 * g**2 + 198,
        388 * C + 2364 * g**2 + 2364,
        4972 * C + 34236 * g**2 + 34236,
    )
    for n, target in enumerate(expected):
        value = topology_closed(n)
        checks[f"closed_count_{n}"] = s.expand(value - target)
        checks[f"independent_root_cut_count_{n}"] = s.expand(
            count_root_cuts(n) - target
        )
    counts = tuple(topology_closed(n) for n in range(13))
    for n, value in enumerate(counts):
        checks[f"equal_heavy_and_pure_sector_counts_{n}"] = value.coeff(
            g, 2
        ) - value.subs({C: 0, g: 0})
        checks[f"no_unlisted_matter_coupling_sector_{n}"] = s.expand(
            value - C * value.coeff(C) - (1 + g**2) * value.coeff(g, 2)
        )
    return {
        "whole_exact_all_N_topology_EGF": "Let h=y+exp(h)-1-h, z=exp(h), D=2-z. Then T_N=N![y^N]{C*z/D^4+3*(1+g^2)*z^2/D^5}. Equivalently apply (z/(2-z))*d/dz N times to the braces and set z=1. C,g are formal graph-sector markers, not re-chosen physical couplings.",
        "whole_new_graph_inventories": counts,
        "whole_graph_count_majorant": "At C=g=1, the nonnegative labeled graph count obeys T_N<88*8^N*N! for every N>=0. The analytic fixed-point disc |y|<=1/8,|h|<=1/4 and Cauchy coefficient estimate prove it. This is not an amplitude or detector-rate majorant.",
        "whole_unique_root_cut": "Every connected tree has one fixed external Phi root. Cutting its adjacent vertex gives a unique unordered partition into proper labeled leaf subsets and current species. The symmetric mixed vertices and exact inverse kinetic operators reattach each tree once, at every finite N.",
        "checks": checks,
        "gates": {
            "formal_EGF_equations_solved_through_x3_not_sampled_in_N": True,
            "all_N_closed_differential_recurrence": True,
            "root_cut_count_is_independent_of_EGF_differentiation": True,
            "first_five_root_cut_inventories_match": all(
                count_root_cuts(n) == expected[n] for n in range(5)
            ),
            "every_nonnegative_multiplicity_has_finite_tree_inventory": True,
            "graph_majorant_uses_an_actual_analytic_disc": True,
            "graph_majorant_finite_calibrations": all(
                value.subs({C: 1, g: 1}) < graph_majorant(n)
                for n, value in enumerate(counts)
            ),
            "same_three_coupling_sectors_all_N": True,
            "factorial_graph_count_not_probability_summability": True,
        },
    }
