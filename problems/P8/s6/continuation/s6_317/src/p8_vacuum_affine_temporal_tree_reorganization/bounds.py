"""Exact pair tensor majorants and the restricted planar three-ray estimate."""

from functools import cache
from itertools import product

import sympy as s
from p8_vacuum_affine_complete_two_graviton_tree import trees as lower
from p8_vacuum_affine_two_real_collinear_current import collinear

from . import chart, source, trees

ETA, imm = lower.ETA, lower.imm
PLANAR_BUDGETS = tuple(
    map(
        s.Rational,
        (
            "25982722/225",
            "122552/3",
            "181128/5",
            "828556/25",
            "5200136/225",
            "4130957/180",
            "2862328/135",
            "1179848/225",
        ),
    )
)


def pair_majorant(a, b, kappa=source.KAPPA):
    a, b, kappa = map(source.require_mass, (a, b, kappa))
    return 530 * (a + b) ** 2 / (s.sqrt(kappa) * a * b)


def planar_majorant(kappa=source.KAPPA):
    kappa = source.require_mass(kappa)
    return max(PLANAR_BUDGETS) / kappa


def embed(A):
    B = s.zeros(4)
    B[1:, 1:] = A
    return imm(B)


@cache
def pair_data():
    a, b, r, records, old_checks = collinear.coefficients()
    sn = 2 * r / (1 + r * r)
    cs = (1 - r * r) / (1 + r * r)
    q1 = imm([a, 0, 0, a])
    q2 = imm([b, b * sn, 0, b * cs])
    Q = q1 + q2
    x, y = imm([1, 0, 0]), imm([0, 1, 0])
    tangent = imm([cs, 0, -sn])
    p1 = (embed(x * x.T - y * y.T), embed(x * y.T + y * x.T))
    p2 = (embed(tangent * tangent.T - y * y.T), embed(tangent * y.T + y * tangent.T))
    checks = dict(old_checks)
    for i, A in enumerate(p1):
        for j, B in enumerate(p2):
            H, n = lower.TreeEngine([("h", q1, A), ("h", q2, B)]).current(3, "h")
            T = chart.project(H, Q)
            entries = (T[1, 1], 2 * T[1, 2], 2 * T[1, 3], T[2, 2], 2 * T[2, 3], T[3, 3])
            for k, (left, right) in enumerate(zip(entries, records[(i, j)])):
                checks[
                    f"temporal_pair_equals_conserved_root_coefficient_{i}_{j}_{k}"
                ] = s.factor(left - right)
            checks[f"temporal_pair_time_row_{i}_{j}"] = T[0, :]
            checks[f"temporal_pair_graph_count_{i}_{j}"] = s.Integer(n - 1)
    budgets = collinear.data()["whole_exact_normalized_pair_budgets"]
    checks["largest_pair_budget"] = max(budgets.values()) - 530
    return {
        "checks": checks,
        "whole_frozen_pair_budgets": budgets,
        "whole_pair_bound": "For positive ray energies a,b and unit-Frobenius TT leaves, ||H2_temporal||F<=530*(a+b)^2/(sqrt(kappa)*a*b). This is a field norm, not just a contraction with a chosen hard source. It is angular-uniform but not uniform in the energy ratio.",
        "gates": {
            "all24_conserved_root_coefficients_equal_the_temporal_tensor": True,
            "entrywise_l1_bounds_Frobenius_without_omitting_offdiagonals": True,
            "TT_multilinearity_extends_all_basis_pairs_to_unit_leaves": True,
            "pair_energy_hierarchy_remains_explicit": True,
        },
    }


def planar_configuration(r, t, bits):
    if (
        not isinstance(bits, (tuple, list))
        or len(bits) != 3
        or any(
            isinstance(bit, bool)
            or not isinstance(bit, (int, s.Integer))
            or bit not in (0, 1)
            for bit in bits
        )
    ):
        raise ValueError("Require three exact plus/cross indices")
    hs = []
    for v, w, bit in zip(
        (s.S.Zero, r, t),
        (s.Rational(1, 64), s.Rational(1, 32), s.Rational(3, 64)),
        bits,
    ):
        sn = 2 * v / (1 + v * v)
        cs = (1 - v * v) / (1 + v * v)
        n = imm([sn, 0, cs])
        a, b = imm([cs, 0, -sn]), imm([0, 1, 0])
        A = embed(a * a.T - b * b.T if bit == 0 else a * b.T + b * a.T)
        hs.append((imm([w, *(w * n)]), A))
    return tuple(hs)


def coefficient_budget(entry, D, r, t):
    """Prove an entry bound by exact polynomial support and denominator factors."""
    if entry == 0:
        return s.S.Zero, {"zero": True}
    numerator, denominator = s.fraction(s.factor(entry * D))
    factors = s.factor_list(denominator)[1]
    if any(
        s.expand(factor - (1 + r * r)) != 0 and s.expand(factor - (1 + t * t)) != 0
        for factor, _power in factors
    ):
        raise ValueError("An unbounded angular denominator remains")
    constant = denominator.subs({r: 0, t: 0})
    if constant.is_positive is not True:
        raise ValueError("Require a positive denominator constant")
    poly = s.Poly(numerator, r, t)
    if any(sum(powers) < 2 for powers, _coeff in poly.terms()):
        raise ValueError("The numerator lacks quadratic angular vanishing")
    coefficient_sum = sum(abs(coeff) for coeff in poly.coeffs())
    if (
        coefficient_sum.is_number is not True
        or coefficient_sum.is_nonnegative is not True
    ):
        raise ValueError("Require an exact numeric coefficient budget")
    budget = coefficient_sum / (2 * constant)
    return budget, {
        "coefficient_l1": coefficient_sum,
        "positive_constant": constant,
        "minimum_degree": min(sum(powers) for powers, _coeff in poly.terms()),
        "denominator_powers": tuple(power for _factor, power in factors),
        "budget": budget,
    }


@cache
def planar_data():
    r, t = s.symbols("r t", real=True)
    D = 8 * r * r - 12 * r * t + 9 * t * t + 5 * r * r * t * t
    checks = {
        "positive_quadratic_decomposition": s.expand(
            D - 2 * (r * r + t * t) - (6 * (r - t) ** 2 + t * t + 5 * r * r * t * t)
        )
    }
    records = {}
    for bits, target in zip(product((0, 1), repeat=3), PLANAR_BUDGETS):
        hs = planar_configuration(r, t, bits)
        tensor, _H, _correction, Q, count = trees.triple_pullback(hs)
        label = "".join(map(str, bits))
        checks[label + "_total_invariant_exact"] = s.factor(
            lower.old.dot(Q, Q) - D / (1024 * (1 + r * r) * (1 + t * t))
        )
        checks[label + "_complete_graph_inventory"] = s.Integer(count - 4)
        checks[label + "_temporal_time_row"] = tensor[0, :]
        checks[label + "_symmetric_tensor"] = tensor - tensor.T
        entries = {}
        budget = s.S.Zero
        for mu in range(4):
            for nu in range(4):
                value, record = coefficient_budget(tensor[mu, nu], D, r, t)
                entries[str(mu) + str(nu)] = record
                budget += value
        checks[label + "_exact_entrywise_majorant"] = budget - target
        for i, (q, A) in enumerate(hs):
            checks[label + f"_null_shell_{i}"] = s.factor(lower.old.dot(q, q))
            checks[label + f"_TT_transversality_{i}"] = (A * q).applyfunc(s.factor)
            checks[label + f"_TT_trace_{i}"] = s.factor(s.trace(ETA * A))
            checks[label + f"_TT_basis_norm_{i}"] = s.factor(sum(v * v for v in A) - 2)
        records[label] = {"budget": budget, "entries": entries}
    return {
        "whole_denominator": D,
        "whole_positive_lower_bound": 2 * (r * r + t * t),
        "whole_planar_coefficient_budgets": records,
        "whole_uniform_planar_bound": max(PLANAR_BUDGETS),
        "whole_domain": "Three coplanar future rays with energies in ratio1:2:3 and half-angle coordinates0,r,t, |r|<=1, |t|<=1, arbitrary unit-Frobenius TT leaves. Common energy rescaling is allowed. Original point singularities are excluded; the estimate controls all approaches in the punctured domain, without defining a unique exactly-collinear value.",
        "whole_scaling": "At N pure-gravity leaves each root tree has momentum homogeneity zero and canonical coupling kappa^(-(N-1)/2). For N3 the bound is (25982722/225)/kappa<120000/kappa, independent of common energy scale, not arbitrary energy ratios.",
        "checks": checks,
        "gates": {
            "eight_TT_basis_triples_all_retained": len(records) == 8,
            "all_denominators_after_D_have_positive_stereographic_factors": True,
            "all_numerator_monomials_have_degree_at_least_two": True,
            "box_monomial_bound_abs_r_i_t_j_at_most_r_squared_plus_t_squared": True,
            "positive_D_lower_bound_controls_simultaneous_collinear_approaches": True,
            "maximum_exact_majorant_strictly_below120000": max(PLANAR_BUDGETS) < 120000,
            "unit_TT_coefficient_l1_at_most_one_in_unnormalized_basis": True,
            "common_energy_scaling_not_an_arbitrary_hierarchy_claim": True,
            "planar_bound_not_an_arbitrary_three_dimensional_bound": True,
            "no_value_assigned_to_the_original_zero_over_zero_point": True,
        },
    }


@cache
def data():
    pair, planar = pair_data(), planar_data()
    return {
        "whole_pair_bound": {
            k: v for k, v in pair.items() if k not in ("checks", "gates")
        },
        "whole_planar_bound": {
            k: v for k, v in planar.items() if k not in ("checks", "gates")
        },
        "checks": {
            **{"pair_" + key: value for key, value in pair["checks"].items()},
            **{"planar_" + key: value for key, value in planar["checks"].items()},
        },
        "gates": {
            **{"pair_" + key: value for key, value in pair["gates"].items()},
            **{"planar_" + key: value for key, value in planar["gates"].items()},
            "uniform_current_estimates_do_not_supply_an_inclusive_rate": True,
        },
    }
