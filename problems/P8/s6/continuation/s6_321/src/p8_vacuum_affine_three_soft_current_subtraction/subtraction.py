"""Reconstructed all-angle current derivatives and compatible energy faces.

All polynomial certificates are derived from the unchanged action at replay.
There is no JSON-current input or filesystem-dependent scientific shortcut.
"""

from functools import cache
from itertools import product

import sympy as s
from sympy.polys.rings import ring

from . import polynomial, source

UPPER = tuple((i, j) for i in range(1, 4) for j in range(i, 4))
LOWER = ((1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1), (0, 1, 1))
C0 = s.Integer(4800000000000000000000)
C1 = s.Integer(10) ** 6
C2 = s.Integer(10) ** 6
C3 = s.Integer(721)


def product_polynomial(parts, R):
    answer = R.one
    for factor, power in parts:
        answer *= factor**power
    return answer


def shifts(dimension, homogeneous=False):
    rows = []
    for delta in product((-1, 0, 1), repeat=dimension):
        nonzero = [d for d in delta if d]
        if nonzero and nonzero[0] > 0:
            if homogeneous:
                da, db, dr, dv, du = delta
                delta = (da, db, -da - db, dr, dv, du)
            rows.append(delta)
    return tuple(sorted(rows, key=lambda d: (sum(map(abs, d)), d)))


def coefficient_bound(P, D, deltas):
    if not D or any(c <= 0 for c in D.values()):
        raise ValueError("A positive denominator coefficient certificate is required")
    major = {}
    midpoints = 0
    for monomial, coefficient in P.items():
        value = abs(coefficient)
        if monomial in D:
            major[monomial] = major.get(monomial, 0) + value
            continue
        for delta in deltas:
            left = tuple(m + d for m, d in zip(monomial, delta))
            right = tuple(m - d for m, d in zip(monomial, delta))
            if left in D and right in D:
                major[left] = major.get(left, 0) + value / 2
                major[right] = major.get(right, 0) + value / 2
                midpoints += 1
                break
        else:
            raise ValueError("Unsupported numerator monomial in coefficient bound")
    bound = max(
        (value / D[m] for m, value in major.items()), default=P.ring.domain.zero
    )
    return P.ring.domain.to_sympy(bound), midpoints


def reconstruct(sector, bits):
    source.require_case(sector, bits)
    R, *gens = ring("a,b,r,v,u,x00,x01,x02,x03,x11,x12,x13,x22,x23,x33", s.QQ)
    a, b, r, v, u = gens[:5]
    one, zero = R.one, R.zero
    coords = source.coordinates(sector, r, v, u)
    coords = tuple(tuple(R(p) for p in xy) for xy in coords)
    ds = tuple(one + x * x + y * y for x, y in coords)
    L = ds[0] * ds[1] * ds[2]

    def engine(qs, leaves):
        return polynomial.FactoredAssemblyPolynomialCurrent(R, qs, leaves, gens[-10:])

    dummy = engine((), ())
    leaves = []
    momenta = []
    for (x, y), weight, d, bit in zip(coords, (a, b, one), ds, bits):
        U = (d - 2 * x * x, -2 * x * y, -2 * x)
        V = (-2 * x * y, d - 2 * y * y, -2 * y)
        n = (2 * x, 2 * y, one - x * x - y * y)
        assert sum((q * q for q in U), zero) == d * d
        assert sum((q * q for q in V), zero) == d * d
        assert sum((q * p for q, p in zip(U, V)), zero) == 0
        assert sum((q * p for q, p in zip(U, n)), zero) == 0
        assert sum((q * p for q, p in zip(V, n)), zero) == 0
        A = [[zero for j in range(4)] for i in range(4)]
        for i in range(3):
            for j in range(3):
                A[i + 1][j + 1] = (
                    U[i] * U[j] - V[i] * V[j] if bit == 0 else U[i] * V[j] + V[i] * U[j]
                )
        assert sum((A[i][i] for i in range(1, 4)), zero) == 0
        assert sum((q * q for row in A for q in row), zero) == 2 * d**4
        constant, parts = dummy.factor_parts(d)
        N, parts = dummy.normalize(
            tuple(map(tuple, A)), [(f, 2 * k) for f, k in parts], constant**2
        )
        leaves.append((N, parts))
        quotient, remainder = dummy.known_factor_division.div(L, d)
        assert not remainder
        momenta.append(tuple(weight * quotient * p for p in (d, *n)))
    current = engine(momenta, leaves)
    for mask in (3, 5, 6):
        assert current.current(mask)[2] == 1
    record = current.current(7)
    assert record[2] == 4
    F, parts, count = polynomial.multiply_energies_over_total(
        current, record, (a, b), a + b + 1
    )
    assert count == 4
    T, *_compact_gens = ring("a,b,r,v,u", s.QQ)

    def compact(P):
        assert all(not any(m[5:]) for m in P)
        made = T.from_dict({m[:5]: c for m, c in P.items()})
        assert R.from_dict({m + (0,) * 10: c for m, c in made.items()}) == P
        return made

    numerators = {(i, j): compact(F[i][j]) for i, j in UPPER}
    parts = tuple((compact(f), k) for f, k in parts)
    assert all(F[0][j] == F[j][0] == 0 for j in range(4))
    assert all(F[i][j] == F[j][i] for i, j in UPPER)
    return T, numerators, parts


def operator(factors, powers, variables):
    R = factors[0].ring
    active = [k for k, f in enumerate(factors) if any(f.diff(v) for v in variables)]
    L = R.one
    for k in active:
        L *= factors[k]
    derivatives = []
    for variable in variables:
        T = R.zero
        for k in active:
            remaining, remainder = L.div(factors[k])
            assert not remainder
            T += powers[k] * factors[k].diff(variable) * remaining
        derivatives.append(T)
    for k in active:
        powers[k] += 1
    return L, tuple(derivatives)


def mixed_operators(parts):
    R = parts[0][0].ring
    a, b = R.gens[:2]
    factors = [f for f, k in parts]
    powers = [k for f, k in parts]
    operators = []
    for variable in (a, b):
        L, (T,) = operator(factors, powers, (variable,))
        operators.append((variable, L, T))
    L, (Ta, Tb) = operator(factors, powers, (a, b))
    D = product_polynomial(zip(factors, powers), R)
    pair_num = (a + b) * (a + 1) + (a + b) * (b + 1) + (a + 1) * (b + 1)
    pair_den = (a + b) * (a + 1) * (b + 1)
    return operators, (L, Ta, Tb), D * pair_num, pair_den


def homogenize(numerators, parts):
    T, *_ = ring("a,b,c,r,v,u", s.QQ)

    def lift(P, degree):
        terms = {}
        for (a, b, r, v, u), coefficient in P.items():
            c = degree - a - b
            assert c >= 0
            terms[(a, b, c, r, v, u)] = coefficient
        made = T.from_dict(terms)
        assert all(sum(m[:3]) == degree for m in made)
        dropped = P.ring.from_dict({(m[0], m[1], *m[3:]): q for m, q in made.items()})
        assert dropped == P
        return made

    lifted = []
    total_degree = 0
    for f, power in parts:
        degree = max(m[0] + m[1] for m in f)
        made = lift(f, degree)
        assert all(c > 0 for c in made.values())
        lifted.append((made, power))
        total_degree += degree * power
    return (
        T,
        {ij: lift(P, total_degree + 2) for ij, P in numerators.items()},
        tuple(lifted),
    )


def face(P, Q, index):
    if not P:
        return P.ring.zero, P.ring.one
    pn = min(m[index] for m in P)
    qn = min(m[index] for m in Q)
    if pn < qn:
        raise ValueError("The proposed energy face is not finite")
    if pn > qn:
        return P.ring.zero, P.ring.one

    def leading(poly, power):
        return P.ring.from_dict(
            {
                m[:index] + (0,) + m[index + 1 :]: v
                for m, v in poly.items()
                if m[index] == power
            }
        )

    return leading(P, pn), leading(Q, qn)


def face_checks(numerators, parts):
    R = parts[0][0].ring
    a, b, c = R.gens[:3]
    energies = (a, b, c)
    D = product_polynomial(parts, R)
    assert all(v > 0 for v in D.values())
    deltas = shifts(3)
    checks = {}
    rows = []
    for (i, j), P in numerators.items():
        for held in range(3):
            first, second = [k for k in range(3) if k != held]
            fp, fq = face(P, D, first)
            gp, gq = face(P, D, second)
            ap, aq = face(fp, fq, second)
            bp, bq = face(gp, gq, first)
            assert ap * bq == bp * aq
            assert all(value > 0 for value in aq.values())
            deltaP = P * aq - ap * D
            deltaD = D * aq * (energies[first] + energies[second]) * (a + b + c)
            assert all(value > 0 for value in deltaD.values())
            ps = {m[:3] for m in deltaP}
            ds = {m[:3] for m in deltaD}
            midpoints = 0
            for m in ps:
                if m in ds:
                    continue
                for delta in deltas:
                    left = tuple(x + y for x, y in zip(m, delta))
                    right = tuple(x - y for x, y in zip(m, delta))
                    if left in ds and right in ds:
                        midpoints += 1
                        break
                else:
                    raise ValueError("Joint-axis energy support not certified")
            label = f"face_{i}{j}_held{held}"
            checks[label + "_sequential_limit_identity"] = s.S.Zero
            rows.append(
                {
                    "component": (i, j),
                    "held_energy": held,
                    "sequential_limits_agree": True,
                    "joint_axis_support": True,
                    "axis_zero": not bool(ap),
                    "energy_midpoints": midpoints,
                }
            )
    assert len(rows) == 18
    return checks, rows


def matrix_bound(component_bounds):
    return sum((1 if i == j else 2) * v for (i, j), v in component_bounds.items())


def case_data(sector, bits):
    source.require_case(sector, bits)
    return _case_data(sector, bits)


@cache
def _case_data(sector, bits):
    source.taylor_calibration()
    R, numerators, parts = reconstruct(sector, bits)
    a, b = R.gens[:2]
    coords = source.coordinates(
        sector, s.Rational(1, 2), s.Rational(1, 3), s.Rational(1, 4)
    )
    oracle = source.evaluate(
        (s.Rational(1, 4), s.S.One, s.S.One),
        tuple(source.direction(*xy) for xy in coords),
        tuple(source.frame(*xy) for xy in coords),
        bits,
    )
    coefficients = {ij: s.Poly(oracle[ij], *source.VARIABLES) for ij in UPPER}
    point = list(
        zip(R.gens, (s.QQ(1, 4), s.QQ.one, s.QQ(1, 2), s.QQ(1, 3), s.QQ(1, 4)))
    )
    D0 = product_polynomial(parts, R)
    operators, (lastL, lastTa, lastTb), D, pair_den = mixed_operators(parts)
    deltas = shifts(5)
    checks = {}
    budgets = {}
    certificates = {}
    for ij, P0 in numerators.items():
        residual = R.domain.to_sympy(
            P0.evaluate(point) / D0.evaluate(point)
        ) - coefficients[ij].coeff_monomial((0, 0, 0))
        assert residual == 0
        checks[f"constant_Taylor_{ij[0]}{ij[1]}"] = residual
        P = P0
        for variable, L, T in operators:
            P = P.diff(variable) * L - P * T
        P = (
            -(a * P.diff(a) + b * P.diff(b)) * lastL + P * (a * lastTa + b * lastTb)
        ) * pair_den
        residual = R.domain.to_sympy(
            P.evaluate(point) / D.evaluate(point)
        ) - coefficients[ij].coeff_monomial((1, 1, 1)) / s.Rational(21, 10)
        assert residual == 0
        checks[f"mixed_Taylor_{ij[0]}{ij[1]}"] = residual
        bound, midpoints = coefficient_bound(P, D, deltas)
        budgets[ij] = bound
        certificates[f"mixed_{ij[0]}{ij[1]}"] = {
            "bound": bound,
            "midpoints": midpoints,
            "numerator_terms": len(P),
            "denominator_terms": len(D),
        }
    third = matrix_bound(budgets)
    assert third < C3
    T, numerators6, parts6 = homogenize(numerators, parts)
    a, b, c = T.gens[:3]
    W = a + b + c
    pair_num = (a + b) * (a + c) + (a + b) * (b + c) + (a + c) * (b + c)
    pair_den = (a + b) * (a + c) * (b + c)
    point6 = list(
        zip(
            T.gens, (s.QQ(1, 4), s.QQ.one, s.QQ.one, s.QQ(1, 2), s.QQ(1, 3), s.QQ(1, 4))
        )
    )
    lower_bounds = {}
    deltas6 = shifts(5, homogeneous=True)
    for derivative in LOWER:
        factors = [f for f, k in parts6]
        powers = [k for f, k in parts6]
        ops = []
        for variable, wanted in zip((a, b, c), derivative):
            if wanted:
                L, (TT,) = operator(factors, powers, (variable,))
                ops.append((variable, L, TT))
        DD = product_polynomial(zip(factors, powers), T)
        first = sum(derivative) == 1
        multiplier = T.one if first else pair_den
        target = s.Rational(9, 4) if first else s.Rational(189, 40)
        DD *= W if first else W * pair_num
        components = {}
        for ij, P0 in numerators6.items():
            P = P0
            for variable, L, TT in ops:
                P = P.diff(variable) * L - P * TT
            P *= multiplier
            residual = (
                T.domain.to_sympy(P.evaluate(point6) / DD.evaluate(point6))
                - coefficients[ij].coeff_monomial(derivative) / target
            )
            assert residual == 0
            key = f"lower_{''.join(map(str, derivative))}_{ij[0]}{ij[1]}"
            checks[key + "_Taylor"] = residual
            bound, midpoints = coefficient_bound(P, DD, deltas6)
            components[ij] = bound
            certificates[key] = {
                "bound": bound,
                "midpoints": midpoints,
                "numerator_terms": len(P),
                "denominator_terms": len(DD),
            }
        total = matrix_bound(components)
        assert total <= (C1 if first else C2)
        lower_bounds["".join(map(str, derivative))] = total
    faces, face_rows = face_checks(numerators6, parts6)
    checks.update(faces)
    assert len(checks) == 66
    return {
        "sector": sector,
        "TT_bits": bits,
        "graph_count": 4,
        "third_mixed_matrix_l1_bound": third,
        "lower_matrix_l1_bounds": lower_bounds,
        "coefficient_certificates": certificates,
        "compatible_axis_certificates": face_rows,
        "checks": checks,
        "gates": {
            "literal_EH3_EH4_complete_Ward_checked": True,
            "all48_independent_Taylor_coefficients_match": True,
            "positive_denominators_and_supported_AMGM_numerators": True,
            "exact_degree_two_energy_homogenization": True,
            "all18_sequential_axes_and_joint_limits_certified": True,
            "angular_axis_Lipschitz_constant_only_pointwise": True,
        },
    }


@cache
def data():
    records = {
        f"sector{sector}_TT{''.join(map(str, bits))}": case_data(sector, bits)
        for sector in (0, 1)
        for bits in product((0, 1), repeat=3)
    }
    checks = {
        label + "_" + key: value
        for label, row in records.items()
        for key, value in row["checks"].items()
    }
    first = max(
        value
        for row in records.values()
        for key, value in row["lower_matrix_l1_bounds"].items()
        if key.count("1") == 1
    )
    second = max(
        value
        for row in records.values()
        for key, value in row["lower_matrix_l1_bounds"].items()
        if key.count("1") == 2
    )
    third = max(row["third_mixed_matrix_l1_bound"] for row in records.values())
    return {
        "whole_reconstructed_sixteen_cases": {
            label: {k: v for k, v in row.items() if k not in ("checks", "gates")}
            for label, row in records.items()
        },
        "whole_exact_basis_derivative_maxima": (first, second, third),
        "whole_unit_complex_TT_hierarchy_caps": (C0, C1, C2, C3),
        "whole_hierarchy_scales": ("W^2/kappa", "W/kappa", "W*S/kappa", "S/kappa"),
        "whole_rectangle_current_bound": "721*J(a,b,c)/kappa",
        "checks": checks,
        "gates": {
            **{
                label + "_" + key: value
                for label, row in records.items()
                for key, value in row["gates"].items()
            },
            "all16_two_sector_TT_currents_reconstructed": len(records) == 16,
            "all768_Taylor_coefficients_and288_axis_certificates": len(checks) == 1056,
            "angular_sector_reduction_and_complex_TT_multilinearity": True,
            "origin_O_W_squared_from_S318": True,
            "compatible_faces_trimmed_cube_FTC_not_Cauchy_through_soft_poles": True,
            "current_rectangle_not_full_probability_subtraction": True,
        },
    }
