"""Exact unchanged EH3/EH4 currents with known polynomial denominator factors.

No global arithmetic replacement. Root coefficients, all known-factor
quotients and every complete source Ward identity are checked exactly.
"""

from itertools import product

import sympy as s
from p8_vacuum_affine_complete_two_graviton_tree.trees import instance_cache
from sympy.polys.rings import ring

SIGNS = (1, -1, -1, -1)


def submasks(mask):
    part = mask
    while True:
        yield part
        if not part:
            break
        part = (part - 1) & mask


class PolynomialCurrent:
    def __init__(self, ring, momenta, leaves):
        self.ring = ring
        self.z = ring.zero
        self.o = ring.one
        self.q = tuple(tuple(q) for q in momenta)
        self.leaves = leaves
        self.memo = {}
        self.amemo = {}
        self.eta = tuple(
            tuple(self.o * SIGNS[i] if i == j else self.z for j in range(4))
            for i in range(4)
        )

    def zeros(self):
        return tuple(tuple(self.z for j in range(4)) for i in range(4))

    def add(self, A, B):
        return tuple(tuple(A[i][j] + B[i][j] for j in range(4)) for i in range(4))

    def scale(self, c, A):
        return tuple(tuple(c * A[i][j] for j in range(4)) for i in range(4))

    def mul(self, A, B):
        return tuple(
            tuple(
                sum(
                    (A[i][k] * B[k][j] for k in range(4) if A[i][k] and B[k][j]), self.z
                )
                for j in range(4)
            )
            for i in range(4)
        )

    def trace(self, A):
        return sum((A[i][i] for i in range(4)), self.z)

    def trace_eta(self, A):
        return sum((SIGNS[i] * A[i][i] for i in range(4)), self.z)

    def factor_parts(self, p):
        constant, parts = s.factor_list(p.as_expr(), *self.ring.symbols)
        return self.ring.domain.convert(constant), [
            (self.ring.from_expr(f), k) for f, k in parts
        ]

    def merge_parts(self, *groups, lcm=False):
        rows = []
        for group in groups:
            for f, k in group:
                for index, (g, j) in enumerate(rows):
                    if f == g:
                        rows[index] = (g, max(j, k) if lcm else j + k)
                        break
                else:
                    rows.append((f, k))
        return rows

    def denominator(self, parts):
        d = self.o
        for f, k in parts:
            d *= f**k
        return d

    def normalize(self, N, parts, constant=1):
        N = self.scale(self.ring.domain.one / self.ring.domain.convert(constant), N)
        parts = self.merge_parts(parts)
        remaining = []
        for f, k in parts:
            removed = 0
            for iteration in range(k):
                rows = []
                valid = True
                for row in N:
                    made = []
                    for p in row:
                        q, remainder = p.div(f)
                        if remainder:
                            valid = False
                            break
                        made.append(q)
                    if not valid:
                        break
                    rows.append(tuple(made))
                if not valid:
                    break
                N = tuple(rows)
                removed += 1
            if removed < k:
                remaining.append((f, k - removed))
        return N, remaining

    def momentum(self, mask):
        return tuple(
            sum((self.q[i][j] for i in range(len(self.q)) if mask >> i & 1), self.z)
            for j in range(4)
        )

    def vertex(self, fields, momenta):
        kernel = self

        class Jet:
            def __init__(self):
                self._memo = {}
                self.fields = fields
                self.full = (1 << len(fields)) - 1
                self.base = []
                for A, p in zip(fields, momenta):
                    self.base.append(
                        {
                            (r, m, n): SIGNS[r]
                            * (
                                SIGNS[m] * p[m] * A[r][n]
                                + SIGNS[n] * p[n] * A[r][m]
                                - SIGNS[r] * p[r] * A[m][n]
                            )
                            / 2
                            for r, m, n in product(range(4), repeat=3)
                        }
                    )

            @instance_cache
            def inverse(self, mask):
                if not mask:
                    return kernel.eta
                total = kernel.zeros()
                for i, A in enumerate(fields):
                    if mask >> i & 1:
                        total = kernel.add(
                            total,
                            kernel.scale(
                                -2,
                                kernel.mul(
                                    kernel.mul(self.inverse(mask ^ (1 << i)), A),
                                    kernel.eta,
                                ),
                            ),
                        )
                return total

            @instance_cache
            def det(self, mask):
                chosen = [A for i, A in enumerate(fields) if mask >> i & 1]
                if not chosen:
                    return kernel.o
                if len(chosen) == 1:
                    return kernel.trace_eta(chosen[0])
                assert len(chosen) == 2
                A, B = chosen
                return kernel.trace_eta(A) * kernel.trace_eta(B) - 2 * kernel.trace(
                    kernel.mul(kernel.mul(kernel.mul(kernel.eta, A), kernel.eta), B)
                )

            @instance_cache
            def density(self, mask):
                total = kernel.zeros()
                for part in submasks(mask):
                    total = kernel.add(
                        total, kernel.scale(self.det(part), self.inverse(mask ^ part))
                    )
                return total

            @instance_cache
            def gamma(self, mask):
                result = {idx: kernel.z for idx in product(range(4), repeat=3)}
                for i, G in enumerate(self.base):
                    if not mask >> i & 1:
                        continue
                    I = self.inverse(mask ^ (1 << i))
                    for r, m, n in result:
                        result[r, m, n] += sum(
                            (
                                2 * I[r][z] * SIGNS[z] * G[z, m, n]
                                for z in range(4)
                                if I[r][z] and G[z, m, n]
                            ),
                            kernel.z,
                        )
                return result

            def bilinear(self, M, G, H):
                ht = [sum((H[z, r, z] for z in range(4)), kernel.z) for r in range(4)]
                total = kernel.z
                for m, n in product(range(4), repeat=2):
                    if not M[m][n]:
                        continue
                    first = sum(
                        (G[r, m, n] * ht[r] for r in range(4) if G[r, m, n] and ht[r]),
                        kernel.z,
                    )
                    second = sum(
                        (
                            G[r, m, z] * H[z, n, r]
                            for r, z in product(range(4), repeat=2)
                            if G[r, m, z] and H[z, n, r]
                        ),
                        kernel.z,
                    )
                    total += M[m][n] * (first - second)
                return total

            def gravity(self):
                total = kernel.z
                for A in submasks(self.full):
                    rest = self.full ^ A
                    if rest.bit_count() < 2:
                        continue
                    for B in submasks(rest):
                        C = rest ^ B
                        if B and C:
                            total += self.bilinear(
                                self.density(A), self.gamma(B), self.gamma(C)
                            )
                return -total / 2

        return Jet().gravity()

    def branch_source(self, fields, momenta):
        rows = [[self.z for j in range(4)] for i in range(4)]
        for i in range(4):
            for j in range(i, 4):
                basis = [[self.z for n in range(4)] for m in range(4)]
                basis[i][j] = basis[j][i] = self.o
                value = self.vertex((tuple(map(tuple, basis)), *fields), momenta)
                if i != j:
                    value /= 2
                rows[i][j] = rows[j][i] = value
        return tuple(map(tuple, rows))

    def amputated(self, mask):
        if mask in self.amemo:
            return self.amemo[mask]
        from p8_vacuum_affine_complete_two_graviton_tree.trees import partitions

        branches = []
        count = 0
        root = tuple(-x for x in self.momentum(mask))
        for row in partitions(mask):
            if len(row) < 2:
                continue
            children = [self.current(part) for part in row]
            fields = tuple(N for N, d, n in children)
            parts = self.merge_parts(*(d for N, d, n in children))
            moms = (root, *(self.momentum(part) for part in row))
            branches.append((self.branch_source(fields, moms), parts))
            count += s.prod(n for N, d, n in children)
        common = self.merge_parts(*(d for N, d in branches), lcm=True)
        D = self.denominator(common)
        N = self.zeros()
        for J, parts in branches:
            quotient, remainder = D.div(self.denominator(parts))
            assert not remainder
            N = self.add(N, self.scale(quotient, J))
        N, parts = self.normalize(N, common)
        Q = self.momentum(mask)
        for j in range(4):
            assert sum((SIGNS[i] * Q[i] * N[i][j] for i in range(4)), self.z) == 0, (
                "complete polynomial Ward",
                mask,
                j,
            )
        self.amemo[mask] = (N, parts, int(count))
        return self.amemo[mask]

    def current(self, mask):
        if mask in self.memo:
            return self.memo[mask]
        if mask & (mask - 1) == 0:
            N, parts = self.leaves[mask.bit_length() - 1]
            return N, parts, 1
        J, parts, count = self.amputated(mask)
        Q = self.momentum(mask)
        square = sum((SIGNS[i] * Q[i] * Q[i] for i in range(4)), self.z)
        T = tuple(
            tuple(
                (Q[0] * Q[0] if i == j else self.z) - Q[i] * Q[j]
                if i > 0 and j > 0
                else self.z
                for j in range(4)
            )
            for i in range(4)
        )
        TR = self.mul(T, J)
        N = self.add(self.scale(-1, self.mul(TR, T)), self.scale(self.trace(TR) / 2, T))
        c1, p1 = self.factor_parts(square)
        c2, p2 = self.factor_parts(Q[0])
        N, parts = self.normalize(
            N, self.merge_parts(parts, p1, [(f, 4 * k) for f, k in p2]), c1 * c2**4
        )
        self.memo[mask] = (N, parts, count)
        return self.memo[mask]

    def sympy_matrix(self, record):
        N, parts, count = record
        D = self.denominator(parts).as_expr()
        return s.ImmutableMatrix([[p.as_expr() / D for p in row] for row in N]), count


class PackedPolynomialCurrent(PolynomialCurrent):
    def __init__(self, ring, momenta, leaves, root_generators):
        super().__init__(ring, momenta, leaves)
        self.root_generators = tuple(root_generators)
        assert len(self.root_generators) == 10
        assert tuple(ring.gens[-10:]) == self.root_generators
        self.root_start = len(ring.gens) - 10

    def branch_source(self, fields, momenta):
        indices = [(i, j) for i in range(4) for j in range(i, 4)]
        root = [[self.z for j in range(4)] for i in range(4)]
        for variable, (i, j) in zip(self.root_generators, indices):
            root[i][j] = root[j][i] = variable
        value = self.vertex((tuple(map(tuple, root)), *fields), momenta)
        entries = [{} for i in range(10)]
        for monomial, coefficient in value.items():
            powers = monomial[self.root_start :]
            assert sum(powers) == 1 and max(powers) == 1, ("nonlinear root", monomial)
            index = powers.index(1)
            target = monomial[: self.root_start] + (0,) * 10
            assert target not in entries[index]
            entries[index][target] = coefficient
        rows = [[self.z for j in range(4)] for i in range(4)]
        for (i, j), coefficients in zip(indices, entries):
            polynomial = self.ring.from_dict(coefficients)
            if i != j:
                polynomial /= 2
            rows[i][j] = rows[j][i] = polynomial
        return tuple(map(tuple, rows))


class DirectProjectionPolynomialCurrent(PackedPolynomialCurrent):
    def current(self, mask):
        if mask in self.memo:
            return self.memo[mask]
        if mask & (mask - 1) == 0:
            N, parts = self.leaves[mask.bit_length() - 1]
            return N, parts, 1
        J, parts, count = self.amputated(mask)
        Q = self.momentum(mask)
        signs = (1, -1, -1, -1)
        square = sum((signs[i] * Q[i] * Q[i] for i in range(4)), self.z)
        trace = self.trace_eta(J)
        S = tuple(
            tuple(
                -signs[i] * signs[j] * J[i][j]
                + (signs[i] * trace / 2 if i == j else self.z)
                for j in range(4)
            )
            for i in range(4)
        )
        rows = [[self.z for j in range(4)] for i in range(4)]
        for i in range(1, 4):
            for j in range(1, 4):
                rows[i][j] = (
                    Q[0] ** 2 * S[i][j]
                    + Q[i] * Q[0] * S[0][j]
                    + Q[j] * Q[0] * S[i][0]
                    + Q[i] * Q[j] * S[0][0]
                )
        c1, p1 = self.factor_parts(square)
        c2, p2 = self.factor_parts(Q[0])
        N, parts = self.normalize(
            tuple(map(tuple, rows)),
            self.merge_parts(parts, p1, [(f, 2 * k) for f, k in p2]),
            c1 * c2**2,
        )
        self.memo[mask] = (N, parts, count)
        return self.memo[mask]


class CoefficientRingDivision:
    def __init__(self, original_ring):
        self.original = original_ring
        self.rings = {}
        self.factors = {}

    def configuration(self, index):
        if index not in self.rings:
            coefficient_ring = self.original.drop(self.original.gens[index])
            outer_ring, _ = ring(
                (self.original.symbols[index],), coefficient_ring.to_domain()
            )
            self.rings[index] = (coefficient_ring, outer_ring)
        return self.rings[index]

    def pack(self, polynomial, index):
        coefficient_ring, outer_ring = self.configuration(index)
        grouped = {}
        for monomial, coefficient in polynomial.items():
            degree = monomial[index]
            rest = monomial[:index] + monomial[index + 1 :]
            grouped.setdefault((degree,), {})[rest] = coefficient
        return outer_ring.from_dict(
            {
                degree: coefficient_ring.from_dict(terms)
                for degree, terms in grouped.items()
            }
        )

    def unpack(self, polynomial, index):
        terms = {}
        for (degree,), coefficient in polynomial.items():
            for rest, value in coefficient.items():
                monomial = rest[:index] + (degree,) + rest[index:]
                assert monomial not in terms
                terms[monomial] = value
        return self.original.from_dict(terms)

    def factor_configuration(self, factor):
        if factor not in self.factors:
            candidates = []
            for index in range(self.original.ngens):
                degree = max(m[index] for m in factor)
                if not degree:
                    continue
                leading = [m for m in factor if m[index] == degree]
                # A nonzero rational leading coefficient is a unit in the
                # polynomial coefficient ring, so no rational-function
                # coefficient domain or hidden GCD reduction is needed.
                if len(leading) == 1 and sum(leading[0]) == degree:
                    candidates.append((degree, index))
            if not candidates:
                self.factors[factor] = None
            else:
                _, index = min(candidates)
                self.factors[factor] = (index, self.pack(factor, index))
        return self.factors[factor]

    def div(self, polynomial, factor):
        if not polynomial:
            return self.original.zero, self.original.zero
        configuration = self.factor_configuration(factor)
        if configuration is None:
            quotient, remainder = polynomial.div(factor)
        else:
            index, packed_factor = configuration
            packed = self.pack(polynomial, index)
            quotient, remainder = packed.div(packed_factor)
            quotient = self.unpack(quotient, index)
            remainder = self.unpack(remainder, index)
        assert polynomial == quotient * factor + remainder
        return quotient, remainder


class CoefficientRingPolynomialCurrent(DirectProjectionPolynomialCurrent):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.known_factor_division = CoefficientRingDivision(self.ring)

    def normalize(self, N, parts, constant=1):
        N = self.scale(self.ring.domain.one / self.ring.domain.convert(constant), N)
        parts = self.merge_parts(parts)
        parts.sort(key=lambda item: (max(sum(m) for m in item[0]), len(item[0])))
        remaining = []
        for factor, power in parts:
            if max(sum(m) for m in factor) > 2:
                remaining.append((factor, power))
                continue
            removed = 0
            for iteration in range(power):
                rows = []
                valid = True
                for row in N:
                    made = []
                    for polynomial in row:
                        quotient, remainder = self.known_factor_division.div(
                            polynomial, factor
                        )
                        if remainder:
                            valid = False
                            break
                        made.append(quotient)
                    if not valid:
                        break
                    rows.append(tuple(made))
                if not valid:
                    break
                N = tuple(rows)
                removed += 1
            if removed < power:
                remaining.append((factor, power - removed))
        return N, remaining


class FactoredAssemblyPolynomialCurrent(CoefficientRingPolynomialCurrent):
    def amputated(self, mask):
        if mask in self.amemo:
            return self.amemo[mask]
        from p8_vacuum_affine_complete_two_graviton_tree.trees import partitions

        branches = []
        count = 0
        root = tuple(-x for x in self.momentum(mask))
        for row in partitions(mask):
            if len(row) < 2:
                continue
            children = [self.current(part) for part in row]
            fields = tuple(N for N, d, n in children)
            parts = self.merge_parts(*(d for N, d, n in children))
            moms = (root, *(self.momentum(part) for part in row))
            branches.append((self.branch_source(fields, moms), parts))
            count += s.prod(n for N, d, n in children)
        common = self.merge_parts(*(d for N, d in branches), lcm=True)
        D = self.denominator(common)
        N = self.zeros()
        for J, parts in branches:
            quotient = self.o
            for factor, power in common:
                branch_power = next((k for f, k in parts if f == factor), 0)
                assert 0 <= branch_power <= power
                quotient *= factor ** (power - branch_power)
            assert quotient * self.denominator(parts) == D
            N = self.add(N, self.scale(quotient, J))
        N, parts = self.normalize(N, common)
        Q = self.momentum(mask)
        signs = (1, -1, -1, -1)
        for j in range(4):
            assert sum((signs[i] * Q[i] * N[i][j] for i in range(4)), self.z) == 0, (
                "complete polynomial Ward",
                mask,
                j,
            )
        self.amemo[mask] = (N, parts, int(count))
        return self.amemo[mask]


def multiply_energies_over_total(engine, record, energies, total):
    """Exact monomial cancellation without repeating irrelevant trial divisions.

    The denominator identity is checked independently of the tensor numerator.
    No coprimality or irreducibility assumption is required for this operation.
    """
    numerator, parts, count = record
    old_denominator = engine.denominator(parts)
    remaining = list(parts)
    multiplier = engine.o
    for energy in energies:
        for index, (factor, power) in enumerate(remaining):
            if factor == energy:
                assert power >= 1
                if power == 1:
                    remaining.pop(index)
                else:
                    remaining[index] = (factor, power - 1)
                break
        else:
            multiplier *= energy
    energy_product = engine.o
    for energy in energies:
        energy_product *= energy
    assert multiplier * old_denominator == energy_product * engine.denominator(
        remaining
    )
    new_numerator = engine.scale(multiplier, numerator)
    return new_numerator, engine.merge_parts(remaining, [(total, 1)]), count
