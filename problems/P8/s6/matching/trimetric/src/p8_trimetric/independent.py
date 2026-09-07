"""Separate exact Fraction matrix jets and coefficientwise determinant replay.

No SymPy or main-engine imports. Only the immutable sparse polynomial
primitive is reused from S6.7; none of its physical verdicts is imported.
"""

from dataclasses import dataclass
from fractions import Fraction as Q
from itertools import combinations, permutations

from p8_composite_modes.independent import Poly


def exact(value):
    if isinstance(value, bool) or not isinstance(value, (int, Q)):
        raise TypeError("Only exact int/Fraction matrix-jet inputs are accepted")
    return Q(value)


@dataclass(frozen=True)
class Dual:
    value: Q
    tangent: Q = Q(0)

    def __post_init__(self):
        object.__setattr__(self, "value", exact(self.value))
        object.__setattr__(self, "tangent", exact(self.tangent))

    def __add__(self, other):
        other = other if isinstance(other, Dual) else Dual(other)
        return Dual(self.value+other.value, self.tangent+other.tangent)

    __radd__ = __add__

    def __neg__(self):
        return Dual(-self.value, -self.tangent)

    def __sub__(self, other):
        return self+-other

    def __rsub__(self, other):
        return other+-self

    def __mul__(self, other):
        other = other if isinstance(other, Dual) else Dual(other)
        return Dual(self.value*other.value, self.value*other.tangent+self.tangent*other.value)

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = other if isinstance(other, Dual) else Dual(other)
        if not other.value:
            raise ZeroDivisionError("A matrix-jet inverse needs a nonzero constant denominator")
        return Dual(self.value/other.value,
                    (self.tangent*other.value-self.value*other.tangent)/other.value**2)

    def __rtruediv__(self, other):
        return Dual(other)/self


def transpose(a):
    return [list(row) for row in zip(*a, strict=True)]


def multiply(a, b):
    return [[sum(a[i][k]*b[k][j] for k in range(len(b))) for j in range(len(b[0]))] for i in range(len(a))]


def add(a, b):
    return [[x+y for x, y in zip(arow, brow, strict=True)] for arow, brow in zip(a, b, strict=True)]


def scale(coefficient, a):
    return [[coefficient*x for x in row] for row in a]


def identity():
    return [[Q(i == j) for j in range(4)] for i in range(4)]


def trace(a):
    return sum(a[i][i] for i in range(len(a)))


def determinant(a):
    total = 0
    for perm in permutations(range(4)):
        sign = (-1)**sum(perm[i] > perm[j] for i in range(4) for j in range(i+1, 4))
        term = sign
        for i, j in enumerate(perm):
            term *= a[i][j]
        total += term
    return total


def inverse(a):
    # Exact Gauss-Jordan inversion over Fraction or the first-jet ring.
    work = [[x if isinstance(x, Dual) else exact(x) for x in row]+unit
            for row, unit in zip(a, identity(), strict=True)]
    for column in range(4):
        candidates = [i for i in range(column, 4) if
                      (work[i][column].value if isinstance(work[i][column], Dual) else work[i][column])]
        if not candidates:
            raise ZeroDivisionError("Singular constant matrix in exact jet")
        pivot = candidates[0]
        work[column], work[pivot] = work[pivot], work[column]
        denominator = work[column][column]
        work[column] = [x/denominator for x in work[column]]
        for row in range(4):
            if row != column:
                coefficient = work[row][column]
                work[row] = [x-coefficient*y for x, y in zip(work[row], work[column], strict=True)]
    return [row[4:] for row in work]


def literal_density(e, v, u, b, pg, pf, bg=0, bf=0):
    source = add(scale(pg, e), scale(pf, v))
    return -2*determinant(u)*(b+trace(multiply(inverse(u), source)))-2*bg*determinant(e)-2*bf*determinant(v)


def jet_matrix(a, row, column):
    return [[Dual(x, int((i, j) == (row, column))) for j, x in enumerate(line)] for i, line in enumerate(a)]


def matrix_jet_checks():
    e = [[Q(x) for x in row] for row in ((2, 1, 0, 0), (0, 3, 1, 0), (0, 0, 2, 1), (1, 0, 0, 2))]
    v = [[Q(x) for x in row] for row in ((1, 0, 1, 0), (0, 2, 0, 1), (0, 0, 3, 0), (0, 0, 0, 1))]
    u = [[Q(x) for x in row] for row in ((2, 1, 0, 0), (0, 1, 1, 0), (0, 0, 3, 1), (0, 0, 0, 2))]
    b, pg, pf, bg, bf = Q(-7), Q(2), Q(3), Q(-1), Q(4)
    inv = inverse(u)
    it = transpose(inv)
    det = determinant(u)
    source = add(scale(pg, e), scale(pf, v))
    scalar = b+trace(multiply(inv, source))
    ee = add(scale(-2*pg*det, it), scale(-2*bg*determinant(e), transpose(inverse(e))))
    ev = add(scale(-2*pf*det, it), scale(-2*bf*determinant(v), transpose(inverse(v))))
    eu = scale(-2*det, add(scale(scalar, it), scale(-1, multiply(multiply(it, transpose(source)), it))))
    checked = 0
    for which, gradient in enumerate((ee, ev, eu)):
        for row in range(4):
            for column in range(4):
                arguments = [e, v, u]
                arguments[which] = jet_matrix(arguments[which], row, column)
                actual = literal_density(*arguments, b, pg, pf, bg, bf)
                if actual.tangent != gradient[row][column]:
                    raise ValueError("Independent full-matrix first-jet Euler identity failed")
                checked += 1
    eta = [[Q((1 if i == 0 else -1) if i == j else 0) for j in range(4)] for i in range(4)]
    k = [[Q(1)], [Q(1)], [Q(0)], [Q(1)]]
    potential = Q(2, 3)
    h = multiply(multiply(transpose(u), eta), u)
    raised = multiply(inverse(h), k)
    kinetic = multiply(transpose(k), raised)[0][0]
    ju = scale(det, add(scale(kinetic/2-potential, it),
                        scale(-1, multiply(multiply(multiply(eta, u), raised), transpose(raised)))))
    for row in range(4):
        for column in range(4):
            uj = jet_matrix(u, row, column)
            hj = multiply(multiply(transpose(uj), eta), uj)
            yj = multiply(multiply(transpose(k), inverse(hj)), k)[0][0]
            actual = determinant(uj)*(yj/2-potential)
            if actual.tangent != ju[row][column]:
                raise ValueError("Independent canonical matter first-jet Euler identity failed")
            checked += 1
    return checked


def polynomial_checks():
    roots = [Poly.variable(name) for name in ("A", "B", "C", "x")]
    plus, minus = Poly(1), Poly(1)
    for root in roots:
        plus *= 1+root
        minus *= 1-root
    e2 = sum((roots[i]*roots[j] for i, j in combinations(range(4), 2)), Poly())
    e4 = roots[0]*roots[1]*roots[2]*roots[3]
    c, r, x = (Poly.variable(name) for name in ("A", "B", "x"))
    beta = [c*r**n for n in range(5)]
    values = {"geometric_middle_minor": beta[1]*beta[3]-beta[2]**2,
              "full_beta4_extended_relative_determinant": -4+2*(plus+minus)-4*(e2+e4),
              "truncated_vs_exact_constant_V_remainder":
                  (Q(2, 3)-x)*(1+x*Q(1, 2))**3-Q(2, 3)+x**2+Q(2, 3)*x**3+x**4*Q(1, 8)}
    if any(not value.is_zero() for value in values.values()):
        raise ValueError("Independent coefficientwise auxiliary-parent identity failed")
    return dict.fromkeys(values, "0")


def checks():
    count = matrix_jet_checks()
    identities = polynomial_checks()
    fixtures = {"original_g_flat_residual": Q(-2), "original_f_flat_residual": Q(-2),
                "beta4_extension_g_flat_residual": -2*Q(1)-2*Q(-1),
                "beta4_extension_f_flat_residual": -2*Q(1)-2*Q(-1),
                "beta4_extension_u_flat_residual": -2*(Q(-6)+8-2),
                "extension_mass_squared_G_F_q_1": Q(2),
                "exact_constant_V_action_at_formal_cancel": Q(2, 3)/(1+Q(1, 3))**3,
                "required_antisymmetric_constraint_coefficient": (Q(1)-Q(4, 3))/2}
    if count != 64 or tuple(fixtures.values()) != (-2, -2, 0, 0, 0, 2, Q(9, 32), Q(-1, 6)):
        raise ValueError("Independent exact domain/control fixture failed")
    return {"coefficientwise_identities": identities,
            "full_matrix_first_jet_directions_checked": count,
            "Fraction_vacuum_source_and_extension_fixtures": {key: str(value) for key, value in fixtures.items()}}
