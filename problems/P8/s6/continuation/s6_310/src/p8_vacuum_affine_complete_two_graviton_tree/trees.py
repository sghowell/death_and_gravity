"""Exact labeled rooted trees from the covariant action, through six points."""

from functools import cache, wraps
from itertools import permutations, product

import sympy as s
from p8_vacuum_affine_minimal_gravity_radiation import vertices as old

ETA = s.ImmutableMatrix(s.diag(1, -1, -1, -1))
ZERO = s.ImmutableMatrix(s.zeros(4))
VECTOR_ZERO = s.ImmutableMatrix(s.zeros(4, 1))
KINDS = ("phi", "H", "h")
imm = s.ImmutableMatrix


def clean(x):
    return s.factor(x)


def tr(A):
    return s.trace(ETA * A)


@cache
def determinant_coefficient(fields):
    count = len(fields)
    if count == 0:
        return s.S.One
    if count == 1:
        return tr(fields[0])
    if count == 2:
        A, B = fields
        return tr(A) * tr(B) - 2 * s.trace(ETA * A * ETA * B)
    if count == 3:
        return sum(
            tr(A) * tr(B) * tr(C) / 6
            - tr(A) * s.trace(ETA * B * ETA * C)
            + s.Rational(4, 3) * s.trace(ETA * A * ETA * B * ETA * C)
            for A, B, C in permutations(fields)
        )
    raise ValueError(
        "Only scalar metric vertices through three gravitons are implemented"
    )


@cache
def inverse_coefficient(fields):
    if not fields:
        return ETA
    result = s.zeros(4)
    for order in permutations(fields):
        value = ETA
        for H in order:
            value = value * H * ETA
        result += (-2) ** len(fields) * value
    return imm(result)


@cache
def density_inverse_coefficient(fields):
    result = s.zeros(4)
    for mask in range(1 << len(fields)):
        chosen = tuple(A for i, A in enumerate(fields) if mask >> i & 1)
        rest = tuple(A for i, A in enumerate(fields) if not (mask >> i & 1))
        result += determinant_coefficient(chosen) * inverse_coefficient(rest)
    return imm(result)


@cache
def scalar_vertex(p, q, fields, mass):
    return clean(
        -(p.T * ETA * density_inverse_coefficient(fields) * ETA * q)[0]
        - mass * determinant_coefficient(fields)
    )


@cache
def quartic_gravity(fields, momenta):
    gs = tuple(old.connection(A, p) for A, p in zip(fields, momenta))
    result = s.S.Zero
    for a, b, c, d in permutations(range(4)):
        A, B = fields[a], fields[b]
        M1 = tr(A) * ETA / 2 - ETA * A * ETA
        M2 = (
            ETA * (tr(A) * tr(B) / 8 - s.trace(ETA * A * ETA * B) / 4)
            - tr(A) * ETA * B * ETA / 2
            + ETA * A * ETA * B * ETA
        )
        result += (
            old.bilinear(M2, gs[c], gs[d])
            + old.bilinear(M1, old.second(B, gs[c]), gs[d])
            + old.bilinear(M1, gs[c], old.second(B, gs[d]))
            + old.bilinear(ETA, old.second(A, old.second(B, gs[c])), gs[d])
            + old.bilinear(ETA, gs[c], old.second(A, old.second(B, gs[d])))
            + old.bilinear(ETA, old.second(A, gs[c]), old.second(B, gs[d]))
        )
    return clean(-8 * result)


@cache
def cubic_gravity(fields, momenta):
    return clean(old.cubic(fields, momenta))


@cache
def partitions(mask):
    if not mask:
        return ((),)
    first = mask & -mask
    rest = mask ^ first
    result = []
    for row in partitions(rest):
        result.append((first, *row))
        for i in range(len(row)):
            result.append((*row[:i], row[i] | first, *row[i + 1 :]))
    # Each unordered partition appears exactly once.
    assert len({frozenset(row) for row in result}) == len(result)
    return tuple(result)


ROOT_BASIS = []
for i in range(4):
    for j in range(i, 4):
        basis = s.zeros(4)
        basis[i, j] = basis[j, i] = 1
        ROOT_BASIS.append((i, j, imm(basis)))


def instance_cache(function):
    """Cache only within the finite tree evaluator that owns each current."""

    @wraps(function)
    def cached(self, *args):
        values = self._memo.setdefault(function.__name__, {})
        if args not in values:
            values[args] = function(self, *args)
        return values[args]

    return cached


class TreeEngine:
    def __init__(self, legs, heavy=128, cubic=0, contact=0, kappa=1):
        self.legs = tuple(
            (kind, imm(p), imm(value) if kind == "h" else s.sympify(value))
            for kind, p, value in legs
        )
        self.heavy, self.cubic, self.contact, self.kappa = map(
            s.sympify, (heavy, cubic, contact, kappa)
        )
        self.full = (1 << len(self.legs)) - 1
        self._memo = {}

    @instance_cache
    def momentum(self, mask):
        return imm(
            sum(
                (leg[1] for i, leg in enumerate(self.legs) if mask >> i & 1),
                VECTOR_ZERO,
            )
        )

    def permitted(self, tags):
        nf, nh, ng = (tags.count(kind) for kind in KINDS)
        if nf == 2 and nh == 0 and 1 <= ng <= 3:
            return True
        if nh == 2 and nf == 0 and 1 <= ng <= 3:
            return True
        if nf == 2 and nh == 1 and 0 <= ng <= 2:
            return self.cubic != 0
        if nf == 4 and nh == 0 and 0 <= ng <= 2:
            return self.contact != 0
        return nf == nh == 0 and ng in (3, 4)

    @instance_cache
    def vertex(self, tags, momenta, values):
        assert self.permitted(tags)
        phi = [i for i, k in enumerate(tags) if k == "phi"]
        heavy = [i for i, k in enumerate(tags) if k == "H"]
        gravitons = [i for i, k in enumerate(tags) if k == "h"]
        hs = tuple(values[i] for i in gravitons)
        ps = tuple(momenta[i] for i in gravitons)
        scalar_product = s.prod(values[i] for i, k in enumerate(tags) if k != "h")
        r = len(hs)
        if len(phi) == 2 and not heavy:
            value = scalar_vertex(momenta[phi[0]], momenta[phi[1]], hs, s.S.One)
        elif len(heavy) == 2 and not phi:
            value = scalar_vertex(momenta[heavy[0]], momenta[heavy[1]], hs, self.heavy)
        elif len(phi) == 2 and len(heavy) == 1:
            value = -self.cubic * determinant_coefficient(hs)
        elif len(phi) == 4:
            value = self.contact * determinant_coefficient(hs)
        else:
            value = (
                cubic_gravity(hs, ps) if r == 3 else quartic_gravity(hs, ps)
            ) * self.kappa
        return clean(scalar_product * value / self.kappa ** s.Rational(r, 2))

    @instance_cache
    def amputated(self, mask, kind):
        result = s.zeros(4) if kind == "h" else s.S.Zero
        graph_count = 0
        root_momentum = imm(-self.momentum(mask))
        for row in partitions(mask):
            if len(row) < 2:
                continue
            for child_tags in product(KINDS, repeat=len(row)):
                tags = (kind, *child_tags)
                if not self.permitted(tags):
                    continue
                children = tuple(
                    self.current(part, tag) for part, tag in zip(row, child_tags)
                )
                counts = [item[1] for item in children]
                if not all(counts):
                    continue
                child_values = tuple(item[0] for item in children)
                moms = (root_momentum, *(self.momentum(part) for part in row))
                graph_count += s.prod(counts)
                if kind == "h":
                    for i, j, basis in ROOT_BASIS:
                        value = self.vertex(tags, moms, (basis, *child_values))
                        result[i, j] += value if i == j else value / 2
                        if i != j:
                            result[j, i] += value / 2
                else:
                    result += self.vertex(tags, moms, (s.S.One, *child_values))
        return (
            imm(result.applyfunc(clean)) if kind == "h" else clean(result),
            int(graph_count),
        )

    @instance_cache
    def current(self, mask, kind):
        if mask & (mask - 1) == 0:
            leg = self.legs[mask.bit_length() - 1]
            return (
                (leg[2], 1)
                if leg[0] == kind
                else (ZERO if kind == "h" else s.S.Zero, 0)
            )
        source, count = self.amputated(mask, kind)
        if count == 0:
            return (ZERO if kind == "h" else s.S.Zero, 0)
        P = self.momentum(mask)
        square = old.dot(P, P)
        if kind == "h":
            if square == 0:
                raise ValueError("A massless internal propagator is on its pole")
            value = -(ETA * source * ETA - ETA * s.trace(ETA * source) / 2) / square
            return imm(value.applyfunc(clean)), count
        mass = s.S.One if kind == "phi" else self.heavy
        if square == mass:
            raise ValueError("A massive internal propagator is on its pole")
        return clean(-source / (square - mass)), count

    def amplitude(self):
        return self.amputated(self.full, "phi")


def amplitude(points, gravitons, **parameters):
    assert all(
        s.simplify(v) == 0
        for v in sum(points, VECTOR_ZERO)
        + sum((q for q, _eps in gravitons), VECTOR_ZERO)
    )
    legs = [("phi", p, 1) for p in points[1:]]
    legs.extend(("h", q, eps) for q, eps in gravitons)
    return TreeEngine(legs, **parameters).amplitude()


def exact_real(value):
    if value is None or isinstance(value, (str, bool, float, s.Float)):
        raise TypeError("Require exact real entries")
    value = s.sympify(value)
    if (
        value.has(s.Float)
        or value.is_number is not True
        or value.is_real is not True
        or value.is_finite is not True
    ):
        raise ValueError("Require finite exact real entries")
    return value


def exact_array(value):
    if isinstance(value, s.MatrixBase):
        for entry in value:
            exact_real(entry)
    elif isinstance(value, (tuple, list)):
        for entry in value:
            if isinstance(entry, (tuple, list)):
                for item in entry:
                    exact_real(item)
            else:
                exact_real(entry)
    else:
        raise TypeError("Require an exact matrix or nested list")
    return imm(value)


def original_amplitude(points, gravitons):
    """Original-parameter physical TT coefficient; no integration or soft cut."""
    from . import source

    if not isinstance(points, (tuple, list)) or len(points) != 4:
        raise ValueError("Require four all-outgoing massive scalar momenta")
    if not isinstance(gravitons, (tuple, list)) or len(gravitons) > 2:
        raise ValueError("Only zero, one or two real gravitons are implemented")
    ps = []
    for p in points:
        p = exact_array(p)
        if p.shape != (4, 1):
            raise ValueError("Require four-vectors")
        for x in p:
            exact_real(x)
        if clean(old.dot(p, p) - 1) != 0:
            raise ValueError("Require mass-one external shells")
        ps.append(p)
    E = -ps[0][0]
    if (
        not s.Rational(5, 4) <= E <= 2
        or ps[1][0] != -E
        or ps[0][1:, 0] + ps[1][1:, 0] != s.zeros(3, 1)
    ):
        raise ValueError("Require the stated incoming center-of-mass domain")
    if any(p[0] <= 0 for p in ps[2:]):
        raise ValueError("Require outgoing positive-energy massive legs")
    hs = []
    for item in gravitons:
        if not isinstance(item, (tuple, list)) or len(item) != 2:
            raise ValueError("Require momentum-polarization pairs")
        q, eps = map(exact_array, item)
        if q.shape != (4, 1) or eps.shape != (4, 4):
            raise ValueError("Require a four-vector and symmetric tensor")
        for x in (*q, *eps):
            exact_real(x)
        if q[0] <= 0 or clean(old.dot(q, q)) != 0:
            raise ValueError("Require future null radiation")
        if (
            eps != eps.T
            or any(eps[0, j] != 0 for j in range(4))
            or tr(eps) != 0
            or eps * q != s.zeros(4, 1)
        ):
            raise ValueError("Require spatial transverse-traceless polarizations")
        hs.append((q, eps))
    if sum((q[0] for q, _ in hs), s.S.Zero) > s.Rational(1, 8):
        raise ValueError("Require total emitted energy at most 1/8")
    if sum(ps, VECTOR_ZERO) + sum((q for q, _ in hs), VECTOR_ZERO) != s.zeros(4, 1):
        raise ValueError("Require exact total momentum conservation")
    if len(hs) == 2 and old.dot(hs[0][0], hs[1][0]) == 0:
        raise ValueError(
            "The point evaluator excludes an exactly collinear propagator pole"
        )
    return amplitude(
        tuple(ps),
        hs,
        heavy=source.HEAVY_MASS2,
        cubic=source.CUBIC,
        contact=source.CONTACT,
        kappa=source.KAPPA,
    )
