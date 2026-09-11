"""Exact ordered matrix words and isotropic sphere contractions."""

from functools import cache

import sympy as s

d, h, h1, z = s.symbols("d H Hdot z", real=True)
f, f2, g = s.symbols("nFn nF2n nGn", real=True)
TF2 = s.Symbol("tr_FF", real=True)


def reduce_word(word):
    out = []
    for letter in word:
        if letter == "P" and out and out[-1] == "P":
            continue
        out.append(letter)
    return tuple(out)


class NC:
    def __init__(self, items=0):
        if isinstance(items, NC):
            self.terms = dict(items.terms)
        elif isinstance(items, dict):
            self.terms = {k: s.expand(v) for k, v in items.items() if v != 0}
        else:
            v = s.sympify(items)
            self.terms = {(): v} if v != 0 else {}

    @classmethod
    def letter(cls, name):
        return cls({(name,): s.Integer(1)})

    def __add__(self, other):
        other = NC(other)
        result = dict(self.terms)
        for key, value in other.terms.items():
            result[key] = s.expand(result.get(key, 0) + value)
        return NC(result)

    __radd__ = __add__

    def __neg__(self):
        return NC({k: -v for k, v in self.terms.items()})

    def __sub__(self, other):
        return self + -NC(other)

    def __rsub__(self, other):
        return NC(other) + -self

    def __mul__(self, other):
        other = NC(other)
        result = {}
        for left, a in self.terms.items():
            for right, b in other.terms.items():
                word = reduce_word(left + right)
                result[word] = result.get(word, 0) + a * b
        return NC(result)

    __rmul__ = __mul__

    def __truediv__(self, value):
        return NC({k: v / value for k, v in self.terms.items()})

    def __pow__(self, n):
        if type(n) is not int or n < 0:
            raise ValueError("Require a nonnegative native ordered power")
        out = NC(1)
        for _ in range(n):
            out = out * self
        return out

    @property
    def T(self):
        return NC({tuple(reversed(k)): v for k, v in self.terms.items()})

    @property
    def H(self):
        return NC({tuple(reversed(k)): s.conjugate(v) for k, v in self.terms.items()})

    def real(self):
        return NC({k: s.re(v) for k, v in self.terms.items()})

    def cancel(self):
        return NC({k: s.cancel(v) for k, v in self.terms.items()})

    def evaluate(self, matrices):
        n = next(iter(matrices.values())).rows
        out = s.zeros(n)
        for word, coefficient in self.terms.items():
            value = s.eye(n)
            for letter in word:
                value = value * matrices[letter]
            out += coefficient * value
        return out.applyfunc(s.expand)


def cyclic(word):
    if not word:
        return ()
    variants = []
    for row in (tuple(word), tuple(reversed(word))):
        variants.extend(row[j:] + row[:j] for j in range(len(row)))
    return min(variants)


@cache
def trace_word(word):
    word = cyclic(word)
    if not word:
        return d
    if word == ("F",):
        return s.Integer(0)
    if word == ("G",):
        return TF2
    return s.Symbol("tr_" + "".join(word), real=True)


def trace(poly):
    return s.expand(
        sum(value * trace_word(word) for word, value in NC(poly).terms.items())
    )


@cache
def pairings(nodes):
    if not nodes:
        return ((),)
    first = nodes[0]
    result = []
    for j in range(1, len(nodes)):
        rest = nodes[1:j] + nodes[j + 1 :]
        for tail in pairings(rest):
            result.append(((first, nodes[j]),) + tail)
    return tuple(result)


def projection_key(words):
    return tuple(sorted(min(tuple(w), tuple(reversed(w))) for w in words if w))


@cache
def projection_average(key):
    if not key:
        return s.Integer(1)
    count = len(key)
    total = 0
    for pairs in pairings(tuple(range(2 * count))):
        delta = {}
        for a, b in pairs:
            delta[a] = b
            delta[b] = a
        unused = set(range(count))
        product = 1
        while unused:
            start = 2 * min(unused)
            current = start
            word = ()
            while True:
                edge = current // 2
                unused.remove(edge)
                word += key[edge] if current % 2 == 0 else tuple(reversed(key[edge]))
                current = delta[current ^ 1]
                if current == start:
                    break
            product *= trace_word(word)
        total += product
    denominator = s.prod(d + 2 * j for j in range(count))
    return s.factor(total / denominator)


def trace_structure(word):
    if "P" not in word:
        return trace_word(word), ()
    first = word.index("P")
    row = word[first + 1 :] + word[:first]
    segments = []
    current = ()
    for letter in row:
        if letter == "P":
            segments.append(current)
            current = ()
        else:
            current += (letter,)
    segments.append(current)
    return s.Integer(1), tuple(segments)


@cache
def angular_trace(poly_terms):
    total = 0
    for word, coefficient in poly_terms:
        scalar, segments = trace_structure(word)
        for powers, c in s.Poly(coefficient, f, f2, g).terms():
            extra = (
                (("F",),) * powers[0]
                + (("F", "F"),) * powers[1]
                + (("G",),) * powers[2]
            )
            total += c * scalar * projection_average(projection_key(segments + extra))
    return s.factor(total)


def average(poly):
    return angular_trace(tuple(sorted(NC(poly).terms.items())))


def radial(poly, denominator_half):
    result = 0
    for (n,), coefficient in s.Poly(poly, z).terms():
        result += coefficient * s.rf(d / 2, n) / s.rf(denominator_half, n)
    return s.factor(result)


@cache
def data():
    T2 = trace_word(("F", "F"))
    T3 = trace_word(("F", "F", "F"))
    T4 = trace_word(("F", "F", "F", "F"))
    TFG = trace_word(("F", "G"))
    TFFG = trace_word(("F", "F", "G"))
    D2 = d * (d + 2)
    D3 = D2 * (d + 4)
    D4 = D3 * (d + 6)
    cases = {
        "tracefree_first": ((("F",),), 0),
        "unimodular_second_jet": ((("G",),), T2 / d),
        "ordered_square": ((("F", "F"),), T2 / d),
        "ordered_cube": ((("F", "F", "F"),), T3 / d),
        "ordered_fourth": ((("F", "F", "F", "F"),), T4 / d),
        "tracefree_pair": ((("F",), ("F",)), 2 * T2 / D2),
        "mixed_first_second": ((("F",), ("G",)), 2 * TFG / D2),
        "first_and_square": ((("F",), ("F", "F")), 2 * T3 / D2),
        "square_and_second": ((("F", "F"), ("G",)), (T2**2 + 2 * TFFG) / D2),
        "two_squares": ((("F", "F"), ("F", "F")), (T2**2 + 2 * T4) / D2),
        "three_tracefree": ((("F",), ("F",), ("F",)), 8 * T3 / D3),
        "two_tracefree_and_second": (
            (("F",), ("F",), ("G",)),
            (2 * T2**2 + 8 * TFFG) / D3,
        ),
        "two_tracefree_and_square": (
            (("F",), ("F",), ("F", "F")),
            (2 * T2**2 + 8 * T4) / D3,
        ),
        "four_tracefree": (
            (("F",), ("F",), ("F",), ("F",)),
            (12 * T2**2 + 48 * T4) / D4,
        ),
        "noncommuting_mixed_square": (
            (("F", "G"), ("G", "F")),
            (
                TFG**2
                + trace_word(("F", "G", "F", "G"))
                + trace_word(("F", "F", "G", "G"))
            )
            / D2,
        ),
    }
    checks = {
        name: s.factor(projection_average(projection_key(words)) - expected)
        for name, (words, expected) in cases.items()
    }
    for n in range(1, 5):
        checks[f"all_endpoint_pairing_count_{n}"] = len(
            pairings(tuple(range(2 * n)))
        ) - s.factorial2(2 * n - 1)
    for n in range(5):
        checks[f"physical_dimension_inverse_frequency_radial_ratio_{n}"] = s.rf(
            s.Rational(3, 2), n
        ) / s.rf(s.Rational(1, 2), n) - (2 * n + 1)
        checks[f"physical_dimension_inverse_cubic_radial_ratio_{n}"] = (
            s.rf(s.Rational(3, 2), n) / s.rf(s.Rational(3, 2), n) - 1
        )
    return {
        "algebra": "Matrix words preserve order; only adjacent rank-one projectors P=nn^t satisfy P^2=P. Cyclic trace rotations and word reversal use symmetric F,G, not commutation.",
        "sphere_contraction": "Each product of r quadratic projections is averaged by all pairings of its2r endpoints divided by d(d+2)...(d+2r-2). Each contraction loop is a trace of its oriented matrix-edge words.",
        "unimodular_constraints": "At E=I, F=E' and G=E'' obey trF=0 and trG=trF^2. Ordered projections n^t FG n and n^t GF n agree as scalars, but trFGFG and trFFGG need not agree.",
        "radial_identity": "Integral z^n/omega^(2alpha) divided by its n0 value equals (d/2)_n/(alpha)_n after meromorphic continuation. alpha=1/2 at order2 and3/2 at order4; retain the d dependence before the dimension limit.",
        "explicit_sphere_moments": {
            name: expected for name, (_, expected) in cases.items()
        },
        "checks": checks,
        "gates": {
            "fourth_order_pairings_complete": len(pairings(tuple(range(8)))) == 105,
            "noncommuting_trace_classes_distinct": trace_word(("F", "G", "F", "G"))
            != trace_word(("F", "F", "G", "G")),
            "projector_not_identity": trace_word(()) == d
            and trace_structure(("P",)) == (1, ((),)),
            "no_dimension_three_angular_shortcut": d
            in projection_average((("F",), ("F",))).free_symbols,
        },
    }
