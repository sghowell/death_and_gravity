"""Euclidean FLRW orthonormal-frame derivatives and bilinear contractions."""
from collections import Counter
from functools import cache, lru_cache
from itertools import combinations, pairwise, product

import sympy as sp

H = sp.symbols("H0:8", real=True)
k = sp.Symbol("physical_spatial_momentum", real=True)
mass = sp.Symbol("constant_reference_mass", positive=True)


class ModeTensor:
    def __init__(self, name, sign):
        if type(sign) is not int or (name, sign) not in (("plus", 1), ("minus", -1)):
            raise ValueError("Require one of the two specified Fourier legs")
        self.sign = sign
        self.temporal = sp.symbols(name+"_Yt0:5", real=True)
        self.spatial = sp.symbols(name+"_Ys0:5", real=True)
        self.flows = {**dict(pairwise(H)), k: -H[0]*k,
                      **dict(pairwise(self.temporal)), **dict(pairwise(self.spatial))}
        self._covariant = lru_cache(maxsize=None, typed=True)(self._covariant)

    def time_derivative(self, value):
        return sum(sp.diff(value, var)*self.flows[var] for var in value.free_symbols if var in self.flows)

    def covariant(self, sequence, i, j):
        if type(sequence) is not tuple or len(sequence) > 4:
            raise ValueError("Require native four-dimensional indices and derivative order <=4")
        if any(type(index) is not int or not 0 <= index < 4 for index in sequence+(i, j)):
            raise ValueError("Require native four-dimensional indices and derivative order <=4")
        return self._covariant(sequence, i, j)

    def _covariant(self, sequence, i, j):
        if not sequence:
            return self.temporal[0] if i == j == 0 else self.spatial[0] if i == j else sp.Integer(0)
        a, rest = sequence[0], sequence[1:]
        previous = self.covariant(rest, i, j)
        if a == 0:
            value = self.time_derivative(previous)
        elif a == 1:
            value = self.sign*sp.I*k*previous
        else:
            value = sp.Integer(0)
        if a != 0:
            inner = rest+(i, j)
            for position, index in enumerate(inner):
                if index == 0:
                    replacement, factor = a, -H[0]
                elif index == a:
                    replacement, factor = 0, H[0]
                else:
                    continue
                changed = inner[:position]+(replacement,)+inner[position+1:]
                value += factor*self.covariant(changed[:-2], changed[-2], changed[-1])
        return sp.expand(value)


plus, minus = ModeTensor("plus", 1), ModeTensor("minus", -1)


def ricci(i, j):
    if i != j:
        return sp.Integer(0)
    return -3*(H[1]+H[0]**2) if i == 0 else -H[1]-3*H[0]**2


def riemann(a, b, c, d):
    if a == b or c == d:
        return sp.Integer(0)
    sectional = -(H[1]+H[0]**2) if 0 in (a, b) else -H[0]**2
    if a == c and b == d:
        return sectional
    if a == d and b == c:
        return -sectional
    return sp.Integer(0)


def scalar_curvature():
    return -6*H[1]-12*H[0]**2


def swap_legs(expression):
    mapping = {k: -k}
    for left, right in ((plus.temporal, minus.temporal), (plus.spatial, minus.spatial)):
        mapping.update(dict(zip(left, right)))
        mapping.update(dict(zip(right, left)))
    return expression.xreplace(mapping)


def contraction(factors):
    if type(factors) is not tuple or not factors:
        raise ValueError("Require an explicit immutable tensor contraction")
    arities = {"R": 0, "Ric": 2, "Riem": 4, "Y": 2}
    for factor in factors:
        if type(factor) is not tuple or len(factor) != 3:
            raise ValueError("Require kind, ordered derivatives and tensor indices")
        kind, sequence, indices = factor
        if kind not in arities or type(sequence) is not str or type(indices) is not str:
            raise ValueError("Unknown tensor or invalid Einstein labels")
        if len(indices) != arities[kind] or (kind != "Y" and sequence) or len(sequence) > 4:
            raise ValueError("Invalid tensor rank or derivative order")
        if any(label not in "abcdefghijklmnopqrstuvwxyz" for label in sequence+indices):
            raise ValueError("Require lowercase single-letter Einstein labels")
    counts = Counter("".join(sequence+indices for _, sequence, indices in factors))
    if any(number != 2 for number in counts.values()) or sum(kind == "Y" for kind, _, _ in factors) != 2:
        raise ValueError("Require exactly two mass insertions and paired Euclidean indices")
    return _contraction(factors)


@cache
def _contraction(factors):
    counts = Counter("".join(sequence+indices for _, sequence, indices in factors))
    labels = tuple(counts)
    assignments = ({label: index for label, index in zip(labels, indices)}
                   for indices in product(range(4), repeat=len(labels)))
    pieces = []
    for indices in assignments:
        value, leg = sp.Integer(1), 0
        for kind, sequence, slots in factors:
            evaluated = tuple(indices[label] for label in slots)
            if kind == "R":
                component = scalar_curvature()
            elif kind == "Ric":
                component = ricci(*evaluated)
            elif kind == "Riem":
                component = riemann(*evaluated)
            elif kind == "Y":
                mode = plus if leg == 0 else minus
                component = mode.covariant(tuple(indices[label] for label in sequence), *evaluated)
                leg += 1
            else:
                raise ValueError("Unknown tensor in a quadratic mass-insertion invariant")
            if component == 0:
                value = sp.Integer(0)
                break
            value *= component
        if value != 0:
            pieces.append(value)
    unsymmetrized = sp.expand(sum(pieces))
    return sp.expand((unsymmetrized+swap_legs(unsymmetrized))/2)


def checks():
    out = {"Euclidean_FLRW_Ricci_from_Riemann_"+str(i)+str(j): sp.expand(
        sum(riemann(r, i, r, j) for r in range(4))-ricci(i, j)) for i in range(4) for j in range(4)}
    trace0 = plus.temporal[0]+3*plus.spatial[0]
    actual = sum(plus.covariant((a, a), i, i) for a in range(4) for i in range(4))
    expected = plus.time_derivative(plus.time_derivative(trace0))+3*H[0]*plus.time_derivative(trace0)-k**2*trace0
    out["covariant_trace_Laplacian_in_noncoordinate_frame"] = sp.expand(actual-expected)
    out["spatial_covariant_mass_derivative_retains_connection"] = sp.expand(
        plus.covariant((2,), 0, 2)-H[0]*(plus.temporal[0]-plus.spatial[0]))
    for a, b in combinations(range(4), 2):
        out["covariant_mass_commutator_"+str(a)+str(b)] = sp.ImmutableMatrix([[sp.expand(
            plus.covariant((a, b), i, j)-plus.covariant((b, a), i, j)
            -sum(riemann(a, b, i, r)*plus.covariant((), r, j)
                 +riemann(a, b, j, r)*plus.covariant((), i, r) for r in range(4)))
            for j in range(4)] for i in range(4)])
    return out
