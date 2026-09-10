"""Exact dimension-symbolic Dirac traces and two-loop quadratic slope tensors."""

from functools import cache
from itertools import product

import sympy as s

x, y, z, A, B, C, d, e = s.symbols("x y z A B C d epsilon")
I = s.I


@cache
def pairings(word):
    if not word:
        return ((1, ()),)
    result = []
    for j in range(1, len(word)):
        rest = word[1:j] + word[j + 1 :]
        for sign, pairs in pairings(rest):
            result.append(((-1) ** (j + 1) * sign, ((word[0], word[j]),) + pairs))
    return tuple(result)


def contract(pairs):
    pairs = list(pairs)
    result = 1
    for index in ("mu", "nu"):
        hits = [i for i, p in enumerate(pairs) if index in p]
        if not hits:
            continue
        if len(hits) == 1:
            assert pairs[hits[0]] == (index, index)
            result *= d
            pairs.pop(hits[0])
            continue
        assert len(hits) == 2, (index, pairs)
        i, j = hits
        other = lambda p, current=index: p[1] if p[0] == current else p[0]
        new = (other(pairs[i]), other(pairs[j]))
        pairs.pop(j)
        pairs.pop(i)
        pairs.append(new)
    for a, b in pairs:
        result *= x if a == b == "k" else y if a == b == "l" else z
        assert a in ("k", "l") and b in ("k", "l")
    return result


@cache
def traceword(word):
    if len(word) % 2:
        return 0
    return 4 * s.expand(sum(sign * contract(pairs) for sign, pairs in pairings(word)))


def trace(*factors):
    total = 0
    for choices in product(*(list(f.items()) for f in factors)):
        word = sum((w for w, c in choices), ())
        coeff = s.prod(c for w, c in choices)
        total += coeff * traceword(word)
    return s.expand(total)


Sk = {(): 1 / A, ("k",): -I / A}
Sl = {(): 1 / B, ("l",): -I / B}
Tk = {
    (): -1 / A**2 + 4 * x / (d * A**3),
    ("k",): I / A**2 - 4 * I * x / (d * A**3) + 2 * I / (d * A**2),
}
Tl = {
    (): -1 / B**2 + 4 * y / (d * B**3),
    ("l",): I / B**2 - 4 * I * y / (d * B**3) + 2 * I / (d * B**2),
}
mu = {("mu",): 1}
nu = {("nu",): 1}
scalar_self = 2 * trace(Sk, Sl, Sk, Tk)
scalar_vertex = (
    trace(Sk, Tk, Sl, Sl)
    + trace(Sk, Sk, Tl, Sl)
    - trace(Sk, Sk, nu, Sk, Sl, nu, Sl, Sl) / d
)
gauge_self = 2 * trace(Sk, mu, Sl, mu, Sk, Tk)
gauge_vertex = (
    trace(mu, Sk, Tk, mu, Sl, Sl)
    + trace(mu, Sk, Sk, mu, Tl, Sl)
    - trace(mu, Sk, Sk, nu, Sk, mu, Sl, nu, Sl, Sl) / d
)
