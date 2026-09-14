"""Exact finite exterior jet identities; not a spatial truncation regulator."""

from itertools import product

import sympy as s


class G:
    def __init__(self, terms=None):
        self.terms = {m: s.expand(a) for m, a in (terms or {}).items() if a != 0}
        self.terms = {m: a for m, a in self.terms.items() if a != 0}

    @staticmethod
    def cast(a):
        return a if isinstance(a, G) else G({(): s.sympify(a)})

    def __add__(self, b):
        b = G.cast(b)
        out = dict(self.terms)
        for m, a in b.terms.items():
            out[m] = out.get(m, 0) + a
        return G(out)

    __radd__ = __add__

    def __neg__(self):
        return G({m: -a for m, a in self.terms.items()})

    def __sub__(self, b):
        return self + -G.cast(b)

    def __rsub__(self, b):
        return G.cast(b) + -self

    def __mul__(self, b):
        b = G.cast(b)
        out = {}
        for m, a in self.terms.items():
            for n, z in b.terms.items():
                if set(m) & set(n):
                    continue
                sign = (-1) ** sum(i > j for i in m for j in n)
                key = tuple(sorted(m + n))
                out[key] = out.get(key, 0) + sign * a * z
        return G(out)

    __rmul__ = __mul__

    def __truediv__(self, b):
        return self * s.sympify(b) ** -1

    def __bool__(self):
        return bool(self.terms)

    def count(self):
        return len(self.terms)

    def data(self):
        return {str(m): str(a) for m, a in self.terms.items()}


class Jets:
    def __init__(self):
        self.symbols = {}
        self.reverse = {}
        self.rules = {}
        self.odd_rules = {}
        self._sg_cache = {}

    def field(self, name, indices=()):
        key = (name, tuple(sorted(indices)))
        if key not in self.symbols:
            symbol = s.Symbol(name + "".join("_" + str(i) for i in key[1]), real=True)
            assert symbol not in self.reverse
            self.symbols[key] = symbol
            self.reverse[symbol] = key
        return G.cast(self.symbols[key])

    def ghost(self, i, indices=()):
        if i == 0:
            return G()
        return G({((i, tuple(sorted(indices))),): s.S.One})

    def d(self, a, mu):
        a = G.cast(a)
        out = G()
        for mon, coef in a.terms.items():
            for sym in coef.free_symbols:
                if sym in self.reverse:
                    name, indices = self.reverse[sym]
                    out += G({mon: s.diff(coef, sym)}) * self.field(
                        name, indices + (mu,)
                    )
            for pos, (i, indices) in enumerate(mon):
                out += (
                    G({mon[:pos]: coef})
                    * self.ghost(i, indices + (mu,))
                    * G({mon[pos + 1 :]: s.S.One})
                )
        return out

    def derivatives(self, a, indices):
        for mu in indices:
            a = self.d(a, mu)
        return a

    def sg(self, i, indices=()):
        key = (i, indices, self.odd_rules.get(i))
        if key in self._sg_cache:
            return self._sg_cache[key]
        if i in self.odd_rules:
            out = self.odd_rules[i]
        elif i in (1, 2, 3):
            out = sum((self.ghost(j) * self.ghost(i, (j,)) for j in range(1, 4)), G())
        else:
            raise KeyError("No BRST rule for odd generator " + str(i))
        result = self.derivatives(out, indices)
        self._sg_cache[key] = result
        return result

    def brst(self, a):
        a = G.cast(a)
        out = G()
        for mon, coef in a.terms.items():
            for sym in coef.free_symbols:
                if sym in self.reverse:
                    name, indices = self.reverse[sym]
                    out += self.derivatives(self.rules[name], indices) * G(
                        {mon: s.diff(coef, sym)}
                    )
            for pos, (i, indices) in enumerate(mon):
                out += (
                    (-1) ** pos
                    * G({mon[:pos]: coef})
                    * self.sg(i, indices)
                    * G({mon[pos + 1 :]: s.S.One})
                )
        return out

    def transport(self, a):
        return sum((self.ghost(j) * self.d(a, j) for j in range(1, 4)), G())

    def divergence_ghost(self):
        return sum((self.ghost(j, (j,)) for j in range(1, 4)), G())

    def add_tensor(self, prefix, rank, variance, weight=0, symmetric=False):
        components = {}
        for idx in product(range(1, 4), repeat=rank):
            key = tuple(sorted(idx)) if symmetric else idx
            name = prefix + "".join(str(i) for i in key)
            components[idx] = self.field(name)
        for idx, a in components.items():
            rule = self.transport(a) + weight * self.divergence_ghost() * a
            for slot, typ in enumerate(variance):
                for j in range(1, 4):
                    changed = idx[:slot] + (j,) + idx[slot + 1 :]
                    rule += (
                        typ
                        * self.ghost(
                            j if typ == 1 else idx[slot],
                            (idx[slot] if typ == 1 else j,),
                        )
                        * components[changed]
                    )
            key = tuple(sorted(idx)) if symmetric else idx
            self.rules[prefix + "".join(str(i) for i in key)] = rule
        return components

    def add_affine(self, compensated=False):
        comp = {
            idx: self.field("Gamma" + "".join(map(str, idx)))
            for idx in product(range(4), repeat=3)
        }
        for (a, b, c), value in comp.items():
            rule = self.transport(value) + self.ghost(a, (b, c))
            for d in range(4):
                rule += (
                    -self.ghost(a, (d,)) * comp[d, b, c]
                    + self.ghost(d, (b,)) * comp[a, d, c]
                    + self.ghost(d, (c,)) * comp[a, b, d]
                )
            if compensated and a == b:
                rule -= self.d(self.divergence_ghost(), c) / 4
            self.rules["Gamma" + str(a) + str(b) + str(c)] = rule
        return comp
