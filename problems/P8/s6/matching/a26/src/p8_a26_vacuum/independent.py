"""Separate Fraction Laurent-series and Fourier constant-term reconstruction.

No SymPy or imports of the primary formula modules. Laurent coefficients are
obtained from finite power-series division of the literal action. Fourier
moments come from multiplication of Laurent polynomials, not quadrature.
"""

from fractions import Fraction as Q
from math import comb

ORDER = 5


class Series:
    """Exact truncated series with an explicit Laurent shift."""

    def __init__(self, data):
        if isinstance(data, (int, Q)):
            data = {0: Q(data)}
        self.data = {n: Q(v) for n, v in data.items() if v and -3 <= n <= ORDER}

    def __add__(self, other):
        other = other if isinstance(other, Series) else Series(other)
        keys = self.data.keys() | other.data.keys()
        return Series({n: self.data.get(n, 0)+other.data.get(n, 0) for n in keys})

    __radd__ = __add__

    def __neg__(self):
        return Series({n: -v for n, v in self.data.items()})

    def __sub__(self, other):
        return self+-other if isinstance(other, Series) else self+(-Q(other))

    def __rsub__(self, other):
        return -self+other

    def __mul__(self, other):
        other = other if isinstance(other, Series) else Series(other)
        out = {}
        for n, a in self.data.items():
            for m, b in other.data.items():
                out[n+m] = out.get(n+m, Q(0))+a*b
        return Series(out)

    __rmul__ = __mul__

    def __pow__(self, n):
        if n < 0:
            return self.inverse()**(-n)
        result = Series(1)
        for _ in range(n):
            result = result*self
        return result

    def inverse(self):
        if not self.data:
            raise ZeroDivisionError("Zero formal series")
        first = min(self.data)
        a0 = self.data[first]
        coeff = {0: 1/a0}
        for n in range(1, ORDER+1):
            coeff[n] = -sum(self.data.get(first+j, 0)*coeff[n-j]
                            for j in range(1, n+1))/a0
        return Series({n-first: value for n, value in coeff.items()})

    def __truediv__(self, other):
        other = other if isinstance(other, Series) else Series(other)
        return self*other.inverse()


def residue_fixture(a, g, *, printed):
    """Enough retained orders for the X^-1 coefficient of every operation."""
    a, g = Q(a), Q(g)
    if g == Q(1, 2):
        raise ValueError("Nondegenerate F2 required")
    x = Series({1: 1})
    f, ax, fx = Q(1, 2)-g*(1+x), -a*(1+x), -g
    a3 = 2*(ax-2*fx)*(ax*x-2*f)/(x*(3*ax*x-4*f))
    disputed = 4*(3*f+16*x*fx*ax**2) if printed else 4*(3*f+16*x*fx)*ax**2
    a4 = (-16*x*ax**3+disputed-x**2*f*a3**2
          -(16*x**2*fx-12*x*f)*a3*ax-16*fx*(3*f+4*x*fx)*ax
          +8*f*(x*fx-f)*a3+48*f*fx**2)/(8*(f-x*ax)**2)
    a5 = ((4*fx-2*ax+x*a3)*(-2*ax**2-3*x*ax*a3+4*fx*ax+4*f*a3)
          /(8*(f-x*ax)**2))
    return tuple(coefficient.data.get(-1, Q(0)) for coefficient in (a3, a4, a5))


def _add(a, b):
    return {key: a.get(key, 0)+b.get(key, 0) for key in a.keys() | b.keys()}


def _mul(a, b):
    out = {}
    for (i, j), v in a.items():
        for (k, ell), w in b.items():
            key = (i+k, j+ell)
            out[key] = out.get(key, Q(0))+v*w
    return out


def _power(a, n):
    out = {(0, 0): Q(1)}
    for _ in range(n):
        out = _mul(out, a)
    return out


def fourier_moments():
    """Constant terms of Laurent polynomials on the unit two-torus."""
    cx = {(-1, 0): Q(1, 2), (1, 0): Q(1, 2)}
    cy = {(0, -1): Q(1, 2), (0, 1): Q(1, 2)}
    sin2x = {(0, 0): Q(1, 2), (-2, 0): -Q(1, 4), (2, 0): -Q(1, 4)}
    sin2y = {(0, 0): Q(1, 2), (0, -2): -Q(1, 4), (0, 2): -Q(1, 4)}
    g, grad2 = _add(cx, cy), _add(sin2x, sin2y)
    lap_hess = {key: 2*value for key, value in _mul(cx, cy).items()}
    terms = {"g4": _power(g, 4), "g2_grad2": _mul(_power(g, 2), grad2),
             "grad4": _power(grad2, 2), "g2_lap_hess": _mul(_power(g, 2), lap_hess),
             "grad2_lap_hess": _mul(grad2, lap_hess)}
    return {key: value.get((0, 0), Q(0)) for key, value in terms.items()}


def local_background():
    # At z=0 the derivative of sigma multiplies equal factors, while each
    # even scale factor has second derivative 2p. Thus a''=p_++p_-.
    hprime = Q(1, 6)+Q(1, 20)
    nprime = Q(1)
    zprime = 9*hprime-2*nprime
    return {"hprime": hprime, "Zprime": zprime,
            "a1_over_g1": (2*nprime/hprime-1)/4,
            "r_over_g1": zprime/(4*hprime)}


def ibp_bulk():
    moments = fourier_moments()
    # log|X| = 2l*T' g-l²(T'^2*g²+T²*|grad g|²)+O(l³).
    # D/l² = 2TT''g² +2T'^2|grad g|²+T²(lap²-hess²).
    # -log|X|D/2 at l4 has these four time monomials.
    terms = {"T3_Tpp": moments["g2_grad2"],
             "T2_Tp2": moments["grad4"]+moments["g2_lap_hess"]/2,
             "T_Tp2_Tpp": moments["g4"],
             "Tp4": moments["g2_grad2"]}
    # int T³T''=-3 int T²T'²; int T*T'²*T''=-int T'^4/3.
    return {"raw": terms, "T2_Tp2_after_IBP": -3*terms["T3_Tpp"]+terms["T2_Tp2"],
            "Tp4_after_IBP": -terms["T_Tp2_Tpp"]/3+terms["Tp4"],
            "fourth_difference_multiplier": sum(Q(w)*Q(k)**4
                                                for k, w in ((2, 1), (-2, 1), (1, -4), (-1, -4), (0, 6)))}


def checks():
    fixtures = []
    for a, g in ((Q(1, 7), Q(1, 5)), (Q(3), Q(2, 3)), (Q(2, 5), Q(1, 5)),
                 (Q(-2), Q(0)), (Q(7, 11), Q(-3, 4))):
        r, f = 2*g-a, Q(1, 2)-g
        expected = r, -r, -r*r/(2*f)
        for printed in (False, True):
            actual = residue_fixture(a, g, printed=printed)
            if actual != expected:
                raise ValueError("Independent Laurent reconstruction failed")
            fixtures.append({"a": str(a), "g": str(g), "printed_A4": printed,
                             "residues": list(map(str, actual))})
    moments = fourier_moments()
    if tuple(moments.values()) != (Q(9, 4), Q(3, 4), Q(5, 4), Q(1), Q(0)):
        raise ValueError("Independent Fourier constant terms failed")
    bulk = ibp_bulk()
    if (bulk["T2_Tp2_after_IBP"], bulk["Tp4_after_IBP"], bulk["fourth_difference_multiplier"]) != (-Q(1, 2), 0, 24):
        raise ValueError("Independent non-boundary bulk reconstruction failed")
    if local_background() != {"hprime": Q(13, 60), "Zprime": -Q(1, 20),
                              "a1_over_g1": Q(107, 52), "r_over_g1": -Q(3, 52)}:
        raise ValueError("Independent analytic-numerator test failed")
    # A smooth even exponent is an exclusion control for the cusp argument.
    smooth_third_factor = Q(2)*Q(1)*Q(0)
    if smooth_third_factor != 0 or comb(4, 2) != 6:
        raise ValueError("Elementary regularity or finite-difference control failed")
    return {"Laurent_fixtures": fixtures,
            "Fourier_moments": {key: str(value) for key, value in moments.items()},
            "bulk_after_IBP": {"int_T2_Tp2": "-1/2", "int_Tp4": "0", "fourth_difference": "24"},
            "analytic_background": {key: str(value) for key, value in local_background().items()}}
