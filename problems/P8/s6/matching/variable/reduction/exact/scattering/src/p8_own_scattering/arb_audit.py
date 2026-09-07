"""Independent univariate Arb jets of the implicit stationary branch.

No primary Fraction jet, symbolic stationary module or interval engine is
imported here. Formal coefficients are factorial-normalized Taylor
coefficients. The implicit unknown is solved one coefficient at a time,
not by importing partial derivatives from the other implementation.
"""
from fractions import Fraction

from flint import arb, ctx


def rational(value):
    if type(value) is not int and not isinstance(value, Fraction):
        raise TypeError("An exact int/Fraction is required")
    value = Fraction(value)
    return arb(value.numerator)/value.denominator


def ball(lower, upper):
    if type(lower) is not int and not isinstance(lower, Fraction):
        raise TypeError("Exact rational endpoints are required")
    if type(upper) is not int and not isinstance(upper, Fraction):
        raise TypeError("Exact rational endpoints are required")
    lo, hi = Fraction(lower), Fraction(upper)
    if lo > hi:
        raise ValueError("Reversed ball endpoints")
    return arb(rational((lo+hi)/2), rational((hi-lo)/2))


class Jet:
    """Fixed-order real-ball Taylor algebra, with no silent order changes."""

    def __init__(self, coefficients):
        self.a = tuple(value if isinstance(value, arb) else rational(value) for value in coefficients)
        if not self.a:
            raise ValueError("At least a constant coefficient is required")

    @property
    def order(self):
        return len(self.a)-1

    def lift(self, other):
        if isinstance(other, Jet):
            if other.order != self.order:
                raise ValueError("Different Taylor orders require explicit truncation")
            return other
        return Jet((other, *([0]*self.order)))

    def __add__(self, other):
        other = self.lift(other)
        return Jet(tuple(a+b for a, b in zip(self.a, other.a, strict=True)))

    __radd__ = __add__

    def __neg__(self):
        return Jet(tuple(-value for value in self.a))

    def __sub__(self, other):
        return self+-self.lift(other)

    def __rsub__(self, other):
        return self.lift(other)+-self

    def __mul__(self, other):
        other = self.lift(other)
        return Jet(tuple(sum((self.a[j]*other.a[n-j] for j in range(n+1)), arb(0))
                         for n in range(self.order+1)))

    __rmul__ = __mul__

    def inverse(self):
        if not (self.a[0] > 0 or self.a[0] < 0):
            raise ZeroDivisionError("The Taylor constant must exclude zero")
        values = [1/self.a[0]]
        for n in range(1, self.order+1):
            values.append(-values[0]*sum((self.a[j]*values[n-j] for j in range(1, n+1)), arb(0)))
        return Jet(values)

    def __truediv__(self, other):
        return self*self.lift(other).inverse()

    def __rtruediv__(self, other):
        return self.lift(other)*self.inverse()

    def __pow__(self, exponent):
        if type(exponent) is not int:
            raise TypeError("An integer Taylor power is required")
        if exponent < 0:
            return self.inverse()**(-exponent)
        result = self.lift(1)
        for _ in range(exponent):
            result *= self
        return result

    def root3(self):
        if not self.a[0] > 0:
            raise ValueError("A positive constant is required for the real cube-root jet")
        values = [self.a[0].root(3)]
        for n in range(1, self.order+1):
            known = sum((values[i]*values[j]*values[n-i-j]
                         for i in range(n) for j in range(n)
                         if 0 <= n-i-j < n), arb(0))
            values.append((self.a[n]-known)/(3*values[0]**2))
        return Jet(values)

    def truncate(self, order):
        if type(order) is not int:
            raise TypeError("An integer Taylor order is required")
        if not 0 <= order <= self.order:
            raise ValueError("Truncation cannot create unavailable coefficients")
        return Jet(self.a[:order+1])

    def derivative(self):
        if self.order == 0:
            raise ValueError("No derivative coefficient remains")
        return Jet(tuple(n*self.a[n] for n in range(1, self.order+1)))


def formulas(v, zeta, delta):
    c, d = 2+delta, 1+v
    jj = c*d**4*(1-7*v)+12*v
    jv = c*d**3*(-3-35*v)+12
    gap = c-2/d**4
    q = 32*(1-v)/(c*d**14)
    ratio = d**8*jj/(8*c*(1-v))
    log_r = 8/d+jv/jj+1/(1-v)
    hh = gap*(jv/jj-6/d)-8/d**5
    one_minus = 1-v*gap*zeta
    b = (one_minus/ratio).root3()
    ff = 2*(v*zeta*hh-log_r)/3
    fixed = 3*b*ff**2/(2*q)
    return {"b": b, "F": ff, "T": fixed, "Q": q, "D": gap, "H": hh, "one_minus": one_minus}


def domain():
    return {"v": (Fraction(0), Fraction(1, 10000)),
            "delta": (Fraction(0), Fraction(1, 100)),
            "zeta": (Fraction(11), Fraction(13))}


def enclosure(precision=256, box=None):
    if type(precision) is not int:
        raise TypeError("Precision must be an actual integer")
    if precision < 128:
        raise ValueError("At least 128 bits are required")
    reference = domain()
    if box is None:
        box = reference
    if not isinstance(box, dict) or set(box) != set(reference):
        raise ValueError("The box must specify v, delta and zeta")
    for name, pair in box.items():
        if not isinstance(pair, tuple) or len(pair) != 2:
            raise TypeError("Each box axis is a pair of rational endpoints")
        for endpoint in pair:
            if type(endpoint) is not int and not isinstance(endpoint, Fraction):
                raise TypeError("Only exact endpoints are admitted")
        if not reference[name][0] <= pair[0] <= pair[1] <= reference[name][1]:
            raise ValueError("The requested box lies outside the proved domain")
    with ctx.workprec(precision):
        v0, z0, delta = (ball(*box[name]) for name in ("v", "zeta", "delta"))
        v = Jet((v0, 1, 0, 0))
        z = [z0, arb(0), arb(0), arb(0)]
        base = formulas(v, Jet(z), delta)
        F0, H0, gap0, one_minus = (base[name].a[0] for name in ("F", "H", "D", "one_minus"))
        partial_z = base["T"].a[0]*(-v0*gap0/(3*one_minus)+2*(2*v0*H0/3)/F0)
        divisor = 1-partial_z
        if not divisor > 0:
            raise ValueError("The implicit derivative denominator must be positive")
        for n in (1, 2, 3):
            known = formulas(v, Jet(z), delta)["T"].a[n]
            z[n] = known/divisor
        final = formulas(v, Jet(z), delta)
        b, ff, q = (final[name] for name in ("b", "F", "Q"))
        bprime = b.derivative()
        lapse = 2*bprime/ff.truncate(2)
        mass_numerator = 2*q.truncate(2)*lapse/(b**2).truncate(2)
        kinetic = (b**3).truncate(2)/lapse
        ell = kinetic.a[1]/kinetic.a[0]
        ell_v = 2*kinetic.a[2]/kinetic.a[0]-ell**2
        pump = ell+2*v0*ell_v+v0*ell**2
        return {"precision_bits": precision, "zeta": tuple(z), "T_zeta": partial_z,
                "b": b.a, "F": ff.a, "N": lapse.a, "k": kinetic.a,
                "A": mass_numerator.a, "A_v": mass_numerator.a[1],
                "ell": ell, "ell_v": ell_v, "pump": pump}


def cover_boxes():
    """An exact eight-box cover, not a collection of sample points."""
    boxes = [{}]
    for name, (lo, hi) in domain().items():
        middle = (lo+hi)/2
        boxes = [{**box, name: interval} for box in boxes for interval in ((lo, middle), (middle, hi))]
    return boxes


def report():
    records = []
    for box in cover_boxes():
        data = enclosure(box=box)
        checks = {"mass_numerator_derivative_below129": bool(abs(data["A_v"]) < 129),
                  "canonical_pump_below22": bool(abs(data["pump"]) < 22),
                  "implicit_denominator_positive": bool(1-data["T_zeta"] > 0),
                  "scale_derivative_negative": bool(data["b"][1] < 0),
                  "kinetic_positive": bool(data["k"][0] > 0)}
        if not all(checks.values()):
            raise ValueError("An independent Arb derivative-cover comparison failed")
        records.append({"box": {name: [str(lo), str(hi)] for name, (lo, hi) in box.items()}, "checks": checks})
    return {"precision_bits": 256, "exact_partition_box_count": 8,
            "whole_domain_not_sampled_grid": True, "independent_univariate_implicit_recurrence": True,
            "boxes": records}
