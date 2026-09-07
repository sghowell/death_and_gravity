"""Independent Fraction intervals, bivariate jets, and norm reconstruction.

No SymPy or primary-module imports. The c derivatives are computed by
truncated Taylor arithmetic applied directly to the normalization formulas,
not by reading a primary derivative or rounded sampling. Positive square
roots are enclosed using integer arithmetic. All box subdivisions cover the
entire closed box, including u=0 and c=2 for the ANALYTIC coefficients only.
"""

from fractions import Fraction as Q
from functools import cache
from math import factorial, isqrt


def rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, Q)):
        raise TypeError("Independent arithmetic requires int or Fraction")
    return Q(value)


def upper_round(value):
    """Outward exact rounding only, to keep the replay manifest readable."""
    scaled = value*10**6
    return Q(-(-scaled.numerator//scaled.denominator), 10**6)


class Interval:
    def __init__(self, lower, upper=None):
        self.lo = rational(lower)
        self.hi = self.lo if upper is None else rational(upper)
        if self.lo > self.hi:
            raise ValueError("Interval endpoints are reversed")

    def __add__(self, other):
        other = other if isinstance(other, Interval) else Interval(other)
        return Interval(self.lo+other.lo, self.hi+other.hi)

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        return self+-other if isinstance(other, Interval) else self+(-rational(other))

    def __rsub__(self, other):
        return -self+other

    def __mul__(self, other):
        other = other if isinstance(other, Interval) else Interval(other)
        products = [a*b for a in (self.lo, self.hi) for b in (other.lo, other.hi)]
        return Interval(min(products), max(products))

    __rmul__ = __mul__

    def inverse(self):
        if self.lo <= 0 <= self.hi:
            raise ValueError("A reciprocal interval must exclude zero")
        return Interval(1/self.hi, 1/self.lo)

    def __truediv__(self, other):
        other = other if isinstance(other, Interval) else Interval(other)
        return self*other.inverse()

    def __rtruediv__(self, other):
        return self.inverse()*other

    def __pow__(self, exponent):
        exponent = rational(exponent)
        if exponent.denominator == 2:
            return self.sqrt()**exponent.numerator
        if exponent.denominator != 1:
            raise ValueError("Only integer and half-integer interval powers are needed")
        n = exponent.numerator
        if n < 0:
            return self.inverse()**(-n)
        if n == 0:
            return Interval(1)
        values = [self.lo**n, self.hi**n]
        lower = 0 if n % 2 == 0 and self.lo <= 0 <= self.hi else min(values)
        return Interval(lower, max(values))

    def sqrt(self):
        if self.lo < 0:
            raise ValueError("A real square-root interval must be nonnegative")

        def bracket(value):
            # A fixed rational refinement avoids a coarse root enclosure for
            # integers. floor(sqrt(n/d)*10^12) is exact via integer square root.
            scale = 10**12
            floor = isqrt(value.numerator*value.denominator*scale**2)
            denom = value.denominator*scale
            lower = Q(floor, denom)
            upper = lower if floor**2 == value.numerator*value.denominator*scale**2 else Q(floor+1, denom)
            return lower, upper

        return Interval(bracket(self.lo)[0], bracket(self.hi)[1])

    def absolute_upper(self):
        return max(abs(self.lo), abs(self.hi))

    def pair(self):
        return self.lo, self.hi


class Jet:
    """Taylor coefficients u^i*c^j, 0<=i,j<=1, over exact intervals."""

    def __init__(self, value, *, variable=None):
        if isinstance(value, dict):
            self.data = value
        else:
            self.data = {(0, 0): value if isinstance(value, Interval) else Interval(value)}
            if variable is not None:
                if variable not in ((1, 0), (0, 1)):
                    raise ValueError("The only two independent jet variables are u,c")
                self.data[variable] = Interval(1)

    def get(self, index):
        return self.data.get(index, Interval(0))

    def __add__(self, other):
        other = other if isinstance(other, Jet) else Jet(other)
        return Jet({key: self.get(key)+other.get(key) for key in self.data.keys() | other.data.keys()})

    __radd__ = __add__

    def __neg__(self):
        return Jet({key: -value for key, value in self.data.items()})

    def __sub__(self, other):
        return self+-other if isinstance(other, Jet) else self+(-Jet(other))

    def __rsub__(self, other):
        return -self+other

    def __mul__(self, other):
        other = other if isinstance(other, Jet) else Jet(other)
        result = {}
        for (i, j), value in self.data.items():
            for (k, ell), item in other.data.items():
                key = i+k, j+ell
                if key[0] <= 1 and key[1] <= 1:
                    result[key] = result.get(key, Interval(0))+value*item
        return Jet(result)

    __rmul__ = __mul__

    def __pow__(self, exponent):
        exponent = rational(exponent)
        if exponent.denominator == 1 and exponent >= 0:
            out = Jet(1)
            for _ in range(exponent.numerator):
                out = out*self
            return out
        base = self.get((0, 0))
        nilpotent = Jet({key: value for key, value in self.data.items() if key != (0, 0)})
        # nilpotent^3=0 at these bidegrees. Binomial coefficients are exact.
        return (Jet(base**exponent)+nilpotent*(base**(exponent-1))*exponent
                +nilpotent*nilpotent*(base**(exponent-2))*exponent*(exponent-1)*Q(1, 2))

    def __truediv__(self, other):
        other = other if isinstance(other, Jet) else Jet(other)
        return self*other**-1

    def __rtruediv__(self, other):
        return self**-1*other

    def u_derivative(self):
        # Only the u^0 coefficients are used after this operation; they include
        # the exact mixed u,c coefficient. No nonexistent u^2 jet is inferred.
        return Jet({(0, j): self.get((1, j)) for j in (0, 1)})


def profiles(u_box, c_box, k_box):
    u, c, k = Jet(u_box, variable=(1, 0)), Jet(c_box, variable=(0, 1)), Jet(k_box)
    d = 1+u*u
    denom = c*d**12+8
    w1, w2 = c*d**12/denom, 8/denom
    theta = 6*u*(c*d**12-8)/(d*denom)
    omega = -24*(2*c)**Q(1, 2)*u*d**5/denom
    dcross = 48*(2*c)**Q(1, 2)*d**5/denom
    ns, nr = theta.u_derivative()+theta*theta, -theta.u_derivative()+theta*theta
    q2, cf2 = k*k/d**4, c*c*d**8/4
    dg = (8*c)**Q(1, 2)*d**6/denom*(cf2-1)
    result = {
        "A": q2*(w1+cf2*w2)-ns,
        "C": q2*dg-2*omega.u_derivative()-2*omega*theta,
        "d_cross": dcross, "f_cross": -dcross,
        "E": q2*dg-2*omega*theta,
        "b_analytic": q2*(w2+cf2*w1)-nr-4*omega*omega,
        "j_light": c**Q(1, 2)*d**9/denom**Q(1, 2),
        "j_heavy": -2*Jet(2)**Q(1, 2)*d**3/denom**Q(1, 2),
    }
    return result


@cache
def analytic_boxes():
    """Closed cover of |u|<=.1, 2<=c<=2.01, 1<=kbar<=2.

    This is an interval/jet proof, not a mesh of point evaluations. Subdivision
    only reduces dependency overestimation. Pole-containing B is absent.
    """
    count = 40
    keys = ("A", "C", "d_cross", "f_cross", "E", "b_analytic", "j_light", "j_heavy")
    derivatives = dict.fromkeys(keys, Q(0))
    sources = dict.fromkeys(("j_light", "j_heavy"), Q(0))
    for index in range(count):
        lo, hi = Q(-1, 10)+Q(index, 5*count), Q(-1, 10)+Q(index+1, 5*count)
        data = profiles(Interval(lo, hi), Interval(2, Q(201, 100)), Interval(1, 2))
        for key in keys:
            derivatives[key] = max(derivatives[key], data[key].get((0, 1)).absolute_upper())
        for key in sources:
            sources[key] = max(sources[key], data[key].get((0, 0)).absolute_upper())
    for key in keys[:6]:
        if derivatives[key] >= 8:
            raise ValueError(f"The independently differentiated analytic box exceeded 8: {key}")
    for key in sources:
        if sources[key] >= 1 or derivatives[key] >= 1:
            raise ValueError(f"The physical source projection or derivative exceeded 1: {key}")
    return {"u_subintervals": count, "closed_u_box": (Q(-1, 10), Q(1, 10)),
            "closed_c_box": (Q(2), Q(201, 100)), "closed_k_box": (Q(1), Q(2)),
            "analytic_c_derivative_abs_upper": {key: upper_round(value) for key, value in derivatives.items()},
            "physical_source_projection_abs_upper": {key: upper_round(value) for key, value in sources.items()}}


def norm_constants():
    mass = Q(126)+80*(Q(41, 10)/32+Q(121, 10)/64)
    outer = (29+(35+14*Q(13, 4))/3+Q(200, 3)/10+Q(160, 3)/300+Q(14, 3)/30)
    difference = 1+8*Q(1, 10)**3+(12+8*Q(13, 4)+Q(8, 3))/3000+Q(8, 3)/30000+Q(8, 3)/10000
    inner = 29+(35+14*Q(13, 4))*Q(3, 8)+Q(200, 3)/8+Q(1280, 3)*Q(3, 512)+Q(14, 3)*Q(3, 64)
    middle = 128*2*4*2*2+Q(1, 8)*4*2*2+1152*2*2+Q(1, 8)*2*2+128*2
    return {"mass_remainder_upper": mass, "outer_remainder_over_radius_upper": outer,
            "outer_difference_over_delta_rminus2_upper": difference,
            "inner_remainder_over_radius_upper": inner,
            "exp_two_thirds_upper": Q(41, 21), "middle_telescoping_coefficient": middle,
            "full_telescoping_coefficient": Q(4*9000+48+2),
            "source_coefficient_before_enlargement": Q(8+160000*16),
            "nontrivial_adapted_transfer_error_upper": Q(1, 250)}


def checks():
    constants = norm_constants()
    caps = {"mass_remainder_upper": 152, "outer_remainder_over_radius_upper": 64,
            "outer_difference_over_delta_rminus2_upper": 2,
            "inner_remainder_over_radius_upper": 72, "exp_two_thirds_upper": 2,
            "middle_telescoping_coefficient": 9000, "full_telescoping_coefficient": 40000,
            "source_coefficient_before_enlargement": 3000000}
    if any(constants[key] >= cap for key, cap in caps.items()):
        raise ValueError("An independent strict rational norm comparison failed")
    # e^9>sum_0^8 9^n/n!>2002 implies 0<csch(pi*mu)<1/1000,
    # using pi>3,mu>3 and exp(-9)<1. This is not a decimal Gamma fit.
    lower_exp9 = sum(Q(9)**n/factorial(n) for n in range(9))
    if lower_exp9 <= 2002:
        raise ValueError("The elementary nonzero-reflection enclosure failed")
    return {"norms": constants, "analytic_interval_jet_cover": analytic_boxes(),
            "exp9_lower_taylor": lower_exp9,
            "strict_reflection_modulus_upper": Q(1, 1000)}
