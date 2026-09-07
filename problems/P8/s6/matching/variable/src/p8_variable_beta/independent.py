"""Independent Fraction/Taylor/coefficient engine; no SymPy or primary imports.

Only the immutable sparse polynomial primitive is reused. Taylor jets are
implemented here, including their derivative and inversion combinatorics.
"""

from dataclasses import dataclass
from fractions import Fraction as Q

from p8_composite_modes.independent import NAMES, Poly


def exact(value):
    if isinstance(value, bool) or not isinstance(value, (int, Q)):
        raise TypeError("Exact int/Fraction arithmetic required")
    return Q(value)


@dataclass(frozen=True)
class Jet:
    """Taylor coefficients through order two, not divided derivative guesses."""
    value: Q
    first: Q = Q(0)
    second_coefficient: Q = Q(0)

    def __post_init__(self):
        for name in ("value", "first", "second_coefficient"):
            object.__setattr__(self, name, exact(getattr(self, name)))

    def __add__(self, other):
        other = other if isinstance(other, Jet) else Jet(other)
        return Jet(self.value+other.value, self.first+other.first, self.second_coefficient+other.second_coefficient)

    __radd__ = __add__

    def __neg__(self):
        return Jet(-self.value, -self.first, -self.second_coefficient)

    def __sub__(self, other):
        return self+-other

    def __rsub__(self, other):
        return other+-self

    def __mul__(self, other):
        other = other if isinstance(other, Jet) else Jet(other)
        return Jet(self.value*other.value, self.first*other.value+self.value*other.first,
                   self.second_coefficient*other.value+self.first*other.first+self.value*other.second_coefficient)

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = other if isinstance(other, Jet) else Jet(other)
        if not other.value:
            raise ZeroDivisionError("An invertible constant Taylor coefficient is required")
        reciprocal = Jet(1/other.value, -other.first/other.value**2,
                         other.first**2/other.value**3-other.second_coefficient/other.value**2)
        return self*reciprocal

    def __rtruediv__(self, other):
        return Jet(other)/self

    def __pow__(self, power):
        if type(power) is not int or power < 0:
            raise TypeError("An exact nonnegative integer power is required")
        value = Jet(1)
        for _ in range(power):
            value *= self
        return value


def family(u, c):
    """Independent raw arithmetic construction of the actual profiles."""
    d = 1+u*u
    a, b = d*d, 2/(d*d)
    y = b/a
    h, hp = 4*u/d, 4*(1-u*u)/(d*d)
    chi = 1/(10*a**3)
    n = 2*(y**3/c-1)*hp
    k = n-chi*chi
    b1 = y**3*hp/(c*(c-y))
    b4 = 3*h*h/(2*c*c)-b1/y**3
    b0 = (3*h*h-n/2)/2-3*b1*y
    return {"a": a, "b": b, "y": y, "h": h, "h_u": hp, "h_f": -h/c,
            "D_f_h_f": -hp/(c*c), "chi_speed": chi, "nbar": n, "kbar": k,
            "b0": b0, "b1": b1, "b4": b4, "P": 2*b1}


def background_fixtures():
    cases = [(Q(0), Q(1)), (Q(1, 10), Q(1)), (Q(0), Q(4)),
             (Q(1, 10), Q(4)), (Q(-1, 10), Q(5, 2)),
             (Q(1, 20), Q(3)), (Q(0), Q(201, 100))]
    result = []
    for u, c in cases:
        jets = family(Jet(u, 1), c)
        f = {key: value.value for key, value in jets.items()}
        y, h, hp, n, P = (f[key] for key in ("y", "h", "h_u", "nbar", "P"))
        U, V = f["b0"]+3*f["b1"]*y, f["b1"]+f["b4"]*y**3
        beta_force = 2*(jets["b0"].first+(c+3*y)*jets["b1"].first+c*y**3*jets["b4"].first)
        residuals = {"g_lapse": 3*h*h-n/2-2*U, "f_lapse": 3*h*h/(c*c)-2*V/y**3,
                     "g_null": -2*hp-n-(y-c)*P, "f_null": 2*hp/(c*c)-(c-y)*P/(c*y**3),
                     "chi": jets["chi_speed"].first+3*h*f["chi_speed"],
                     "clock_times_speed": jets["kbar"].first/2+3*h*f["kbar"]+beta_force,
                     "sourced_Bianchi": 3*P*(-y*h/c-h)-2*(jets["b1"].first+y**3*jets["b4"].first)}
        if any(residuals.values()) or f["kbar"] <= 0:
            raise ValueError("Independent actual metric/clock/matter fixture failed")
        result.append({"u": u, "c": c, "values": f, "actual_Euler_residuals": residuals})
    return result


def canonical_fixtures():
    result = []
    for c in map(Q, (1, Q(5, 2), 3, 4, Q(201, 100))):
        f = family(Jet(0, 1), c)
        y, a, b = f["y"], f["a"], f["b"]
        K1, K2 = a**3/8, b**3/(8*c)
        total, relative = K1+K2, K1*K2/(K1+K2)
        w1, w2 = K1/total, K2/total
        mass = y*f["P"]*(1+c/y**3)
        Ns = total.second_coefficient/total.value-total.first**2/(4*total.value**2)
        Nr = relative.second_coefficient/relative.value-relative.first**2/(4*relative.value**2)
        omega_prime_squared = w2.second_coefficient**2/(w1.value*w2.value)
        cl, ch = w1+w2*c*c/y**2, w2+w1*c*c/y**2
        D_squared = w1*w2*(c*c/y**2-1)**2
        variation = 2*mass.second_coefficient/mass.value**2
        if variation != -c*(7*c*c+146*c-240)/(8*(c+8)**2):
            raise ValueError("Independent second-jet mass variation failed")
        result.append({"c": c, "K1": K1.value, "K2": K2.value,
                       "mass_squared": mass.value, "N_sum": Ns, "N_relative": Nr,
                       "omega_prime_squared": omega_prime_squared, "c_light_squared": cl.value,
                       "c_heavy_squared": ch.value, "D_squared": D_squared.value,
                       "mass_variation_ratio": variation})
    return result


def multiply_series(a, b):
    return [sum((a[j]*b[n-j] for j in range(n+1)), Poly()) for n in range(3)]


def elementary_series(roots):
    values = [[Poly(1), Poly(), Poly()]]
    for root in roots:
        new = [[Poly(), Poly(), Poly()] for _ in range(len(values)+1)]
        for index, coefficient in enumerate(values):
            product = multiply_series(coefficient, root)
            for order in range(3):
                new[index][order] += coefficient[order]
                new[index+1][order] += product[order]
        values = new
    return values


def derivative(polynomial, name):
    index = NAMES.index(name)
    terms = {}
    for monomial, coefficient in polynomial.terms.items():
        if monomial[index]:
            reduced = list(monomial)
            reduced[index] -= 1
            terms[tuple(reduced)] = coefficient*monomial[index]
    return Poly(terms)


def identities():
    c, y, eps, x = (Poly.variable(name) for name in ("c", "y", "u", "x"))
    beta = [Poly.variable(f"beta{i}") for i in range(5)]
    roots = [[c, Poly(), Poly()], [y, y*Q(1, 2), y*Q(1, 8)],
             [y, -y*Q(1, 2), y*Q(1, 8)], [y, Poly(), Poly()]]
    elementary = elementary_series(roots)
    tt = -2*sum((beta[i]*elementary[i][2] for i in range(5)), Poly())
    expected = -y*(beta[1]+beta[2]*(c+y)+beta[3]*c*y)*Q(1, 2)
    # Cancel the exact epsilon² factor coefficientwise in c*d⁴-2.
    cc, d = 2+eps**2, 1+eps**2*x**2
    raw = cc*d**4-2
    index = NAMES.index("u")
    quotient = {}
    for monomial, coefficient in raw.terms.items():
        if monomial[index] < 2:
            raise ValueError("The removable inner denominator factor is absent")
        powers = list(monomial)
        powers[index] -= 2
        quotient[tuple(powers)] = coefficient
    removed = Poly(quotient)
    leading = Poly({monomial: coefficient for monomial, coefficient in removed.terms.items() if monomial[index] == 0})
    exact, algebraic = (1+8*x*x)*Q(1, 96), (1+8*x*x)*Q(1, 80)
    return {"literal_norm_two_TT_elementary_coefficient": tt-expected,
            "exact_inner_denominator_division": eps**2*removed-raw,
            "limiting_inner_denominator": leading-(1+8*x*x),
            "exact_limiting_forced_particular": (1+8*x*x)*derivative(derivative(exact, "x"), "x")+80*exact-(1+8*x*x),
            "algebraic_omission_residual_one_fifth": (1+8*x*x)*derivative(derivative(algebraic, "x"), "x")+80*algebraic-(1+8*x*x)*Q(6, 5),
            "forced_particular_ratio_five_sixths": exact-algebraic*Q(5, 6)}


def checks():
    polys = identities()
    if any(not value.is_zero() for value in polys.values()):
        raise ValueError("An independent variable-beta polynomial identity failed")
    return {"coefficientwise_identities": dict.fromkeys(polys, "0"),
            "actual_background_Fraction_fixtures": background_fixtures(),
            "canonical_second_jet_Fraction_fixtures": canonical_fixtures()}
