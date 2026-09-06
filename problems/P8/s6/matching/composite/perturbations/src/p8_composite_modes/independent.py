"""Independent coefficientwise Fraction arithmetic, without SymPy/model imports."""

from fractions import Fraction

NAMES = ("c", "y", "alpha", "beta", "m4", "p", "rho", "beta0", "beta1", "beta2", "beta3", "beta4",
         "A", "B", "C", "x", "z", "u", "v", "accel", "root")
ZERO = (0,)*len(NAMES)


class Poly:
    """Small sparse exact polynomial ring; no floating inputs or inversions."""

    def __init__(self, value=0):
        if isinstance(value, Poly):
            self.terms = dict(value.terms)
        elif isinstance(value, dict):
            if any(len(key) != len(NAMES) or any(not isinstance(power, int) or power < 0 for power in key)
                   or not isinstance(coefficient, (int, Fraction)) for key, coefficient in value.items()):
                raise TypeError("Exact nonnegative polynomial monomials required")
            self.terms = {key: Fraction(coefficient) for key, coefficient in value.items() if coefficient}
        elif isinstance(value, (int, Fraction)) and not isinstance(value, bool):
            self.terms = {ZERO: Fraction(value)} if value else {}
        else:
            raise TypeError("Polynomial coefficients must be exact int/Fraction")

    @classmethod
    def variable(cls, name):
        powers = list(ZERO)
        powers[NAMES.index(name)] = 1
        return cls({tuple(powers): 1})

    def __add__(self, other):
        result = dict(self.terms)
        for key, coefficient in Poly(other).terms.items():
            result[key] = result.get(key, 0)+coefficient
        return Poly(result)

    __radd__ = __add__

    def __neg__(self):
        return Poly({key: -coefficient for key, coefficient in self.terms.items()})

    def __sub__(self, other):
        return self+-Poly(other)

    def __rsub__(self, other):
        return Poly(other)+-self

    def __mul__(self, other):
        result = {}
        for first, a in self.terms.items():
            for second, b in Poly(other).terms.items():
                key = tuple(x+y for x, y in zip(first, second, strict=True))
                result[key] = result.get(key, 0)+a*b
        return Poly(result)

    __rmul__ = __mul__

    def __pow__(self, exponent):
        if not isinstance(exponent, int) or exponent < 0:
            raise TypeError("Only nonnegative integer powers")
        result = Poly(1)
        for _ in range(exponent):
            result = result*self
        return result

    def reduce_root(self):
        """Quotient by root²=7/3; normal form remains coefficientwise exact."""
        index = NAMES.index("root")
        result = {}
        for key, coefficient in self.terms.items():
            reduced = list(key)
            quotient, reduced[index] = divmod(reduced[index], 2)
            reduced = tuple(reduced)
            result[reduced] = result.get(reduced, 0)+coefficient*Fraction(7, 3)**quotient
        return Poly(result)

    def evaluate(self, values):
        total = Fraction()
        for powers, coefficient in self.terms.items():
            term = coefficient
            for name, power in zip(NAMES, powers, strict=True):
                if power:
                    value = values[name]
                    if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
                        raise TypeError("Evaluation inputs must be exact int/Fraction")
                    term *= Fraction(value)**power
            total += term
        return total

    def is_zero(self):
        return not self.terms


def multiply_series(first, second):
    """Truncated order-two series; coefficient ring is the independent Poly."""
    return [sum((first[j]*second[n-j] for j in range(n+1)), Poly()) for n in range(3)]


def elementary_series(roots):
    result = [[Poly(1), Poly(), Poly()]]
    for root in roots:
        next_result = [[Poly(), Poly(), Poly()] for _ in range(len(result)+1)]
        for degree, coefficient in enumerate(result):
            product = multiply_series(coefficient, root)
            for order in range(3):
                next_result[degree][order] += coefficient[order]
                next_result[degree+1][order] += product[order]
        result = next_result
    return result


def identities():
    variables = {name: Poly.variable(name) for name in NAMES}
    c, y, alpha, beta, m4, p, rho = (variables[name] for name in NAMES[:7])
    bs = [variables[f"beta{j}"] for j in range(5)]
    r, s, ell = alpha+beta*y, alpha+beta*c, c+y
    roots = [[c, Poly(), Poly()], [y, y*Fraction(1, 2), y*Fraction(1, 8)],
             [y, -y*Fraction(1, 2), y*Fraction(1, 8)], [y, Poly(), Poly()]]
    elementary = elementary_series(roots)
    tt_interaction = -m4*sum((bs[j]*elementary[j][2] for j in range(5)), Poly())
    volume = [Poly(1), Poly(), Poly()]
    for root in roots[1:]:
        volume = multiply_series(volume, [alpha+beta*root[0], beta*root[1], beta*root[2]])
    tt_matter = s*p*volume[2]
    mu = y*(m4*(bs[1]+bs[2]*(c+y)+bs[3]*c*y)-alpha*beta*r*s*p)

    # Independent quotient-rule calculation for the composite matter
    # action numerator N=r²[s² n W-(rho-p)D²], denominator 2D.
    # First w-derivatives are represented after multiplication by H=2(c+y)².
    null, d0, w0 = rho+p, r*s, r**2
    dh, wh = -alpha*beta*ell, 2*alpha*beta*y
    numerator0 = r**2*(s**2*null*w0-(rho-p)*d0**2)
    numerator_h = r**2*(s**2*null*wh-2*(rho-p)*d0*dh)
    matter_cleared = numerator_h*d0-numerator0*dh
    expected_cleared = 2*d0**2*alpha*beta*(r**2*rho*ell+r*s*y*null)
    # Shift potential: coefficient of L in (1+yz)^2(1+Lz+cy z²).
    coefficient_l = [Poly(), Poly(1), 2*y, y**2, Poly()]
    pot_derivative = -m4*sum((bs[j]*coefficient_l[j] for j in range(5)), Poly())
    P = m4*(bs[1]+2*bs[2]*y+bs[3]*y**2)

    aa, bb, cc, x, z, u, v = (variables[name] for name in ("A", "B", "C", "x", "z", "u", "v"))
    determinant = aa*bb+aa*cc+bb*cc
    un = aa*(bb+cc)*x+bb*cc*z
    vn = aa*cc*x+bb*(aa+cc)*z
    reduced_numerator = aa*(x*determinant-un)**2+bb*(z*determinant-vn)**2+cc*(un-vn)**2
    du, dv = u*determinant-un, v*determinant-vn
    completion = (aa+cc)*du**2-2*cc*du*dv+(bb+cc)*dv**2
    lagrangian = aa*(x-u)**2+bb*(z-v)**2+cc*(u-v)**2
    accel, root = variables["accel"], variables["root"]
    yp, cp = -root, -(5+12*accel)*root*Fraction(1, 7)
    jet = yp*(yp-cp)*Fraction(1, 2)
    return {"TT_source_polynomial": tt_interaction+tt_matter+mu*Fraction(1, 4),
            "transverse_shift_matter_quotient_rule": matter_cleared-expected_cleared,
            "transverse_shift_potential_polynomial": pot_derivative+P,
            "first_shift_stationarity": (aa+cc)*un-cc*vn-aa*x*determinant,
            "second_shift_stationarity": (bb+cc)*vn-cc*un-bb*z*determinant,
            "stationary_vector_action": reduced_numerator-aa*bb*cc*determinant*(x-z)**2,
            "full_completed_square": lagrangian*determinant**2-aa*bb*cc*determinant*(x-z)**2-completion,
            "quadratic_root_generic_acceleration": (jet-(1-6*accel)*Fraction(1, 3)).reduce_root()}


def checks():
    values = identities()
    if any(not value.is_zero() for value in values.values()):
        raise ValueError("An independent Fraction polynomial identity failed")
    # Exact positive spring fixture; independent source-normalization control.
    aa, bb, cc = Fraction(1, 2), Fraction(1, 2), Fraction(1, 2)
    harmonic = 1/(1/aa+1/bb+1/cc)
    literal_K, literal_U = harmonic/8, Fraction(1, 16)
    printed_K, printed_U = Fraction(1, 24), Fraction(1, 16)
    frequency_difference = printed_U/printed_K-literal_U/literal_K
    free_mu, cd_mu, rescaled_mu = (Fraction(1, 3)-2*a for a in (Fraction(7, 36), Fraction(4), Fraction(4, 25)))
    if not (literal_K > 0 and frequency_difference == -Fraction(3, 2)
            and free_mu == -Fraction(1, 18) and cd_mu == -Fraction(23, 3)
            and rescaled_mu == Fraction(1, 75)):
        raise ValueError("An independent normalization or scale control failed")
    return {"coefficientwise_identities": dict.fromkeys(values, "0"),
            "exact_spring_K_at_mu_zero": str(literal_K),
            "printed_4_12_minus_literal_frequency_fixture": str(frequency_difference),
            "free_mu_u2": str(free_mu), "frozen_CD_mu_u2": str(cd_mu),
            "m_tau_5_mu_u2": str(rescaled_mu)}
