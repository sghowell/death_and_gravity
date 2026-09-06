"""Independent exact Fraction Laurent-polynomial variation; no SymPy.

All identities are coefficientwise in free jet and coupling variables.
Negative powers occur only for nonzero scale factors and lapses. This is
not interpolation, a sampled sign check, or a physical negative-lapse model.
"""

from fractions import Fraction
from functools import cache

NAMES = ("a", "b", "Ng", "Nf", "ad", "bd", "Ngd", "Nfd", "add", "bdd",
         "G", "F", "Kg", "Kf", "Vg", "Vf", "m4", "beta0", "beta1", "beta2", "beta3", "beta4")
ZERO_POWER = (0,)*len(NAMES)
TIME_JETS = {"a": "ad", "b": "bd", "Ng": "Ngd", "Nf": "Nfd", "ad": "add", "bd": "bdd"}


def exact(value):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
        raise TypeError("Only exact int/Fraction coefficients are accepted")
    return Fraction(value)


class Laurent:
    def __init__(self, terms=None):
        prepared = {tuple(power): exact(value) for power, value in (terms or {}).items()}
        self.terms = {power: value for power, value in prepared.items() if value}
        if any(len(power) != len(NAMES) or any(isinstance(x, bool) or not isinstance(x, int) for x in power) for power in self.terms):
            raise ValueError("Invalid exponent vector")

    @classmethod
    def constant(cls, value):
        return cls({ZERO_POWER: exact(value)})

    @classmethod
    def variable(cls, name):
        power = list(ZERO_POWER)
        power[NAMES.index(name)] = 1
        return cls({tuple(power): Fraction(1)})

    def __add__(self, other):
        other = other if isinstance(other, Laurent) else Laurent.constant(other)
        values = dict(self.terms)
        for power, coefficient in other.terms.items():
            values[power] = values.get(power, Fraction(0))+coefficient
        return Laurent(values)

    __radd__ = __add__

    def __neg__(self):
        return Laurent({power: -coefficient for power, coefficient in self.terms.items()})

    def __sub__(self, other):
        other = other if isinstance(other, Laurent) else Laurent.constant(other)
        return self+-other

    def __rsub__(self, other):
        return -self+other

    def __mul__(self, other):
        other = other if isinstance(other, Laurent) else Laurent.constant(other)
        values = {}
        for left, lc in self.terms.items():
            for right, rc in other.terms.items():
                power = tuple(a+b for a, b in zip(left, right, strict=True))
                values[power] = values.get(power, Fraction(0))+lc*rc
        return Laurent(values)

    __rmul__ = __mul__

    def __pow__(self, exponent):
        if isinstance(exponent, bool) or not isinstance(exponent, int):
            raise TypeError("Exponent must be an exact integer")
        if exponent < 0:
            if len(self.terms) != 1:
                raise ValueError("Only nonzero monomials can be inverted")
            power, coefficient = next(iter(self.terms.items()))
            return Laurent({tuple(exponent*x for x in power): coefficient**exponent})
        result = Laurent.constant(1)
        for _ in range(exponent):
            result = result*self
        return result

    def __truediv__(self, other):
        other = other if isinstance(other, Laurent) else Laurent.constant(other)
        return self*other**-1

    def derivative(self, name):
        index = NAMES.index(name)
        values = {}
        for power, coefficient in self.terms.items():
            if power[index]:
                reduced = list(power)
                reduced[index] -= 1
                values[tuple(reduced)] = coefficient*power[index]
        return Laurent(values)

    def dt(self):
        return sum((self.derivative(name)*Laurent.variable(next_name) for name, next_name in TIME_JETS.items()),
                   Laurent.constant(0))

    def evaluate(self, values):
        prepared = {name: exact(value) for name, value in values.items()}
        return sum((coefficient*product(prepared[name]**exponent for name, exponent in zip(NAMES, power, strict=True) if exponent)
                    for power, coefficient in self.terms.items()), Fraction(0))


def product(values):
    result = Fraction(1)
    for value in values:
        result *= value
    return result


@cache
def derive():
    v = {name: Laurent.variable(name) for name in NAMES}
    a, b, ng, nf = (v[name] for name in ("a", "b", "Ng", "Nf"))
    ad, bd, g, f, m4 = (v[name] for name in ("ad", "bd", "G", "F", "m4"))
    beta0, beta1, beta2, beta3, beta4 = (v["beta"+str(index)] for index in range(5))
    y, c = b/a, nf/ng
    # Generate all e_n from the product over the four eigenvalues.
    elementary = [Laurent.constant(1)]
    for eigenvalue in (c, y, y, y):
        updated = [Laurent.constant(0) for _ in range(len(elementary)+1)]
        for index, value in enumerate(elementary):
            updated[index] += value
            updated[index+1] += eigenvalue*value
        elementary = updated
    potential = -m4*ng*a**3*sum((v["beta"+str(index)]*value for index, value in enumerate(elementary)), Laurent.constant(0))
    lagrangian = (-3*g*a*ad**2/ng-3*f*b*bd**2/nf+potential
                  +a**3*v["Kg"]/(2*ng)-ng*a**3*v["Vg"]+b**3*v["Kf"]/(2*nf)-nf*b**3*v["Vf"])
    rho_g = v["Kg"]/(2*ng**2)+v["Vg"]
    p_g = v["Kg"]/(2*ng**2)-v["Vg"]
    rho_f = v["Kf"]/(2*nf**2)+v["Vf"]
    p_f = v["Kf"]/(2*nf**2)-v["Vf"]
    rho_gv = m4*(beta0+3*beta1*y+3*beta2*y**2+beta3*y**3)
    p_gv = -m4*(beta0+beta1*(2*y+c)+beta2*(y**2+2*c*y)+beta3*c*y**2)
    rho_fv = m4*(beta4+3*beta3*y**-1+3*beta2*y**-2+beta1*y**-3)
    p_fv = -m4*(beta4+beta3*(2*y**-1+c**-1)+beta2*(y**-2+2*c**-1*y**-1)
                +beta1*c**-1*y**-2)
    hg, hf = ad/(ng*a), bd/(nf*b)
    actual = {"EL_Ng": lagrangian.derivative("Ng")/a**3,
              "EL_Nf": lagrangian.derivative("Nf")/b**3,
              "EL_a": (lagrangian.derivative("a")-lagrangian.derivative("ad").dt())/(3*ng*a**2),
              "EL_b": (lagrangian.derivative("b")-lagrangian.derivative("bd").dt())/(3*nf*b**2)}
    expected = {"EL_Ng": 3*g*hg**2-rho_g-rho_gv,
                "EL_Nf": 3*f*hf**2-rho_f-rho_fv,
                "EL_a": g*(2*hg.dt()/ng+3*hg**2)+p_g+p_gv,
                "EL_b": f*(2*hf.dt()/nf+3*hf**2)+p_f+p_fv}
    polynomial = m4*(beta1+2*beta2*y+beta3*y**2)
    numerator = ng*bd-nf*ad
    checks = {key: actual[key]-expected[key] for key in actual}
    checks.update({"potential_lapse_g_density": -potential.derivative("Ng")/a**3-rho_gv,
                   "potential_lapse_f_density": -potential.derivative("Nf")/b**3-rho_fv,
                   "potential_scale_g_pressure": potential.derivative("a")/(3*ng*a**2)-p_gv,
                   "potential_scale_f_pressure": potential.derivative("b")/(3*nf*b**2)-p_fv,
                   "g_interaction_null_factor": rho_gv+p_gv-(y-c)*polynomial,
                   "f_interaction_null_factor": rho_fv+p_fv-(c-y)*polynomial/(c*y**3),
                   "weighted_interaction_cancellation": rho_gv+p_gv+c*y**3*(rho_fv+p_fv),
                   "undivided_g_Bianchi": rho_gv.dt()/ng+3*hg*(rho_gv+p_gv)-3*polynomial*numerator/(ng**2*a),
                   "undivided_f_Bianchi": rho_fv.dt()/nf+3*hf*(rho_fv+p_fv)+3*polynomial*numerator/(ng*nf*a*c*y**3)})
    wrong_polynomial = m4*(beta1+beta2*y+beta3*y**2)
    controls = {"omitting_lapse_ratio_in_weight": rho_gv+p_gv+y**3*(rho_fv+p_fv),
                "wrong_beta2_Bianchi_coefficient": rho_gv+p_gv-(y-c)*wrong_polynomial}
    return {"L": lagrangian, "potential": potential, "actual": actual, "checks": checks, "controls": controls}


def exact_checks():
    data = derive()
    if any(value.terms for value in data["checks"].values()):
        raise ValueError("A Fraction Laurent-polynomial coefficient identity failed")
    if any(not value.terms for value in data["controls"].values()):
        raise ValueError("A Fraction polynomial omission control failed to fire")
    return {"arithmetic": "exact Fraction Laurent coefficients; arbitrary independent couplings and jets",
            "zero_coefficient_residuals": dict.fromkeys(data["checks"], 0),
            "full_Euler_equation_term_counts": {key: len(value.terms) for key, value in data["actual"].items()},
            "omission_control_nonzero_term_counts": {key: len(value.terms) for key, value in data["controls"].items()}}
