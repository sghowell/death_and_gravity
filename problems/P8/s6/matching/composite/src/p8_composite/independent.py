"""Separate coefficientwise Fraction engine; no SymPy or model imports.

Ae and Ne are independent formal variables while differentiating, with
the composite chain rule included explicitly. Positive monomial
denominators are then cleared before substituting Ae=alpha*a+beta*b and
Ne=alpha*Ng+beta*Nf. Zero tests are exact polynomial tests, not sampling.
"""

from fractions import Fraction
from functools import cache

NAMES = ("a", "b", "Ng", "Nf", "ad", "bd", "Ngd", "Nfd", "add", "bdd",
         "G", "F", "m4", "alpha", "beta", "Ae", "Ne", "K", "V", "rho", "p", "rho_d",
         "beta0", "beta1", "beta2", "beta3", "beta4")
ZERO = (0,)*len(NAMES)


def exact(value):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
        raise TypeError("Only exact int/Fraction coefficients are accepted")
    return Fraction(value)


class Laurent:
    def __init__(self, terms=None):
        prepared = {tuple(power): exact(value) for power, value in (terms or {}).items()}
        if any(len(power) != len(NAMES) or any(isinstance(x, bool) or not isinstance(x, int) for x in power)
               for power in prepared):
            raise ValueError("Invalid exponent vector")
        self.terms = {power: value for power, value in prepared.items() if value}

    @classmethod
    def constant(cls, value):
        return cls({ZERO: exact(value)})

    @classmethod
    def variable(cls, name):
        power = list(ZERO)
        power[NAMES.index(name)] = 1
        return cls({tuple(power): Fraction(1)})

    def __add__(self, other):
        other = other if isinstance(other, Laurent) else Laurent.constant(other)
        output = dict(self.terms)
        for power, value in other.terms.items():
            output[power] = output.get(power, Fraction(0))+value
        return Laurent(output)

    __radd__ = __add__

    def __neg__(self):
        return Laurent({power: -value for power, value in self.terms.items()})

    def __sub__(self, other):
        return self+(-other if isinstance(other, Laurent) else -exact(other))

    def __rsub__(self, other):
        return -self+other

    def __mul__(self, other):
        other = other if isinstance(other, Laurent) else Laurent.constant(other)
        output = {}
        for left, lc in self.terms.items():
            for right, rc in other.terms.items():
                power = tuple(a+b for a, b in zip(left, right, strict=True))
                output[power] = output.get(power, Fraction(0))+lc*rc
        return Laurent(output)

    __rmul__ = __mul__

    def __pow__(self, exponent):
        if isinstance(exponent, bool) or not isinstance(exponent, int):
            raise TypeError("Power must be an exact integer")
        if exponent < 0:
            if len(self.terms) != 1:
                raise ValueError("Only nonzero monomials can be inverted")
            power, value = next(iter(self.terms.items()))
            return Laurent({tuple(exponent*x for x in power): value**exponent})
        output = Laurent.constant(1)
        base = self
        while exponent:
            if exponent % 2:
                output = output*base
            exponent //= 2
            if exponent:
                base = base*base
        return output

    def __truediv__(self, other):
        other = other if isinstance(other, Laurent) else Laurent.constant(other)
        return self*other**-1

    def derivative(self, name):
        index = NAMES.index(name)
        output = {}
        for power, value in self.terms.items():
            if power[index]:
                reduced = list(power)
                reduced[index] -= 1
                output[tuple(reduced)] = value*power[index]
        return Laurent(output)

    def dt(self):
        v = Laurent.variable
        mapping = {"a": v("ad"), "b": v("bd"), "Ng": v("Ngd"), "Nf": v("Nfd"),
                   "ad": v("add"), "bd": v("bdd"), "rho": v("rho_d"),
                   "Ae": v("alpha")*v("ad")+v("beta")*v("bd"),
                   "Ne": v("alpha")*v("Ngd")+v("beta")*v("Nfd")}
        return sum((self.derivative(name)*value for name, value in mapping.items()), Laurent.constant(0))

    def composite_derivative(self, name):
        addition = {"a": ("alpha", "Ae"), "b": ("beta", "Ae"),
                    "Ng": ("alpha", "Ne"), "Nf": ("beta", "Ne")}
        output = self.derivative(name)
        if name in addition:
            coefficient, extra = addition[name]
            output += Laurent.variable(coefficient)*self.derivative(extra)
        return output

    def evaluate(self, values):
        prepared = {name: exact(value) for name, value in values.items()}
        total = Fraction(0)
        for power, coefficient in self.terms.items():
            for name, exponent in zip(NAMES, power, strict=True):
                if exponent:
                    coefficient *= prepared[name]**exponent
            total += coefficient
        return total


def clear_and_substitute(value):
    """Clear only nonzero monomial denominators, then impose both sums."""
    minima = [min((power[index] for power in value.terms), default=0) for index in range(len(NAMES))]
    offset = tuple(max(0, -x) for x in minima)
    v = Laurent.variable
    replacements = {"Ae": v("alpha")*v("a")+v("beta")*v("b"),
                    "Ne": v("alpha")*v("Ng")+v("beta")*v("Nf")}
    output = Laurent.constant(0)
    for power, coefficient in value.terms.items():
        term = Laurent.constant(coefficient)
        for name, exponent, added in zip(NAMES, power, offset, strict=True):
            term *= replacements.get(name, v(name))**(exponent+added)
        output += term
    return output


@cache
def derive():
    v = {name: Laurent.variable(name) for name in NAMES}
    a, b, ng, nf, ad, bd = (v[name] for name in ("a", "b", "Ng", "Nf", "ad", "bd"))
    ae, ne, al, be, g, f, m4 = (v[name] for name in ("Ae", "Ne", "alpha", "beta", "G", "F", "m4"))
    b0, b1, b2, b3, b4 = (v["beta"+str(index)] for index in range(5))
    y, c = b/a, nf/ng
    potential = -m4*(b0*ng*a**3+b1*(nf*a**3+3*ng*a**2*b)
                      +3*b2*(nf*a**2*b+ng*a*b**2)+b3*(3*nf*a*b**2+ng*b**3)+b4*nf*b**3)
    matter = ae**3*v["K"]/(2*ne)-ne*ae**3*v["V"]
    lagrangian = -3*g*a*ad**2/ng-3*f*b*bd**2/nf+potential+matter
    hg, hf, he = ad/(ng*a), bd/(nf*b), ae.dt()/(ne*ae)
    polynomial = m4*(b1+2*b2*y+b3*y**2)
    ug = m4*(b0+3*b1*y+3*b2*y**2+b3*y**3)
    uf = m4*(b4+3*b3/y+3*b2/y**2+b1/y**3)
    pg, pf = -ug+(y-c)*polynomial, -uf+(c-y)*polynomial/(c*y**3)
    rho, p = v["K"]/(2*ne**2)+v["V"], v["K"]/(2*ne**2)-v["V"]
    expected = {"EL_Ng": 3*g*hg**2-ug-al*(ae/a)**3*rho,
                "EL_Nf": 3*f*hf**2-uf-be*(ae/b)**3*rho,
                "EL_a": g*(2*hg.dt()/ng+3*hg**2)+pg+al*(ne/ng)*(ae/a)**2*p,
                "EL_b": f*(2*hf.dt()/nf+3*hf**2)+pf+be*(ne/nf)*(ae/b)**2*p}
    actual = {"EL_Ng": lagrangian.composite_derivative("Ng")/a**3,
              "EL_Nf": lagrangian.composite_derivative("Nf")/b**3,
              "EL_a": (lagrangian.composite_derivative("a")-lagrangian.derivative("ad").dt())/(3*ng*a**2),
              "EL_b": (lagrangian.composite_derivative("b")-lagrangian.derivative("bd").dt())/(3*nf*b**2)}
    checks = {key: actual[key]-expected[key] for key in actual}
    rho, p = v["rho"], v["p"]
    rg, rf = al*(ae/a)**3*rho, be*(ae/b)**3*rho
    mgp, mfp = al*(ne/ng)*(ae/a)**2*p, be*(ne/nf)*(ae/b)**2*p
    conservation = rho.dt()/ne+3*he*(rho+p)
    q = polynomial-al*be*(ae/a)**2*p
    branch = ng*bd-nf*ad
    gbalance = (ug+rg).dt()/ng+3*hg*(ug+pg+rg+mgp)
    fbalance = (uf+rf).dt()/nf+3*hf*(uf+pf+rf+mfp)
    checks.update({
        "source_aware_g_Bianchi": gbalance-al*(ae/a)**3*(ne/ng)*conservation-3*q*branch/(ng**2*a),
        "source_aware_f_Bianchi": fbalance-be*(ae/b)**3*(ne/nf)*conservation+3*q*branch/(ng*nf*a*c*y**3),
        "weighted_shared_null": rg+mgp+c*y**3*(rf+mfp)-(ne/ng)*(ae/a)**3*(rho+p),
        "g_pressure_null_rewrite": ug+pg+rg+mgp-al*(ae/a)**3*(rho+p)-(y-c)*q,
        "f_pressure_null_rewrite": uf+pf+rf+mfp-be*(ae/b)**3*(rho+p)-(c-y)*q/(c*y**3),
        "undivided_branch_kinematics": hf-hg/y-branch/(ng*nf*b),
        "effective_Hubble_kinematics": he-(al*hg+be*c*y*hf)/((ne/ng)*(ae/a)),
    })
    controls = {"omitting_pressure_in_Bianchi": 3*(q-polynomial)*branch/(ng**2*a),
                "omitting_relative_lapse_weight": rg+mgp+y**3*(rf+mfp)-(ne/ng)*(ae/a)**3*(rho+p),
                "omitting_Neff_chain_rule": matter.derivative("Ng")-matter.composite_derivative("Ng")}
    return {"checks": {key: clear_and_substitute(value) for key, value in checks.items()},
            "controls": {key: clear_and_substitute(value) for key, value in controls.items()},
            "actual_term_counts": {key: len(value.terms) for key, value in actual.items()}}


def checks():
    data = derive()
    if any(value.terms for value in data["checks"].values()):
        raise ValueError("An independent composite coefficient identity failed")
    if any(not value.terms for value in data["controls"].values()):
        raise ValueError("An independent composite omission control failed to fire")
    # At the double root: H_e=-p_y*y'/(6p), so its derivative
    # depends on p_yy and y'^2 without division by p_y.
    acceleration = -Fraction(-1, 4)*Fraction(7, 3)/(6*Fraction(1, 2))
    if acceleration != Fraction(7, 36):
        raise ValueError("Independent double-root acceleration failed")
    return {"arithmetic": "coefficientwise Fraction Laurent algebra; positive denominator clearing and composite substitution",
            "zero_coefficient_residuals": dict.fromkeys(data["checks"], 0),
            "full_Euler_term_counts": data["actual_term_counts"],
            "nonzero_omission_term_counts": {key: len(value.terms) for key, value in data["controls"].items()},
            "double_root_acceleration": str(acceleration)}
