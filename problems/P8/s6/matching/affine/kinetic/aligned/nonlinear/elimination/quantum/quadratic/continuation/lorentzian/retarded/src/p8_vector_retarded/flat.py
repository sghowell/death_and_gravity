"""Independent flat static-energy and frequency-pole checks."""
from functools import cache

import sympy as sp
from p8_vector_bubble import pole as frozen

m, w = sp.symbols("mass omega", positive=True)
alpha, beta, n, frequency2 = sp.symbols("alpha beta source external_frequency_squared", real=True)
z = 1-m**2/w**2


def order(value):
    if type(value) is not int or value not in (0, 1, 2):
        raise ValueError("Require native frequency-squared order 0, 1 or 2")
    return value


def coefficient(value):
    value = order(value)
    bubble = m**4*beta**2/(4*w**3)+w*(beta+alpha*z)**2/8
    return sp.expand(bubble/(4*w**2)**value-(alpha**2*z*w/2 if value == 0 else 0))


def radial_residue(power):
    if type(power) is not int or power not in (-3, -1, 1):
        raise ValueError("Require one of the three divergent odd frequency powers")
    # Coefficients of 1/epsilon, omitting pi^-2. Derived from the Gamma radial integral.
    return {1: -m**4/32, -1: -m**2/8, -3: sp.Rational(1, 4)}[power]


def integrated_residue(value):
    value = order(value)
    integrand = coefficient(value)
    return sp.factor(sum(integrand.coeff(w, power)*radial_residue(power) for power in (1, -1, -3)))


def frequency_kernel():
    """Analytic expression with frequency2=(Omega+i0)^2 in retarded use."""
    return (m**4*beta**2/w+w**3*(beta+alpha*z)**2/2)/(4*w**2-frequency2)-alpha**2*z*w/2


@cache
def checks():
    out = {}
    for j in (0, 1, 2):
        actual = sp.diff(frequency_kernel(), frequency2, j).subs(frequency2, 0)/sp.factorial(j)
        out["exact_sine_transform_frequency_coefficient_"+str(2*j)] = sp.factor(actual-coefficient(j))
    q = w**2-m**2
    frequencies = {"T": sp.sqrt(q+m**2*(1+beta*n)),
                   "L": sp.sqrt((1+beta*n)*(q/(1+alpha*n)+m**2))}
    bubble = {"T": m**4*beta**2/(8*w**3), "L": w*(beta+alpha*z)**2/8}
    contact = {"T": 0, "L": -alpha**2*z*w/2}
    for sector in ("T", "L"):
        susceptibility = -sp.diff(frequencies[sector]/2, n, 2).subs(n, 0)
        out[sector+"_Abel_static_response_equals_vacuum_energy_derivative"] = sp.factor(
            bubble[sector]+contact[sector]-susceptibility)
    eps = sp.Symbol("epsilon")
    for power in (1, -1, -3):
        exponent = -sp.Rational(power, 2)
        expression = m**(3-2*eps-2*exponent)*sp.gamma(exponent-sp.Rational(3, 2)+eps)/(
            (4*sp.pi)**(sp.Rational(3, 2)-eps)*sp.gamma(exponent))
        out["radial_Gamma_residue_"+str(power)] = sp.simplify(
            sp.limit(eps*expression, eps, 0)*sp.pi**2-radial_residue(power))
    variables = {frozen.m2: m**2, frozen.P2: -frequency2,
                 frozen.PA: -alpha*frequency2, frozen.PA2: -alpha**2*frequency2,
                 frozen.TA: alpha+3*beta, frozen.TA2: alpha**2+3*beta**2}
    for number, name in enumerate(("constant", "second", "fourth")):
        target = frozen.parameter_integral()[name].subs(variables).expand().coeff(frequency2, number)/32
        out["retarded_frequency_pole_matches_frozen_Feynman_"+str(2*number)] = sp.factor(integrated_residue(number)-target)
    X, time = sp.symbols("X time", real=True)
    positive = sp.Symbol("positive_frequency", positive=True)
    C = X*sp.exp(-2*sp.I*positive*time)
    out["Kubo_sign"] = sp.simplify(sp.expand_complex(sp.I*(C-sp.conjugate(C)))-2*X*sp.sin(2*positive*time))
    damping = sp.Symbol("Abel_damping", positive=True)
    # Laplace sine integral obtained as the imaginary part of (damping-2iw)^-1.
    out["Abel_static_sine_integral"] = sp.simplify(
        sp.limit(sp.im(1/(damping-2*sp.I*w)), damping, 0, dir="+")-1/(2*w))
    return out


@cache
def controls():
    return {"omitted_longitudinal_contact_per_momentum": alpha**2*z*w/2,
            "omitted_contact_UV_residue_times_pi_squared":
                sp.factor(alpha**2*(radial_residue(1)-m**2*radial_residue(-1))/2),
            "static_residue_times_pi_squared": integrated_residue(0),
            "second_frequency_residue_times_pi_squared": integrated_residue(1),
            "fourth_frequency_residue_times_pi_squared": integrated_residue(2)}
