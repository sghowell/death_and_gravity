"""Flat zero-past two-current response, with three exact subtractions."""
from functools import cache

import sympy as sp
from p8_vector_metric_local import canonical, jets, variation

z, h = sp.symbols("momentum_fraction positive_h", positive=True)
alpha, beta = sp.symbols("mass_lapse_first_a mass_lapse_first_b", real=True)
alpha2, beta2 = sp.symbols("mass_lapse_second_a mass_lapse_second_b", real=True)
w, m, s, y = sp.symbols("omega mass laplace radial_fraction", positive=True)
N, Z = sp.symbols("physical_lapse physical_logscale", real=True)


def sector(value):
    if type(value) is not str or value not in ("T", "L"):
        raise ValueError("Require native sector T or L")
    return value


def order(value):
    if type(value) is not int or value not in (1, 2, 3):
        raise ValueError("Require native radial moment order 1, 2 or 3")
    return value


def pair(value):
    value = sector(value)
    return sp.ImmutableMatrix([beta*(1-z), 2*(1-z)] if value == "T"
                              else [beta+alpha*z, 2*(1+z)])


@cache
def matrix():
    return sp.ImmutableMatrix(sp.expand(2*pair("T")*pair("T").T+pair("L")*pair("L").T))


def actual(value):
    return value.subs({alpha: 4/(9*h), beta: 28/(81*h)}, simultaneous=True).applyfunc(sp.factor)


@cache
def coefficients():
    return tuple(sp.ImmutableMatrix(matrix().applyfunc(lambda entry, j=j: sp.expand(entry).coeff(z, j)))
                 for j in range(3))


def moment(number, denominator):
    """Integral_0^1 y^(2n)/(d-y^2) dy, d>1; symbolic formula only."""
    number = order(number)
    if isinstance(denominator, (bool, str, float)) or denominator in (sp.true, sp.false):
        raise ValueError("Require an exact symbolic denominator, not a converted scalar")
    d = sp.sympify(denominator)
    if d.has(sp.Float) or d in (sp.oo, -sp.oo, sp.zoo, sp.nan):
        raise ValueError("Require a finite exact symbolic denominator")
    if d.is_number and (d.is_real is not True or d <= 1):
        raise ValueError("Require d>1 for a numerical denominator")
    base = sp.atanh(1/sp.sqrt(d))/sp.sqrt(d)
    return d**number*base-sum(d**k/sp.Integer(2*(number-k)-1) for k in range(number))


def subtracted_normalized(denominator):
    """64*pi^2/s^4 times the remainder after Taylor orders 0,2,4."""
    return sp.ImmutableMatrix(-sum((coefficient*moment(j+1, denominator)
                                    for j, coefficient in enumerate(coefficients())), sp.zeros(2))/4)


@cache
def finite_asymptotic():
    return sp.ImmutableMatrix(sum((coefficient*weight for coefficient, weight in zip(
        coefficients(), (sp.Rational(1, 4), sp.Rational(1, 3), sp.Rational(23, 60)))), sp.zeros(2)))


def contact(value):
    """Expectation of -H_,ij at fixed physical canonical coordinates."""
    return _contact(sector(value))


@cache
def _contact(value):
    q = w**2-m**2
    am = 1+alpha*(N-1)+alpha2*(N-1)**2/2
    bm = 1+beta*(N-1)+beta2*(N-1)**2/2
    qp = q*sp.exp(-2*Z)
    if value == "T":
        g2, frequency2 = sp.exp(Z)/N, N**2*(qp+m**2*bm)
    else:
        g2 = sp.exp(3*Z)*m**2*am*qp/(N*(qp+m**2*am))
        frequency2 = N**2*bm*(qp/am+m**2)
    point = {N: 1, Z: 0}
    g0 = g2.subs(point)
    expectation = w*(g0/g2+g2*frequency2/(g0*w**2))/4
    vertices = sp.ImmutableMatrix([[-sp.factor(sp.diff(expectation, left, right).subs(point))
                                   for right in (N, Z)] for left in (N, Z)])
    vacuum = -sp.sqrt(frequency2)/2
    static = sp.ImmutableMatrix([[sp.factor(sp.diff(vacuum, left, right).subs(point))
                                 for right in (N, Z)] for left in (N, Z)])
    pairs = pair(value).subs(z, 1-m**2/w**2)
    bubble_static = w*pairs*pairs.T/8
    return {"contact": vertices, "static_vacuum_hessian": static,
            "static_check": sp.ImmutableMatrix((vertices+bubble_static-static).applyfunc(sp.factor))}


@cache
def checks():
    out = {}
    mapping = {jets.D: 3, jets.z: z, jets.alpha[0]: alpha, jets.beta[0]: beta}
    for value in ("T", "L"):
        weights = canonical.data(value)["weights"]
        independent = sp.ImmutableMatrix([(weights[name][0]-weights[name][1]).subs(mapping)
                                          for name in ("N", "Z")])
        out[value+"_pair_from_physical_Hamiltonian"] = sp.ImmutableMatrix(
            (independent-pair(value)).applyfunc(sp.factor))
        out[value+"_static_response_includes_contact"] = contact(value)["static_check"]
    denominator = sp.Symbol("radial_denominator", positive=True)
    exact = 1/(4*w**2+s**2)
    taylor = 1/(4*w**2)-s**2/(4*w**2)**2+s**4/(4*w**2)**3
    out["three_subtracted_retarded_denominator"] = sp.factor(
        exact-taylor+s**6/((4*w**2)**3*(4*w**2+s**2)))
    for number in (1, 2, 3):
        polynomial = sum(denominator**k*y**(2*(number-k)-2) for k in range(number))
        out["radial_polynomial_division_"+str(number)] = sp.factor(
            y**(2*number)/(denominator-y**2)
            -denominator**number/(denominator-y**2)+polynomial)
        primitive = denominator**number*sp.atanh(y/sp.sqrt(denominator))/sp.sqrt(denominator)-sum(
            denominator**k*y**(2*(number-k)-1)/sp.Integer(2*(number-k)-1) for k in range(number))
        out["independent_radial_primitive_"+str(number)] = sp.factor(
            sp.diff(primitive, y)-y**(2*number)/(denominator-y**2))
    # p=m*y/sqrt(1-y^2), omega=m/sqrt(1-y^2).
    p = m*y/sp.sqrt(1-y**2)
    jacobian = p**2*sp.diff(p, y)
    omega = m/sp.sqrt(1-y**2)
    original = -s**6*jacobian*omega**3/(2*64*omega**6*(4*omega**2+s**2))/(2*sp.pi**2)
    transformed = -s**6*y**2/(256*sp.pi**2*(4*m**2+s**2*(1-y**2)))
    out["radial_measure_and_Kubo_normalization"] = sp.factor(original-transformed)
    source = sp.ImmutableMatrix([alpha+beta, 4])
    out["rank_one_high_frequency_pair"] = sp.ImmutableMatrix((matrix().subs(z, 1)-source*source.T).applyfunc(sp.factor))
    out["highest_log_normalization_from_radial_pole"] = sp.ImmutableMatrix(
        (matrix().subs(z, 1)/8-2*(source/4)*(source/4).T).applyfunc(sp.factor))
    return out


@cache
def adiabatic_checks():
    out = {}
    frozen = {field: 0 for row in (jets.H, jets.alpha[1:], jets.beta[1:],
                                   jets.alpha2, jets.beta2) for field in row}
    frozen.update({jets.D: 3, jets.z: z, jets.alpha[0]: alpha, jets.beta[0]: beta})
    for value in ("T", "L"):
        for j in (1, 2):
            for left, output in enumerate(("N", "Z")):
                coefficient = variation.data(value)["coefficients"][output][j]
                for right, row in enumerate((jets.n, jets.v)):
                    actual = sp.expand(coefficient).coeff(row[2*j]).subs(frozen, simultaneous=True)
                    expected = (-1)**j*pair(value)[left]*pair(value)[right]/(2*4**j)
                    out[value+output+str(right)+"_flat_adiabatic_equals_frequency_Taylor_"+str(2*j)] = sp.factor(actual-expected)
    return out
