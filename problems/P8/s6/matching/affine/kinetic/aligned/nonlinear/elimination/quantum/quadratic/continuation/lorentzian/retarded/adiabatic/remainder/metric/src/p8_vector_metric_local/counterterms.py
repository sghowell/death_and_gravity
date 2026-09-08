"""The fixed covariant subtraction action, varied twice before D=3."""
from functools import cache

import sympy as sp
from p8_vector_clock_matching import counterterms as first
from p8_vector_lorentzian_pole import response as mass_pole
from p8_vector_quadratic import tensors

from . import jets

mass = tensors.mass
e = sp.Symbol("metric_jet_parameter", real=True)


def truncate(value):
    return sp.Poly(sp.expand(value), e).terms()


def keep(value, degree=2):
    return sp.Add(*(coefficient*e**powers[0] for powers, coefficient in truncate(value) if powers[0] <= degree))


@cache
def geometry():
    D, H = jets.D, jets.H
    invN = 1-e*jets.n[0]+e**2*jets.n[0]**2
    proper = [keep((H[0]+e*jets.v[1])*invN)]
    for _ in range(3):
        proper.append(keep(invN*jets.time(proper[-1])))
    hp, acc = proper[0], keep(proper[1]+proper[0]**2)
    expansion2 = keep(hp**2)
    R = keep(2*D*acc+D*(D-1)*expansion2)
    Ricci2 = keep(D**2*acc**2+D*(acc+(D-1)*expansion2)**2)
    Riemann2 = keep(4*D*acc**2+2*D*(D-1)*expansion2**2)
    a4V = keep(-R**2/8+sp.Rational(29, 60)*Ricci2-Riemann2/15)
    a4sc = keep(R**2/72-Ricci2/180+Riemann2/180)
    scalar = first.scalar_curvature_tensor()
    R0 = scalar["R_D"]
    cnn = -R0/12+sp.Rational(5, 6)*D*(H[1]+H[0]**2)
    csp = R0/12-sp.Rational(5, 6)*(H[1]+D*H[0]**2)
    ga = {1: -2*cnn, 2: -scalar["gradient_normal"]/2}
    gb = {1: 6*csp, 2: 3*scalar["gradient_spatial_average"]/2}

    def first_change(value):
        return sp.expand(sum(sp.diff(value, H[j])*proper[j].coeff(e, 1) for j in range(4)))

    return {"R": R, "a4V": a4V, "a4sc": a4sc, "ga": ga, "gb": gb,
            "delta_ga": {j: first_change(ga[j]) for j in (1, 2)},
            "delta_gb": {j: first_change(gb[j]) for j in (1, 2)},
            "proper_H_jets": proper}


def density(order):
    return _density(jets.order(order))


@cache
def _density(order):
    D, n, v = jets.D, jets.n[0], jets.v[0]
    a, b, a2, b2 = jets.alpha[0], jets.beta[0], jets.alpha2[0], jets.beta2[0]
    volume1, volume2 = n+D*v, D*n*v+D**2*v**2/2
    if order == 0:
        c1 = (3*a+9*b)*n/2
        c2 = ((3*a2+9*b2)/4+(3*a**2+6*a*b+15*b**2)/8)*n**2
        return sp.expand(c2+volume1*c1+3*volume2)
    data = geometry()
    invariant = data["R"] if order == 1 else 2*data["a4V"]
    pure = keep((1+e*volume1+e**2*volume2)*invariant).coeff(e, 2)
    linear = -(a*n*data["delta_ga"][order]+b*n*data["delta_gb"][order]
               +(a2*n**2/2+volume1*a*n)*data["ga"][order]
               +(b2*n**2/2+volume1*b*n)*data["gb"][order])
    return sp.expand(pure+linear)


def generic_operator(output, order):
    output, order = jets.output(output), jets.order(order)
    return _generic_operator(output, order)


@cache
def _generic_operator(output, order):
    row = jets.n if output == "N" else jets.v
    value, result = density(order), sp.Integer(0)
    for j, field in enumerate(row):
        term = sp.diff(value, field)
        if term == 0:
            continue
        for _ in range(j):
            term = -jets.weighted(term)
        result += term
    return sp.expand(result)


def operator(output, order, dimension_jet=0):
    output, order = jets.output(output), jets.order(order)
    if type(dimension_jet) is not int or dimension_jet not in (0, 1):
        raise ValueError("Require native dimensional jet zero or one")
    return _operator(output, order, dimension_jet)


@cache
def _operator(output, order, dimension_jet):
    main = sp.diff(generic_operator(output, order), jets.D, dimension_jet).subs(jets.D, 3)
    result = mass**(4-2*order)*jets.actual(main)
    if order and output == "N":
        mapping = {field: jets.n[j] for j, field in enumerate(mass_pole.n)}
        result += 2*mass_pole.operator(2*order, dimension_jet).subs(mapping).subs(tensors.k, 0)
    return jets.linear_clean(result)


@cache
def controls():
    data = geometry()
    v, n, H = jets.v, jets.n, jets.H
    expected = [v[1]-H[0]*n[0],
                v[2]-2*H[1]*n[0]-H[0]*n[1]]
    out = {"proper_H_variation_"+str(j): sp.expand(data["proper_H_jets"][j].coeff(e, 1)-expected[j])
           for j in range(2)}
    # Constant-invariant divergences are omitted only as whole compact-support
    # Euler-null actions. The mass-weighted box R in ga4/gb4 remains present.
    out["pure_scalar_coefficients_are_fixed_before_dimension_limit"] = sp.expand(
        geometry()["a4V"].subs(e, 0)
        -(-geometry()["R"].subs(e, 0)**2/8
          +sp.Rational(29, 60)*(jets.D**2*(H[1]+H[0]**2)**2
                              +jets.D*(H[1]+jets.D*H[0]**2)**2)
          -(4*jets.D*(H[1]+H[0]**2)**2+2*jets.D*(jets.D-1)*H[0]**4)/15))
    return out
