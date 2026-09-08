"""Independent physical-signature frame contractions for the local UV pole."""
from collections import Counter
from functools import cache

import sympy as sp
from p8_vector_quadratic import invariants, tensors
from p8_vector_quadratic_matching import contractions as euclidean
from p8_vector_quadratic_matching.geometry import dimension


def eta(index):
    if type(index) is not int or not 0 <= index < 6:
        raise ValueError("Require a native symbolic frame-class index")
    return -1 if index == 0 else 1


def covariant(mode, sequence, i, j):
    if mode not in (tensors.plus, tensors.minus) or type(sequence) is not tuple or len(sequence) > 2:
        raise ValueError("Require a frozen Fourier leg and at most two ordered derivatives")
    if any(type(index) is not int or not 0 <= index < 6 for index in sequence+(i, j)):
        raise ValueError("Require native symbolic frame-class indices")
    return _covariant(mode, sequence, i, j)


@cache
def _covariant(mode, sequence, i, j):
    if not sequence:
        return -mode.temporal[0] if i == j == 0 else mode.spatial[0] if i == j else sp.Integer(0)
    a, rest = sequence[0], sequence[1:]
    previous = covariant(mode, rest, i, j)
    if a == 0:
        value = mode.time_derivative(previous)
    elif a == 1:
        value = mode.sign*sp.I*tensors.k*previous
    else:
        value = sp.Integer(0)
    if a != 0:
        inner = rest+(i, j)
        for position, index in enumerate(inner):
            if index == 0:
                replacement = a
            elif index == a:
                replacement = 0
            else:
                continue
            changed = inner[:position]+(replacement,)+inner[position+1:]
            value -= tensors.H[0]*covariant(mode, changed[:-2], changed[-2], changed[-1])
    return sp.expand(value)


def ricci(i, j):
    if any(type(index) is not int or not 0 <= index < 6 for index in (i, j)):
        raise ValueError("Require native symbolic frame-class indices")
    if i != j:
        return sp.Integer(0)
    h0, h1 = tensors.H[:2]
    return -dimension*(h1+h0**2) if i == 0 else h1+dimension*h0**2


def riemann(a, b, c, d):
    if any(type(index) is not int or not 0 <= index < 6 for index in (a, b, c, d)):
        raise ValueError("Require native symbolic frame-class indices")
    if a == b or c == d:
        return sp.Integer(0)
    sectional = -(tensors.H[1]+tensors.H[0]**2) if 0 in (a, b) else tensors.H[0]**2
    if a == c and b == d:
        return sectional
    if a == d and b == c:
        return -sectional
    return sp.Integer(0)


@cache
def contraction(factors):
    if factors not in tuple(value for _, value in invariants.ZERO+invariants.SECOND):
        raise ValueError("Only the frozen zero- and two-derivative bases are supported")
    labels = tuple(Counter("".join(sequence+indices for _, sequence, indices in factors)))
    output = sp.Integer(0)
    for mapping, multiplicity in euclidean.assignments(labels):
        value, leg = multiplicity*sp.prod(eta(index) for index in mapping.values()), 0
        for kind, derivatives, indices in factors:
            slots = tuple(mapping[label] for label in indices)
            if kind == "R":
                component = 2*dimension*tensors.H[1]+dimension*(dimension+1)*tensors.H[0]**2
            elif kind == "Ric":
                component = ricci(*slots)
            elif kind == "Riem":
                component = riemann(*slots)
            elif kind == "Y":
                component = covariant(tensors.plus if leg == 0 else tensors.minus,
                                      tuple(mapping[label] for label in derivatives), *slots)
                leg += 1
            else:
                raise ValueError("Unknown covariant tensor")
            if component == 0:
                value = 0
                break
            value *= component
        output += value
    return sp.expand((output+tensors.swap_legs(output))/2)


@cache
def second_loop_pole():
    # Lorentzian loop density is minus the covariantly continued Euclidean
    # density. This is not the sign of the bare-action counterterm.
    return -sp.expand(sum(coefficient*contraction(factors) for coefficient, factors in invariants.SECOND)/(48*tensors.mass**2))


@cache
def zero_loop_control():
    return -sp.expand(sum(coefficient*contraction(factors) for coefficient, factors in invariants.ZERO)).subs(dimension, 3)
