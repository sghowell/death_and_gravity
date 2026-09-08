"""Exact symbolic-dimensional Einstein sums by transverse index partitions."""
from collections import Counter
from functools import cache

import sympy as sp
from p8_vector_quadratic import invariants, tensors

from .geometry import dimension


def assignments(labels):
    if type(labels) is not tuple or not 1 <= len(labels) <= 4 or len(set(labels)) != len(labels):
        raise ValueError("Require one to four distinct immutable Einstein labels")
    if any(type(label) is not str or len(label) != 1 or label not in "abcdefghijklmnopqrstuvwxyz" for label in labels):
        raise ValueError("Require lowercase single-letter Einstein labels")
    def generate(position, mapping, anonymous_count):
        if position == len(labels):
            multiplicity = sp.prod(dimension-1-j for j in range(anonymous_count))
            yield dict(mapping), multiplicity
            return
        label = labels[position]
        for index in range(anonymous_count+3):
            mapping[label] = index
            yield from generate(position+1, mapping, anonymous_count+(index == anonymous_count+2))
        del mapping[label]
    yield from generate(0, {}, 0)


def mode_component(mode, sequence, i, j):
    # The second-order basis has rank <=4. O(D-1) reflection symmetry
    # forces every anonymous transverse index to occur an even number of
    # times. At most two distinct anonymous indices can remain; their
    # component is exactly represented in the frozen four-dimensional frame.
    if mode not in (tensors.plus, tensors.minus) or type(sequence) is not tuple or len(sequence) > 2:
        raise ValueError("Require a frozen mode and derivative order at most two")
    indices = sequence+(i, j)
    if any(type(index) is not int or not 0 <= index < 6 for index in indices):
        raise ValueError("Require native frame-class indices between zero and five")
    counts = Counter(index for index in indices if index >= 2)
    if any(count % 2 for count in counts.values()):
        return sp.Integer(0)
    if len(sequence) > 2 or len(counts) > 2:
        raise ValueError("Only the certified rank-four second-derivative basis is supported")
    relabel = {index: place+2 for place, index in enumerate(sorted(counts))}
    mapped = tuple(relabel.get(index, index) for index in indices)
    return mode.covariant(mapped[:-2], mapped[-2], mapped[-1])


def ricci(i, j):
    if any(type(index) is not int or not 0 <= index < 6 for index in (i, j)):
        raise ValueError("Require native frame-class indices between zero and five")
    if i != j:
        return sp.Integer(0)
    h0, h1 = tensors.H[:2]
    return -dimension*(h1+h0**2) if i == 0 else -h1-dimension*h0**2


@cache
def contraction(factors):
    if factors not in tuple(value for _, value in invariants.SECOND):
        raise ValueError("Only the frozen second-derivative basis is supported")
    counts = Counter("".join(sequence+indices for _, sequence, indices in factors))
    if any(count != 2 for count in counts.values()):
        raise ValueError("Require paired Einstein indices")
    output = sp.Integer(0)
    for mapping, multiplicity in assignments(tuple(counts)):
        value, leg = multiplicity, 0
        for kind, derivatives, indices in factors:
            slots = tuple(mapping[label] for label in indices)
            if kind == "R":
                component = -2*dimension*tensors.H[1]-dimension*(dimension+1)*tensors.H[0]**2
            elif kind == "Ric":
                component = ricci(*slots)
            elif kind == "Riem":
                component = tensors.riemann(*slots)
            elif kind == "Y":
                component = mode_component(tensors.plus if leg == 0 else tensors.minus,
                                           tuple(mapping[label] for label in derivatives), *slots)
                leg += 1
            else:
                raise ValueError("Unknown continued tensor kind")
            if component == 0:
                value = sp.Integer(0)
                break
            value *= component
        output += value
    return sp.expand((output+tensors.swap_legs(output))/2)


@cache
def second_pole():
    return sp.expand(sum(coefficient*contraction(factors) for coefficient, factors in invariants.SECOND)/(48*tensors.mass**2))


@cache
def checks():
    out = {"exact_index_partition_count_rank_"+str(rank):
           sp.expand(sum(multiplicity for _, multiplicity in assignments(tuple("abcdef"[:rank])))-(dimension+1)**rank)
           for rank in range(1, 5)}
    out["continued_second_pole_returns_four_dimensions"] = sp.expand(
        second_pole().subs(dimension, 3)-invariants.evaluate(2)/(48*tensors.mass**2))
    for index, (_, factors) in enumerate(invariants.SECOND):
        out["continued_individual_second_invariant_"+str(index)] = sp.expand(
            contraction(factors).subs(dimension, 3)-tensors.contraction(factors))
    return out
