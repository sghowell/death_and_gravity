"""Literal local vertex words and complete covariant metric variations."""

from functools import cache
from itertools import permutations

import sympy as s
from p8_vacuum_affine_heavy_parent_one_loop import loops

SIX_ASSIGNMENTS = tuple(permutations(range(6)))
FOUR_ASSIGNMENTS = tuple(permutations(range(4)))
POSITIONS = tuple((i, j) for i in range(4) for j in range(i, 4))
GVARS = s.symbols("g0:10")
HVARS = s.symbols("h0:10")
BASIS = ("phi4", "phi2Y", "Y2", "Gal", "phi_vHv", "phi_box_Y", "phi2_Hdiff", "Y_Hdiff")


def symbolic_grams():
    gram, polarization = s.zeros(5), s.zeros(5)
    for (i, j), g, h in zip(POSITIONS, GVARS, HVARS):
        gram[i, j] = gram[j, i] = g
        polarization[i, j] = polarization[j, i] = h
    return gram, polarization


def quartic_word(name, p, gram):
    a, b, c, d = p
    if name == "phi4":
        return s.S.One
    if name == "phi2Y":
        return -gram[c, d]
    if name == "Y2":
        return gram[a, b] * gram[c, d]
    if name == "Gal":
        return loops.galileon_vertex_word(p, gram)
    if name == "phi_vHv":
        return gram[b, c] * gram[c, d]
    if name == "phi_box_Y":
        return gram[b, b] * gram[c, d]
    diff = gram[c, c] * gram[d, d] - gram[c, d] ** 2
    if name == "phi2_Hdiff":
        return diff
    if name == "Y_Hdiff":
        return -gram[a, b] * diff
    raise ValueError("Unknown local quartic basis word")


@cache
def quartic_polynomials():
    gram, _ = symbolic_grams()
    return {
        name: s.expand(s.Add(*(quartic_word(name, p, gram) for p in FOUR_ASSIGNMENTS)))
        for name in BASIS
    }


def quartic_vertex(name, gram):
    if name not in BASIS:
        raise ValueError("Unknown quartic vertex")
    return quartic_polynomials()[name].subs(
        {g: gram[i, j] for (i, j), g in zip(POSITIONS, GVARS)}
    )


def quartic_connection(name, p, gram, polarization, ak, connection_trace):
    a, b, c, d = p
    if name == "Gal":
        zb = (
            ak[a] * polarization[b, c]
            + ak[c] * polarization[a, b]
            - ak[b] * polarization[a, c]
        )
        zc = (
            ak[b] * polarization[c, d]
            + ak[d] * polarization[b, c]
            - ak[c] * polarization[b, d]
        )
        return (
            zb * (gram[d, d] - gram[c, d])
            - gram[a, b] * zc
            + gram[a, b] * gram[b, c] * connection_trace[d]
        )
    if name == "phi_vHv":
        return -(
            ak[b] * polarization[c, d]
            + ak[d] * polarization[b, c]
            - ak[c] * polarization[b, d]
        )
    if name == "phi_box_Y":
        return -connection_trace[b] * gram[c, d]
    if name in ("phi2_Hdiff", "Y_Hdiff"):
        value = (
            2 * (ak[c] + ak[d]) * polarization[c, d]
            - ak[c] * polarization[d, d]
            - ak[d] * polarization[c, c]
        )
        value -= connection_trace[c] * gram[d, d] + connection_trace[d] * gram[c, c]
        return value if name == "phi2_Hdiff" else -gram[a, b] * value
    return s.S.Zero


def quartic_contact(name, gram, polarization, ak, trace, k_e_p):
    """Canonical metric contact, including volume and all connection terms."""
    if name not in BASIS:
        raise ValueError("Unknown quartic contact")
    polynomial = quartic_polynomials()[name]
    sub = {g: gram[i, j] for (i, j), g in zip(POSITIONS, GVARS)}
    sub.update({h: polarization[i, j] for (i, j), h in zip(POSITIONS, HVARS)})
    metric = -2 * sum(h * s.diff(polynomial, g) for g, h in zip(GVARS, HVARS))
    connection_trace = tuple(2 * k_e_p[i] - ak[i] * trace for i in range(4))
    connection = s.Add(
        *(
            quartic_connection(name, p, gram, polarization, ak, connection_trace)
            for p in FOUR_ASSIGNMENTS
        )
    )
    return s.expand((metric + trace * polynomial).subs(sub) + connection)


def shifted_vertex(name, gram, ak, index):
    if type(index) is not int or not 0 <= index < 4:
        raise ValueError("Require one of four scalar emission legs")
    shifted = s.Matrix(gram)
    for j in range(4):
        shifted[index, j] += ak[j]
        shifted[j, index] += ak[j]
    return s.expand(quartic_vertex(name, shifted))


def six_tt_variation(
    name,
    p,
    gram,
    polarization,
    ak,
    external_connection_only=False,
    connection_only=False,
):
    """Literal six-field variation for null TT; volume and trace vanish there."""
    a, b, c, d, e, f = p
    if name == "phi6":
        return s.S.Zero
    if name == "phi4Y":
        return s.S.Zero if connection_only else 2 * polarization[e, f]
    if name == "phi2Y2":
        return (
            s.S.Zero
            if connection_only
            else -2
            * (polarization[c, d] * gram[e, f] + gram[c, d] * polarization[e, f])
        )
    if name == "Y3":
        return (
            s.S.Zero
            if connection_only
            else 2
            * (
                polarization[a, b] * gram[c, d] * gram[e, f]
                + gram[a, b] * polarization[c, d] * gram[e, f]
                + gram[a, b] * gram[c, d] * polarization[e, f]
            )
        )
    zb = (
        ak[a] * polarization[b, c]
        + ak[c] * polarization[a, b]
        - ak[b] * polarization[a, c]
    )
    zc = (
        ak[b] * polarization[c, d]
        + ak[d] * polarization[b, c]
        - ak[c] * polarization[b, d]
    )
    connection = (
        zb * (gram[d, d] - gram[c, d])
        if not external_connection_only or b < 4
        else s.S.Zero
    )
    connection += (
        -gram[a, b] * zc if not external_connection_only or c < 4 else s.S.Zero
    )
    metric = 2 * (
        polarization[a, b] * gram[b, c] * gram[d, d]
        + gram[a, b] * polarization[b, c] * gram[d, d]
        + gram[a, b] * gram[b, c] * polarization[d, d]
        - polarization[a, b] * gram[b, c] * gram[c, d]
        - gram[a, b] * polarization[b, c] * gram[c, d]
        - gram[a, b] * gram[b, c] * polarization[c, d]
    )
    if connection_only:
        metric = s.S.Zero
    value = metric + connection
    if name == "phi2_L3_minus_L4":
        return value
    if name != "Y_L3_minus_L4":
        raise ValueError("Unknown degree-six scalar word")
    return -gram[e, f] * value + (
        0
        if connection_only
        else 2 * polarization[e, f] * loops.galileon_vertex_word(p, gram)
    )


def generic_radiative_data():
    a = s.symbols("hard_a0:5")
    b = s.symbols("TT_b0:6")
    pairs = tuple((i, j) for i in range(4) for j in range(i + 1, 4))
    gram = s.eye(5)
    gram[4, 4] = 0
    for (i, j), value in zip(pairs, (*a, -2 - s.Add(*a))):
        gram[i, j] = gram[j, i] = value
    for i in range(4):
        gram[i, 4] = gram[4, i] = -sum(gram[i, j] for j in range(4))
    polarization = s.zeros(5)
    for (i, j), value in zip(pairs, b):
        polarization[i, j] = polarization[j, i] = value
    for i in range(4):
        polarization[i, i] = -sum(polarization[i, j] for j in range(4) if j != i)
    return gram, polarization
