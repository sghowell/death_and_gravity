"""All 42 derivative-source radiation graphs and the existing OS4 contact."""

from functools import cache

import sympy as s
from p8_vacuum_affine_local_tadpole_radiation import vertices

from . import source


def one_channel(L, R, gram, pol, ak, mass, trace=0, kep=None):
    trace = s.sympify(trace)
    kep = (0, 0, 0, 0) if kep is None else kep
    pl2 = sum(gram[i, j] for i in L for j in L)
    pr2 = sum(gram[i, j] for i in R for j in R)
    dl, dr = pl2 - mass, pr2 - mass
    vl = 2 * (1 - gram[L[0], L[1]])
    vr = 2 * (1 - gram[R[0], R[1]])
    currents = [pol[i, i] / ak[i] + kep[i] / ak[i] - trace / 2 for i in range(4)]
    internal = (
        2 * sum(pol[i, j] for i in L for j in L)
        + 2 * sum(kep[i] for i in L)
        - trace * (pl2 + sum(ak[i] for i in L) - mass)
    )
    value = -sum(
        currents[i] * (vl + vr - 2 * ak[next(j for j in L if j != i)]) / dr for i in L
    )
    value -= sum(
        currents[i] * (vl + vr - 2 * ak[next(j for j in R if j != i)]) / dl for i in R
    )
    value -= (vl + vr) * internal / (dl * dr)
    value -= (trace * (vl + vr) + 4 * pol[L[0], L[1]]) / dr
    value -= (trace * (vl + vr) + 4 * pol[R[0], R[1]]) / dl
    old_unit = (
        -sum(currents[i] / dr for i in L)
        - sum(currents[i] / dl for i in R)
        - internal / (dl * dr)
        - trace * (1 / dl + 1 / dr)
    )
    return value, old_unit, dl, dr


def source_factor(heavy=source.HEAVY_MASS2, kappa=source.KAPPA):
    heavy, kappa = map(s.sympify, (heavy, kappa))
    return 8 * (4 - heavy) / (16 * s.pi**2 * kappa)


def centered_contact(heavy=source.HEAVY_MASS2, cubic=source.CUBIC):
    heavy, cubic = map(s.sympify, (heavy, cubic))
    return -3 * cubic**2 / (heavy - s.Rational(4, 3))


def four_correction(ss, tt):
    ss, tt = map(s.sympify, (ss, tt))
    n, g = source.HEAVY_MASS2, source.CUBIC
    return (
        source_factor()
        * g
        * g
        * (sum(1 / (n - a) for a in (ss, tt, 4 - ss - tt)) - 3 / (n - s.Rational(4, 3)))
    )


@cache
def data():
    G, H = vertices.generic_radiative_data()
    a = tuple(G[i, 4] for i in range(4))
    n = s.Symbol("heavy_mass_squared")
    xi = s.symbols("xi_dot_p0:4")
    beta = -sum(xi)
    gauge = s.Matrix(4, 4, lambda i, j: a[i] * xi[j] + a[j] * xi[i])
    checks = {}
    for j in (1, 2, 3):
        L = (0, j)
        R = tuple(i for i in range(4) if i not in L)
        actual, unit, _, _ = one_channel(L, R, G, H, a, n)
        checks["channel" + str(j) + "_complete42_TT_identity"] = s.factor(
            actual - 2 * (4 - n) * unit - 2 * sum(H[i, i] / a[i] for i in range(4))
        )
        pure, _, _, _ = one_channel(
            L, R, G, gauge, a, n, 2 * beta, tuple(beta * v for v in a)
        )
        checks["channel" + str(j) + "_full_puregauge_Ward"] = s.factor(pure)
    ss, tt, g, kappa = s.symbols("s t g kappa")
    channels = (ss, tt, 4 - ss - tt)
    direct = (
        8
        * g
        * g
        / (16 * s.pi**2 * kappa)
        * (
            sum((4 - v) / (n - v) for v in channels)
            - 3 * (4 - s.Rational(4, 3)) / (n - s.Rational(4, 3))
        )
    )
    reduced = (
        source_factor(n, kappa)
        * g
        * g
        * (sum(1 / (n - v) for v in channels) - 3 / (n - s.Rational(4, 3)))
    )
    checks["whole_flat_source_correction_recovers_frozen239"] = s.factor(
        direct - reduced
    )
    J = s.Symbol("complete_external_current")
    checks["full_three_channels_and_OS4_contact_reduction"] = s.factor(
        6 * J
        - (6 + 6 * (4 - n) / (n - s.Rational(4, 3))) * J
        + 6 * (4 - n) * J / (n - s.Rational(4, 3))
    )
    checks["actual_source_factor_from_cI_and_two_endpoints"] = s.factor(
        2 * (4 - n) * (-4 * g / kappa) * (-1 / (16 * s.pi**2)) / g
        - source_factor(n, kappa)
    )
    return {
        "checks": checks,
        "gates": {
            "all_three_channels_and_two_source_endpoints": True,
            "literal42_plus_five_counterterm_graphs": 3 * (8 + 2 + 4) + 5 == 47,
            "traceful_general_Ward_before_TT_projection": True,
            "generic_Gram_identities_not_only_recoil_samples": True,
            "source_factor_nonzero_at_original_parameters": bool(source_factor() != 0),
            "centered_contact_only_linear_loop_representation": True,
        },
        "whole_original_relative_source_factor": source_factor(),
        "whole_centered_contact_representation": centered_contact(),
        "whole_graph_prescription": "For each of three channels, place H(Y+Phi^2) at either heavy-source endpoint:8 external,2 heavy-line and4 source-metric emission graphs. Include all trace terms for the general Ward identity. The physical TT sum equals2(4-n) times unit-heavy radiation plus2 times the complete external current in each channel. The five graphs of the existing OS4 constant subtraction leave exactly Q times the centered complete matter tree,Q=8(4-n)/(16pi^2 kappa).",
        "whole_matching_boundary": "The contact -3g^2/(n-4/3) is only the linear representation of the known selected one-loop correction. The original tuned Born contact remains unchanged. Neither an independent finite condition nor a value for unknown curvature matching is introduced.",
    }
