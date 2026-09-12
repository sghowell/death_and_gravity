from functools import cache

import mpmath as mp
import sympy as s
from p8_vacuum_affine_ordered_scalar_symbol import geometry, jets
from p8_vacuum_affine_spatial_symbol import sectors
from p8_vector_clock_matching import continuation
from p8_vector_state import wkb


@cache
def full_dimensional_WKB_coefficient(kind, order):
    if order == 0:
        return s.Integer(1)
    u, z = wkb.u, wkb.z
    U = continuation.coefficients(kind)["U"]
    lam = wkb.background()["lambda"]

    def time_slope(value):
        return s.factor(
            s.diff(value, u) + wkb.background()["z_prime"] * s.diff(value, z)
        )

    previous = [full_dimensional_WKB_coefficient(kind, n) for n in range(order)]
    inverse = [s.Integer(1)]
    for n in range(1, order):
        inverse.append(
            s.factor(-sum(previous[j] * inverse[n - j] for j in range(1, n + 1)))
        )
    DS = [s.Integer(0)] + [
        time_slope(previous[n]) - 2 * n * lam * previous[n] for n in range(1, order)
    ]
    rate = [lam] + [
        s.factor(sum(DS[j] * inverse[n - j] for j in range(1, n + 1)))
        for n in range(1, order)
    ]
    value = (
        -sum(previous[j] * previous[order - j] for j in range(1, order))
        - time_slope(rate[-1]) / 2
        + (order - 1) * lam * rate[-1]
        + sum(rate[j] * rate[order - 1 - j] for j in range(order)) / 4
        - (U if order == 1 else 0)
    )
    return s.factor(value / 2)


@cache
def dimensional_coefficient_function(kind):
    expressions = []
    for order in range(1, 5):
        value = full_dimensional_WKB_coefficient(kind, order)
        expressions.extend((value, s.diff(value, wkb.u), s.diff(value, wkb.z)))
    return s.lambdify(
        (wkb.u, wkb.z, continuation.local.dimension), expressions, "mpmath", cse=True
    )


@cache
def amplitude_function():
    a, mass = s.symbols("a mass")
    point = {"a": a}
    args = [a, mass]
    for key in ("kt", "kl", "lt", "ll"):
        point[key] = s.symbols("f_" + key + " p_" + key + " o_" + key)
        args.extend(point[key])
    expressions = dict(sectors.scalar_amplitudes(point, mass))
    _fk, pk, ok = point["kl"]
    _fl, pl, ol = point["ll"]
    expressions["C"] = -pk * pl / (2 * a * a * ok * ol)
    return s.lambdify(args, list(expressions.values()), "mpmath", cse=True)


def scalar_values(point, mass):
    args = [point["a"], mass]
    for key in ("kt", "kl", "lt", "ll"):
        args.extend(point[key])
    return amplitude_function()(*args)


def scaled_point(u, x, k, l, phase_sign=-1, orders=4, dimension=3):
    a = (1 + u * u) ** 2
    H = 4 * u / (1 + u * u)
    point = {"a": a}
    frequencies = {}
    for prefix, v in (("k", k), ("l", l)):
        omega = mp.sqrt((v.T * v)[0] / a**2 + 1000**2 * x * x)
        z = 1 - 1000**2 * x * x / omega**2
        domega = -H * z * omega
        dz = -2 * H * z * (1 - z)
        for suffix, kind in (("t", "transverse"), ("l", "longitudinal")):
            coeff = dimensional_coefficient_function(kind)(u, z, dimension)
            W, dW = omega, domega
            for q in range(1, orders + 1):
                P, du, dzP = coeff[3 * (q - 1) : 3 * q]
                W += P * x ** (2 * q) * omega ** (1 - 2 * q)
                dW += x ** (2 * q) * (
                    (du + dzP * dz) * omega ** (1 - 2 * q)
                    + (1 - 2 * q) * P * omega ** (-2 * q) * domega
                )
            f = 1 / mp.sqrt(2 * W)
            rate = (
                (dimension - 2) * H / 2
                if suffix == "t"
                else ((dimension - 2) / 2 + z) * H
            )
            p = (phase_sign * 1j * W - x * (dW / (2 * W) + rate)) * f
            point[prefix + suffix] = (f, p, omega)
            frequencies[prefix + suffix] = W
    return point, frequencies


def four_inverse_phases(freq):
    return [
        1 / (freq[a] + freq[b])
        for a, b in (("kt", "lt"), ("kt", "ll"), ("kl", "lt"), ("kl", "ll"))
    ]


GEOS = ("00", "01", "10", "11", "TL", "LT", "LL", "LC", "CL", "CC")
PAIRING = (
    (0, 0, 0),
    (0, 1, 0),
    (1, 0, 0),
    (1, 1, 0),
    (2, 2, 1),
    (3, 3, 2),
    (4, 4, 3),
    (4, 5, 3),
    (5, 4, 3),
    (5, 5, 3),
)


def mode_endpoint_rows(x, time, angle, transfer, dimension, samples=20):
    n = mp.matrix([mp.sqrt(1 - angle * angle), 0, angle])
    k = n
    ell = -n + mp.matrix([0, 0, x * transfer])
    detector, df = scaled_point(time, x, k, ell, 1, 4, dimension)
    detector_amps = scalar_values(detector, 1000 * x)
    gd = four_inverse_phases(df)
    clock_radius = mp.mpf("1e-5")
    phases = [mp.exp(2j * mp.pi * i / samples) for i in range(samples)]
    source_amps = []
    source_g = []
    for phase in phases:
        source, freq = scaled_point(
            time + clock_radius * phase, x, k, ell, -1, 4, dimension
        )
        source_amps.append(scalar_values(source, 1000 * x))
        source_g.append(four_inverse_phases(freq))
    ampjet = [
        [
            mp.fsum(source_amps[i][a] * phases[i] ** (-r) for i in range(samples))
            / samples
            / clock_radius**r
            for r in range(5)
        ]
        for a in range(6)
    ]
    gjet = [
        [
            mp.fsum(source_g[i][a] * phases[i] ** (-r) for i in range(samples))
            / samples
            / clock_radius**r
            for r in range(5)
        ]
        for a in range(4)
    ]
    out = {}
    for geo, (left, right, sector) in zip(GEOS, PAIRING):
        g = gjet[sector]
        for r in range(5):
            row = [
                ampjet[right][d - r] / mp.factorial(r) if d >= r else mp.mpc(0)
                for d in range(5)
            ]
            for j in range(5):
                if r <= j:
                    out[geo, j, r] = (
                        1j * (-1j) ** j * gd[sector] * detector_amps[left] * row[0]
                    )
                row = [
                    (d + 1) * mp.fsum(g[h] * row[d + 1 - h] for h in range(d + 2))
                    for d in range(len(row) - 1)
                ]
    return out


@cache
def complete_dimensional_endpoint_function():
    rows = jets.endpoint_products()
    keys = [(key, degree) for key, row in rows.items() for degree in range(len(row))]
    fn = s.lambdify(
        (geometry.d, jets.t, jets.u, jets.p, jets.m),
        [rows[key][degree] for key, degree in keys],
        "mpmath",
        cse=True,
    )
    return keys, fn
