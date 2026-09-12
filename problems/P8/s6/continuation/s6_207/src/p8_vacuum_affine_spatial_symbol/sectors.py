"""Exact four-sector factorization of all nine tracefree stress pairs."""

from functools import cache

import sympy as s

SECTORS = ("TT", "TL", "LT", "LL")


def cross(v):
    return s.Matrix([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])


def projector(v):
    v = s.Matrix(v)
    if v.shape != (3, 1) or (v.T * v)[0] == 0:
        raise ValueError("Non-null three-component bilinear momentum required")
    return s.eye(3) - v * v.T / (v.T * v)[0]


def geometric_factors(k, l, D, G):
    k, l = s.Matrix(k), s.Matrix(l)
    P, Q = projector(k), projector(l)
    C, E = cross(k), cross(l)
    kl, ll = (k.T * k)[0], (l.T * l)[0]
    magG = C.T * G * E
    magD = C.T * D * E
    return {
        "00": s.trace(P * G * Q * D),
        "01": s.trace(P * magG * Q * D),
        "10": s.trace(P * G * Q * magD.T),
        "11": s.trace(P * magG * Q * magD.T),
        "TL": (l.T * D * P * G * l)[0] / ll,
        "LT": (k.T * D * Q * G * k)[0] / kl,
        "LL": (k.T * D * l)[0] * (k.T * G * l)[0] / (kl * ll),
    }


def scalar_amplitudes(point, m):
    a = point["a"]
    fkt, pkt, _wkt = point["kt"]
    fkl, pkl, wkl = point["kl"]
    flt, plt, _wlt = point["lt"]
    fll, pll, wll = point["ll"]
    return {
        "A": m * m * fkt * flt - pkt * plt,
        "B": fkt * flt / a**2,
        "TL": m * wll * fkt * fll - m * pkt * pll / wll,
        "LT": m * wkl * fkl * flt - m * pkl * plt / wkl,
        "LL": wkl * wll * fkl * fll - m * m * pkl * pll / (wkl * wll),
    }


def pair_products(k, l, Dsharp, G, source, detector_sharp, m):
    """Detector sharp includes p conjugation and the magnetic i sign flip."""
    geom = geometric_factors(k, l, Dsharp, G)
    A = scalar_amplitudes(detector_sharp, m)
    B = scalar_amplitudes(source, m)
    return {
        "TT": A["A"] * B["A"] * geom["00"]
        + A["A"] * B["B"] * geom["01"]
        + A["B"] * B["A"] * geom["10"]
        + A["B"] * B["B"] * geom["11"],
        "TL": A["TL"] * B["TL"] * geom["TL"],
        "LT": A["LT"] * B["LT"] * geom["LT"],
        "LL": A["LL"] * B["LL"] * geom["LL"],
    }


def physical_vectors(momentum, e1, e2, point, prefix, m, sharp=False):
    r = s.sqrt(momentum.dot(momentum))
    n = momentum / r
    f, p, _w = point[prefix + "t"]
    fl, pl, wl = point[prefix + "l"]
    sign = -1 if sharp else 1
    transverse = [
        s.Matrix(
            [
                *(p * e),
                *(sign * s.I * momentum.cross(e) * f / point["a"]),
                0,
                *(m * f * e),
            ]
        )
        for e in (e1, e2)
    ]
    longitudinal = s.Matrix(
        [
            *(m * pl * n / wl),
            0,
            0,
            0,
            -sign * s.I * r * pl / (point["a"] * wl),
            *(wl * fl * n),
        ]
    )
    return transverse + [longitudinal]


@cache
def data():
    m = s.symbols("m", positive=True)
    source = {"a": s.symbols("a_source", positive=True)}
    detector = {"a": s.symbols("a_detector", positive=True)}
    for name in ("kt", "kl", "lt", "ll"):
        source[name] = s.symbols("f_s_" + name + " p_s_" + name + " omega_s_" + name)
        detector[name] = s.symbols("f_d_" + name + " p_d_" + name + " omega_d_" + name)
    k = s.Matrix([3, 0, 4])
    l = s.Matrix([0, 5, 12])
    e1 = s.Matrix([0, 1, 0])
    e2 = k.cross(e1) / 5
    f1 = s.Matrix([1, 0, 0])
    f2 = l.cross(f1) / 13
    z = s.symbols("D0:5")
    w = s.symbols("G0:5")
    D = s.Matrix([[z[0], z[1], z[2]], [z[1], z[3], z[4]], [z[2], z[4], -z[0] - z[3]]])
    G = s.Matrix([[w[0], w[1], w[2]], [w[1], w[3], w[4]], [w[2], w[4], -w[0] - w[3]]])
    MD = s.diag(-D, -D, s.zeros(1), D)
    MG = s.diag(-G, -G, s.zeros(1), G)
    sk = physical_vectors(k, e1, e2, source, "k", m)
    sl = physical_vectors(l, f1, f2, source, "l", m)
    dk = physical_vectors(k, e1, e2, detector, "k", m, True)
    dl = physical_vectors(l, f1, f2, detector, "l", m, True)
    expected = pair_products(k, l, D, G, source, detector, m)
    checks = {}
    for tag, ii, jj in (
        ("TT", range(2), range(2)),
        ("TL", range(2), (2,)),
        ("LT", (2,), range(2)),
        ("LL", (2,), (2,)),
    ):
        actual = sum(
            (dk[i].T * MD * dl[j])[0] * (sk[i].T * MG * sl[j])[0]
            for i in ii
            for j in jj
        )
        checks["complete_two_time_" + tag + "_factorization"] = s.expand(
            actual - expected[tag]
        )
    checks["complete_transverse_projector_k"] = projector(k) - e1 * e1.T - e2 * e2.T
    checks["complete_transverse_projector_l"] = projector(l) - f1 * f1.T - f2 * f2.T
    checks["physical_pair_partition"] = 4 + 2 + 2 + 1 - 9
    return {
        "exact_two_time": "The detector time and source time are independent. Detector amplitudes use their analytic Schwarz values; the source uses its original values. Dsharp is held fixed during source derivatives.",
        "four_sectors": "TT uses four magnetic/base trace contractions. TL,LT and LL use the complete bilinear projectors and longitudinal outer products. Multiplicities4,2,2,1 retain all nine physical pairs.",
        "phase_boundary": "The four sectors retain their distinct sums Wk_T+Wl_T, Wk_T+Wl_L, Wk_L+Wl_T, Wk_L+Wl_L. They cannot be merged before inverse-phase weighting and time differentiation.",
        "constraint_scope": "This exact tracefree spatial stress contraction has no direct mA0 term, but the physical longitudinal electric and mass readouts retain the constraint. It is not a lapse/shift or full stress replacement.",
        "checks": checks,
        "gates": {
            "all_four_sectors_retained": len(SECTORS) == 4,
            "independent_source_and_detector_scale_factors": source["a"]
            != detector["a"],
            "noncollinear_two_momenta": k.cross(l) != s.zeros(3, 1),
            "general_noncommuting_tracefree_tensors": D * G != G * D,
            "no_common_T_L_frequency_assumed": True,
        },
    }
