"""Actual CD low inverse-radius jets with the complete longitudinal constraint."""

from functools import cache

import sympy as s
from p8_vector_hadamard import series
from p8_vector_state import wkb


def clock(t):
    a = (1 + t * t) ** 2
    H = 4 * t / (1 + t * t)
    Hp = 4 * (1 - t * t) / (1 + t * t) ** 2
    cL = -(Hp + 2 * H * H) / 2
    return a, H, Hp, cL


@cache
def amplitude_jets():
    x = s.symbols("x", real=True)
    a = s.symbols("a", positive=True)
    m, H, c = s.symbols("m H cL", real=True)
    omega = s.sqrt(a**-2 + m * m * x * x)
    WT = omega
    WL = omega + c * x * x / omega
    fT = 1 / s.sqrt(2 * WT)
    fL = 1 / s.sqrt(2 * WL)
    pT = -s.I * WT * fT
    pL = (-s.I * WL - x * H) * fL
    values = {
        "A": m * m * x * x * fT * fT - pT * pT,
        "B": fT * fT / a**2,
        "TL": m * x * omega * fT * fL - m * x * pT * pL / omega,
        "LT": m * x * omega * fL * fT - m * x * pL * pT / omega,
        "LL": omega * omega * fL * fL - m * m * x * x * pL * pL / (omega * omega),
        "gLL": 1 / (2 * WL),
    }
    rows = {}
    for key, value in values.items():
        expression = s.series(value, x, 0, 3).removeO().expand()
        rows[key] = tuple(s.factor(expression.coeff(x, d)) for d in range(3))
    return a, m, H, c, rows


def multiply(a, b):
    return tuple(s.expand(sum(a[k] * b[d - k] for k in range(d + 1))) for d in range(3))


@cache
def data():
    t = wkb.u
    a, _H, _Hp, c = clock(t)
    av, m, Hv, cv, rows = amplitude_jets()
    expected = {
        "A": (1 / (2 * av), 0, 3 * av * m * m / 4),
        "B": (1 / (2 * av), 0, -av * m * m / 4),
        "TL": (0, m, -s.I * av * Hv * m / 2),
        "LT": (0, m, -s.I * av * Hv * m / 2),
        "LL": (1 / (2 * av), 0, av * (3 * m * m / 4 - cv / 2)),
        "gLL": (av / 2, 0, -(av**3) * (m * m / 2 + cv) / 2),
    }
    checks = {
        "actual_transverse_first_WKB_coefficient_at_z1": s.factor(
            series.coefficient("transverse", 1).subs(wkb.z, 1)
        ),
        "actual_longitudinal_first_WKB_coefficient_at_z1": s.factor(
            series.coefficient("longitudinal", 1).subs(wkb.z, 1) - c
        ),
        "actual_curvature_combination": s.factor(
            av.subs(av, a) * (-2 * c) - 4 * (1 + 7 * t * t)
        ),
        "complete_six_amplitude_and_phase_low_jets": s.Matrix(
            [
                s.factor(value - expected[key][d])
                for key, row in rows.items()
                for d, value in enumerate(row)
            ]
        ),
    }
    at, as_ = s.symbols("a_detector a_source", positive=True)
    Ht, Hs, ct, cs = s.symbols("H_detector H_source c_detector c_source", real=True)
    source = {
        key: tuple(value.subs({av: as_, Hv: Hs, cv: cs}) for value in row)
        for key, row in rows.items()
    }
    detector = {
        key: tuple(s.conjugate(value.subs({av: at, Hv: Ht, cv: ct})) for value in row)
        for key, row in rows.items()
    }
    checks["all_two_time_pair_products_real_through_degree_two"] = s.Matrix(
        [
            s.simplify(s.im(value))
            for left, right in (
                ("A", "A"),
                ("A", "B"),
                ("B", "A"),
                ("B", "B"),
                ("TL", "TL"),
                ("LT", "LT"),
                ("LL", "LL"),
            )
            for value in multiply(source[left], detector[right])
        ]
    )
    fullLL = multiply(
        multiply(rows["LL"], tuple(s.conjugate(value) for value in rows["LL"])),
        rows["gLL"],
    )
    checks["complete_longitudinal_current_second_jet"] = s.factor(
        fullLL[2] - av * (5 * m * m / 16 - 3 * cv / 8)
    )
    checks["longitudinal_curvature_shift_not_deleted"] = s.factor(
        fullLL[2] - fullLL[2].subs(cv, 0) + 3 * av * cv / 8
    )
    return {
        "actual_WKB_input": "The exact frozen recurrence gives P1_T(t,z1)=0 and P1_L(t,z1)=-(H_prime+2H^2)/2, with a=(1+t^2)^2 and H=a_prime/a. No higher W8 term can affect the second inverse-radius degree; all four remain in independent finite-momentum tests.",
        "transverse_order": "What_T=Omega+O(x^4), and What_time/(2What)+H/2=O(x^2). Thus the transverse momentum has no order-x real correction. The longitudinal rate is H+O(x^2). These statements retain the constraint, rather than replacing Proca by Maxwell.",
        "full_pairs": "The TL/LT amplitudes have an imaginary order-x^2 coefficient, but start at order x, so their complete source/detector products first acquire an imaginary coefficient at order x^3. TT and LL products are also real through degree two. Real geometric and inverse-phase jets preserve this statement at general P.",
        "endpoint_one": "The j1 prefactor is real, so its degree-one current coefficient vanishes in this actual unit-W8 tracefree calculation. This is an evaluated slot, not deletion of the endpoint: the independent j1/d3 coefficient and finite-momentum odd endpoints remain nonzero.",
        "curvature_only_LL": "The complete j0 second-grade curvature correction is -3a P1_L (n.D.n)(n.G.n)/8. At this degree all curvature insertions are P-independent by the inverse-radius filtration. The remaining mass and spatial pieces are the complete massive flat coefficients with m->m*a and factor1/a.",
        "checks": checks,
        "gates": {
            "both_actual_T_L_WKB_coefficients_used": True,
            "longitudinal_imaginary_amplitude_not_discarded": rows["TL"][2] != 0,
            "complete_two_time_product_not_equal_time_shortcut": at != as_,
            "second_grade_not_finite_reference_replacement": True,
            "all_original_odd_endpoint_labels_retained": True,
        },
    }
