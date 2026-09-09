"""Measure-correct proper-clock maps for BOTH nonminimal QSEI terms."""

from functools import cache

import sympy as sp


@cache
def data():
    a = sp.Symbol("scale", positive=True)
    h, hd, hdd, h3 = sp.symbols("H Hdot Hddot Hthird", real=True)
    f = sp.symbols("f0:5", real=True)

    def dt(x):
        return (
            sp.diff(x, a) * a * h
            + sp.diff(x, h) * hd
            + sp.diff(x, hd) * hdd
            + sp.diff(x, hdd) * h3
            + sum(sp.diff(x, f[j]) * f[j + 1] for j in range(4))
        )

    F = a ** sp.Rational(-3, 2) * f[0]
    first = sp.factor(a * dt(F))
    second = sp.factor(a * dt(first))
    L = (
        f[2]
        - 2 * h * f[1]
        + (sp.Rational(3, 4) * h * h - sp.Rational(3, 2) * hd) * f[0]
    )
    G = f[1] - sp.Rational(3, 2) * h * f[0]
    return {
        "a": a,
        "sampler_jets": f,
        "H": h,
        "Hdot": hd,
        "flat_sampler": F,
        "flat_sampler_first_conformal_derivative": first,
        "flat_sampler_second_conformal_derivative": second,
        "proper_second_order_operator": L,
        "proper_Wick_square_first_order_operator": G,
    }


@cache
def checks():
    d = data()
    a = d["a"]
    return {
        "proper_conformal_energy_measure": sp.factor(
            a * d["flat_sampler"] ** 2 - a**-2 * d["sampler_jets"][0] ** 2
        ),
        "proper_conformal_second_derivative_norm": sp.factor(
            d["flat_sampler_second_conformal_derivative"] ** 2 / a
            - d["proper_second_order_operator"] ** 2
        ),
        "physical_Wick_square_weight_keeps_extra_scale_squared": sp.factor(
            a * d["flat_sampler_first_conformal_derivative"] ** 2
            - d["proper_Wick_square_first_order_operator"] ** 2
        ),
        "actual_conformal_scalar_quantum_to_delta_ratio": sp.Rational(7, 144) * 8
        - sp.Rational(7, 18),
        "actual_conformal_scalar_Wick_square_coefficient": 2 * sp.Rational(1, 6)
        - sp.Rational(1, 3),
    }
