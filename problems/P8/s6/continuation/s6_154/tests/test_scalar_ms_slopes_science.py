"""Independent regulated sunset, kernel-derivative and OS remainder checks."""

from functools import cache

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_scalar_ms_slopes import (
    audit,
    bounds,
    calibration,
    cauchy,
    conversion,
    sunset,
)
from p8_vacuum_two_loop_light_pole import graphs


@pytest.mark.parametrize("name", tuple(audit.residuals()))
def test_exact_identity(name):
    assert audit.residuals()[name] == 0


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[v[0] for v in audit.bad_cases()]
)
def test_invalid_input(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@cache
def gauss(n):
    with mp.workdps(40):
        nodes, weights = mp.gauss_quadrature(n, "legendre")
        return tuple(((nodes[k] + 1) / 2, weights[k] / 2) for k in range(n))


@cache
def samples():
    with mp.workdps(36):
        out = []
        for u, wu in gauss(48):
            r = u**8
            for t, wt in gauss(32):
                a, v = 1 + r * (1 + t), 1 + t + r * t
                U = r * v / a**2
                P = r * t / (a * v)
                jac = 6 * 8 * u**7 * wu * wt
                out.append((U, P, jac * t / v**3, jac * r / a**3))
        return tuple(out)


@pytest.mark.parametrize("e", (mp.mpf(5) / 8, mp.mpf(3) / 4))
@pytest.mark.parametrize("z", (0, s.Rational(1, 2), 1))
def test_full_dimensional_sunset_derivative_and_inverse_sign(e, z):
    with mp.workdps(32):
        z = mp.mpf(str(z)) if z != s.Rational(1, 2) else mp.mpf(".5")
        ell = mp.mpf(2)

        # Raw Pi is positive times Gamma(-1+2e); no OS projection.
        def raw(u):
            total = sum(
                jac * U ** (-2 + e) * (1 - u * P) ** (1 - 2 * e)
                for U, P, weight, jac in samples()
            )
            return mp.exp(2 * mp.euler * e + 2 * ell * e) * mp.gamma(-1 + 2 * e) * total

        actual = mp.diff(raw, z)
        derivative = (
            mp.exp(2 * mp.euler * e + 2 * ell * e)
            * mp.gamma(2 * e)
            * sum(
                weight * U**e * (1 - z * P) ** (-2 * e)
                for U, P, weight, jac in samples()
            )
        )
        assert actual > 0
        assert abs(actual - derivative) < mp.mpf("1e-25")


@cache
def independent_moments():
    with mp.workdps(30):

        def integral(which):
            def outer(t):
                def inner(r):
                    v, a = 1 + t + r * t, 1 + r * (1 + t)
                    U, P = r * v / a**2, r * t / (a * v)
                    factor = mp.log(U) if which == "T1" else -mp.log1p(-P)
                    return 6 * t / v**3 * factor

                return mp.quad(inner, [0, 1])

            return mp.quad(outer, [0, 1])

        return integral("T1"), integral("C")


@pytest.mark.parametrize("ell", (0, 2, -1))
def test_complex_regulator_finite_part_without_fitting_MS_reference(ell):
    with mp.workdps(30):
        T1, C = independent_moments()
        expected = mp.mpf(ell) / 2 + T1 / 2 + C
        pole, finite = 0, 0
        for k in range(16):
            e = mp.mpf(".025") * mp.exp(2j * mp.pi * k / 16)
            integral = sum(
                weight * mp.exp(e * (mp.log(U) - 2 * mp.log1p(-P)))
                for U, P, weight, jac in samples()
            )
            regular = (
                mp.exp(2 * mp.euler * e + 2 * ell * e)
                * mp.gamma(1 + 2 * e)
                * integral
                / 2
            )
            pole += regular / 16
            finite += regular / e / 16
        assert abs(pole - mp.mpf(1) / 4) < mp.mpf("1e-16")
        assert abs(finite - expected) < mp.mpf("1e-15")


def test_local_integral_enclosures_use_positive_moments_not_numeric_fit():
    with mp.workdps(28):
        T1, C = independent_moments()
        assert -6 < T1 < 0
        assert 0 < C < mp.mpf(1) / 16
        assert abs(sum(w for U, P, w, jac in samples()) - mp.mpf(".5")) < mp.mpf(
            "1e-17"
        )
        for U, P, w, jac in samples():
            assert 0 < U <= mp.mpf(1) / 3
            assert 0 < P <= mp.mpf(1) / 9
            assert w > 0


@pytest.mark.parametrize("M", (33, 40, 1000, 1000000))
@pytest.mark.parametrize("z", (1, mp.mpc(1, ".5"), mp.mpc("1.5", ".5")))
def test_full_mixed_and_squared_heavy_OS_remainders(M, z):
    with mp.workdps(32):
        M, z = mp.mpf(M), mp.mpc(z)
        v = z - 1
        grid = [0, 1 / M, 1]
        D = lambda x: x * M + (1 - x) ** 2
        b = lambda x: x * (1 - x) / D(x)
        mixed = mp.quad(lambda x: -mp.log1p(-b(x) * v) - b(x) * v, grid)
        squared = mp.quad(
            lambda x: (
                x / (D(x) - x * (1 - x) * v)
                - x / D(x)
                - v * x * x * (1 - x) / D(x) ** 2
            ),
            grid,
        )
        geometric = mp.quad(
            lambda x: x * b(x) ** 2 * v**2 / (D(x) * (1 - b(x) * v)), grid
        )
        assert abs(squared - geometric) < mp.mpf("1e-27")
        assert abs(mixed) <= abs(v) ** 2 / (6 * M**2 * (1 - 1 / M)) + mp.mpf("1e-29")
        assert abs(squared) <= abs(v) ** 2 / (3 * M**3 * (1 - 1 / M)) + mp.mpf("1e-29")


@pytest.mark.parametrize("M", (33, 40, 1000))
def test_mass_derivative_fixes_squared_heavy_slope(M):
    with mp.workdps(32):
        M = mp.mpf(M)

        def first(m):
            return mp.quad(
                lambda x: x * (1 - x) / (x * m + (1 - x) ** 2), [0, 1 / m, 1]
            )

        squared = mp.quad(
            lambda x: x * x * (1 - x) / (x * M + (1 - x) ** 2) ** 2, [0, 1 / M, 1]
        )
        assert abs(-mp.diff(first, M) - squared) < mp.mpf("1e-27")
        assert 0 < first(M) < 1 / (2 * M)
        assert 0 < squared < 1 / (2 * M**2)


@pytest.mark.parametrize("kind,choice", graphs.cases())
def test_every_parent_graph_remains_in_complete_partition(kind, choice):
    label = kind + "_" + graphs.group(kind, choice)
    assert label in cauchy.data()["same_full_graph_partition"]
    graph = graphs.refine(kind, choice)
    assert len(graph["edges"]) - len(graph["vertices"]) + 1 == 2
    if kind == "sunset" and choice == (0, 0):
        ends = [n for n in graph["vertices"] if graph["external"][n]]
        for core in graphs.uv(graph):
            if core[1] == 1:
                vertices = {n for j in core[0] for n in graph["edges"][j]}
                assert set(ends) <= vertices


@pytest.mark.parametrize("alpha1", (1, -2, s.Rational(3, 7)))
def test_nested_slope_has_no_pole_product_but_outer_mass_can(alpha1):
    e, a0, ell, p = s.symbols("e alpha0 ell mass_pole")
    alpha = a0 + e * alpha1
    assert s.limit(alpha**2, e, 0) == a0**2
    assert s.limit(alpha * (p / e + ell) - a0 * p / e, e, 0) == a0 * ell + alpha1 * p


def test_early_outer_OS_would_erase_the_required_MS_slope():
    z, a, b, c = s.symbols("z a b c")
    f = a + b * (z - 1) + c * (z - 1) ** 2
    OS = f - f.subs(z, 1) - (z - 1) * s.diff(f, z).subs(z, 1)
    assert s.diff(f, z).subs(z, 1) == b
    assert s.diff(OS, z).subs(z, 1) == 0
    assert sunset.data()["UV_slope_pole"] != 0


def test_actual_complete_bound_does_not_normalize_MS_field_automatically():
    d = calibration.data()
    assert all(d["bounds"].values())
    assert (
        len(d["actual_scalar_MS_slope_enclosure"]["finite_MS_slope_group_uppers"]) == 7
    )
    assert len(audit.frontier()) == 9
    assert len(audit.matching()) == 12
    assert sum(r["status"] != "OPEN" for r in audit.matching()) == 5
    assert len(audit.residuals()) == 55
    assert len(audit.gates()) == 23
    assert audit.controls()["rejected_inputs"] == 77
    assert conversion.data()["MS_outer_OS_quadratic_conversion_upper"] != 0
    assert bounds.enclosure(1, 1, 40)["finite_MS_slope_absolute_upper"] > 0
