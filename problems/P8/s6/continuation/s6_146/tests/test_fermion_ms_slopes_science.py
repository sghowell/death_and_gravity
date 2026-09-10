"""Independent dimensional, matrix, Gamma, mass-difference and scope checks."""

import copy

import mpmath as mp
import pytest
import sympy as sp
from p8_vacuum_fermion_ms_slopes import (
    audit,
    calibration,
    masters,
    reference,
    remainder,
)
from p8_vacuum_fermion_ms_slopes import tensors as t


@pytest.mark.parametrize("name,value", list(audit.residuals().items()))
def test_exact_identity(name, value):
    assert value == 0, name


@pytest.mark.parametrize("name,value", list(audit.gates().items()))
def test_proof_gate(name, value):
    assert value is True, name


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_invalid_input(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def gamma_matrices():
    sigma = (
        mp.matrix([[0, 1], [1, 0]]),
        mp.matrix([[0, -1j], [1j, 0]]),
        mp.matrix([[1, 0], [0, -1]]),
    )
    g = []
    for s in sigma:
        g.append(
            mp.matrix(
                [
                    [0, 0, -1j * s[0, 0], -1j * s[0, 1]],
                    [0, 0, -1j * s[1, 0], -1j * s[1, 1]],
                    [1j * s[0, 0], 1j * s[0, 1], 0, 0],
                    [1j * s[1, 0], 1j * s[1, 1], 0, 0],
                ]
            )
        )
    g.append(mp.matrix([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]]))
    g.append(g[0] * g[1] * g[2] * g[3])
    return g


@pytest.mark.parametrize("dim", (3, 4, 5))
@pytest.mark.parametrize("config", (0, 1))
@pytest.mark.parametrize(
    "kind", ("scalar_self", "scalar_vertex", "gauge_self", "gauge_vertex")
)
def test_explicit_matrix_Laplacian(dim, config, kind):
    with mp.workdps(45):
        gam = gamma_matrices()
        trace = lambda a: sum(a[i, i] for i in range(4))
        kvals = [sp.Rational(j + 1 + config, 5) for j in range(dim)]
        lvals = [sp.Rational((-1) ** j * (j + 2), 7 + config) for j in range(dim)]
        k = mp.matrix([mp.mpf(int(v.p)) / int(v.q) for v in kvals])
        l = mp.matrix([mp.mpf(int(v.p)) / int(v.q) for v in lvals])

        def S(v):
            return (
                mp.eye(4) - 1j * sum((gam[j] * v[j] for j in range(dim)), mp.zeros(4))
            ) / (1 + sum(a * a for a in v))

        Sk, Sl = S(k), S(l)
        xx = sum(v * v for v in kvals)
        yy = sum(v * v for v in lvals)
        zz = sum(a * b for a, b in zip(kvals, lvals))
        funcs = {
            "scalar_self": lambda p: 2 * trace(Sk * Sl * Sk * S(k + p)),
            "scalar_vertex": lambda p: trace(Sk * S(k + p) * S(l + p) * Sl),
            "gauge_self": lambda p: (
                2 * sum(trace(Sk * g * Sl * g * Sk * S(k + p)) for g in gam[:dim])
            ),
            "gauge_vertex": lambda p: sum(
                trace(g * Sk * S(k + p) * g * S(l + p) * Sl) for g in gam[:dim]
            ),
        }
        result = 0
        for j in range(dim):
            direction = mp.zeros(dim, 1)
            direction[j] = 1
            result += mp.diff(lambda v, axis=direction: funcs[kind](v * axis), 0, 2) / (
                2 * dim
            )
        expr = getattr(t, kind).subs(
            {t.x: xx, t.y: yy, t.z: zz, t.A: 1 + xx, t.B: 1 + yy, t.d: dim}
        )
        exact = mp.mpf(str(sp.N(expr, 45)))
        assert abs(result - exact) < mp.mpf("1e-38")


@pytest.mark.parametrize(
    "a,b,c", ((1, 1, 1), (1, 2, 1), (2, 2, 1), (3, 1, 1), (2, 3, 0), (3, 3, -1))
)
def test_master_Gamma_formula_against_exact_shift_reduction(a, b, c):
    with mp.workdps(45):
        e = mp.mpc(".23", ".11")
        d = 4 - 2 * e
        direct = (
            mp.gamma(d / 2 - c)
            * mp.gamma(a + b + c - d)
            * mp.gamma(a + c - d / 2)
            * mp.gamma(b + c - d / 2)
            / (
                mp.gamma(a)
                * mp.gamma(b)
                * mp.gamma(d / 2)
                * mp.gamma(a + b + 2 * c - d)
            )
        )
        expr = masters.master(a, b, c).subs(
            t.e, sp.Rational(23, 100) + sp.I * sp.Rational(11, 100)
        )
        re, im = sp.N(expr, 45).as_real_imag()
        reduced = mp.mpc(str(re), str(im)) * mp.gamma(e) ** 2
        assert abs(direct - reduced) < mp.mpf("1e-36")


@pytest.mark.parametrize(
    "a,b,c,e", ((1, 1, 1, ".75"), (1, 2, 1, ".6"), (2, 2, 1, ".6"), (3, 1, 1, ".6"))
)
def test_independent_convergent_Schwinger_integrals(a, b, c, e):
    with mp.workdps(55):
        e = mp.mpf(e)
        d = 4 - 2 * e
        first = mp.quad(lambda v: v ** (c - 1) / (1 + v) ** (d / 2), [0, 1, mp.inf])
        second = mp.quad(
            lambda u: u ** (a + b + c - d - 1) * mp.exp(-u), [0, 1, mp.inf]
        )
        third = mp.quad(
            lambda z: z ** (a + c - d / 2 - 1) * (1 - z) ** (b + c - d / 2 - 1),
            [0, 0.5, 1],
        )
        numerical = first * second * third / (mp.gamma(a) * mp.gamma(b) * mp.gamma(c))
        exact = (
            mp.gamma(d / 2 - c)
            * mp.gamma(a + b + c - d)
            * mp.gamma(a + c - d / 2)
            * mp.gamma(b + c - d / 2)
            / (
                mp.gamma(a)
                * mp.gamma(b)
                * mp.gamma(d / 2)
                * mp.gamma(a + b + 2 * c - d)
            )
        )
        assert abs(numerical / exact - 1) < mp.mpf("1e-12")


@pytest.mark.parametrize("a,b", ((2, 2), (2, 3), (3, 1)))
def test_zero_chord_power_reduces_to_independent_tadpoles(a, b):
    e = t.e
    expected = (
        masters.shifted_gamma_ratio(e, a - 2)
        * masters.shifted_gamma_ratio(e, b - 2)
        / (sp.factorial(a - 1) * sp.factorial(b - 1))
    )
    assert sp.factor(masters.master(a, b, 0) - expected) == 0


@pytest.mark.parametrize(
    "sector,pole2,pole1,finite",
    (
        ("scalar", 3, sp.Rational(-5, 2), sp.Rational(49, 12)),
        ("gauge", -6, 5, sp.Rational(-37, 6)),
    ),
)
def test_independent_complex_MS_Laurent_extraction(sector, pole2, pole1, finite):
    with mp.workdps(65):
        values = []
        for j in range(32):
            e = mp.mpf(".0001") * mp.exp(2j * mp.pi * j / 32)
            A = mp.exp(mp.euler * e) * mp.gamma(1 + e)
            if sector == "scalar":
                R = (
                    2
                    * (2 * e - 3)
                    * (6 * e**3 - 14 * e**2 - 4 * e + 3)
                    / (3 * (e - 2) * (e - 1) * (2 * e + 1))
                )
                ct = 6 * (1 - e) * (1 - 2 * e / 3)
            else:
                R = (
                    4
                    * (2 * e - 3)
                    * (6 * e**3 - 14 * e**2 - e + 3)
                    / (3 * (e - 2) * (2 * e + 1))
                )
                ct = -12 * (1 - e) * (1 - 2 * e / 3)
            values.append(
                (A * A * R + A * ct - pole2) / e**2
                - (
                    mp.mpf(str(pole1.p)) / int(pole1.q)
                    if isinstance(pole1, sp.Rational)
                    else pole1
                )
                / e
            )
        target = mp.mpf(str(finite.p)) / int(finite.q)
        assert abs(sum(values) / 32 - target) < mp.mpf("1e-48")


def test_independent_fermion_field_and_mass_CT_variation():
    h, z, eta, nu, e, m, mu = sp.symbols("h z eta nu epsilon m mu", positive=True)
    original = (mu / m) ** (2 * e)
    changed = (
        (1 + h * nu) ** 2
        / (1 + h * z) ** 2
        * (mu / (m * (1 + h * eta) / (1 + h * z))) ** (2 * e)
    )
    expected = 2 * (nu - z) * original + (eta - z) * m * sp.diff(original, m)
    assert sp.simplify(sp.diff(changed, h).subs(h, 0) - expected) == 0
    assert sp.simplify(m * sp.diff(original, m) + 2 * e * original) == 0
    # Resetting mu=m before differentiating would incorrectly give zero.
    assert sp.diff(original.subs(mu, m), m) == 0


@pytest.mark.parametrize("a,b", remainder.data()["exponents"])
def test_fractional_sunset_constants_by_independent_integrals(a, b):
    with mp.workdps(55):
        a, b = mp.mpf(int(a.p)) / int(a.q), mp.mpf(int(b.p)) / int(b.q)
        g = mp.mpf("1.75")
        first = mp.quad(lambda v: v ** (g - 1) / (1 + v) ** 2, [0, 1, mp.inf])
        second = mp.quad(lambda v: v ** (a + b + g - 5) * mp.exp(-v), [0, 1, mp.inf])
        third = mp.quad(
            lambda u: u ** (a + g - 3) * (1 - u) ** (b + g - 3), [0, 0.5, 1]
        )
        result = first * second * third / (mp.gamma(a) * mp.gamma(b) * mp.gamma(g))
        exact = (
            mp.gamma(2 - g)
            * mp.gamma(a + b + g - 4)
            * mp.gamma(a + g - 2)
            * mp.gamma(b + g - 2)
            / (mp.gamma(a) * mp.gamma(b) * mp.gamma(a + b + 2 * g - 4))
        )
        assert abs(result / exact - 1) < mp.mpf("1e-12")
        assert 0 < exact < 80000


@pytest.mark.parametrize("x", ("1e-50", ".001", ".3", "1", "100", "1e50"))
def test_full_chord_difference_and_fractional_majorant(x):
    with mp.workdps(100):
        x = mp.mpf(x)
        difference = 1 / (x + 1) - 1 / x
        assert abs(difference + 1 / (x * (x + 1))) < mp.mpf("1e-45")
        assert abs(difference) <= x ** (-mp.mpf("1.75"))


@pytest.mark.parametrize("kind", ("self", "vertex"))
@pytest.mark.parametrize("case", ("comparable", "first_large", "second_large"))
def test_mass_difference_soft_coefficient_pointwise(kind, case):
    with mp.workdps(45):
        m = sp.Integer(720)
        kv = [sp.Integer(v) for v in (10, 20, 30, 40)]
        lv = [sp.Integer(v) for v in (50, -30, 10, 20)]
        if case == "first_large":
            kv = [v * 1000 for v in kv]
        if case == "second_large":
            lv = [v * 1000 for v in lv]
        xx = sum(v * v for v in kv)
        yy = sum(v * v for v in lv)
        zz = sum(a * b for a, b in zip(kv, lv))
        cc = xx + yy - 2 * zz
        Sq = m * m + xx
        Sl = m * m + yy
        # The tensor modules use momenta divided by m; undo that scaling
        # for the coefficient of the physical p^2 in the four propagators.
        expr = t.scalar_self / 2 if kind == "self" else t.scalar_vertex
        coefficient = (
            expr.subs(
                {
                    t.x: xx / m**2,
                    t.y: yy / m**2,
                    t.z: zz / m**2,
                    t.A: Sq / m**2,
                    t.B: Sl / m**2,
                    t.d: 4,
                }
            )
            / m**6
        )
        actual = abs(coefficient / (cc * (cc + 1)))
        a, b = (3, 1) if kind == "self" else (2, 2)
        bound = (
            64
            * 360**2
            / (
                Sq ** sp.Rational(a, 2)
                * Sl ** sp.Rational(b, 2)
                * min(Sq, Sl)
                * cc
                * (cc + 1)
            )
        )
        assert 0 < actual < bound


@pytest.mark.parametrize("location", (4, 8, 30, 4 + 2j))
def test_zero_soft_to_on_shell_derivative_Cauchy_bound(location):
    with mp.workdps(40):
        loc = mp.mpc(location)
        T = lambda s: s * s / (loc - s)
        E = 9 / (abs(loc - 1) - 2)
        derivative = lambda s: (2 * loc * s - s * s) / (loc - s) ** 2
        assert derivative(0) == 0
        assert abs(mp.diff(T, 0)) < mp.mpf("1e-60")
        assert abs(mp.diff(T, 1) - derivative(1)) < mp.mpf("1e-35")
        assert abs(derivative(1)) <= E / 2


@pytest.mark.parametrize("i", range(9))
def test_frontier_row_deletion_rejected(i):
    rows = audit.frontier()
    del rows[i]
    with pytest.raises(ValueError):
        audit.validate_frontier(rows)


@pytest.mark.parametrize("field", ("finite_MS_mass", "canonical_pole", "full_matching"))
def test_scope_promotion_rejected(field):
    rows = copy.deepcopy(audit.frontier())
    next(r for r in rows if r["id"] == audit.TARGETS[0])[field] = "COMPLETE"
    with pytest.raises(ValueError):
        audit.validate_frontier(rows)


def test_no_float_reconstruction_and_actual_scope():
    for value in (
        reference.data()["scalar_raw_rational_Gamma_factor"],
        reference.data()["gauge_raw_rational_Gamma_factor_with_vertex_sign"],
    ):
        assert not value.atoms(sp.Float)
    assert "finite mass" in calibration.data()["enclosure"]["scope"]
    assert (
        calibration.enclosure(720, 0, 1, 144, 0)[
            "primitive_on_shell_MS_slope_absolute_upper"
        ]
        == 0
    )
    assert sum(r["status"] == "UNEVALUATED" for r in audit.frontier()) == 2
