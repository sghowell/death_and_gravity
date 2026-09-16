"""Independent full-D finite limit, physical contour and compact-domain tests."""

import itertools
from functools import cache

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_minimal_gravity_finite import (
    assembly,
    audit,
    bounds,
    masters,
    source,
)

ROWS, GATES = audit.residuals(), audit.gates()
a, b, e = s.symbols("a b e")
row = source.dimensional.master_coefficients(a, b, 1, e)
functions = {
    name: s.lambdify((a, b, e), value, "mpmath") for name, value in row.items()
}
T = s.lambdify((a, b, e), source.poles.newton_shape(a, b, 1, 4 + 2 * e), "mpmath")
R = s.lambdify((a, b, e), source.poles.completion(a, b, 1, 4 + 2 * e), "mpmath")
master_rows = [masters.symbols("m" + str(i)) for i in range(3)]
ell_symbol = s.Symbol("ell")
formal = assembly.known_finite(a, b, 1, ell_symbol, master_rows)
arguments = (
    a,
    b,
    ell_symbol,
    *(row[key] for row in master_rows for key in ("J", "Q", "L", "C", "la")),
)
evaluate = s.lambdify(arguments, formal, "mpmath")
CASES = (("6.25", "-1"), ("9", "-2"), ("16", "-.001"))


def lower_log(z):
    z = mp.mpc(z)
    return mp.log(-z.real) - mp.j * mp.pi if z.imag == 0 and z.real < 0 else mp.log(z)


def numeric_triangle(channel):
    def integrand(x):
        z = abs(channel) * x * (1 - x)
        if channel < 0:
            return -(mp.log(z) + mp.pi / mp.sqrt(z)) / (2 * (1 + z))
        h = z - 1
        real = (
            -(1 - h / 2 + h * h / 3)
            if abs(h) < mp.mpf("1e-25")
            else mp.log(z) / (1 - z)
        )
        return -(real + mp.j * mp.pi / (mp.sqrt(z) * (1 + mp.sqrt(z)))) / 2

    return mp.quad(integrand, [0, mp.mpf(".25"), mp.mpf(".5"), mp.mpf(".75"), 1])


def numeric_masters(channel):
    ell = mp.log(4 * mp.pi) - mp.euler
    la = ell - lower_log(-channel)
    if channel < 0:
        J = mp.quad(lambda x: 1 / (1 - channel * x * (1 - x)), [0, mp.mpf(".5"), 1])
        Q = mp.quad(
            lambda x: mp.log(1 - channel * x * (1 - x)) / (1 - channel * x * (1 - x)),
            [0, mp.mpf(".5"), 1],
        )
        L = mp.quad(lambda x: mp.log(1 - channel * x * (1 - x)), [0, mp.mpf(".5"), 1])
    else:

        def state(r):
            x = r + mp.j * r * (1 - r) * (1 - 2 * r) / 4
            xp = 1 + mp.j * (1 - 6 * r + 6 * r * r) / 4
            return xp, 1 - channel * x * (1 - x)

        chart = [0, mp.mpf(".25"), mp.mpf(".5"), mp.mpf(".75"), 1]
        J = mp.quad(lambda r: state(r)[0] / state(r)[1], chart)
        Q = mp.quad(lambda r: state(r)[0] * lower_log(state(r)[1]) / state(r)[1], chart)
        L = mp.quad(lambda r: state(r)[0] * lower_log(state(r)[1]), chart)
    C = numeric_triangle(channel)
    return {
        "J": J,
        "Q": Q,
        "L": L,
        "C": C,
        "la": la,
        "cm": (ell * J - Q) / 2,
        "b00": 2 + la,
        "bmm": ell - L,
    }


def full_D_direct(vals, ms, eps):
    aa, bb, _ = vals
    ell = mp.log(4 * mp.pi) - mp.euler
    vs = [x * x - 4 * x + 2 for x in vals]
    vd = [v + 2 * eps / (1 + eps) for v in vs]
    value = mp.mpc(0)
    for i, x in enumerate(vals):
        for j, _ in enumerate(vals):
            if i != j:
                value += vd[j] ** 2 * ms[j]["J"] / x * (-1 / eps + ms[i]["la"])
    for i, j in ((0, 1), (1, 0), (2, 0)):
        x, y = vals[i], vals[j]
        z = 4 - x - y
        m = ms[i]
        value += functions["C00mu"](x, y, eps) * m["C"]
        value += functions["C0mumu"](x, y, eps) * (-m["J"] / (2 * eps) + m["cm"])
        value += functions["B00"](x, y, eps) * (-1 / eps + m["b00"])
        value += functions["Bmm"](x, y, eps) * (-1 / eps + m["bmm"])
        Bs = -1 / eps + ell + 2
        q = x - 4
        value -= -(12 * Bs - 8) * y * z / q**2 + (5 * Bs - 4) * y * z / q
    ratio = (3 + 2 * eps) / ((1 + eps) * (1 + 2 * eps))
    value += (-1 / eps + ell) * (2 * ratio * T(aa, bb, eps) - R(aa, bb, eps))
    uv = mp.mpf(203) / 40 * sum(x * x for x in vals) - mp.mpf(169) / 3
    soft = sum(vs[i] * ms[i]["J"] for i in range(3)) / 2 - 1
    return value + uv / eps + 2 * T(aa, bb, eps) * soft / eps


def evaluate_rows(vals, ms):
    ell = mp.log(4 * mp.pi) - mp.euler
    return evaluate(
        vals[0],
        vals[1],
        ell,
        *(m[key] for m in ms for key in ("J", "Q", "L", "C", "la")),
    )


@cache
def calibrated(aa, bb):
    with mp.workdps(80):
        vals = (mp.mpf(aa), mp.mpf(bb))
        vals = (*vals, 4 - sum(vals))
        ms = [numeric_masters(x) for x in vals]
        known = evaluate_rows(vals, ms)
        errors = [
            abs(full_D_direct(vals, ms, eps) - known)
            for eps in (mp.mpf("1e-8"), mp.mpf("1e-16"), mp.mpf("1e-24"))
        ]
        return vals, ms, known, errors


@pytest.mark.parametrize("name", tuple(ROWS))
def test_all_exact_residuals(name):
    value = ROWS[name]
    assert all(
        v == 0 for v in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    assert not s.sympify(value).has(s.Float)


@pytest.mark.parametrize("name", tuple(GATES))
def test_all_written_proof_gates(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_all_original_scope_rejections(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("aa,bb", CASES)
def test_full_D_independent_finite_Laurent_limit(aa, bb):
    _, _, _, errors = calibrated(aa, bb)
    assert errors[2] < errors[1] < errors[0]
    assert errors[2] < mp.mpf("1e-12")


@pytest.mark.parametrize("aa,bb", CASES)
def test_whole_physical_amplitude_all_six_crossings(aa, bb):
    vals, ms, known, _ = calibrated(aa, bb)
    with mp.workdps(80):
        for perm in itertools.permutations(range(3)):
            got = evaluate_rows(tuple(vals[i] for i in perm), [ms[i] for i in perm])
            assert abs(got - known) < mp.mpf("1e-60")


@pytest.mark.parametrize("aa,bb", CASES)
def test_physical_master_normalization_and_known_compact_bound(aa, bb):
    vals, ms, known, _ = calibrated(aa, bb)
    with mp.workdps(80):
        beta = mp.sqrt(1 - 4 / vals[0])
        exactJ = (-4 * mp.atanh(beta) + 2 * mp.j * mp.pi) / (vals[0] * beta)
        assert abs(ms[0]["J"] - exactJ) < mp.mpf("1e-50")
        # Independent massless cut dictionary from the crossed spacelike kernel.
        q = vals[0] - 4
        M = mp.atanh(mp.sqrt(q / (4 + q))) / mp.sqrt(q * (4 + q))
        assert abs(ms[0]["C"].imag + 2 * mp.pi * M) < mp.mpf("1e-35")
        delta = min(-vals[1], -vals[2], mp.mpf(1))
        assert abs(known) < mp.mpf(10) ** 11 * (1 + abs(mp.log(delta))) / delta**2


@pytest.mark.parametrize("channel", ("-12", "-1", "-1e-204", "6.25", "9", "16"))
def test_all_master_majorants_with_extreme_spacelike_window(channel):
    with mp.workdps(80):
        value = mp.mpf(channel)
        row = numeric_masters(value)
        delta = min(abs(value), mp.mpf(1))
        L = abs(mp.log(delta))
        assert abs(row["J"]) <= 200 and abs(row["Q"]) <= 1800 and abs(row["L"]) <= 18
        assert abs(row["cm"]) <= 1200 and abs(row["bmm"]) <= 21
        assert abs(row["la"]) <= 10 + L and abs(row["b00"]) <= 12 + L
        assert abs(row["C"]) < 10 * (1 + L) / mp.sqrt(delta)
        if value > 0:
            assert abs(row["C"]) < 5


def test_independent_threshold_triangle_value_and_sheet():
    with mp.workdps(80):
        got = numeric_triangle(mp.mpf(4))
        assert abs(got - (mp.log(2) - mp.j * mp.pi / 2)) < mp.mpf("1e-35")


@pytest.mark.parametrize("energy", (s.Rational(25, 4), s.Integer(9), s.Integer(16)))
def test_exact_contour_gap_and_lower_sheet_calibration(energy):
    for i in range(65):
        r = s.Rational(i, 64)
        x = masters.contour(r)
        den = s.expand_complex(1 - energy * x * (1 - x))
        norm = s.expand_complex(den * s.conjugate(den))
        assert norm >= s.Rational(1, 10000) and norm < 64
        assert s.im(den) <= 0


@pytest.mark.parametrize("mass", (s.Rational(1, 4), s.Integer(3), s.Integer(10)))
def test_full_symbolic_mass_homogeneity(mass):
    aa, bb = s.Rational(25, 4), -s.Integer(1)
    scaled = [
        {
            key: value / mass if key in ("J", "Q", "C") else value
            for key, value in row.items()
        }
        for row in master_rows
    ]
    base = assembly.known_finite(aa, bb, 1, ell_symbol, master_rows)
    got = assembly.known_finite(mass * aa, mass * bb, mass, ell_symbol, scaled)
    assert s.factor(got - mass**2 * base) == 0


@pytest.mark.parametrize(
    "value",
    (
        True,
        False,
        0,
        -1,
        2,
        1.0,
        s.Float(1),
        None,
        "1",
        s.oo,
        -s.oo,
        s.I,
        s.nan,
        s.Symbol("delta"),
    ),
)
def test_invalid_window_rejected(value):
    with pytest.raises((TypeError, ValueError)):
        bounds.require_delta(value)


@pytest.mark.parametrize(
    "energy,transfer,delta",
    (
        (6, -1, 1),
        (17, -1, 1),
        (9, 0, 1),
        (9, -6, 1),
        (9, -s.Rational(1, 2), 1),
        (9.0, -1, 1),
        (9, "-1", 1),
        (s.Symbol("s"), -1, 1),
    ),
)
def test_invalid_physical_domain_rejected(energy, transfer, delta):
    with pytest.raises((TypeError, ValueError)):
        bounds.require_physical(energy, transfer, delta)


@pytest.mark.parametrize(
    "energy,transfer,delta",
    (
        (s.Rational(25, 4), -1, 1),
        (16, -6, 1),
        (16, -s.Rational(1, 10**204), s.Rational(1, 10**204)),
        (9, -2, 1),
    ),
)
def test_valid_closed_window_domain(energy, transfer, delta):
    a, b, c, d = bounds.require_physical(energy, transfer, delta)
    assert a + b + c == 4 and d == delta


def test_original_parameter_known_bound_and_symbolic_anchor_boundary():
    assert sum(bounds.exact_budget().values()) == 10316388000
    upper = bounds.known_amplitude_bound(bounds.COMPARISON_WINDOW, source.KAPPA)
    # log(10)<3 is used analytically in the written proof; this is a calibration.
    assert s.N(upper / s.Integer(10) ** (-1180), 40) < 1
    data = assembly.data()
    matched = data["whole_matched_finite_representative"]
    assert matched.has(assembly.ALPHA, assembly.BETA, assembly.DK)
    assert "symbolic" in audit.observable()["established"]


def test_whole_original_frontier_and_historical_qualification_unchanged():
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.matching()) == 158 and len(audit.frontier()) == 9
    assert audit.qualifications() == audit.previous.qualifications()
    assert len(audit.qualifications()) == 6
    assert any(r["status"].startswith("REJECTED_") for r in audit.matching())
