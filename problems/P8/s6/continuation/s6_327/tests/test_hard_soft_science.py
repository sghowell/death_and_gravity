"""Exact helicity, integral, original-parameter and subtraction-domain tests."""

from itertools import product

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_known_hard_soft_subtraction import (
    audit,
    integrals,
    interference,
    source,
)

ROWS, GATES = audit.residuals(), audit.gates()


@pytest.mark.parametrize("name", tuple(ROWS))
def test_every_exact_residual(name):
    value = ROWS[name]
    assert all(
        v == 0 for v in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    assert not s.sympify(value).has(s.Float)


@pytest.mark.parametrize("name", tuple(GATES))
def test_every_analytic_proof_gate(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_rejected_inputs_and_scope(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize(
    "bad",
    (
        None,
        True,
        False,
        "1",
        1.0,
        s.Float(1),
        [],
        {},
        s.Matrix([1]),
        s.I,
        s.oo,
        s.nan,
        s.Symbol("x"),
    ),
)
def test_exact_resolution_api_rejects_inexact_or_unresolved_values(bad):
    for call in (
        integrals.interference_kernel_upper,
        integrals.known_hard_soft_upper,
        integrals.logarithmic_addition_upper,
    ):
        with pytest.raises((TypeError, ValueError)):
            call(bad)


@pytest.mark.parametrize("bad", (0, -1, s.Rational(1, 7), 1))
def test_resolution_domain(bad):
    for call in (
        integrals.interference_kernel_upper,
        integrals.known_hard_soft_upper,
        integrals.logarithmic_addition_upper,
    ):
        with pytest.raises(ValueError):
            call(bad)


@pytest.mark.parametrize(
    "x", (s.Rational(1, 10**30), s.Rational(1, 128), s.Rational(1, 8))
)
def test_original_parameter_bounds(x):
    kap, L = s.Integer(10) ** 800, -s.log(x)
    assert (
        integrals.interference_kernel_upper(x)
        == s.Min(s.Integer(10) ** 18, 2 * s.Integer(10) ** 20 * x) / kap
    )
    assert (
        integrals.known_hard_soft_upper(x)
        == s.Min(s.Integer(10) ** 30, 2 * s.Integer(10) ** 32 * x) / kap**2
    )
    assert (
        integrals.logarithmic_addition_upper(x)
        == s.Integer(10) ** 18 * x * (1 + L) / kap**2
        + 3000 * x * x * (L * L + L + s.Rational(1, 2)) / kap**3
    )
    assert source.HEAVY_MASS2 == s.Integer(10) ** 200 / 512 + 2


@pytest.mark.parametrize(
    "angle,phase",
    tuple(
        product((0, s.Rational(1, 3), s.Rational(3, 4)), (1, s.I, (3 + 4 * s.I) / 5))
    ),
)
def test_unpolarized_phase_identity_for_independent_unitary_frames(angle, phase):
    angle = s.sympify(angle)
    c, d = (1 - angle**2) / (1 + angle**2), 2 * angle / (1 + angle**2)
    U = s.diag(phase, s.conjugate(phase)) * s.Matrix([[c, d], [-d, c]])
    F, S = s.Matrix([2, -3]), s.Matrix([5, 7])
    z = ((U * F).conjugate().T * (U * S))[0]
    assert not z.has(s.Float)
    assert s.expand_complex(z) == -11
    assert s.simplify(2 * s.re(z * (3 + 11 * s.I))) == -66


@pytest.mark.parametrize(
    "row",
    interference.original_calibration()["records"],
    ids=("original_0", "original_2"),
)
def test_complete_original47_counterexample_to_per_helicity_phase_deletion(row):
    _, fp, fc, sp, sc, single = row
    assert all(s.im(v) == 0 for v in (fp, fc, sp, sc))
    assert single != 0
    assert s.factor(single - (fp * sc - fc * sp) / 2) == 0
    fh, sh = (fp + s.I * fc) / s.sqrt(2), (sp + s.I * sc) / s.sqrt(2)
    assert s.factor(s.re(s.conjugate(fh) * sh * s.I) + single) == 0


@pytest.mark.parametrize("power", (0, 1))
def test_log_weighted_primitives_and_zero_endpoints(power):
    w, x = s.symbols("w x", positive=True)
    result = x ** (power + 1) * (s.log(1 / x) + s.Rational(1, power + 1)) / (power + 1)
    assert s.expand(s.diff(result, x) - x**power * s.log(1 / x)) == 0
    assert s.limit(result, x, 0, dir="+") == 0
    assert s.expand(s.integrate(w**power * s.log(1 / w), (w, 0, x)) - result) == 0


@pytest.mark.parametrize(
    "y,x",
    tuple(
        product(
            (s.Rational(1, 10**30), s.Rational(1, 64), s.Rational(1, 4), s.S.One),
            (s.Rational(1, 10**20), s.Rational(1, 128), s.Rational(1, 8)),
        )
    ),
)
def test_independent_low_high_integrated_kernel_majorant(y, x):
    with mp.workdps(100):
        yy, xx = mp.mpf(str(s.N(y, 110))), mp.mpf(str(s.N(x, 110)))
        aa, uu = yy / 192, min(xx, yy / 192)
        value = (320 * 10**8 + 102400 * yy**2) * uu
        if xx > aa:
            logratio = mp.log(xx / aa)
            third = yy * (
                logratio - mp.log((1 + (xx / yy) ** 2) / (1 + mp.mpf(1) / 192**2)) / 2
            )
            value += (
                4 * 160**2 * yy**2 * logratio
                + 320 * 330000 * yy * (xx - aa)
                + 320 * 10**16 * third
            )
        value /= 4 * mp.pi**2
        assert 0 < value < min(mp.mpf(10) ** 18, 2 * mp.mpf(10) ** 20 * xx)


@pytest.mark.parametrize(
    "y,x",
    tuple(
        product(
            (s.Rational(1, 10**30), s.Rational(1, 64), s.Rational(1, 4), s.S.One),
            (s.Rational(1, 10**20), s.Rational(1, 128), s.Rational(1, 8)),
        )
    ),
)
def test_independent_log_weighted_triangle_majorant(y, x):
    with mp.workdps(100):
        yy, xx = mp.mpf(str(s.N(y, 110))), mp.mpf(str(s.N(x, 110)))
        uu = min(xx, yy / 192)
        primitive0 = lambda z: z * (1 - mp.log(z))
        primitive1 = lambda z: z * z * (-mp.log(z) + mp.mpf(1) / 2) / 2
        upper = 64 * primitive0(xx) + mp.mpf(10) ** 8 / yy * primitive1(uu)
        if xx > uu:
            upper += 330000 * (primitive1(xx) - primitive1(uu)) + mp.mpf(10) ** 16 * (
                primitive0(xx) - primitive0(uu)
            )
        assert 0 < upper < 2 * mp.mpf(10) ** 16 * primitive0(xx)


def test_real_part_bound_cannot_bound_full_hard_modulus_or_square():
    h = s.Rational(1, 10**100) + s.I * s.Integer(10) ** 100
    assert s.re(h) < 1 and s.expand(h * s.conjugate(h)) > s.Integer(10) ** 200
    assert "modulus or square" in audit.observable()["not_established"]


def test_scope_is_not_full_radiative_loop_or_positive_measure():
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 183
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.qualifications()) == 6
    assert (
        "Finite radiative hard-loop remainder" in audit.observable()["not_established"]
    )
    assert "positive normalized" in audit.observable()["not_established"]
