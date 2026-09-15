"""Independent complete quartic graphs, full-D cuts and finite error control."""

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_whole_quartic_gravity_sector import (
    audit,
    cuts,
    forward,
    proper,
    source,
)

ROWS, GATES = audit.residuals(), audit.gates()


@pytest.mark.parametrize("name", tuple(ROWS))
def test_every_exact_identity(name):
    value = ROWS[name]
    assert all(
        x == 0 for x in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    assert not value.atoms(s.Float)


@pytest.mark.parametrize("name", tuple(GATES))
def test_every_written_proof_gate(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[row[0] for row in audit.bad_cases()]
)
def test_every_unsupported_scope_input(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize(
    "mm,scale2,res2", ((1, 1, 1), (2, 3, "0.5"), ("1.5", 2, "1.5"))
)
def test_independent_entire_D_forward_moments_and_finite_coefficient(mm, scale2, res2):
    with mp.workdps(55):
        mm, scale2, res2 = map(mp.mpf, (mm, scale2, res2))

        def T(aa, ee):
            return mp.quad(
                lambda xx: (
                    (((aa - 2 * mm) ** 2 - 2 * mm * mm / (1 + ee)) / 2)
                    * (mm - aa * xx * (1 - xx)) ** (ee - 1)
                    + (aa - 2 * mm) * (mm - aa * xx * (1 - xx)) ** ee
                ),
                [0, 1],
            )

        def E(aa, ee):
            return mp.quad(
                lambda xx: (
                    (aa + 2 * mm / (1 + ee))
                    * xx
                    * (1 - xx)
                    * (mm - aa * xx * (1 - xx)) ** ee
                ),
                [0, 1],
            )

        def Tsecond(ee):
            return mp.quad(
                lambda xx: (
                    (1 - 2 * ee * xx * (1 - xx))
                    * (mm * (1 - 2 * xx * (1 - xx))) ** (ee - 1)
                    - mm
                    * mm
                    * (ee - 1)
                    * (ee - 2)
                    / (1 + ee)
                    * (xx * (1 - xx)) ** 2
                    * (mm * (1 - 2 * xx * (1 - xx))) ** (ee - 3)
                ),
                [0, 1],
            )

        def Esecond(ee):
            return mp.quad(
                lambda xx: (
                    -2
                    * ee
                    * (xx * (1 - xx)) ** 2
                    * (mm * (1 - 2 * xx * (1 - xx))) ** (ee - 1)
                    + 2
                    * mm
                    * (2 + ee)
                    / (1 + ee)
                    * ee
                    * (ee - 1)
                    * (xx * (1 - xx)) ** 3
                    * (mm * (1 - 2 * xx * (1 - xx))) ** (ee - 2)
                ),
                [0, 1],
            )

        for ee in (mp.mpf("0.0625"), mp.mpf("0.03125")):
            assert abs(
                mp.diff(lambda aa, ee=ee: T(aa, ee), 2 * mm, 2) - Tsecond(ee)
            ) < mp.mpf("1e-45")
            assert abs(
                mp.diff(lambda aa, ee=ee: E(aa, ee), 2 * mm, 2) - Esecond(ee)
            ) < mp.mpf("1e-45")
        B2 = 3 * mp.pi / (8 * mm)
        numericT0 = Tsecond(0)
        numericT1 = mp.diff(Tsecond, 0)
        analyticT1 = (
            3 * mp.pi * mp.log(2 * mm) / 8
            - 3 * mp.catalan / 2
            + mp.mpf("1.25")
            - mp.pi / 4
        ) / mm
        assert abs(numericT0 - B2) < mp.mpf("1e-45")
        assert abs(numericT1 - analyticT1) < mp.mpf("1e-45")

        def full_subtracted(ee):
            raw = (
                mp.gamma(-ee)
                * (4 * mp.pi * scale2) ** (-ee)
                * (-2 * Tsecond(ee) + Esecond(ee))
                / (16 * mp.pi**2)
            )
            soft = B2 * (res2 / scale2) ** ee / (8 * mp.pi**2 * ee)
            return raw - soft

        expected = (
            3 * mp.pi * (mp.euler + mp.log(mm / (2 * mp.pi * res2)) - 1) / 8
            - 3 * mp.catalan / 2
            + mp.mpf("1.75")
        ) / (8 * mp.pi**2 * mm)
        errors = [
            abs(full_subtracted(ee) - expected)
            for ee in (mp.mpf("1e-5"), mp.mpf("1e-6"))
        ]
        assert errors[1] < mp.mpf("0.11") * errors[0], errors
        assert errors[0] < mp.mpf("1e-4")
        library = forward.finite_b20(
            s.Rational(str(mm)), 1, 1, s.Rational(str(res2))
        ).evalf(45)
        assert abs(mp.mpf(str(library)) - expected) < mp.mpf("1e-40")


@pytest.mark.parametrize(
    "m,energy,scale,epsilon",
    (
        (1, 9, 1, "0.0625"),
        (1, 9, 1, "0.125"),
        (2, 17, 3, "0.0625"),
        (2, 17, 3, "0.125"),
    ),
)
def test_independent_whole_tree_moment_phase_cut(m, energy, scale, epsilon):
    with mp.workdps(60):
        epsilon = mp.mpf(epsilon)
        m, energy, scale = map(mp.mpf, (m, energy, scale))
        beta = mp.sqrt(1 - 4 * m / energy)
        q = energy - 4 * m

        def moment(alpha):
            return (
                (energy / 4) ** alpha
                * beta ** (2 * alpha + 1)
                / 2
                * (
                    mp.betainc(alpha + 1, -alpha - mp.mpf("0.5"), 0, 1 - beta * beta)
                    + mp.exp(-mp.j * mp.pi * alpha) * mp.beta(mp.mpf("0.5"), alpha + 1)
                )
            )

        M = moment(epsilon)
        J = moment(epsilon - 1)
        Y = (m * M - moment(epsilon + 1)) / energy
        Ve = (energy - 2 * m) ** 2 - 2 * m * m / (1 + epsilon)
        pref = mp.gamma(-epsilon) * (4 * mp.pi * scale) ** (-epsilon)
        loop = mp.im(
            pref
            * (
                -Ve * J
                - 2 * (energy - 2 * m) * M
                + (energy + 2 * m / (1 + epsilon)) * Y
            )
        ) / (16 * mp.pi**2)
        phase = (
            beta
            / (8 * mp.pi)
            * (4 * mp.pi * scale) ** (-epsilon)
            * energy**epsilon
            * beta ** (2 * epsilon)
            * mp.gamma(1 + epsilon)
            / mp.gamma(2 + 2 * epsilon)
        )
        assert abs(mp.im(pref * M) / (16 * mp.pi**2) - phase / 2) < mp.mpf("1e-50")
        angular_norm = mp.quad(lambda x: (1 - x * x) ** epsilon, [0, 1])
        angular_pole = 1 / (2 * epsilon) + mp.quad(
            lambda x: (1 - x * x) ** epsilon / (1 + x), [0, 1]
        )
        angular_quad = mp.quad(lambda x: x * x * (1 - x * x) ** epsilon, [0, 1])
        p = 4 * Ve / q
        aa = -7 * energy / 4 + 4 * m + 2 * m * m / (energy * (1 + epsilon))
        bb = -q * q / (4 * energy)
        sewn = (
            phase
            / 2
            * (p * angular_pole + aa * angular_norm + bb * angular_quad)
            / angular_norm
        )
        err = abs(loop - sewn)
        assert err < mp.mpf("1e-48"), err


@pytest.mark.parametrize("dimension", (4, 5, 6))
@pytest.mark.parametrize("energy", (4, 5, 6))
def test_independent_whole_two_endpoint_projector_components(dimension, energy):
    eta = s.diag(1, *([-1] * (dimension - 1)))
    q = s.Matrix([0, 2, *([0] * (dimension - 2))])
    P = s.Matrix([energy, 0, 1, *([0] * (dimension - 3))])
    dot = lambda v, w: (v.T * eta * w)[0]
    a = dot(q, q)
    mass = dot(P, P) + a / 4
    Q = q * q.T - eta * a
    tree = 2 * P * P.T - Q / 2
    contraction = sum(
        eta[i, i] * eta[j, j] * Q[i, j] * tree[i, j]
        for i in range(dimension)
        for j in range(dimension)
    )
    contraction -= s.trace(eta * Q) * s.trace(eta * tree) / (dimension - 2)
    assert (tree * eta * q) == s.zeros(dimension, 1)
    assert s.factor(2 * contraction / a - a - 4 * mass / (dimension - 2)) == 0
    assert mass > 0


@pytest.mark.parametrize(
    "bad",
    (True, False, 0, -1, s.Rational(1, 8), 2, 1.0, None, "1", s.Symbol("resolution")),
)
def test_reject_outside_original_resolution_domain(bad):
    with pytest.raises((TypeError, ValueError)):
        forward.require_resolution(bad)


@pytest.mark.parametrize(
    "bad", (True, False, 0, -1, s.Rational(1, 4), 1.0, None, "epsilon")
)
def test_reject_outside_regulator_bound_domain(bad):
    with pytest.raises((TypeError, ValueError)):
        forward.regulator_remainder(bad, -1, 1, 1)


@pytest.mark.parametrize(
    "bad", (True, False, 1.0, s.Float(1), None, "quartic", s.oo, s.I)
)
def test_reject_inexact_or_unsupported_bound_quartic(bad):
    with pytest.raises((TypeError, ValueError)):
        forward.finite_magnitude_bound(bad, 1, 1)


@pytest.mark.parametrize(
    "args",
    (
        (1, 4, s.Rational(1, 16)),
        (1, 3, s.Rational(1, 16)),
        (0, 9, s.Rational(1, 16)),
        (1, 9, 0),
        (1, 9, s.Rational(1, 4)),
        (True, 9, s.Rational(1, 16)),
        (1.0, 9, s.Rational(1, 16)),
        (1, 9, "epsilon"),
    ),
)
def test_reject_outside_entire_scalar_cut_domain(args):
    with pytest.raises((TypeError, ValueError)):
        cuts.require_cut_domain(*args)


@pytest.mark.parametrize("resolution", (s.Rational(1, 4), s.Rational(1, 2), s.S.One))
def test_public_bound_and_resolution_flow(resolution):
    assert forward.regulator_remainder(
        s.Rational(1, 16), s.Rational(-1, 2), 1, resolution
    ) == s.Rational(1, 4)
    assert forward.finite_magnitude_bound(
        s.Rational(-1, 2), 1, resolution
    ) == s.Rational(1, 8)
    r = s.Symbol("positive_resolution_squared", positive=True)
    value = forward.finite_b20(1, -1, 1, r)
    expected = forward.soft_b20(1) / (8 * s.pi**2 * r)
    assert s.simplify(s.diff(value, r) - expected) == 0
    assert s.simplify(value.subs(r, resolution)).is_real


@pytest.mark.parametrize("epsilon", (s.Rational(1, 16), s.Rational(1, 8)))
def test_complete_public_center_integrands_against_independent_differentiation(epsilon):
    a = s.Symbol("channel")
    x = proper.X
    y = x * (1 - x)
    A = 1 - a * y
    V = (a - 2) ** 2 - 2 / (1 + epsilon)
    rawT = V * A ** (epsilon - 1) / 2 + (a - 2) * A**epsilon
    rawE = (a + 2 / (1 + epsilon)) * y * A**epsilon
    actualT, actualE = forward.center_integrands(epsilon, 1)
    # The denominator is positive on the entire integration interval.
    for at in (s.S.Zero, s.Rational(1, 3), s.Rational(1, 2), s.S.One):
        assert (
            s.simplify(s.diff(rawT, a, 2).subs({a: 2, x: at}) - actualT.subs(x, at))
            == 0
        )
        assert (
            s.simplify(s.diff(rawE, a, 2).subs({a: 2, x: at}) - actualE.subs(x, at))
            == 0
        )


def test_whole_zero_moment_and_weighted_boundary_not_arbitrary_constant():
    assert proper.moment(0, proper.S) == 1
    assert proper.weighted_moment(0, proper.S) == s.Rational(1, 6)
    e = s.Symbol("epsilon", positive=True)
    assert proper.moment(e, 0, 2) == 2**e
    assert proper.weighted_moment(e, 0, 2) == 2**e / 6
    assert proper.E_coefficient(0) == 5 * source.MU / 3


def test_original_source_bounds_and_historical_scope():
    assert source.data()["checks"] is not source.previous.data()["checks"]
    assert len(source.previous.data()["checks"]) == 109
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 149
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.qualifications()) == 6
    assert any(row["status"].startswith("REJECTED_") for row in audit.matching())
    assert "Coulomb" in audit.observable()["not_established"]
    assert "Regge" in audit.observable()["not_established"]
    assert abs(source.CONTACT) < 4 * source.CUBIC**2 / source.HEAVY_MASS2
    assert source.CUBIC**2 / (source.KAPPA * source.HEAVY_MASS2) < s.Rational(
        1, 10**1005
    )


def test_entire_scientific_payload_recursively_exact():
    def exact(value):
        if isinstance(value, float):
            pytest.fail("Python float in scientific payload")
        if isinstance(value, (s.Basic, s.MatrixBase)):
            assert not value.has(s.Float)
        elif isinstance(value, dict):
            for item in value.values():
                exact(item)
        elif isinstance(value, (tuple, list)):
            for item in value:
                exact(item)

    for packet in audit.packets().values():
        exact(packet)
