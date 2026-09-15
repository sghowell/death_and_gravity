"""Independent full tensors, physical-sheet moments and finite paired functional."""

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_heavy_resonance_soft_pole_pairing import (
    audit,
    masters,
    pairing,
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
def test_every_proof_gate(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[row[0] for row in audit.bad_cases()]
)
def test_every_unsupported_input(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("dimension", (4, 5, 6))
@pytest.mark.parametrize("seed", (1, 2, 3))
def test_independent_complete_pair_and_contact_tensors(dimension, seed):
    eta = s.diag(1, *([-1] * (dimension - 1)))
    p = s.Matrix([seed + 5, 1, 2, *([0] * (dimension - 3))])
    r = s.Matrix([seed + 4, 2, -1, *([0] * (dimension - 3))])
    loop = s.Matrix(
        [
            s.Rational(seed, 3),
            s.Rational(1, 5),
            -s.Rational(2, 7),
            *([0] * (dimension - 3)),
        ]
    )
    dot = lambda v, w: (v.T * eta * w)[0]
    mi, mj, z = dot(p, p), dot(r, r), dot(p, r)
    Ti = p * (p + loop).T + (p + loop) * p.T - eta * (dot(p, p + loop) - mi)
    Tj = r * (r - loop).T + (r - loop) * r.T - eta * (dot(r, r - loop) - mj)
    literal = sum(
        eta[i, i] * eta[j, j] * Ti[i, j] * Tj[i, j]
        for i in range(dimension)
        for j in range(dimension)
    )
    literal -= s.trace(eta * Ti) * s.trace(eta * Tj) / (dimension - 2)
    N0 = 4 * z * z - 4 * mi * mj / (dimension - 2)
    k2 = dot(loop, loop)
    di, dj = k2 + 2 * dot(p, loop), k2 - 2 * dot(r, loop)
    assert literal == N0 + 2 * z * k2 - 2 * z * (di + dj)
    assert (
        s.factor(
            literal / (k2 * di * dj)
            - N0 / (k2 * di * dj)
            - 2 * z / (di * dj)
            + 2 * z / (k2 * di)
            + 2 * z / (k2 * dj)
        )
        == 0
    )
    contact = s.trace(eta * Ti) * (1 - s.Rational(dimension, dimension - 2))
    assert contact == -4 * mi / (dimension - 2) + di - k2
    assert literal != N0


@pytest.mark.parametrize("mu,n", ((1, 9), (2, 17), (1, 100)))
@pytest.mark.parametrize("angle", (s.S.Zero, s.Rational(1, 2), s.Rational(3, 4)))
def test_independent_full_three_leg_current_and_TT(mu, n, angle):
    mu, n = map(s.Integer, (mu, n))
    eta = s.diag(1, -1, -1, -1)
    E = s.sqrt(n) / 2
    a = s.sqrt(n - 4 * mu) * s.sqrt(1 - angle * angle) / 2
    b = s.sqrt(n - 4 * mu) * angle / 2
    p = s.Matrix([E, a, 0, b])
    r = s.Matrix([E, -a, 0, -b])
    P = p + r
    q = s.Matrix([1, 0, 0, 1])
    dot = lambda v, w: (v.T * eta * w)[0]
    current = p * p.T / dot(p, q) + r * r.T / dot(r, q) - P * P.T / dot(P, q)
    assert (current * eta * q).applyfunc(s.simplify) == s.zeros(4, 1)
    literal = (
        sum(
            eta[i, i] * eta[j, j] * current[i, j] ** 2
            for i in range(4)
            for j in range(4)
        )
        - s.trace(eta * current) ** 2 / 2
    )
    expected = (
        (n - 4 * mu) ** 2
        / (2 * n)
        * (1 - angle * angle) ** 2
        / (1 - (1 - 4 * mu / n) * angle * angle) ** 2
    )
    assert s.simplify(literal - expected) == 0
    assert literal > 0


@pytest.mark.parametrize("mu,nval", ((1, 9), (2, 17), (1, 100)))
def test_independent_complete_physical_root_moments_and_first_derivative(mu, nval):
    with mp.workdps(65):
        mu, nval = map(mp.mpf, (mu, nval))
        beta = mp.sqrt(1 - 4 * mu / nval)
        c = 1 - beta * beta

        def moment(power):
            return (
                (nval / 4) ** power
                * beta ** (2 * power + 1)
                / 2
                * (
                    mp.betainc(power + 1, -power - mp.mpf("0.5"), 0, c)
                    + mp.exp(-mp.j * mp.pi * power) * mp.beta(mp.mpf("0.5"), power + 1)
                )
            )

        I0 = mp.atanh(beta) / beta
        J0 = -4 * I0 / nval + mp.j * 2 * mp.pi / (nval * beta)
        V0 = 2 * mp.log(2)
        V1 = 2 * mp.log(2) ** 2 - mp.pi**2 / 6
        U1 = 2 * mp.quad(
            lambda z: mp.log(1 - z * z) / (1 + z) - 2 * mp.log(z) / (1 - z * z),
            [beta, 1],
        )
        H0 = -2 * mp.atanh(beta) + mp.j * mp.pi
        H1 = (mp.log(c) ** 2 + mp.pi**2) / 2 + U1 - V1 + mp.j * mp.pi * V0
        J1 = 2 / (nval * beta) * (H1 + mp.log(nval * beta * beta / 4) * H0)
        h = mp.mpf("1e-12")
        numericJ1 = (moment(-1 + h) - moment(-1 - h)) / (2 * h)
        numericJ0 = (moment(-1 + h) + moment(-1 - h)) / 2
        errJ0 = abs(numericJ0 - J0)
        errJ1 = abs(numericJ1 - J1)
        assert errJ0 < mp.mpf("1e-20") and errJ1 < mp.mpf("1e-20"), (
            mu,
            nval,
            errJ0,
            errJ1,
        )
        bubblelog = mp.log(mu) - 2 + 2 * beta * mp.atanh(beta) - mp.j * mp.pi * beta
        numericlog = mp.diff(moment, 0)
        assert abs(moment(0) - 1) < mp.mpf("1e-50")
        assert abs(numericlog - bubblelog) < mp.mpf("1e-50")
        # Independently integrate both physical segments in a nonsingular example.
        for power in (mp.mpf("0.5"), mp.mpf("1.5")):
            literal = (nval / 4) ** power * (
                mp.exp(-mp.j * mp.pi * power)
                * mp.quad(
                    lambda w, power=power: (beta * beta - w * w) ** power, [0, beta]
                )
                + mp.quad(
                    lambda w, power=power: (w * w - beta * beta) ** power, [beta, 1]
                )
            )
            assert abs(literal - moment(power)) < mp.mpf("1e-50")


@pytest.mark.parametrize("mu,n,nu2", ((1, 9, 1), (2, 17, 3)))
def test_independent_whole_D_real_virtual_finite_window_pair(mu, n, nu2):
    with mp.workdps(55):
        mu, n, nu2 = map(mp.mpf, (mu, n, nu2))
        beta = mp.sqrt(1 - 4 * mu / n)
        b = beta * beta
        c = 1 - b
        I0 = mp.atanh(beta) / beta
        F0 = (3 - b + (b - 1) * (b + 3) * I0) / (2 * b * b)
        R3 = (n - 4 * mu) ** 2 * F0 / (2 * n)
        zll = (n - 2 * mu) / 2
        zhl = -n / 2

        def light_moment(a):
            return (
                (n / 4) ** a
                * beta ** (2 * a + 1)
                / 2
                * (
                    mp.betainc(a + 1, -a - mp.mpf("0.5"), 0, c)
                    + mp.exp(-mp.j * mp.pi * a) * mp.beta(mp.mpf("0.5"), a + 1)
                )
            )

        def heavy_moment(a):
            return mp.quad(lambda v: (n * v + mu * (1 - v) ** 2) ** a, [0, 1])

        def B(e):
            aLL = 2 * zll * zll - mu * mu / (1 + e)
            aHL = 2 * zhl * zhl - n * mu / (1 + e)
            return (
                aLL * light_moment(e - 1)
                + 2 * zll * light_moment(e)
                + 2 * (aHL * heavy_moment(e - 1) + 2 * zhl * heavy_moment(e))
                + (2 * mu ** (1 + e) + n ** (1 + e)) / (2 * (1 + e))
            )

        Jll0 = -4 * I0 / n + mp.j * 2 * mp.pi / (n * beta)
        U1 = 2 * mp.quad(
            lambda z: mp.log(1 - z * z) / (1 + z) - 2 * mp.log(z) / (1 - z * z),
            [beta, 1],
        )
        H0 = -2 * mp.atanh(beta) + mp.j * mp.pi
        H1 = (
            (mp.log(c) ** 2 + mp.pi**2) / 2
            + U1
            - (2 * mp.log(2) ** 2 - mp.pi**2 / 6)
            + mp.j * mp.pi * 2 * mp.log(2)
        )
        Jll1 = 2 / (n * beta) * (H1 + mp.log(n * b / 4) * H0)
        Lll = mp.log(mu) - 2 + 2 * beta * mp.atanh(beta) - mp.j * mp.pi * beta
        A = lambda v: n * v + mu * (1 - v) ** 2
        Jhl0 = 2 * I0 / n
        Jhl1 = mp.quad(lambda v: mp.log(A(v)) / A(v), [0, 1])
        Lhl = mp.quad(lambda v: mp.log(A(v)), [0, 1])
        B1 = (
            mu * mu * Jll0
            + (2 * zll * zll - mu * mu) * Jll1
            + 2 * zll * Lll
            + 2 * (n * mu * Jhl0 + (2 * zhl * zhl - n * mu) * Jhl1 + 2 * zhl * Lhl)
            + (2 * mu * (mp.log(mu) - 1) + n * (mp.log(n) - 1)) / 2
        )
        V = n * n - 4 * mu * n + 2 * mu * mu
        B0 = -R3 + mp.j * mp.pi * V / (n * beta)
        h = mp.mpf("1e-11")
        assert abs((B(h) + B(-h)) / 2 - B0) < mp.mpf("1e-16")
        assert abs((B(h) - B(-h)) / (2 * h) - B1) < mp.mpf("1e-16")

        def K(s, e):
            b = 1 - 4 * mu / s
            N = (
                mp.beta(e + 3, mp.mpf("0.5"))
                / 2
                * mp.hyp2f1(2, mp.mpf("0.5"), e + mp.mpf("3.5"), b)
            )
            return (
                (s - 4 * mu) ** 2
                / (4 * mp.pi * s)
                * (1 + 2 * e)
                / (2 + 2 * e)
                * (16 * mp.pi * nu2 * s) ** (-e)
                / mp.gamma(1 + e)
                * N
            )

        W = lambda s: 1 / (s - 2 * mu) ** 3
        L = mp.mpf("0.25")

        def H(y, e):
            return n ** (2 * e) * K(n * (1 + y), e) * W(n * (1 + y))

        def regular(e):
            h0 = H(0, e)
            derivative = mp.diff(lambda y: H(y, e), 0)

            def integrand(y):
                if y < mp.mpf("1e-25"):
                    return derivative * y ** (2 * e)
                return y ** (-1 + 2 * e) * (H(y, e) - h0)

            return mp.quad(integrand, [0, L])

        angularlog = mp.quad(
            lambda x: (1 - x * x) ** 2 * mp.log(1 - x * x) / (1 - b * x * x) ** 2,
            [0, 1],
        )
        K0 = K(n, 0)
        finite = (
            K0 * W(n) * mp.log(L)
            + K0 * W(n) / 2 * (1 + mp.log(n / 4) + angularlog / F0 + mp.re(B1) / R3)
            + regular(0)
        )

        def paired(e):
            realcut = H(0, e) * L ** (2 * e) / (2 * e) + regular(e)
            C = mp.gamma(1 - e) * (4 * mp.pi * nu2) ** (-e) * B(e)
            virtual = mp.re(C) * W(n) / (8 * mp.pi * e)
            return realcut + virtual

        errors = [abs(paired(e) - finite) for e in (mp.mpf("1e-5"), mp.mpf("1e-6"))]
        assert errors[1] < mp.mpf("0.11") * errors[0], errors
        assert errors[0] < mp.mpf("1e-3") * max(abs(finite), mp.mpf("1e-12")), errors


@pytest.mark.parametrize("mu,n", ((1, 9), (2, 17), (1, 100)))
def test_public_master_first_derivative_and_cubic_coefficient(mu, n):
    with mp.workdps(55):
        beta = mp.sqrt(1 - 4 * mp.mpf(mu) / n)

        def M(a):
            return (
                (mp.mpf(n) / 4) ** a
                * beta ** (2 * a + 1)
                / 2
                * (
                    mp.betainc(a + 1, -a - mp.mpf("0.5"), 0, 1 - beta * beta)
                    + mp.exp(-mp.j * mp.pi * a) * mp.beta(mp.mpf("0.5"), a + 1)
                )
            )

        h = mp.mpf("1e-12")
        expected = (M(-1 + h) - M(-1 - h)) / (2 * h)
        value = masters.light_cusp1(mu, n).evalf(45)
        actual = mp.mpc(str(s.re(value)), str(s.im(value)))
        assert abs(actual - expected) < mp.mpf("1e-20")
    assert s.simplify(proper.complete_B(0, mu, n) - proper.B0(mu, n)) == 0
    assert not proper.complete_B(s.Rational(1, 16), mu, n).has(s.Float)


@pytest.mark.parametrize(
    "bad", (True, False, 0, 1, None, "stable_H", "massless", s.Symbol("pair"))
)
def test_unsupported_pair_species_rejected(bad):
    with pytest.raises((TypeError, ValueError)):
        masters.raw_pair_masters(bad)


@pytest.mark.parametrize(
    "args",
    (
        (1, 4, s.Rational(1, 4)),
        (1, 3, s.Rational(1, 4)),
        (1, 9, 1),
        (0, 9, s.Rational(1, 4)),
        (1, 9, 0),
        (1.0, 9, s.Rational(1, 4)),
        (True, 9, s.Rational(1, 4)),
    ),
)
def test_invalid_paired_window_domains_rejected(args):
    with pytest.raises((TypeError, ValueError)):
        pairing.require_window(*args)


@pytest.mark.parametrize("dimension", (s.Integer(4), s.Integer(5), s.Integer(6)))
@pytest.mark.parametrize("mass", (s.S.One, s.Rational(3, 2)))
def test_independent_full_A0_B0_LSZ_collapse_before_raw_pole(dimension, mass):
    e = (dimension - 4) / 2
    B0 = s.Symbol("complete_B0")
    A0 = mass * (1 + 2 * e) * B0 / (1 + e)
    residue = mass * (3 + 2 * e) * B0 / (1 + e)
    pair_onshell = -2 * mass * B0
    contact = 4 * mass * B0 / (dimension - 2) + A0
    assert (
        s.factor(
            pair_onshell
            + contact
            - residue / 2
            + mass * (1 + 2 * e) * B0 / (2 * (1 + e))
        )
        == 0
    )


def test_known_finite_pair_does_not_close_physical_IR_or_original_scope():
    assert len(audit.frontier()) == 9
    assert len(audit.matching()) == 148
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.qualifications()) == 6
    assert any(row["status"].startswith("REJECTED_") for row in audit.matching())
    assert "Coulomb" in audit.observable()["not_established"]
    assert "H-metric" in audit.observable()["not_established"]
    assert source.data()["checks"] is not source.previous.data()["checks"]
    assert len(source.previous.data()["checks"]) == 104


def test_every_scientific_payload_is_recursively_exact():
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
