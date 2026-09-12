"""Independent scalar Legendre, fixed-profile, metric-flow and infrared controls."""

import mpmath as mp
import pytest
import sympy as s
from p8_affine import verify as serializer
from p8_vacuum_affine_reduced_scalar_hamiltonian import (
    audit,
    estimates,
    pullback,
    scalar,
)

PACKETS = (
    scalar.spatial_data,
    scalar.retuning_data,
    scalar.hamiltonian_data,
    pullback.ward_data,
    pullback.clock_data,
    estimates.projection_data,
    estimates.phase_data,
)


@pytest.mark.parametrize("name,value", list(audit.residuals().items()))
def test_integrated_exact_residual(name, value):
    row = list(value) if isinstance(value, s.MatrixBase) else [value]
    assert all(s.factor(v) == 0 for v in row), name


@pytest.mark.parametrize("name,value", list(audit.gates().items()))
def test_integrated_proof_gate(name, value):
    assert value is True, name


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_every_unsupported_input_and_false_closure_rejected(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


def test_actual_scope_and_mathematical_compact_band_acceptance():
    for t in (-s.Rational(1, 2), 0, s.Rational(1, 2)):
        assert audit.require_scope(t)[0] == t
    for p in (0, s.Rational(1, 100), s.Integer(10) ** 6):
        assert audit.require_band(p) == p
    assert audit.validate_scope(audit.frontier(), audit.matching())
    assert len(audit.frontier()) == 9
    assert audit.matching()[:-1] == audit.previous.matching()


def test_actual_phase_comparison_has_only_nonnegative_momentum_powers():
    M0, M2 = estimates.phase_matrices()
    assert M0.shape == M2.shape == (3, 4)
    p = s.Symbol("p", real=True)
    combined = M0 + p * p * M2
    assert combined[2, :] == s.Matrix([[0, 0, s.Rational(1, 2), 0]])
    assert all(s.diff(v, p, 3) == 0 for v in combined)
    assert all(not s.denom(s.factor(v)).has(p) for v in combined)


@pytest.mark.parametrize("packet", PACKETS)
def test_exact_packet_and_exact_serialization(packet):
    d = packet()
    for key, value in d["checks"].items():
        row = list(value) if isinstance(value, s.MatrixBase) else [value]
        assert all(s.factor(v) == 0 for v in row), key
    assert all(bool(v) for v in d["gates"].values())
    serializer.serialize(
        {key: value for key, value in d.items() if key not in ("checks", "gates")}
    )


def number(value):
    value = s.N(value, mp.mp.dps)
    return mp.mpc(str(s.re(value)), str(s.im(value)))


@pytest.mark.parametrize("theta", ("-1.7", "0", "1.9"))
@pytest.mark.parametrize("qval", ("0.000000000001", "2", "10000"))
def test_independent_full_constraint_saddle_and_legendre(theta, qval):
    with mp.workdps(80):
        th, qq = mp.mpf(theta), mp.mpf(qval)
        JJ, dJ, ell, w, E = map(mp.mpf, ("1.3", "0.07", "0.08", "-0.031", "-0.23"))
        A, T = mp.mpf("-0.011"), mp.mpf("0.019")
        v, sigma, pv, ps = map(mp.mpf, ("0.13", "-0.17", "0.29", "0.23"))

        def L(vd, sd, n, b):
            # Independently entered complete coefficient-sector ADM saddle.
            return (
                -3 * vd * vd
                + (JJ + w * w / 2 - 3 * th * th) * n * n
                + 6 * th * n * vd
                + sd * sd / 2
                + w * n * sd
                - 3 * ell * vd * sigma
                + 2 * b * (vd - th * n)
                + ell * b * sigma
                + qq * v * v
                + 2 * E * qq * n * v
                - qq * sigma * sigma / 2
                + dJ * n * n
                + 3 * T * n * v
                + mp.mpf("4.5") * A * v * v
            )

        def eq(values):
            return mp.matrix(
                [
                    mp.diff(
                        L, tuple(values), tuple(1 if j == i else 0 for j in range(4))
                    )
                    - ([pv, ps, 0, 0][i])
                    for i in range(4)
                ]
            )

        zero = mp.matrix([0, 0, 0, 0])
        rhs = -eq(zero)
        M = mp.matrix(4, 4)
        for j in range(4):
            unit = mp.matrix([1 if i == j else 0 for i in range(4)])
            M[:, j] = eq(unit) - eq(zero)
        sol = mp.lu_solve(M, rhs)
        assert mp.norm(eq(sol)) < mp.mpf("1e-65")
        substitutions = {
            scalar.t: 0,
            scalar.Theta: s.Rational(theta),
            scalar.q: s.Rational(qval),
            scalar.J0: s.Rational(13, 10),
            scalar.dJ: s.Rational(7, 100),
            scalar.w: s.Rational(-31, 1000),
            scalar.E: s.Rational(-23, 100),
            scalar.A: s.Rational(-11, 1000),
            scalar.Tc: s.Rational(19, 1000),
            scalar.v: s.Rational(13, 100),
            scalar.sigma: s.Rational(-17, 100),
            scalar.pv: s.Rational(29, 100),
            scalar.ps: s.Rational(23, 100),
        }
        # ell is the actual fixed background expression; replace it before t.
        h = (
            scalar.reduced_hamiltonian()
            .subs(scalar.ell, s.Rational(2, 25))
            .subs(substitutions)
        )
        nformula = (
            (scalar.lapse_numerator() / (2 * (scalar.J0 + scalar.dJ)))
            .subs(scalar.ell, s.Rational(2, 25))
            .subs(substitutions)
        )
        actual = pv * sol[0] + ps * sol[1] - L(*sol)
        assert abs(actual - number(h)) < mp.mpf("1e-60")
        assert abs(sol[2] - number(nformula)) < mp.mpf("1e-65")
        assert abs(sol[3] - pv / 2) < mp.mpf("1e-65")
        if theta == "0":
            assert abs(sol[0] + ell * sigma / 2) < mp.mpf("1e-65")


@pytest.mark.parametrize("time", ("-0.5", "0", "0.5"))
def test_full_order_1024_switch_retuning_metric_volume_hessian(time):
    with mp.workdps(85):
        tt = mp.mpf(time)
        delta = 1 / (2 * (1 + tt * tt) ** 3)
        A, B, PV = map(mp.mpf, ("-0.031", "0.017", "0.043"))
        nv, vv = mp.mpf("0.21"), mp.mpf("-0.13")

        def full(e):
            N = 1 + e * nv
            X = N**-2
            R = 1 + 2 * delta * (X - 1)
            switch = X**1024 / (X**1024 + (1 - X) ** 1024)
            F = -PV + switch * (A + PV + B * (X - 1))
            return N * mp.exp(3 * e * vv) * R ** (-mp.mpf(3) / 4) * F

        observed = mp.diff(full, 0, 2) / 2
        wanted = scalar.retuning_data()["quadratic_addition"].subs(
            {
                s.Symbol("delta", real=True): s.Rational(str(delta)),
                scalar.A: s.Rational(-31, 1000),
                s.Symbol("B", real=True): s.Rational(17, 1000),
                scalar.n: s.Rational(21, 100),
                scalar.v: s.Rational(-13, 100),
            }
        )
        assert abs(observed - number(wanted)) < mp.mpf("1e-70")
        dropped = (
            scalar.retuning_data()["DeltaJ"].subs(
                {
                    s.Symbol("delta", real=True): s.Rational(str(delta)),
                    scalar.A: s.Rational(-31, 1000),
                    s.Symbol("B", real=True): s.Rational(17, 1000),
                }
            )
            * s.Rational(21, 100) ** 2
        )
        assert abs(observed - number(dropped)) > mp.mpf("1e-4")


def literal_metric(n, beta, Q, scale):
    h = scale * scale * mp.expm(Q)
    g = mp.zeros(4)
    g[0, 0] = n * n - (beta.T * h * beta)[0]
    g[0, 1:] = -beta.T * h
    g[1:, 0] = -h * beta
    g[1:, 1:] = -h
    return g


def pair(E, G):
    return mp.fsum(E[i, j] * G[i, j] for i in range(4) for j in range(4))


def sympy_ward_scalar_value(expr, vals):
    sub = {}
    for name, value in vals.items():
        if name in ("Aprime", "Bprime"):
            sub[s.diff(s.Function(name[0])(pullback.t), pullback.t)] = s.Rational(
                str(value)
            )
        else:
            sub[s.Function(name)(pullback.t)] = s.Rational(str(value))
    return number(expr.subs(sub, simultaneous=True).subs(pullback.t, s.Rational(1, 5)))


@pytest.mark.parametrize("momentum", ("0.000000000000001", "0.7", "37"))
def test_independent_density_coordinate_push_and_full_metric_chart(momentum):
    with mp.workdps(105):
        p = mp.mpf(momentum)
        t = mp.mpf("0.2")
        a = (1 + t * t) ** 2
        H = 4 * t / (1 + t * t)
        nD, zD, bD, eta, c, nd, cd = map(
            mp.mpf, ("0.13", "-0.11", "0.17", "0.19", "-0.23", "0.29", "0.31")
        )
        Af = lambda time: 2 + time + 3 * time * time
        Bf = lambda time: -1 + 2 * time - time * time
        density = lambda time: mp.diag([Af(time), Bf(time), Bf(time), Bf(time)])
        J = mp.zeros(4)
        J[0, 0], J[0, 1] = nd, 1j * p * eta
        J[1, 0], J[1, 1] = -1j * cd / p, c

        def pushed(e):
            X = mp.eye(4) + e * J
            return mp.det(X) * (X**-1) * density(t + e * eta) * (X.T**-1)

        dE = mp.matrix(4, 4)
        for i in range(4):
            for j in range(4):
                dE[i, j] = mp.diff(lambda e, i=i, j=j: pushed(e)[i, j], 0)
        betaD = mp.matrix([1j * bD / p, 0, 0])
        betaG = mp.matrix([-1j * cd / p - 1j * p * eta / a**2, 0, 0])
        QD = 2 * zD * mp.eye(3)
        QG = 2 * H * eta * mp.eye(3) + mp.diag([2 * c, 0, 0])

        def metric(e, h):
            return literal_metric(
                1 + e * nD + h * nd, e * betaD + h * betaG, e * QD + h * QG, a
            )

        first, contact = mp.zeros(4), mp.zeros(4)
        for i in range(4):
            for j in range(4):
                first[i, j] = mp.diff(lambda e, i=i, j=j: metric(e, 0)[i, j], 0)
                contact[i, j] = mp.diff(
                    lambda e, h, i=i, j=j: metric(e, h)[i, j], (0, 0), (1, 1)
                )
        observed = pair(dE, first) + pair(density(t), contact)
        vals = {
            "nD": nD,
            "zD": zD,
            "bD": bD,
            "etaG": eta,
            "cG": c,
            "nG": nd,
            "A": Af(t),
            "B": Bf(t),
            "Aprime": mp.diff(Af, t),
            "Bprime": mp.diff(Bf, t),
        }
        wanted = sympy_ward_scalar_value(
            pullback.ward_data()["scalar_ordered_source"], vals
        )
        assert abs(observed - wanted) < mp.mpf("1e-65")
        assert abs(pair(dE, first) - observed) > mp.mpf("1e-4")
        if p < mp.mpf("1e-10"):
            assert abs(pair(dE, first)) > mp.mpf("1e20")
            assert abs(observed) < 100


@pytest.mark.parametrize("time", ("-0.5", "0", "0.5"))
def test_full_clock_composition_second_contact(time):
    with mp.workdps(75):
        t = mp.mpf(time)
        a = (1 + t * t) ** 2
        d = 1 / (2 * (1 + t * t) ** 3)
        nD, nG, vD, vG = map(mp.mpf, ("0.13", "-0.17", "0.19", "0.23"))
        EE = mp.diag([-2, 3, 3, 3])

        def full(e, h, linear=False):
            N = 1 + e * nD + h * nG
            clock = (
                2 * d * (e * nD + h * nG)
                if linear
                else -mp.log(1 + 2 * d * (N**-2 - 1)) / 2
            )
            Q = (2 * (e * vD + h * vG) + clock) * mp.eye(3)
            return pair(EE, literal_metric(N, mp.zeros(3, 1), Q, a))

        difference = mp.diff(lambda e, h: full(e, h), (0, 0), (1, 1)) - mp.diff(
            lambda e, h: full(e, h, True), (0, 0), (1, 1)
        )
        expected = -9 * a * a * (8 * d * d - 6 * d) * nD * nG
        assert abs(difference - expected) < mp.mpf("1e-60")
        assert abs(difference) > mp.mpf("1e-3")


@pytest.mark.parametrize(
    "side,end", (("retarded", -s.Rational(1, 2)), ("advanced", s.Rational(1, 2)))
)
def test_both_prepared_endpoints_and_spatial_projector(side, end):
    t = pullback.t
    eta = (t - end) ** 4
    c = (t - end) ** 5
    n = s.diff(eta, t)
    # Construct the input from a nonzero prepared gauge fixture.
    b = s.diff(c, t) + pullback.p**2 * eta / pullback.a**2
    assert s.factor(pullback.k.primitive(n, side) - eta) == 0
    # Integrate c' after exact cancellation, avoiding an unrelated rational integral.
    assert (
        s.factor(
            pullback.k.primitive(
                s.cancel(b - pullback.p**2 * eta / pullback.a**2), side
            )
            - c
        )
        == 0
    )
    assert eta.subs(t, end) == 0 and c.subs(t, end) == 0
    q = pullback.projected(0, eta, c)
    assert all(s.factor(v) == 0 for v in q.subs(t, end))


@pytest.mark.parametrize("direction", ((1, 0, 0), (1, 2, 3), (-2, 4, 1)))
def test_bounded_projector_arbitrarily_small_nonzero_momentum(direction):
    vec = s.Matrix(direction)
    expected = vec * vec.T / (vec.dot(vec))
    for radius in (s.Rational(1, 10**40), s.Rational(1, 3), s.Integer(10) ** 20):
        P = radius * vec
        actual = P * P.T / (P.dot(P))
        assert actual == expected
        assert s.trace(actual.T * actual) == 1
    assert (expected * expected) == expected


def test_regular_bounce_hamiltonian_is_not_claimed_coercive():
    h = scalar.reduced_hamiltonian().subs(
        {
            scalar.Theta: 0,
            scalar.q: 1,
            scalar.v: 0,
            scalar.ps: 0,
            scalar.sigma: 1,
            scalar.pv: 20,
            scalar.J0: 1,
            scalar.dJ: 0,
            scalar.Tc: 0,
            scalar.A: 0,
            scalar.w: 0,
            scalar.t: 0,
        }
    )
    assert h == -s.Rational(203, 400)
    assert h < 0
    assert "not a coercive" in scalar.hamiltonian_data()["inverse_boundary"]


def test_no_global_zero_mode_or_high_jet_or_inverse_overclaim():
    assert "physical homogeneous" in pullback.clock_data()["zero_transfer"]
    assert "no unproved numerical13-jet" in estimates.phase_data()["phase_input_maps"]
    assert (
        "not a chosen physical EFT cutoff"
        in estimates.phase_data()["compact_band_scope"]
    )
    assert (
        "not identified with the literal sharp-band"
        in estimates.projection_data()["regulator_boundary"]
    )
    assert "nonlocal" in scalar.hamiltonian_data()["inverse_boundary"]


@pytest.mark.parametrize("time", ("-0.5", "0", "0.5"))
@pytest.mark.parametrize("momentum", ("0.000000001", "3"))
def test_literal_three_by_three_extrinsic_curvature_spatial_mean(time, momentum):
    with mp.workdps(75):
        tt, pp = mp.mpf(time), mp.mpf(momentum)
        aa = (1 + tt * tt) ** 2
        HH = 4 * tt / (1 + tt * tt)
        dd = 1 / (2 * (1 + tt * tt) ** 3)
        charge = 1 / (10 * (1 + tt * tt) ** 6)
        vv, nn, bb, vd, ss, sd = map(
            mp.mpf, ("0.13", "-0.17", "0.19", "0.23", "-0.29", "0.31")
        )
        bjets = [
            number(v.subs(scalar.t, s.Rational(time))).real
            for v in scalar.tree.clock_coefficients()["coefficient_clock_N_jets"]["b"]
        ]

        def integrand(angle):
            co, si = mp.cos(angle), mp.sin(angle)

            def density(e, spatial):
                N = 1 + e * nn * co
                RR = 1 + 2 * dd * (N**-2 - 1)
                volume = mp.exp(3 * e * vv * co)
                beta = e * bb * si / pp if spatial else 0
                betax = e * bb * co if spatial else 0
                vx = -e * pp * vv * si if spatial else 0
                vxx = -e * pp * pp * vv * co if spatial else 0
                sx = -e * pp * ss * si if spatial else 0
                rate = HH + e * vd * co - beta * vx
                KK = mp.diag([(rate - betax) / N, rate / N, rate / N])
                trace = mp.fsum(KK[i, i] for i in range(3))
                square = mp.fsum(
                    KK[i, j] * KK[j, i] for i in range(3) for j in range(3)
                )
                bcoef = mp.fsum(
                    bjets[j] * (N - 1) ** j / mp.factorial(j) for j in range(3)
                )
                grav = (
                    N
                    * volume
                    * (
                        RR ** mp.mpf("0.25") * (square - trace * trace) / 2
                        + bcoef * trace
                    )
                )
                R3 = mp.exp(-2 * e * vv * co) * (-4 * vxx - 2 * vx * vx) / aa**2
                curvature = N * volume * RR ** mp.mpf("0.75") * R3 / 2
                matter = (
                    volume
                    * RR ** (-mp.mpf("0.75"))
                    * (charge + e * sd * co - beta * sx) ** 2
                    / (2 * N)
                )
                gradient = (
                    -N
                    * mp.exp(e * vv * co)
                    * RR ** (-mp.mpf("0.25"))
                    * sx
                    * sx
                    / (2 * aa**2)
                )
                return grav + curvature + matter + gradient

            return mp.diff(lambda e: density(e, True) - density(e, False), 0, 2) / 2

        observed = mp.quad(integrand, [0, mp.pi, 2 * mp.pi]) / (2 * mp.pi)
        theta = HH - tt / (1 + tt * tt) ** 4
        qq = pp * pp / aa**2
        expected = (
            qq * vv * vv
            + 2 * (1 - 3 * dd) * qq * nn * vv
            - qq * ss * ss / 2
            + 2 * bb * (vd - theta * nn)
            + charge * bb * ss
        ) / 2
        assert abs(observed - expected) < mp.mpf("1e-58")


@pytest.mark.parametrize("momentum", ("0.000000000000001", "0.7", "37"))
def test_independent_synchronous_detector_coordinate_pull_and_chart(momentum):
    with mp.workdps(105):
        pp, tt = mp.mpf(momentum), mp.mpf("0.2")
        scale = lambda u: (1 + u * u) ** 2
        hubble = lambda u: 4 * u / (1 + u * u)
        zz = lambda u: 1 + 2 * u + u * u
        etaG = lambda u: -1 + u + 3 * u * u
        cG = lambda u: 2 - u + u * u
        etaD, cD, nd, cd = map(mp.mpf, ("0.13", "-0.17", "0.19", "0.23"))
        aa = scale(tt)
        HH = hubble(tt)
        EE = mp.diag([2, -3, -3, -3])

        def QQ(u):
            return 2 * (zz(u) - hubble(u) * etaG(u)) * mp.eye(3) - 2 * cG(u) * mp.diag(
                [1, 0, 0]
            )

        def hG(u, x):
            out = mp.zeros(4)
            out[1:, 1:] = -(scale(u) ** 2) * QQ(u) * mp.exp(1j * pp * x)
            return out

        jac = mp.zeros(4)
        jac[0, 0], jac[0, 1] = nd, -1j * pp * etaD
        jac[1, 0], jac[1, 1] = 1j * cd / pp, cD
        chi = 1j * cD / pp

        def pulled(e):
            X = mp.eye(4) + e * jac
            return pair(EE, X.T * hG(tt + e * etaD, e * chi) * X)

        lie = mp.diff(pulled, 0)
        beta = mp.matrix([1j * cd / pp + 1j * pp * etaD / aa**2, 0, 0])
        QD = 2 * HH * etaD * mp.eye(3) + 2 * cD * mp.diag([1, 0, 0])

        def metric(e, h):
            return pair(
                EE, literal_metric(1 + e * nd, e * beta, e * QD + h * QQ(tt), aa)
            )

        contact = mp.diff(metric, (0, 0), (1, 1))
        observed = -lie + contact
        trace = lambda u: mp.fsum(QQ(u)[i, i] for i in range(3))
        expected = -3 * aa * aa * (etaD * mp.diff(trace, tt) - cD * trace(tt))
        assert abs(observed - expected) < mp.mpf("1e-65")
        assert abs(contact) > mp.mpf("1e-3")
