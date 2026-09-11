"""Independent complete projectors, dimensional integrals and metric variations."""

import json
import runpy
from functools import cache

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_curved_dirac_state import audit as previous
from p8_vacuum_curved_dirac_state import verify as state_verify
from p8_vacuum_curved_dirac_stress import audit, reference, stress, subtraction, verify


@pytest.mark.parametrize("name", list(audit.residuals()))
def test_named_identity(name):
    assert audit.residuals()[name] == 0


@pytest.mark.parametrize("name", list(audit.gates()))
def test_proof_gate(name):
    assert audit.gates()[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_unsupported_inputs(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@cache
def frame_helper():
    return runpy.run_path(
        str(state_verify.ROOT / "tests" / "test_curved_dirac_state_science.py")
    )["exact_frames"]


def rotation(axis, angle):
    c, z = mp.cos(angle), mp.sin(angle)
    if axis == 0:
        return mp.matrix([[1, 0, 0], [0, c, -z], [0, z, c]])
    return mp.matrix([[c, 0, z], [0, 1, 0], [-z, 0, c]])


def projector_vector(angles):
    vec = mp.matrix([0, 0, 1])
    for j in reversed(range(len(angles))):
        vec = rotation(j % 2, angles[j]) * vec
    return vec


@pytest.mark.parametrize("A", ["0.001", "0.05", "0.24"])
@pytest.mark.parametrize("r", ["0.001", "0.1", "0.49"])
@pytest.mark.parametrize("phase", [0, 1, 2])
def test_twenty_finite_angle_projector_bounds(A, r, phase):
    with mp.workdps(90):
        A, r = mp.mpf(A), mp.mpf(r)
        angles = [
            mp.atan(A * r**j * (-1 if (j + phase) % 3 == 0 else 1)) for j in range(20)
        ]
        vec = projector_vector(angles)
        q0 = mp.tan(angles[0])
        b = mp.tan(angles[1]) * (1 + q0 * q0) ** mp.mpf("1.5")
        assert abs(sum(v * v for v in vec) - 1) < mp.mpf("1e-80")
        assert abs(vec[2] - mp.cos(angles[0])) <= 4 * A * A * r * r
        # Here b can exceed Ar slightly; the physical q1 relation is retained.
        assert abs(vec[0] - b) <= 4 * A * r**3 + 8 * A**3 * r


@pytest.mark.parametrize("time", ["-3", "0", "0.1", "4"])
@pytest.mark.parametrize("amplitude", ["0", "1000", "-1000"])
def test_full_twenty_frame_physical_projector_and_subtractions(time, amplitude):
    with mp.workdps(500):
        t, p, m, d, tau = map(mp.mpf, (time, "7e8", "1e8", amplitude, "0.2"))
        pairs = frame_helper()(t, p, m, d, tau, 20)
        angles = [mp.atan(g / e) for e, g in pairs[:20]]
        vec = projector_vector(angles)
        f = lambda z: (z / tau) / (1 + (z / tau) ** 8) ** mp.mpf("0.125")
        mass = lambda z: m + d * f(z)
        qfun = lambda z: p / (1 + z * z) ** 2
        omega = lambda z: mp.sqrt(qfun(z) ** 2 + mass(z) ** 2)
        q = qfun(t)
        M = mass(t)
        w = omega(t)
        qzero = lambda z: (
            (qfun(z) * mp.diff(mass, z) - mass(z) * mp.diff(qfun, z))
            / (2 * omega(z) ** 3)
        )
        q0 = qzero(t)
        b = mp.diff(qzero, t) / (2 * w)
        L = mp.sqrt(t * t + tau * tau)
        E = mp.sqrt(q * q + (mp.mpf(".99") * m) ** 2)
        B = abs(d) / tau * (1 + (t / tau) ** 8) ** (-mp.mpf(9) / 8) + m * L / (
            1 + t * t
        )
        A = 16 * q * B / E**3
        r = mp.mpf(163840) / (L * E)
        assert A < mp.mpf(".25") and r < mp.mpf(".5")
        assert abs(q0) <= A and abs(b) <= A * r
        assert abs(
            pairs[1][1] / pairs[1][0] - b / (1 + q0 * q0) ** mp.mpf("1.5")
        ) < mp.mpf("1e-470")
        rho = -2 * w * vec[2]
        P = -2 * (q * q * vec[2] + q * M * vec[0]) / (3 * w)
        rho_sub = -2 * w + w * q0 * q0
        P_sub = (
            -2 * q * q / (3 * w) + q * q * q0 * q0 / (3 * w) - 2 * q * M * b / (3 * w)
        )
        diagonal = 16 * E * A * A * r * r + mp.mpf("1.5") * E * A**4
        assert abs(rho - rho_sub) <= diagonal
        assert (
            abs(P - P_sub) <= diagonal / 3 + 4 * m * (4 * A * r**3 + 8 * A**3 * r) / 3
        )
        assert abs(vec[2] - mp.cos(angles[0])) <= 4 * A * A * r * r
        assert abs(vec[0] - b) <= 4 * A * r**3 + 8 * A**3 * r
        assert abs(mp.norm(vec) - 1) < mp.mpf("1e-480")


def test_actual_scale_tiny_angle_coherence_is_not_rounded_to_zero():
    with mp.workdps(1200):
        A = mp.mpf("1e-102")
        r = mp.mpf("1e-94")
        angles = [mp.atan(A * r**j) for j in range(20)]
        vec = projector_vector(angles)
        b = A * r * (1 + A * A) ** mp.mpf("1.5")
        dz = vec[2] - mp.cos(angles[0])
        dx = vec[0] - b
        assert dz != 0 and dx != 0
        assert abs(dz) <= 4 * A * A * r * r
        assert abs(dx) <= 4 * A * r**3 + 8 * A**3 * r


@pytest.mark.parametrize("mass", ["0.5", "1", "3", "10"])
def test_four_complete_physical_radial_integrals(mass):
    with mp.workdps(75):
        m = mp.mpf(mass)
        cases = (
            (4, mp.mpf("3.5"), 1 / (5 * m * m)),
            (6, mp.mpf("5.5"), 2 / (63 * m**4)),
            (3, 3, 1 / (4 * m * m)),
            (5, 5, 1 / (24 * m**4)),
        )
        for n, power, target in cases:
            value = mp.quad(
                lambda q, n=n, power=power: q**n / (q * q + m * m) ** power,
                [0, m, mp.inf],
            )
            assert abs(value - target) < mp.mpf("1e-65")


def test_pressure_UV_coefficient_from_fixed_comoving_derivative():
    q, M, Md, Mdd, H, Hd = s.symbols("q M Mdot Mddot H Hdot", real=True)
    w = s.sqrt(q * q + M * M)
    K = Md + H * M
    q0 = q * K / (2 * w**3)
    dot = lambda f: (
        s.diff(f, q) * (-H * q)
        + s.diff(f, M) * Md
        + s.diff(f, Md) * Mdd
        + s.diff(f, H) * Hd
    )
    P2 = q * q * q0 * q0 / (3 * w) - q * M * dot(q0) / (3 * w * w)
    A = K * K / 4 - s.Rational(3, 2) * H * M * K
    B = -M * (Mdd + (Hd - H * H) * M) / 2
    C = s.Rational(3, 2) * M * M * K * Md
    expected = (A * q**4 / w**7 + B * q * q / w**5 + C * q * q / w**7) / 3
    assert s.simplify(P2 - expected) == 0


@pytest.mark.parametrize("epsilon", ["0.125", "0.25", "0.375"])
def test_full_dimensional_radial_measure_at_nonzero_regulator(epsilon):
    with mp.workdps(75):
        e = mp.mpf(epsilon)
        d = 3 - 2 * e
        m = mp.mpf("1.7")
        for radial, denominator in (
            (1, mp.mpf("2.5")),
            (2, mp.mpf("3.5")),
            (1, mp.mpf("3.5")),
        ):
            exponent = d + 2 * radial - 2 * denominator

            def integrand(u, exponent=exponent, denominator=denominator, radial=radial):
                if u >= 0:
                    return mp.exp(exponent * u) / (1 + mp.exp(-2 * u)) ** denominator
                return mp.exp((d + 2 * radial) * u) / (1 + mp.exp(2 * u)) ** denominator

            measure = 2 / (4 * mp.pi) ** (d / 2) / mp.gamma(d / 2)
            direct = (
                measure
                * m**exponent
                * mp.quad(integrand, [-mp.inf, 0, 1, 5, 50, 500, mp.inf])
            )
            gamma = (
                m**exponent
                / (4 * mp.pi) ** (d / 2)
                * mp.gamma(radial + d / 2)
                * mp.gamma(denominator - radial - d / 2)
                / (mp.gamma(d / 2) * mp.gamma(denominator))
            )
            assert abs(direct - gamma) < mp.mpf("1e-63")


@pytest.mark.parametrize("mass", ["0.997", "1", "1.003", "2"])
@pytest.mark.parametrize("hubble", ["-0.3", "0", "0.2"])
def test_nonzero_regulator_finite_stress_matches_covariant_action(mass, hubble):
    with mp.workdps(95):
        M, H = mp.mpf(mass), mp.mpf(hubble)
        v, acc, Hd, mu = map(mp.mpf, ("0.07", "-0.03", "0.4", "1"))
        Q = 16 * mp.pi**2
        ell = mp.log(M * M / (mu * mu))
        e = mp.mpf("1e-27")
        d = 3 - 2 * e
        common = mp.exp((mp.euler - ell) * e) * mp.gamma(1 + e)
        J0 = 4 / Q * (1 - 2 * e / 3) * common / e
        J1 = 4 / Q * (1 - 2 * e / 3) * (1 - 2 * e / 5) * common / e
        J2 = 8 / (5 * Q * M * M) * (1 - 2 * e / 3) * common
        K = v + H * M
        A = K * K / 4 - 3 * H * M * K / 2
        B = -M * (acc + (Hd - H * H) * M) / 2
        C = 3 * M * M * K * v / 2
        bare_rho = K * K * J0 / 4
        bare_P = (A * J1 + B * J0 + C * J2) / d
        ct_rho = (-v * v - d * (d - 1) * M * M * H * H / 6 - 2 * d * H * M * v / 3) / (
            Q * e
        )
        ct_P = (
            -v * v / 3
            + 2 * M * acc / 3
            + M * M * ((d - 1) * Hd / 3 + d * (d - 1) * H * H / 6)
            + 2 * (d - 1) * H * M * v / 3
        ) / (Q * e)
        massfun = lambda t: M + v * t + acc * t * t / 2
        F = lambda t: (
            -(massfun(t) ** 2) * (mp.log(massfun(t) ** 2 / (mu * mu)) - 1) / (6 * Q)
        )
        kin = -v * v * (ell + mp.mpf(2) / 3) / Q
        rho = kin + 6 * F(0) * H * H + 6 * H * mp.diff(F, 0)
        P = (
            kin
            - 2 * F(0) * (2 * Hd + 3 * H * H)
            - 2 * mp.diff(F, 0, 2)
            - 4 * H * mp.diff(F, 0)
        )
        assert abs(bare_rho + ct_rho - rho) < mp.mpf("1e-23")
        assert abs(bare_P + ct_P - P) < mp.mpf("1e-23")
        # Setting d=3 in the counterterm pressure loses a nonzero finite term.
        frozen_ct = (
            -v * v / 3
            + 2 * M * acc / 3
            + M * M * (2 * Hd / 3 + H * H)
            + 4 * H * M * v / 3
        ) / (Q * e)
        missing = (2 * M * M * Hd / 3 + 5 * M * M * H * H / 3 + 4 * H * M * v / 3) / Q
        assert abs((frozen_ct - ct_P) - missing) < mp.mpf("1e-23")
        assert abs(missing) > mp.mpf("1e-5")


def test_dimensional_FR_and_kinetic_counterterms_by_independent_metric_variation():
    # Reduced FR action after IBP in logarithmic scale factor b=log(a).
    dim = s.Symbol("d", positive=True)
    scale, lapse = s.symbols("b N", real=True)
    v, acc, Ndot, F, Fdot, Fddot = s.symbols("v acc Ndot F Fdot Fddot", real=True)
    D = lambda f: (
        s.diff(f, scale) * v
        + s.diff(f, v) * acc
        + s.diff(f, lapse) * Ndot
        + s.diff(f, F) * Fdot
        + s.diff(f, Fdot) * Fddot
    )
    L = s.exp(dim * scale) * (dim * (dim - 1) * F * v * v + 2 * dim * Fdot * v) / lapse
    rho = s.simplify(-s.diff(L, lapse).subs(lapse, 1) / s.exp(dim * scale))
    P = s.simplify(
        (s.diff(L, scale) - D(s.diff(L, v))).subs({lapse: 1, Ndot: 0})
        / (dim * s.exp(dim * scale))
    )
    assert s.simplify(rho - dim * (dim - 1) * F * v * v - 2 * dim * Fdot * v) == 0
    assert (
        s.simplify(
            P
            + 2 * (dim - 1) * F * acc
            + dim * (dim - 1) * F * v * v
            + 2 * Fddot
            + 2 * (dim - 1) * v * Fdot
        )
        == 0
    )
    M, Md, Mdd = s.symbols("M Md Mdd", real=True)
    Fsub = {F: -M * M / 6, Fdot: -M * Md / 3, Fddot: -(Md * Md + M * Mdd) / 3}
    rhoct = s.expand(rho.subs(Fsub) - Md * Md)
    Pct = s.expand(P.subs(Fsub) - Md * Md)
    assert (
        s.expand(
            rhoct
            + Md * Md
            + dim * (dim - 1) * M * M * v * v / 6
            + 2 * dim * v * M * Md / 3
        )
        == 0
    )
    assert (
        s.expand(
            Pct
            + Md * Md / 3
            - 2 * M * Mdd / 3
            - M * M * ((dim - 1) * acc / 3 + dim * (dim - 1) * v * v / 6)
            - 2 * (dim - 1) * v * M * Md / 3
        )
        == 0
    )


def test_full_Euler_lapse_and_scale_variations_before_dimension_limit():
    dim = s.Symbol("d", positive=True)
    b, N = s.symbols("b N", real=True)
    v, acc, jerk, Ndot, Nddot = s.symbols("v acc jerk Ndot Nddot", real=True)
    C = dim * (dim - 1) * (dim - 2)
    L = (
        C
        * s.exp(dim * b)
        * (4 * v * v * acc / N**3 - 4 * v**3 * Ndot / N**4 + (dim + 1) * v**4 / N**3)
    )
    D = lambda f: (
        s.diff(f, b) * v
        + s.diff(f, v) * acc
        + s.diff(f, acc) * jerk
        + s.diff(f, N) * Ndot
        + s.diff(f, Ndot) * Nddot
    )
    rho = s.simplify(
        -(s.diff(L, N) - D(s.diff(L, Ndot))).subs({N: 1, Ndot: 0, Nddot: 0})
        / s.exp(dim * b)
    )
    P = s.simplify(
        (s.diff(L, b) - D(s.diff(L, v)) + D(D(s.diff(L, acc)))).subs(
            {N: 1, Ndot: 0, Nddot: 0}
        )
        / (dim * s.exp(dim * b))
    )
    assert s.factor(rho + C * (dim - 3) * v**4) == 0
    assert s.factor(P - C * (dim - 3) * (v**4 + 4 * v * v * acc / dim)) == 0
    eps = s.Symbol("eps")
    rho_finite = s.limit(s.Rational(11, 720) * rho.subs(dim, 3 - 2 * eps) / eps, eps, 0)
    P_finite = s.limit(s.Rational(11, 720) * P.subs(dim, 3 - 2 * eps) / eps, eps, 0)
    assert s.factor(rho_finite - s.Rational(11, 60) * v**4) == 0
    assert (
        s.factor(
            P_finite + s.Rational(11, 60) * (v**4 + s.Rational(4, 3) * v * v * acc)
        )
        == 0
    )
    assert rho.subs(dim, 3) == 0 and rho_finite != 0


def test_curvature_basis_has_a_real_finite_R_squared_difference():
    dim = s.Symbol("D")
    R2, Ric2, Riem2 = s.symbols("R2 Ric2 Riem2")
    literal = (5 * R2 - 8 * Ric2 - 7 * Riem2) / 360
    E = Riem2 - 4 * Ric2 + R2
    C = Riem2 - 4 * Ric2 / (dim - 2) + 2 * R2 / ((dim - 1) * (dim - 2))
    eps = s.Symbol("eps")
    difference = s.limit(
        (literal - s.Rational(11, 360) * E + C / 20).subs(dim, 4 - 2 * eps) / (2 * eps),
        eps,
        0,
    )
    assert s.simplify(difference + Ric2 / 20 - R2 / 72) == 0
    # On conformally flat four-dimensional metrics: C^2=0.
    assert (
        s.simplify(
            difference
            - E / 40
            + R2 / 360
            - (Riem2 - 2 * Ric2 + R2 / 3) * (-s.Rational(1, 40))
        )
        == 0
    )
    assert difference.subs({R2: 36, Ric2: 12}) != 0


def test_finite_R_squared_choice_changes_actual_bounce_energy():
    b, N, v, acc, jerk, Ndot, Nddot = s.symbols("b N v acc jerk Ndot Nddot", real=True)
    R = -6 * (acc / N**2 - v * Ndot / N**3 + 2 * v * v / N**2)
    L = s.exp(3 * b) * N * R * R
    D = lambda f: (
        s.diff(f, b) * v
        + s.diff(f, v) * acc
        + s.diff(f, acc) * jerk
        + s.diff(f, N) * Ndot
        + s.diff(f, Ndot) * Nddot
    )
    rho = s.simplify(-(s.diff(L, N) - D(s.diff(L, Ndot))) / s.exp(3 * b))
    assert rho.subs({N: 1, Ndot: 0, Nddot: 0, v: 0, acc: 4}) == 576
    assert -s.Rational(576, 360) == -s.Rational(8, 5)


@pytest.mark.parametrize("x", ["-10", "-1", "-0.1", "0", "1e-180", "0.1", "1", "10"])
def test_actual_scale_complete_local_profile_and_Newton_reference(x):
    with mp.workdps(850):
        x = mp.mpf(x)
        m = mp.mpf(10) ** 200
        d = 3 * mp.mpf(10) ** 197
        tau = mp.mpf(10) ** -100
        t = tau * x
        H = 4 * t / (1 + t * t)
        Hd = 4 * (1 - t * t) / (1 + t * t) ** 2
        shape = lambda z: z / (1 + z**8) ** mp.mpf("0.125")
        Q = 16 * mp.pi**2
        bound = reference.enclosures(10**200, 3 * 10**197, s.Rational(1, 10**100))
        limit = lambda key: mp.mpf(str(s.N(bound[key], 840)))
        pair_rho = mp.mpf(0)
        pair_P = mp.mpf(0)
        for sign in (-1, 1):
            M = m + sign * d * shape(x)
            Md = sign * d / tau * mp.diff(shape, x)
            Mdd = sign * d / (tau * tau) * mp.diff(shape, x, 2)
            ell = mp.log(M * M / (m * m))
            D = M * M * (ell - 1) + m * m
            Dt = 2 * M * Md * ell
            Dtt = 2 * (Md * Md + M * Mdd) * ell + 4 * Md * Md
            rho = (-Md * Md * (ell + mp.mpf(2) / 3) - D * H * H - H * Dt) / Q
            P = (
                -Md * Md * (ell + mp.mpf(2) / 3)
                + D * (2 * Hd + 3 * H * H) / 3
                + Dtt / 3
                + 2 * H * Dt / 3
            ) / Q
            assert abs(ell) <= limit("absolute_logarithm")
            assert abs(D) <= limit("absolute_D_function")
            assert abs(Dt) <= limit("absolute_D_time_derivative")
            assert abs(Dtt) <= limit("absolute_D_second_time_derivative")
            if x != 0:
                assert D != 0
            pair_rho += rho
            pair_P += P
        assert abs(3 * pair_rho) <= limit("local_two_derivative_energy")
        assert abs(3 * pair_P) <= limit("local_two_derivative_pressure")
        if x == 0:
            assert pair_rho + pair_P == 0
            assert pair_rho < 0 and pair_P > 0


def test_reference_zero_removes_only_local_second_order_not_curved_species():
    m = 10**200
    tau = s.Rational(1, 10**100)
    value = stress.enclosures(m, 0, tau)
    assert value["complete_potential_upper"] == 0
    assert (
        value["complete_local_two_derivative_components"]["local_two_derivative_energy"]
        == 0
    )
    assert (
        value["complete_local_two_derivative_components"][
            "local_two_derivative_pressure"
        ]
        == 0
    )
    assert (
        value["complete_remainder_components"]["complete_subtracted_energy_remainder"]
        > 0
    )
    assert value["finite_Euler_energy_upper"] > 0
    assert value["complete_absolute_curved_energy_upper"] > 0
    # The finite reference coefficient is nonzero despite its zero-field subtraction.
    d = subtraction.data()
    symbols = d["symbols"]
    point = {
        symbols["M"]: symbols["m"],
        symbols["Mdot"]: 0,
        symbols["Mddot"]: 0,
        symbols["ell"]: 0,
    }
    assert s.simplify(d["finite_MS_rho_second_order"].subs(point)) != 0
    assert s.simplify(d["Newton_referenced_rho_second_order"].subs(point)) == 0


def test_actual_twenty_frame_projector_precision_refinement():
    outputs = []
    for precision in (4800, 5200):
        with mp.workdps(precision):
            m = mp.mpf(10) ** 200
            d = 3 * mp.mpf(10) ** 197
            tau = mp.mpf(10) ** -100
            t = mp.mpf(".1") * tau
            p = 7 * m
            pairs = frame_helper()(t, p, m, d, tau, 20)
            vec = projector_vector([mp.atan(g / e) for e, g in pairs[:20]])
            q0 = pairs[0][1] / pairs[0][0]
            b = pairs[1][1] / pairs[1][0] * (1 + q0 * q0) ** mp.mpf("1.5")
            dz = vec[2] - 1 / mp.sqrt(1 + q0 * q0)
            dx = vec[0] - b
            q = p / (1 + t * t) ** 2
            E = mp.sqrt(q * q + (mp.mpf(".99") * m) ** 2)
            L = mp.sqrt(t * t + tau * tau)
            B = d / tau * (1 + (t / tau) ** 8) ** (-mp.mpf(9) / 8) + m * L / (1 + t * t)
            A = 16 * q * B / E**3
            r = mp.mpf(163840) / (L * E)
            assert pairs[-1][1] != 0 and pairs[-1][0] - pairs[-2][0] != 0
            assert dz != 0 and dx != 0
            assert abs(dz) <= 4 * A * A * r * r
            assert abs(dx) <= 4 * A * r**3 + 8 * A**3 * r
            outputs.append((+dz, +dx))
    with mp.workdps(5200):
        for lo, hi in zip(*outputs):
            assert abs((lo - hi) / hi) < mp.mpf("1e-100")


def test_complete_exact_sum_scope_and_serializer_boundaries():
    assert audit.frontier() == previous.frontier()
    assert audit.matching()[:-1] == previous.matching()
    assert (
        "NOT_FULL_INTERACTING_PARENT_OR_BACKGROUND_RESPONSE"
        in audit.matching()[-1]["status"]
    )
    d = stress.data()
    actual = d["actual_uniform_exact_rational_enclosures"]
    assert 0 < actual["complete_absolute_curved_energy_upper"] < 10**789
    assert 0 < actual["complete_absolute_curved_pressure_upper"] < 10**789
    assert d["active_copies"] == 6 and d["all_curved_copies"] == 42
    assert "controlled C2 background" in d["boundary"]
    assert "zero bounce density" in d["boundary"]
    for mod in audit.MODULES:
        json.dumps(verify.serialize(verify.payload(mod.data())))
    assert all(type(v) is str for v in d["decimal_diagnostics_only"].values())
    assert all(
        type(v) is str for v in reference.data()["decimal_diagnostics_only"].values()
    )
    with pytest.raises(ValueError, match="Inexact or nonfinite"):
        verify.serialize(s.Float(1))
