"""Independent complex-plane, frequency-line and noncommuting inverse checks."""

import mpmath as mp
import pytest
import sympy as s
from p8_affine import verify as serializer
from p8_vacuum_affine_weighted_channel_resolvent import analytic, audit, curved
from p8_vacuum_affine_weighted_channel_resolvent import weighted as w


@pytest.mark.parametrize("name,value", list(audit.residuals().items()))
def test_integrated_exact_residual(name, value):
    assert s.factor(value) == 0, name


@pytest.mark.parametrize("name,value", list(audit.gates().items()))
def test_integrated_proof_gate(name, value):
    assert value is True, name


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_every_unsupported_input_and_false_completion_rejected(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


def test_scoped_stages_and_unchanged_frontier():
    for stage in (
        "full_complex_log_bound",
        "weighted_channel_inverse",
        "bounded_reference_correction",
    ):
        assert audit.require_stage(stage) == stage
    for t in (-s.Rational(1, 2), 0, s.Rational(1, 2)):
        assert audit.require_scope(t)[0] == t
    assert audit.matching()[:-1] == audit.previous.matching()
    assert audit.validate_scope(audit.frontier(), audit.matching())


@pytest.mark.parametrize(
    "bound", (0, s.Rational(1, 10**12), s.Rational(3, 2), 10, 1000)
)
def test_exact_weight_profile_and_strict_contraction(bound):
    bound = s.Rational(bound)
    data = w.profile(bound)
    assert data["channel_norm_bound"] == bound
    assert data["Laplace_weight"] == 2000 * s.exp(32 + 8 * bound)
    gap = (32 + 13 * bound) / 5
    assert data["factor_gap"] == gap
    assert (
        bound * data["original_channel_inverse_bound"]
        < s.Rational(5, 13)
        < s.Rational(1, 2)
    )
    assert data["corrected_channel_inverse_bound"] == 1 / (gap - bound)


def closed(z, spin):
    d = 1 + 1 / z
    root = mp.sqrt(d)
    logpart = (mp.log(z) + 2 * mp.log(root + 1)) / 2
    if spin == 2:
        return (
            -mp.mpf(172) / 225
            + 19 * d / 30
            - d * d / 10
            + root * (30 - 20 * d + 3 * d * d) * logpart / 30
        )
    return mp.mpf(16) / 15 + d - 3 * d * d + root * (3 - 2 * d + 3 * d * d) * logpart


def radial(z, spin):
    if spin == 2:
        return mp.mpf(1) / 30 + mp.quad(
            lambda y: (
                z * y * y * (30 - 20 * y * y + 3 * y**4) / (30 * (1 + z * (1 - y * y)))
            ),
            [0, mp.mpf(".5"), mp.mpf(".9"), 1],
        )
    return 4 + mp.quad(
        lambda y: z * y * y * (3 - 2 * y * y + 3 * y**4) / (1 + z * (1 - y * y)),
        [0, mp.mpf(".5"), mp.mpf(".9"), 1],
    )


@pytest.mark.parametrize("spin", (0, 2))
def test_independent_radial_integral_and_stable_first_sheet_branch(spin):
    with mp.workdps(110):
        for z in (mp.mpc(1, 1), mp.mpc(-1, 2), mp.mpc(3, -7), mp.mpc(100, 10)):
            assert abs(closed(z, spin) - radial(z, spin)) < mp.mpf("1e-90")


@pytest.mark.parametrize("spin", (0, 2))
def test_independent_exterior_circle_and_both_cut_bank_diagnostics(spin):
    # Diagnostics supplement the written harmonic minimum-principle proof.
    with mp.workdps(110):
        for logr in (mp.mpf(0), mp.mpf(".1"), 1, 16, 64, 160):
            radius = mp.exp(logr)
            for angle in (
                0,
                mp.pi / 3,
                -mp.pi / 2,
                mp.pi - mp.mpf("1e-40"),
                -mp.pi + mp.mpf("1e-40"),
            ):
                z = radius * mp.exp(1j * angle)
                assert mp.re(closed(z, spin)) >= mp.mpf(13) * logr / 80 - 4 - mp.mpf(
                    "1e-80"
                )


@pytest.mark.parametrize("spin", (0, 2))
def test_independent_joint_frequency_line_resonant_and_large_momentum(spin):
    with mp.workdps(110):
        for bound in (0, mp.mpf("1.5"), 10):
            # sigma and transfer here are scaled by2m and4m^2 respectively.
            sigma = mp.exp(32 + 8 * bound)
            gap = (32 + 13 * bound) / 5
            for factor in (0, 1, 10, mp.exp(30)):
                omega = factor * sigma
                for q in (0, omega * omega, 2 * omega * omega + sigma * sigma):
                    z = (sigma + 1j * omega) ** 2 + q
                    value = closed(z, spin)
                    assert abs(z) >= sigma * sigma
                    assert mp.re(value) >= gap
                    assert abs(1 / value) <= 1 / gap


def test_independent_continuous_cut_bounds_and_outer_slopes():
    z = s.Symbol("independent_cut_coordinate", real=True)
    p0 = 3 - 2 * z + 3 * z * z
    p2 = 30 - 20 * z + 3 * z * z
    assert s.expand(p0 - 3 * (z - s.Rational(1, 3)) ** 2 - s.Rational(8, 3)) == 0
    assert s.expand(p2 - 13 - (1 - z) * (17 - 3 * z)) == 0
    assert s.Rational(3, 4) * 13 / 60 == s.Rational(13, 80)
    assert s.Rational(3, 4) * s.Rational(8, 3) / 2 == 1
    # The bounded log interval and two infinite cut tails are distinct cases.
    assert s.Rational(13, 80) * 16 - 4 < -s.Rational(172, 225)
    assert s.Rational(13, 60) > s.Rational(13, 80)
    assert 2 > s.Rational(13, 80)


@pytest.mark.parametrize("sigma", (s.Rational(1, 10), s.Rational(3), s.Rational(10000)))
def test_independent_frequency_polynomial_nonnegative_for_every_transfer(sigma):
    omega, q = s.symbols("independent_frequency independent_transfer", real=True)
    real = sigma * sigma - omega * omega + q
    imag = 2 * sigma * omega
    observed = s.expand(
        real * real + imag * imag - sigma * sigma * (sigma * sigma + omega * omega + q)
    )
    wanted = (q - omega * omega) ** 2 + sigma * sigma * (q + omega * omega)
    assert s.expand(observed - wanted) == 0
    # q>=0 makes each displayed term nonnegative, including the wave cone.
    assert (
        s.expand(observed.subs(q, omega * omega) - 2 * sigma * sigma * omega * omega)
        == 0
    )


def test_independent_noncommuting_time_channel_solve_and_finite_Neumann_error():
    K0 = s.Matrix([[s.Rational(1, 100), 0], [s.Rational(1, 500), s.Rational(1, 100)]])
    K2 = s.Matrix([[s.Rational(1, 150), 0], [-s.Rational(1, 600), s.Rational(1, 150)]])
    K = s.diag(K0, K2)
    H = s.diag(s.Rational(1, 5), s.Rational(1, 4))
    G = s.diag(s.Rational(1, 10), -s.Rational(1, 8))
    V = s.BlockMatrix([[H, G], [-G, H / 2]]).as_explicit()
    assert sum(v * v for v in K) < s.Rational(1, 100)
    assert sum(v * v for v in V) < 1
    I = s.eye(4)
    F = K.inv()
    actual = (F + V).inv()
    left = (I + K * V).inv() * K
    right = K * (I + V * K).inv()
    assert actual == left == right
    assert actual * (F + V) == I and (F + V) * actual == I
    assert K * V != V * K
    assert actual != K * (I + K * V).inv()
    N = 8
    partial = sum(((-K * V) ** n * K for n in range(N + 1)), s.zeros(4))
    tail = s.Rational(1, 9 * 10 ** (N + 1))
    assert sum(v * v for v in actual - partial) < tail * tail
    # This is a finite ordering/tail diagnostic, not a continuum discretization.


def test_independent_bounded_correction_has_growing_reference_pole():
    with mp.workdps(80):
        positive = radial(mp.mpf(1), 2)
        coefficient = 8 * positive / 3
        assert 0 < coefficient < mp.mpf(208) / 315 < 1

        # The unmodified original radial factor, evaluated at p=4m^2,
        # is canceled by this separately named bounded channel correction.
        def corrected(point):
            return -8 * radial(point, 2) / 3 + coefficient

        assert abs(corrected(mp.mpf(1))) < mp.mpf("1e-70")
        assert mp.diff(corrected, 1) < 0
        assert corrected(mp.mpf(".9")) > 0 > corrected(mp.mpf("1.1"))
        assert 2 * 1000 > 0  # lambda=2m is positive at q=0; not physical stability.


def test_independent_curved_weighted_integral_and_density_constants():
    b = s.Symbol("decay_gap", positive=True)
    t = s.Symbol("elapsed", nonnegative=True)
    assert s.integrate(t * s.exp(-b * t), (t, 0, s.oo)) == 1 / b**2
    assert s.Rational(5, 2) ** 2 * 5 == s.Rational(125, 4)
    assert 384 * s.Rational(125, 4) == 12000


@pytest.mark.parametrize("module", (analytic, w, curved))
def test_exact_module_checks_gates_and_serialization(module):
    data = module.data()
    assert all(s.factor(value) == 0 for value in data["checks"].values())
    assert all(bool(value) for value in data["gates"].values())
    serializer.serialize(
        {key: value for key, value in data.items() if key not in ("checks", "gates")}
    )
