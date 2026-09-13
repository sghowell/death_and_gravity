"""Independent regulated Gaussian, diagram, logarithm and first-loop scope checks."""

import itertools

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_heavy_scalar_one_loop import audit, gaussian, width
from p8_vacuum_affine_heavy_scalar_one_loop import self_energy as se
from p8_vacuum_affine_heavy_scalar_tree_matching import audit as previous
from p8_vacuum_affine_heavy_scalar_tree_matching import model


def mpq(value):
    value = s.Rational(value)
    return mp.mpf(int(value.p)) / int(value.q)


@pytest.mark.parametrize("name,value", list(audit.residuals().items()))
def test_every_exact_scalar_entry(name, value):
    assert s.cancel(value) == 0, name


@pytest.mark.parametrize("name,value", list(audit.gates().items()))
def test_every_proof_gate(name, value):
    assert value is True, name


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_every_unsupported_input_rejected(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


@pytest.mark.parametrize("sites,mass2", [(2, 9), (3, 9), (3, 100)])
def test_independent_full_nonscalar_gaussian_hessian_and_schur(sites, mass2):
    phi = s.Matrix(s.symbols("p:" + str(sites)))
    H = s.Matrix(s.symbols("h:" + str(sites)))
    lap = s.zeros(sites)
    for i in range(sites - 1):
        lap[i, i] += 1
        lap[i + 1, i + 1] += 1
        lap[i, i + 1] -= 1
        lap[i + 1, i] -= 1
    KH = mass2 * s.eye(sites) + lap
    KL = s.eye(sites) + lap
    g = s.Rational(2, 7)
    D = mass2 - 2
    C = -g * g * (s.Rational(3, D) - s.Rational(2, D * D))
    f = phi.applyfunc(lambda y: y * y)
    G = KH.inv()
    full = (
        (phi.T * KL * phi)[0] / 2
        + (H.T * KH * H)[0] / 2
        - g * (H.T * f)[0] / 2
        - C * sum(y**4 for y in phi) / 24
    )
    center = g * G * f / 2
    eliminated = full.subs(dict(zip(H, center)), simultaneous=True)
    candidate = (
        (phi.T * KL * phi)[0] / 2
        - C * sum(y**4 for y in phi) / 24
        - g * g * (f.T * G * f)[0] / 8
    )
    assert s.expand(eliminated - candidate) == 0
    expected = (
        KL
        - C * s.diag(*f) / 2
        - g * g * s.diag(*(G * f)) / 2
        - g * g * s.diag(*phi) * G * s.diag(*phi)
    )
    assert (s.hessian(eliminated, tuple(phi)) - expected).applyfunc(
        s.cancel
    ) == s.zeros(sites)
    blocks = s.hessian(full, tuple(phi) + tuple(H)).subs(
        dict(zip(H, center)), simultaneous=True
    )
    assert (
        blocks[:sites, :sites]
        - blocks[:sites, sites:] * G * blocks[sites:, :sites]
        - expected
    ).applyfunc(s.cancel) == s.zeros(sites)
    assert (G * s.ones(sites, 1) - s.ones(sites, 1) / mass2) == s.zeros(sites, 1)
    vector = s.Matrix([(-1) ** i * (i + 1) for i in range(sites)])
    assert (vector.T * (s.eye(sites) / mass2 - G) * vector)[0] >= 0


@pytest.mark.parametrize("phi", [-2, s.Rational(1, 2), 3])
def test_actual_one_dimensional_gaussian_integral_with_source(phi):
    with mp.workdps(60):
        n = mp.mpf(9)
        g = mp.mpf(2) / 7
        p = mpq(phi)
        center = g * p * p / (2 * n)
        direct = mp.quad(
            lambda h: mp.exp(-n * h * h / 2 + g * h * p * p / 2),
            [-mp.inf, center, mp.inf],
        )
        expected = mp.sqrt(2 * mp.pi / n) * mp.exp(g * g * p**4 / (8 * n))
        assert abs(direct / expected - 1) < mp.mpf("1e-55")


def test_independently_enumerated_mixed_and_heavy_tadpole_contractions():
    slots = [(vertex, leg) for vertex in range(2) for leg in range(2)]
    routes = list(itertools.permutations(slots, 2))
    mixed = sum(a[0] != b[0] for a, b in routes)
    reducible = sum(a[0] == b[0] for a, b in routes)
    assert mixed == 8 and reducible == 4
    vertex_and_expansion_denominator = 2**2 * s.factorial(2)
    assert s.Rational(mixed, vertex_and_expansion_denominator) == 1
    assert s.Rational(reducible, vertex_and_expansion_denominator) == s.Rational(1, 2)
    heavy_external_routes = s.factorial(2)
    light_internal_pairings = s.factorial(2)
    assert (
        heavy_external_routes
        * light_internal_pairings
        / vertex_and_expansion_denominator
        == s.Rational(1, 2)
    )
    contact_routes = len(list(itertools.permutations(range(4), 2)))
    assert s.Rational(contact_routes, s.factorial(4)) == s.Rational(1, 2)


def test_full_trace_log_second_and_fourth_field_jets():
    marker = s.Symbol("background_marker")
    K = s.Matrix([[3, -1], [-1, 4]])
    W = s.Matrix([[2, s.Rational(1, 3)], [s.Rational(1, 3), 5]])
    expression = s.log((K + marker * marker * W).det() / K.det()) / 2
    first = s.diff(expression, marker, 2).subs(marker, 0) / 2
    second = s.diff(expression, marker, 4).subs(marker, 0) / 24
    assert s.cancel(first - s.trace(K.inv() * W) / 2) == 0
    assert s.cancel(second + s.trace(K.inv() * W * K.inv() * W) / 4) == 0
    # A light-independent heavy determinant cannot reproduce this nonzero jet.
    assert first != 0 and second != 0


def test_heavy_one_point_counterterm_and_full_source_square():
    g, m2, I0, p, j = s.symbols("g m2 I0 p j")
    shifted = -((g * p * p / 2 - j) ** 2) / (2 * m2)
    induced = s.diff(shifted, j, p, p).subs({p: 0, j: 0})
    assert s.cancel(induced - g / m2) == 0
    onepoint_ct = g * I0 / 2
    reduced_tadpole = -g * g * I0 / (2 * m2)
    assert s.cancel(reduced_tadpole + induced * onepoint_ct) == 0
    assert reduced_tadpole != 0


def test_dimensional_tadpole_and_mixed_bubble_signs():
    g, C, A0, B0, d, P = s.symbols("g C A0 B0 d P")
    bubble = (s.I * g) ** 2 * s.I**2 * (s.I * B0 / (16 * s.pi**2)) / s.I
    tadpole = (s.I * C / 2) * s.I * (s.I * A0 / (16 * s.pi**2)) / s.I
    assert s.cancel(bubble - g * g * B0 / (16 * s.pi**2)) == 0
    assert s.cancel(tadpole + C * A0 / (32 * s.pi**2)) == 0
    first_insertion = (s.I / d) * (s.I * P) * (s.I / d)
    assert s.cancel(first_insertion + s.I * P / d**2) == 0
    assert s.diff(s.I / (d + P), P).subs(P, 0) == -s.I / d**2


@pytest.mark.parametrize("mass2", [9, 100, 10000])
@pytest.mark.parametrize("shift", [s.Rational(-3), 1 + 2 * s.I, s.Rational(1, 20)])
def test_independent_complex_subtracted_bubble_quadrature(mass2, shift):
    with mp.workdps(80):
        real, imag = shift.as_real_imag()
        delta = mp.mpc(mpq(real), mpq(imag))
        if shift == s.Rational(1, 20):
            delta *= mass2
        F = lambda y: (1 - y) ** 2 + y * mass2
        aa = lambda y: y * (1 - y) / F(y)
        B = lambda v: (
            -mp.quad(lambda y: mp.log(y * mass2 + 1 - y - y * (1 - y) * v), [0, 1])
        )
        direct = B(1 + delta) - B(1) - delta * mp.diff(B, 1)
        stable = mp.quad(lambda y: -mp.log1p(-aa(y) * delta) - aa(y) * delta, [0, 1])
        majorant = abs(delta) ** 2 / (6 * mass2**2 * (1 - abs(delta) / mass2))
        assert abs(direct - stable) < mp.mpf("1e-68")
        assert 0 < abs(stable) < majorant
        slope = mp.diff(B, 1)
        assert 0 < slope < mp.mpf(1) / (2 * mass2)
        assert mp.diff(B, 1, 2) > 0


@pytest.mark.parametrize("M", [3, 10, 100])
def test_physical_and_pseudothresholds_independently(M):
    x = s.Symbol("x")
    Q = lambda energy: x * M * M + 1 - x - x * (1 - x) * energy
    assert s.expand(Q((M + 1) ** 2) - ((M + 1) * x - 1) ** 2) == 0
    assert s.expand(Q((M - 1) ** 2) - ((M - 1) * x + 1) ** 2) == 0
    assert Q((M + 1) ** 2).subs(x, s.Rational(1, M + 1)) == 0
    for y in [0, s.Rational(1, 7), s.Rational(1, 2), 1]:
        assert Q((M - 1) ** 2).subs(x, y) > 0


@pytest.mark.parametrize("mass,energy", [(1, 5), (1, 16), (3, 25), (10, 144)])
def test_independent_log_branch_absorptive_integral(mass, energy):
    with mp.workdps(80):
        M = mp.mpf(mass)
        S = mp.mpf(energy)
        disc = (S - M * M - 1) ** 2 - 4 * M * M
        root0 = (S + 1 - M * M - mp.sqrt(disc)) / (2 * S)
        root1 = (S + 1 - M * M + mp.sqrt(disc)) / (2 * S)
        epsilon = mp.mpf("1e-60")
        b0 = -mp.quad(
            lambda y: mp.log(y * M * M + 1 - y - y * (1 - y) * S - mp.j * epsilon),
            [0, root0, root1, 1],
        )
        expected = mp.pi * mp.sqrt(disc) / S
        assert abs(mp.im(b0) - expected) < mp.mpf("1e-50")
        assert mp.im(b0) > 0
        wrong = -b0
        assert mp.im(wrong) < 0


def test_decay_normalization_and_formal_second_sheet_sign():
    g, M, beta = s.symbols("g M beta", positive=True)
    complete_phase = beta / (8 * s.pi)
    rate = g * g * complete_phase / (2 * M * 2)
    absorptive = g * g * beta / (32 * s.pi)
    assert s.cancel(rate - absorptive / M) == 0
    # The same-sign first-sheet upper boundary continues to the outgoing sheet.
    pole_shift = -s.I * absorptive
    assert s.im(pole_shift) < 0
    assert s.im(-pole_shift) > 0


@pytest.mark.parametrize("direction", [1, -1, (1 + 2 * s.I) / 3])
def test_actual_hierarchy_first_loop_disk_and_mass_coefficient(direction):
    with mp.workdps(500):
        n = mpq(model.MASS2)
        g2 = mpq(model.G2)
        C = mpq(model.CONTACT)
        real, imag = (
            direction.as_real_imag()
            if isinstance(direction, s.Basic)
            else (s.Integer(direction), s.Integer(0))
        )
        delta = (mpq(real) + mp.j * mpq(imag)) * mpq(se.RADIUS)
        F = lambda y: (1 - y) ** 2 + y * n
        alpha = lambda y: y * (1 - y) / F(y)
        partitions = [0, 1 / n, 1 / mp.sqrt(n), mp.mpf("0.5"), 1]
        rem = (
            g2
            / (16 * mp.pi**2)
            * mp.quad(
                lambda y: -mp.log1p(-alpha(y) * delta) - alpha(y) * delta, partitions
            )
        )
        relative = rem / delta
        assert 0 < abs(relative) < mpq(se.EPSILON_BOUND)
        assert abs(1 / (1 + relative) - 1) < mp.mpf("1e-209")
        slope = g2 / (16 * mp.pi**2) * mp.quad(alpha, partitions)
        mass = -C / (32 * mp.pi**2) - g2 / (16 * mp.pi**2) * mp.quad(
            lambda y: mp.log(F(y)), partitions
        )
        assert 0 < slope < mp.mpf("1e-207")
        assert abs(mass) < mp.mpf("1e-7")
        # On-shell finite constants are necessary; MSbar alone is not on shell.
        assert mass != 0 and slope != 0


def test_actual_named_rational_margins_without_floating_arithmetic():
    assert 0 < se.RADIUS / model.MASS2 < se.R0 < 1
    assert 0 < se.EPSILON_BOUND / (1 - se.EPSILON_BOUND) < s.Rational(1, 10**209)
    assert 0 < se.SLOPE_BOUND < s.Rational(1, 10**207)
    assert 0 < se.MASS_BOUND < s.Rational(1, 10**7)
    assert (
        s.Rational(1, 10**208)
        < width.WIDTH_LOWER
        < width.WIDTH_UPPER
        < s.Rational(1, 10**207)
    )
    assert 1 - 4 / model.MASS2 > s.Rational(999, 1000) ** 2


@pytest.mark.parametrize(
    "value",
    [0, 1, s.I, 1 + se.RADIUS, 1 + s.I * se.RADIUS, 1 + (3 + 4 * s.I) * se.RADIUS / 5],
)
def test_exact_complex_domain_accepts_boundary_and_anchor(value):
    result = audit.require_observable("V2S-T1", value, 1)
    assert result == ("V2S-T1", s.sympify(value), 1)


@pytest.mark.parametrize(
    "stage",
    [
        "exact_finite_regulator_heavy_integration",
        "complete_nonlocal_one_light_loop_generator",
        "first_light_on_shell_self_energy",
        "one_loop_truncated_complex_disk_comparison",
        "first_heavy_absorptive_decay_coefficient",
    ],
)
def test_only_positive_scopes_are_accepted(stage):
    assert audit.require_stage(stage) == stage


def test_complete_frozen_frontier_and_first_loop_boundary():
    assert audit.frontier() == previous.frontier()
    assert audit.matching()[:-1] == previous.matching()
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 90
    assert len(audit.residuals()) == 58 and audit.scalar_entry_count() == 58
    assert len(audit.gates()) == 25 and len(audit.controls()) == 9
    assert audit.rejected_inputs() == 228
    assert (
        "not an evaluated renormalized four-point"
        in gaussian.data()["complete_formal_one_light_loop"]
    )
    assert (
        "not a proof of exact physical mass or LSZ"
        in se.data()["subtracted_coefficient"]
    )
    assert "one-loop-truncated" in se.data()["actual_propagator_comparison"]
    assert "second sheet" in width.data()["outgoing_pole_scope"]
    assert "metric-dependent" in gaussian.data()["mixed_loop_warning"]
