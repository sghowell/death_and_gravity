"""Independent angular, spectral, hierarchy and whole-Q calibrations."""

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_unequal_mass_regge_residue import (
    audit,
    regge,
    source,
    waves,
)

ROWS, GATES = audit.residuals(), audit.gates()


@pytest.mark.parametrize("name", tuple(ROWS))
def test_every_exact_identity(name):
    value = ROWS[name]
    assert all(
        x == 0 for x in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    assert not s.sympify(value).has(s.Float)


@pytest.mark.parametrize("name", tuple(GATES))
def test_every_written_proof_gate(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_every_scope_guard(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def q2(value):
    if value > 10:
        return (
            2
            * mp.hyp2f1(mp.mpf("1.5"), 2, mp.mpf("3.5"), 1 / value**2)
            / (15 * value**3)
        )
    return (3 * value * value - 1) * mp.atanh(1 / value) / 2 - 3 * value / 2


ANGULAR = (
    ("light", 1, 8, 6),
    ("light", 1, 32, 16),
    ("light", 2, 80, 24),
    ("light", 1, 8, 40),
    ("light", 1, 32, 200),
    ("light", 2, 80, 500),
    ("heavy", 1, 8, 40),
    ("heavy", 1, 32, 200),
    ("heavy", 2, 80, 500),
)


@pytest.mark.parametrize("species,mass,heavy,energy", ANGULAR)
def test_full_independent_Bose_projection(species, mass, heavy, energy):
    with mp.workdps(70):
        if species == "light":
            a = mp.mpf(heavy) + (energy - 4 * mass) / 2
            b = mp.mpf(energy - 4 * mass) / 2
        else:
            a = mp.mpf(energy) / 2 - heavy
            b = mp.sqrt((energy - 4 * mass) * (energy - 4 * heavy)) / 2
        direct = mp.quad(
            lambda z: (1 / (a - b * z) + 1 / (a + b * z)) * (3 * z * z - 1) / 2,
            [-1, 0, 1],
        ) / (32 * mp.pi)
        closed = q2(a / b) / (8 * mp.pi * b)
        assert mp.almosteq(direct, closed, rel_eps=mp.mpf("1e-55"))


def q_integral(value, delta):
    eta = mp.acosh(value)
    d = mp.exp(-2 * eta)
    integrand = lambda w: 2 * (1 - w * w) ** (2 + delta) / mp.sqrt(1 - d * (1 - w * w))
    return mp.exp(-(3 + delta) * eta) * mp.quad(integrand, [0, mp.mpf(".5"), 1])


@pytest.mark.parametrize("value", ("1.0001", "1.1", "2", "10", "1e30"))
@pytest.mark.parametrize("delta", (".001", ".1", "1"))
def test_whole_Q_order_bound(value, delta):
    with mp.workdps(70):
        value, delta = mp.mpf(value), mp.mpf(delta)
        zero = q_integral(value, 0)
        assert mp.almosteq(zero, q2(value), rel_eps=mp.mpf("1e-55"))
        continued = q_integral(value, delta)
        assert mp.almosteq(
            continued, mp.legenq(2 + delta, 0, value, type=3), rel_eps=mp.mpf("1e-55")
        )
        ratio = continued / zero
        assert 0 <= 1 - ratio <= delta * (mp.acosh(value) + mp.mpf(".25"))


@pytest.mark.parametrize("heavy", (8, 32, 512))
def test_full_two_mass_parameter_moments(heavy):
    with mp.workdps(70):
        heavy = mp.mpf(heavy)
        D = lambda z: (1 - z) ** 2 + heavy * z
        light = mp.quad(
            lambda z: heavy**2 * z * z * (1 - z) ** 3 / D(z) ** 2,
            [0, mp.mpf(".1"), mp.mpf(".5"), 1],
        )
        hh = mp.quad(
            lambda z: heavy**2 * z**3 * (1 - z) ** 2 / D(z) ** 2,
            [0, mp.mpf(".1"), mp.mpf(".5"), 1],
        )
        assert mp.mpf(1) / 30 < light + hh < mp.mpf(1) / 3
        assert 0 < mp.mpf(1) / 12 - hh < 2 / (5 * heavy)


def test_original_extreme_hierarchy_without_float_roundoff():
    with mp.workdps(280):
        n = mp.mpf(int(source.HEAVY_MASS2.p)) / int(source.HEAVY_MASS2.q)
        a = 1 / n
        delta = mp.mpf(10) ** -98
        D = lambda z: z + a * (1 - z) ** 2
        nodes = [0, a, mp.sqrt(a), mp.mpf("1e-10"), mp.mpf(".25"), mp.mpf(".5"), 1]
        light = mp.quad(lambda z: z * z * (1 - z) ** 3 / D(z) ** 2, nodes)
        heavy = mp.quad(lambda z: z**3 * (1 - z) ** 2 / D(z) ** 2, nodes)
        assert 0 < mp.mpf(1) / 3 - light - heavy < 7 * delta / 5
        assert abs(heavy / (light + heavy) - mp.mpf(1) / 4) < 2 * delta


def twice_cut(energy, mass, heavy, species):
    active = mass if species == "light" else heavy
    if energy <= 4 * active:
        return mp.mpf(0)
    if species == "light":
        a = heavy + (energy - 4 * mass) / 2
        b = (energy - 4 * mass) / 2
    else:
        a = energy / 2 - heavy
        b = mp.sqrt((energy - 4 * mass) * (energy - 4 * heavy)) / 2
    partial = q2(a / b) / (8 * mp.pi * b)
    return (
        mp.sqrt(1 - 4 * active / energy)
        * (energy - 4 * active)
        / (energy - 4 * mass)
        * partial
    )


@pytest.mark.parametrize("heavy,cap", ((8, 16), (32, 64), (512, 128)))
def test_full_light_window_fraction(heavy, cap):
    with mp.workdps(60):
        n, U = mp.mpf(heavy), mp.mpf(cap)
        low = mp.quad(
            lambda t: twice_cut(t, mp.mpf(1), n, "light") / (mp.pi * t * t),
            [4, (U + 4) / 2, U],
        )
        F = lambda z: (1 - z) ** 2 + n * z
        total = mp.quad(
            lambda z: z * z * (1 - z) ** 2 / F(z) ** 2,
            [0, mp.mpf(".1"), mp.mpf(".5"), 1],
        ) / (48 * mp.pi**2)
        assert 0 < low < (U - 4) / (240 * mp.pi**2 * n**3)
        assert low / total < 6 * U / n


@pytest.mark.parametrize("heavy,cap", ((8, 64), (32, 256), (128, 1024)))
def test_whole_weighted_trajectory_error_majorant(heavy, cap):
    with mp.workdps(50):
        n, U = mp.mpf(heavy), mp.mpf(cap)
        scale = 10 * U
        total = mp.mpf(0)
        for species, active in (("light", mp.mpf(1)), ("heavy", n)):
            lower = 4 * active

            def integrand(t, lower=lower, species=species, active=active):
                if t <= lower:
                    return mp.mpf(0)
                if species == "light":
                    a = n + (t - 4) / 2
                    b = (t - 4) / 2
                else:
                    a = t / 2 - n
                    b = mp.sqrt((t - 4) * (t - 4 * n)) / 2
                loss = (
                    mp.acosh(a / b)
                    + mp.mpf(".25")
                    + mp.log((t - 4) / (t - 4 * active)) / 2
                )
                return twice_cut(t, mp.mpf(1), n, species) * loss / (mp.pi * t * scale)

            total += mp.quad(integrand, [lower, (lower + U) / 2, U])
        q = U - 4
        x = U - 4 * n
        A = 2 * U + 4 * n
        bound = (
            q * q * (mp.log(A / q) + mp.mpf(".75"))
            + x * x * (mp.log(A / x) + mp.log(U / x) / 2 + 1)
        ) / (480 * mp.pi**2 * n**3 * scale)
        assert 0 < total < bound


@pytest.mark.parametrize("species,mass,heavy,energy", ANGULAR)
def test_numeric_domain_guards(species, mass, heavy, energy):
    assert waves.require_domain(species, energy, mass, heavy) == (energy, mass, heavy)
    with pytest.raises(ValueError):
        waves.require_domain(
            species, 4 * (mass if species == "light" else heavy), mass, heavy
        )


@pytest.mark.parametrize(
    "args",
    (
        (32, 100, 0, 0, 1, 8),
        (64, 100, -1, 0, 1, 8),
        (64, 100, 0, -1, 1, 8),
        (64, 100, 0.1, 0, 1, 8),
        (64, 100, 0, 0, 2, 8),
        (64, 0, 0, 0, 1, 8),
    ),
)
def test_conditional_budget_rejects_invalid_inputs(args):
    with pytest.raises((ValueError, TypeError)):
        regge.require_budget(*args)


def test_conditional_budget_is_not_a_matched_number():
    assert regge.require_budget(64, 100, 0, 0, 1, 8) == (64, 100, 0, 0, 1, 8)
    budget = regge.closed_budget()
    assert budget["total"].has(regge.ARC, regge.EPS, regge.SCALE, regge.CAP)
    assert len(audit.frontier()) == 9
    assert len(audit.matching()) == 154
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.qualifications()) == 6


def test_regular_spin2_at_formal_heavy_resonance():
    # The unprojected H pole is P0 and is not declared physically absent.
    assert (
        waves.even_spin2(
            "light", s.Integer(8), s.Integer(1), s.Integer(8), s.Integer(1)
        ).is_finite
        is True
    )
