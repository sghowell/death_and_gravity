"""Independent whole endpoint tensor, pair-cut, bound and retained-scope tests."""

import mpmath as mp
import pytest
import sympy as s
from p8_affine import verify as base
from p8_vacuum_affine_complete_matter_graviton_endpoint import (
    audit,
    cuts,
    endpoint,
    forward,
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
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_every_unsupported_input(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("D", (4, 5, 6))
@pytest.mark.parametrize("choice", (1, 2, 3))
def test_independent_whole_covariant_OS_metric_vertex(D, choice):
    eta = s.diag(1, *([-1] * (D - 1)))
    P = s.Matrix([choice + 3, 0, 1, *([0] * (D - 3))])
    q = s.Matrix([0, choice, *([0] * (D - 2))])
    p, r = P - q / 2, P + q / 2
    dot = lambda a, b: (a.T * eta * b)[0]
    mu = dot(p, p)
    Pi, prime = s.Rational(2, 7), s.Rational(3, 11)
    kinetic = p * r.T + r * p.T - eta * dot(p, r)
    actual = -prime * kinetic + (Pi - mu * prime) * eta
    Q = q * q.T - eta * dot(q, q)
    assert actual == -prime * (2 * P * P.T - Q / 2) + Pi * eta
    assert actual * eta * q == Pi * q
    assert (2 * P * P.T - Q / 2) * eta * q == s.zeros(D, 1)
    assert actual != -prime * (2 * P * P.T - Q / 2)


@pytest.mark.parametrize("D", (4, 5, 6))
@pytest.mark.parametrize("z", (s.Rational(1, 4), s.Rational(2, 3)))
@pytest.mark.parametrize("v", (s.Rational(1, 3), s.Rational(3, 4)))
def test_independent_entire_raw_and_even_triangle_tensor(D, z, v):
    eta = s.diag(1, *([-1] * (D - 1)))
    P = s.Matrix([3, 0, 1, *([0] * (D - 3))])
    q = s.Matrix([0, 2, *([0] * (D - 2))])
    loop = s.Matrix(s.symbols("test_loop0:" + str(D), real=True))
    active = s.Integer(2)
    alpha = (1 - z) * v
    a = z * P + (alpha - 1) * q / 2 - loop
    b = z * P + (alpha + 1) * q / 2 - loop
    raw = a * b.T + b * a.T - eta * ((a.T * eta * b)[0] - active)

    def avg(value):
        poly = s.Poly(s.expand(value), *loop)
        return poly.coeff_monomial(1) + sum(
            eta[i, i] * poly.coeff_monomial(loop[i] ** 2) / D for i in range(D)
        )

    actual = raw.applyfunc(avg)
    qq = (alpha * alpha - 1) / 2
    scalar = (
        (s.Rational(2, D) - 1)
        + active
        - z * z * (P.T * eta * P)[0]
        - (alpha * alpha - 1) * (q.T * eta * q)[0] / 4
    )
    expected = (
        2 * z * z * P * P.T
        + z * alpha * (P * q.T + q * P.T)
        + qq * q * q.T
        + eta * scalar
    )
    assert actual == expected
    assert z * alpha * (P * q.T + q * P.T) != s.zeros(D)
    reflected_a = z * P + (-alpha - 1) * q / 2 - loop
    reflected_b = z * P + (-alpha + 1) * q / 2 - loop
    reflected_raw = (
        reflected_a * reflected_b.T
        + reflected_b * reflected_a.T
        - eta * ((reflected_a.T * eta * reflected_b)[0] - active)
    )
    reflected = reflected_raw.applyfunc(avg)
    assert (actual + reflected) / 2 == 2 * z * z * P * P.T + qq * q * q.T + eta * scalar


@pytest.mark.parametrize(
    "mu,n,energy", ((1, 9, 6), (2, 17, 10), (1, 9, 50), (2, 17, 90))
)
@pytest.mark.parametrize("fraction", (s.S.Zero, s.Rational(1, 2), s.S.One))
def test_independent_exact_parameter_cut_root_and_full_gap(mu, n, energy, fraction):
    mu, n, energy = map(s.Integer, (mu, n, energy))
    active, spectator = (mu, n) if energy < 4 * n else (n, mu)
    v = fraction * s.sqrt(1 - 4 * active / energy)
    A = mu - energy * (1 - v * v) / 4
    B = spectator - active - mu + energy * (1 - v * v) / 2
    C = active - energy * (1 - v * v) / 4
    disc = B * B - 4 * A * C
    z = -2 * C / (B + s.sqrt(disc))
    assert A <= 0 and C <= 0 and B > 0 and disc >= n * n
    assert 0 <= z < 1
    assert s.simplify(A * z * z + B * z + C) == 0
    assert s.simplify(B + 2 * A * z - s.sqrt(disc)) == 0
    literal = endpoint.denominator(active, spectator, energy, z, v, mu)
    assert s.simplify(literal) == 0


@pytest.mark.parametrize(
    "mu,n,energy", ((1, 9, 6), (2, 17, 10), (1, 9, 50), (2, 17, 90))
)
def test_independent_complete_parameter_discontinuity_vs_pair_tensor_sew(mu, n, energy):
    with mp.workdps(60):
        mu, n, energy = map(mp.mpf, (mu, n, energy))
        active, spectator = (mu, n) if energy < 4 * n else (n, mu)
        beta = mp.sqrt(1 - 4 * active / energy)

        def density(v):
            A = mu - energy * (1 - v * v) / 4
            B = spectator - active - mu + energy * (1 - v * v) / 2
            C = active - energy * (1 - v * v) / 4
            root = mp.sqrt(B * B - 4 * A * C)
            z = -2 * C / (B + root)
            return (
                (1 - z) * z * z / root,
                (1 - z) * ((1 - z) ** 2 * v * v - 1) / (2 * root),
            )

        actual = [
            mp.quad(lambda v, index=i: density(v)[index], [0, beta]) / (16 * mp.pi)
            for i in (0, 1)
        ]
        external = energy - 4 * mu
        internal = energy - 4 * active
        d = n + external / 2 if active == mu else energy / 2 - n
        angular = mp.sqrt(external * internal) / 2
        tree = lambda z: 1 / (d - angular * z) + 1 / (d + angular * z)
        a0 = mp.quad(tree, [-1, 0, 1]) / 2
        a2 = 5 * mp.quad(lambda z: tree(z) * (3 * z * z - 1) / 2, [-1, 0, 1]) / 2
        expected = [
            beta * internal * a2 / (160 * mp.pi * external),
            -beta
            * ((energy + 2 * active) * a0 + internal * a2 / 10)
            / (96 * mp.pi * energy),
        ]
        assert max(abs(a - b) for a, b in zip(actual, expected, strict=True)) < mp.mpf(
            "1e-45"
        )
        assert all(
            abs(a + b) > mp.mpf("1e-10") for a, b in zip(actual, expected, strict=True)
        )


@pytest.mark.parametrize("mu", (s.Integer(1), s.Rational(3, 2)))
@pytest.mark.parametrize("ratio", (s.Integer(4), s.Integer(9)))
@pytest.mark.parametrize("z", (s.S.Zero, s.Rational(1, 3), s.S.One))
@pytest.mark.parametrize("v", (s.S.Zero, s.Rational(1, 2), s.S.One))
def test_independent_entire_triangle_density_majorant(mu, ratio, z, v):
    n = mu * ratio
    F_light = mu * (1 - z) ** 2 + n * z
    F_heavy = n * (1 - z) + mu * z * z
    actual = forward.triangle_b20_integrand(mu, n, z, v)
    upper = 3 * ((1 - z) ** 3 / F_light**2 + (1 - z) ** 3 / F_heavy**2)
    assert 0 <= actual <= upper
    A = (1 - z) ** 2 * (1 - v * v) / 4
    for F in (F_light, F_heavy):
        assert F - 2 * mu * A >= F / 2
        assert F + 2 * mu * A <= 3 * F / 2
    assert F_light >= mu * (1 - z) + n * z / 2


@pytest.mark.parametrize("mu", (s.Integer(1), s.Rational(3, 2)))
@pytest.mark.parametrize("x", (s.S.Zero, s.Rational(1, 4), s.Rational(1, 2), s.S.One))
def test_independent_full_bubble_pointwise_majorants(mu, x):
    y = x * (1 - x)
    assert 0 <= y <= s.Rational(1, 4)
    assert -s.log(1 - 2 * y) <= 4 * y
    assert y * y / (mu - 2 * mu * y) <= 2 * y * y / mu
    assert y**3 / (mu - 2 * mu * y) ** 2 <= 4 * y**3 / mu**2


@pytest.mark.parametrize("energy", (40, 50, 70))
def test_entire_above_heavy_cut_retains_both_open_channels(energy):
    total = cuts.total_cut_above_heavy(energy, 1, 9, 1, 0)
    light = cuts.formal_cut("light", energy, 1, 9, 1, 0)
    heavy = cuts.formal_cut("heavy", energy, 1, 9, 1, 0)
    assert all(
        s.simplify(a - b - c) == 0 for a, b, c in zip(total, light, heavy, strict=True)
    )
    assert all(v != 0 for v in light)
    assert any(s.simplify(a - b) != 0 for a, b in zip(total, heavy, strict=True))


@pytest.mark.parametrize("species,threshold", (("light", 4), ("heavy", 36)))
def test_formal_pair_threshold_is_continuous_zero_not_stable_atom_claim(
    species, threshold
):
    assert cuts.formal_cut(species, threshold, 1, 9, 1, 0) == (0, 0)
    assert cuts.exchange_moments(9, 0, 1) == (s.Rational(2, 9), 0)
    assert "exact stable-H" in cuts.data()["physical_domain_scope"]


@pytest.mark.parametrize("coefficient", (1, 2, 3))
def test_independent_complete_crossed_endpoint_with_polynomial_F1_F2(coefficient):
    a, v, mu = s.symbols("channel crossing_v mass", positive=True)
    f1 = lambda x: coefficient * x + x * x / 7 + x**3 / 11
    f2 = lambda x: coefficient + x / 5 + x * x / 13
    energy, crossed = 2 * mu + v, 2 * mu - v

    def channel_contribution(x, y, z):
        N = 2 * mu * mu - 2 * mu * x - y * z
        return -(2 * N + mu * x + x * x / 2) * (f1(x) / x) - (x + 2 * mu) * f2(x)

    actual = channel_contribution(energy, 0, crossed) + channel_contribution(
        crossed, 0, energy
    )
    # Exact t0 continuation of the remaining channel, including its v dependence.
    actual += -2 * (2 * mu * mu - energy * crossed) * coefficient - 2 * mu * f2(0)
    actual = s.diff(actual, v, 2).subs(v, 0) / 2
    expected = s.diff(f1(a), a).subs(a, 2 * mu) - f1(2 * mu) / mu - 2 * coefficient
    expected -= s.diff((a + 2 * mu) * f2(a), a, 2).subs(a, 2 * mu)
    assert s.factor(actual - expected) == 0


@pytest.mark.parametrize("mu,n", ((1, 9), (2, 17), (s.Rational(3, 2), 13)))
def test_constant_curvature_anchor_drops_but_RH_coefficient_does_not(mu, n):
    mu, n = map(s.sympify, (mu, n))
    v, h, c = s.symbols("crossing_v unmatched_h unmatched_constant", real=True)
    energy, crossed = 2 * mu + v, 2 * mu - v
    amplitude = (
        -(energy + 2 * mu) * (c + h / (n - energy))
        - (crossed + 2 * mu) * (c + h / (n - crossed))
        - 2 * mu * (c + h / n)
    )
    actual = s.factor(s.diff(amplitude, v, 2).subs(v, 0) / 2)
    expected = forward.mixing_anchor_b20(h, mu, n, 1)
    assert s.factor(actual - expected) == 0
    assert not actual.has(c)
    assert actual.has(h)
    assert s.diff(actual, h) != 0


def test_entire_F2_retains_independent_unmatched_curved_coordinates():
    value = endpoint.full_f2()
    assert endpoint.ELL in value.free_symbols
    assert endpoint.H in value.free_symbols
    assert s.diff(value, endpoint.ELL) == 1
    assert s.diff(value, endpoint.H) == 1 / (source.N - source.T)
    assert (
        "h is not chosen" in audit.observable()["not_established"]
        or "RH" in audit.observable()["not_established"]
    )


@pytest.mark.parametrize("bad", (True, False, 1.0, s.Float(1), None, "1", s.I, s.oo))
def test_exact_positive_mass_domains_reject_inexact_or_unspecified(bad):
    with pytest.raises((TypeError, ValueError)):
        forward.require_bound_domain(bad, 9)
    with pytest.raises((TypeError, ValueError)):
        cuts.physical_cut("light", 6, bad, 9)


@pytest.mark.parametrize("n", (0, -1, 1, 3, s.Rational(7, 2)))
def test_whole_bound_domain_rejects_heavy_below_required_hierarchy(n):
    with pytest.raises((TypeError, ValueError)):
        forward.require_bound_domain(1, n)


@pytest.mark.parametrize("bad", (True, False, 1.0, s.Float(1), None, "1", s.I, s.oo))
def test_quartic_bound_requires_exact_rational_and_its_own_size_premise(bad):
    with pytest.raises((TypeError, ValueError)):
        forward.known_endpoint_bound(1, 9, 1, bad, 1)
    with pytest.raises((TypeError, ValueError)):
        forward.known_endpoint_bound(1, 9, 1, 1, 1)


@pytest.mark.parametrize(
    "species,energy",
    (("light", 4), ("light", 9), ("light", 10), ("heavy", 36), ("heavy", 20)),
)
def test_endpoint_cut_windows_reject_closed_threshold_or_resonance(species, energy):
    with pytest.raises((TypeError, ValueError)):
        cuts.physical_cut(species, energy, 1, 9, 1, 0)


@pytest.mark.parametrize("bad", ("graviton", "Phi", "", 1, True, None))
def test_unsupported_pair_species_rejected(bad):
    with pytest.raises((TypeError, ValueError)):
        cuts.formal_cut(bad, 6, 1, 9)
    with pytest.raises((TypeError, ValueError)):
        cuts.physical_cut(bad, 6, 1, 9)


def test_actual_known_bound_preserves_all_exact_parameter_margins():
    result = forward.known_endpoint_bound(
        1, source.HEAVY_MASS2, source.CUBIC, source.CONTACT, source.KAPPA
    )
    assert (
        s.factor(
            result
            - 3 * source.CUBIC**2 / (4 * s.pi**2 * source.KAPPA * source.HEAVY_MASS2)
        )
        == 0
    )
    assert 3 * source.CUBIC**2 / (16 * source.KAPPA * source.HEAVY_MASS2) < s.Rational(
        1, 10**1005
    )
    assert all(v > 0 for v in forward.data()["whole_exact_positive_margins"].values())


def test_explicit_upper_Feynman_boundary_includes_heavy_denominator():
    value = endpoint.feynman_boundary_f2(2, mass=1, heavy=9, cubic=1, quartic=0)
    assert isinstance(value, s.Limit)
    assert value.args[3] == s.Symbol("+")
    assert value.args[0].has(s.I)
    assert value.args[0].has(endpoint.H)
    assert not value.has(s.Float)


def test_integer_algebra_apis_preserve_exact_arithmetic():
    values = (
        endpoint.denominator(1, 9, 2, 0, 0, 1),
        endpoint.raw_f2(1, 9, 2, 0, 0, 1),
        endpoint.triangle_OS_integrand(2, 0, 0, 1, 9),
        endpoint.bubble_remainder(2, 1),
        endpoint.f1(2, 1, 9, 1),
        endpoint.full_f2(2, 1, 9, 1, 0),
        *cuts.exchange_moments(10, 1, 1),
        *cuts.physical_cut("light", 6, 1, 9, 1, 0),
        *cuts.total_cut_above_heavy(50, 1, 9, 1, 0),
        forward.triangle_b20_integrand(1, 9, 0, 0),
        *forward.bubble_center_jets(1),
        forward.bubble_b20(1, 9, 1, 0, 1),
        forward.mixing_anchor_b20(1, 1, 9, 1),
        forward.known_endpoint_bound(1, 9, 1, 0, 1),
    )
    assert all(isinstance(v, s.Basic) and not v.has(s.Float) for v in values)


def test_parent_cache_warmth_cannot_change_the_new_packet():
    before = dict(audit.residuals())
    audit.previous.packets()
    assert audit.residuals() == before
    assert source.data()["checks"] is not source.previous.data()["checks"]
    base.certify_residuals(before)


def test_entire_original_scope_and_historical_rejection_retained():
    assert len(audit.frontier()) == 9
    assert len(audit.matching()) == 146
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.qualifications()) == 6
    assert audit.validate_scope(audit.frontier(), audit.matching())
    assert any(r["status"].startswith("REJECTED_") for r in audit.matching())
    assert "Regge" in audit.observable()["not_established"]
    assert "P8" in audit.observable()["not_established"]
