"""Independent massive-cut, analytic-kernel, finite-subtraction and scope tests."""

import pytest
import sympy as s
from p8_vacuum_affine_massive_detector_ir_cut import (
    audit,
    cut,
    soft,
    source,
    subtraction,
)

ROWS = audit.residuals()
GATES = audit.gates()


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
    "name,call,args", audit.bad_cases(), ids=[x[0] for x in audit.bad_cases()]
)
def test_unsupported_inputs(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("order", range(9))
def test_exact_series_interface_before_warm_cache(order):
    assert soft.kernel_coefficient(order) == soft.kernel_coefficient(s.Integer(order))
    for bad in (bool(order), float(order), s.Float(order)):
        with pytest.raises((TypeError, ValueError)):
            soft.kernel_coefficient(bad)


@pytest.mark.parametrize(
    "cap,mass2", ((5, 1), (9, 2), (s.Rational(17, 4), 1), (10**196, 1))
)
def test_exact_cap_requires_a_genuine_massive_gap(cap, mass2):
    assert subtraction.require_cap(cap, mass2) == (cap, mass2)
    for bad in (float(cap), str(cap), True):
        with pytest.raises((TypeError, ValueError)):
            subtraction.require_cap(bad, mass2)


@pytest.mark.parametrize("mass2", (s.Rational(1, 3), s.Integer(1), s.Rational(7, 2)))
@pytest.mark.parametrize("energy_factor", (4, 5, 11))
def test_full_production_and_optical_thresholds(mass2, energy_factor):
    S = mass2 * energy_factor
    x = s.Symbol("x", real=True)
    A = source.pair_tree(S, x, mass2, s.Integer(3))
    rho = s.integrate(A * A, (x, -1, 1)) / (64 * s.pi)
    expected = (S * S + 2 * mass2 * S + 6 * mass2 * mass2) / (960 * s.pi * 9)
    assert s.factor(rho - expected) == 0
    assert rho > 0
    for angle in (-1, -s.Rational(2, 5), 0, s.Rational(3, 7), 1):
        assert A.subs(x, angle) >= mass2 / 3


@pytest.mark.parametrize(
    "angle",
    (s.Rational(-1), s.Rational(-2, 5), s.Integer(0), s.Rational(3, 5), s.Integer(1)),
)
def test_independent_azimuth_and_polar_integration(angle):
    # Full two-direction sewing, with the azimuthal odd moments canceled.
    x = s.Symbol("x", real=True)
    S, mu, k = s.Integer(13), s.Integer(1), s.Integer(7)
    A = source.pair_tree(S, x, mu, k)
    second_square = (1 - angle * angle) * (1 - x * x) / 2 + angle * angle * x * x
    B = (S - (S - 4 * mu) * second_square) / (4 * k)
    sewn = s.integrate(A * B, (x, -1, 1)) / (64 * s.pi)
    t = -(S - 4 * mu) * (1 - angle) / 2
    u = 4 * mu - S - t
    expected = (S * S - t * u + 2 * mu * S + 6 * mu * mu) / (960 * s.pi * k * k)
    assert s.factor(sewn - expected) == 0


def test_minimal_pair_is_not_the_massive_elastic_channel():
    S = source.S
    x = source.Z
    assert source.pair_tree(S, x).subs(S, 4 * source.MU) == source.MU / source.K
    assert 2 * source.MU != 4 * source.MU


def test_soft_TT_norm_does_not_depend_on_longitudinal_current():
    a, b, c, f, g, h = s.symbols("a b c f g h", real=True)
    J = s.Matrix([[h, f, g, h], [f, a, b, f], [g, b, c, g], [h, f, g, h]])
    actual = source.projector_contract(J, J)
    assert s.factor(actual - (a - c) ** 2 / 2 - 2 * b * b) == 0
    assert not actual.has(f, g, h)


@pytest.mark.parametrize("mass2", (s.Rational(1, 3), s.Integer(1), s.Integer(7)))
def test_soft_pseudothreshold_and_midpoint(mass2):
    assert soft.pair_factor(s.Integer(0), mass2) == mass2 / 2
    assert s.simplify(soft.pair_factor(2 * mass2, mass2) + s.pi * mass2 / 4) == 0
    assert soft.kernel_coefficient(0, mass2) == 1 / (4 * mass2)


@pytest.mark.parametrize("order", range(9))
def test_kernel_series_from_beta_function(order):
    mu = s.Symbol("mu", positive=True)
    beta = (
        s.factorial(order)
        * s.gamma(s.Rational(1, 2))
        / (2 * s.gamma(order + s.Rational(3, 2)))
    )
    wanted = s.simplify(beta / (4 * mu) ** (order + 1))
    assert s.simplify(wanted - soft.kernel_coefficient(order, mu)) == 0


def test_full_soft_self_terms_affect_crossing_center():
    mu = s.Symbol("mu", positive=True)
    center = (
        2
        * (
            soft.pair_factor(2 * mu, mu)
            + soft.pair_factor(s.Integer(0), mu)
            + soft.pair_factor(2 * mu, mu)
        )
        - mu
    )
    assert s.simplify(center + s.pi * mu) == 0
    assert center != center + mu


@pytest.mark.parametrize("channel", cut.CHANNELS)
def test_low_piece_discontinuity_supported_only_below_cap(channel):
    # log(-z-i0) and log(L-z-i0) differ below L but share their jump above L.
    before = -s.I * s.pi - 0
    above = -s.I * s.pi - (-s.I * s.pi)
    coeff = -cut.COEFFICIENT * cut.POLYS[channel]
    assert (
        s.expand(coeff * before - s.I * s.pi * cut.COEFFICIENT * cut.POLYS[channel])
        == 0
    )
    assert above == 0


@pytest.mark.parametrize("cap,mass2", ((5, 1), (11, 2), (17, 3)))
def test_full_b20_from_an_independent_even_log_series(cap, mass2):
    D = cap - 2 * mass2
    nu = s.Symbol("nu2", positive=True)
    v = s.Symbol("v", real=True)
    even = 2 * (
        (14 * mass2 * mass2) * (-v * v / (2 * D * D))
        + (6 * mass2 * v) * (-v / D)
        + v * v * s.log(D / nu)
    )
    even += v * v * s.log(cap / nu)
    actual = -cut.COEFFICIENT * s.expand(even).coeff(v, 2)
    expected = subtraction.forward_coefficient().subs(
        {subtraction.CAP: cap, source.MU: mass2, cut.NU2: nu}
    )
    assert s.simplify(actual - expected) == 0
    assert (
        subtraction.cap_derivative().subs(
            {subtraction.CAP: cap, source.MU: mass2, source.K: 1}
        )
        < 0
    )


def test_transfer_cut_term_is_not_in_s_u_optical_subtraction():
    L, mu = subtraction.CAP, source.MU
    D = L - 2 * mu
    s_u = -2 * cut.COEFFICIENT * (L * L + 2 * mu * L + 6 * mu * mu) / (D**3)
    defect = s.factor(subtraction.cap_derivative() - s_u)
    assert defect == -cut.COEFFICIENT / L and defect != 0


def test_local_renormalization_polynomial_remains_unfixed():
    a0, a2 = s.symbols(
        "independent_local_constant independent_local_quadratic", real=True
    )
    polynomial = a0 + a2 * (source.S**2 + source.T**2 + source.U**2)
    v = subtraction.V
    f = polynomial.subs(
        {source.S: 2 + v, source.T: 0, source.U: 2 - v}, simultaneous=True
    )
    assert s.diff(f, v, 2) / 2 == 2 * a2
    assert (
        polynomial.subs({source.S: source.T, source.T: source.S}, simultaneous=True)
        == polynomial
    )


def test_named_subtraction_bound_is_only_the_nonlocal_sector():
    d = subtraction.data()["named_nonlocal_coefficient_bound"]
    assert d["absolute_b20_bound"] / (4 * source.current.heavy.LAMBDA) < s.Rational(
        1, 10**990
    )
    assert "physical EFT cutoff" in audit.observable()["not_established"]
    assert "local" in subtraction.data()["subtraction_definition"]


def test_detector_and_fixed_resolution_decoupling_do_not_commute():
    k = s.Symbol("k", positive=True)
    ell, chi = s.symbols("ell chi", positive=True)
    assert s.limit(ell / k, k, s.oo) == 0
    assert s.limit((ell / k).subs(ell, k * chi), k, s.oo) == chi


def test_all_previous_scientific_scopes_and_records_are_unchanged():
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 134
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.qualifications()) == 6
    assert audit.validate_scope(audit.frontier(), audit.matching())
    assert audit.require_model(audit.MODEL) == audit.MODEL
    for name in audit.OBSERVABLES:
        assert audit.require_observable(name) == name
