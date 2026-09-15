"""Independent covariant insertion, whole Gaussian cut and analytic bound checks."""

import pytest
import sympy as s
from p8_affine import verify as base
from p8_vacuum_affine_covariant_gaussian_four_point_matching import (
    audit,
    forward,
    gaussian,
    local,
    source,
)
from p8_vacuum_affine_gaussian_metric_pole_matching import spectral

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
@pytest.mark.parametrize("spacelike", (False, True))
@pytest.mark.parametrize("choice", (1, 2))
def test_independent_non_rest_frame_conserved_metric_response(D, spacelike, choice):
    eta = s.diag(1, *([-1] * (D - 1)))
    q = s.Matrix(
        [1, choice + 3, 2, *([0] * (D - 3))]
        if spacelike
        else [choice + 3, 1, 2, *([0] * (D - 3))]
    )
    down = eta * q
    q2 = (q.T * eta * q)[0]
    assert q2 < 0 if spacelike else q2 > 0
    projector = s.eye(D) - down * q.T / q2
    seed = s.Matrix(D, D, lambda i, j: s.Rational((i + 1) * (j + 1), 7) + int(i == j))
    stress = projector * seed * projector.T
    assert stress * q == s.zeros(D, 1)
    trace = s.trace(eta * stress)
    bar = stress - eta * trace / (D - 2)
    kappa = s.Integer(3)
    metric = 2 * bar / (kappa * q2)
    div = metric * q
    ricci_old = (
        down * div.T
        + div * down.T
        - q2 * metric
        - down * down.T * s.trace(eta * metric)
    ) / 2
    assert ricci_old == -bar / kappa
    assert s.trace(eta * ricci_old) == 2 * trace / ((D - 2) * kappa)


@pytest.mark.parametrize("D", (4, s.Rational(9, 2), 5))
@pytest.mark.parametrize("column", range(8))
def test_independent_eight_operator_scalar_Fourier_projection(D, column):
    D = s.sympify(D)
    mu, energy, transfer = map(s.Integer, (2, 7, -1))
    crossed = 4 * mu - energy - transfer
    p = s.symbols("incoming_phi0:4")
    gram = s.Matrix(
        [
            [
                mu,
                (energy - 2 * mu) / 2,
                (transfer - 2 * mu) / 2,
                (crossed - 2 * mu) / 2,
            ],
            [
                (energy - 2 * mu) / 2,
                mu,
                (crossed - 2 * mu) / 2,
                (transfer - 2 * mu) / 2,
            ],
            [
                (transfer - 2 * mu) / 2,
                (crossed - 2 * mu) / 2,
                mu,
                (energy - 2 * mu) / 2,
            ],
            [
                (crossed - 2 * mu) / 2,
                (transfer - 2 * mu) / 2,
                (energy - 2 * mu) / 2,
                mu,
            ],
        ]
    )
    X = -(s.Matrix(p).T * gram * s.Matrix(p))[0]
    Y = sum(p) ** 2
    R = -X + D * mu * Y / (D - 2)
    Ric2 = X * X - 2 * mu * X * Y / (D - 2) + D * mu * mu * Y * Y / (D - 2) ** 2
    Weyl2 = 4 * (D - 3) * Ric2 / (D - 2) - D * (D - 3) * R * R / ((D - 1) * (D - 2))
    operators = (
        X * X,
        mu * X * Y,
        mu * mu * Y * Y,
        R * X,
        -X * X + mu * X * Y / (D - 2),
        mu * R * Y,
        R * R,
        Weyl2,
    )
    actual = s.Poly(s.expand(operators[column]), *p).coeff_monomial(s.prod(p))
    a2, a0 = local.projection_matrix(D)[:, column]
    assert actual == a2 * (energy**2 + transfer**2 + crossed**2) + a0 * mu**2


@pytest.mark.parametrize("D", (4, s.Rational(9, 2), 5, 6))
def test_entire_local_nullspace_and_rank(D):
    matrix = local.projection_matrix(D)
    nulls = local.nullspace_columns(D)
    assert matrix.rank() == 2
    assert nulls.rank() == 6
    assert matrix * nulls == s.zeros(2, 6)


@pytest.mark.parametrize("species", ("scalar", "vector"))
@pytest.mark.parametrize("mu", (s.Integer(1), s.Rational(3, 2)))
@pytest.mark.parametrize("n", (s.Integer(2), s.Integer(7)))
def test_independent_full_crossed_amplitude_forward_jet(species, mu, n):
    v, x = s.symbols("independent_crossing_v integration_velocity", real=True)
    energy, transfer = 2 * mu + v, s.S.Zero
    crossed = 4 * mu - energy - transfer
    whole = 0
    for a, b, c in (
        (energy, transfer, crossed),
        (transfer, energy, crossed),
        (crossed, energy, transfer),
    ):
        numerator = 2 * mu * mu - 2 * mu * a - b * c
        trace = (a + 2 * mu) ** 2
        tt = numerator + trace / 6
        w = a * spectral.weight(species, 2, x) / (64 * (4 * n - a * (1 - x * x)))
        r = a * spectral.weight(species, 0, x) / (4608 * (4 * n - a * (1 - x * x)))
        whole += 4 * w * tt + 2 * r * trace
    actual = s.diff(whole, v, 2).subs(v, 0) / 2
    assert s.factor(actual - forward.nonlocal_integrand(species, n, mu, x)) == 0


@pytest.mark.parametrize("species", ("scalar", "vector"))
@pytest.mark.parametrize("mu", (s.Integer(1), s.Rational(3, 2)))
@pytest.mark.parametrize("ratio", (s.Integer(1), s.Rational(3, 2), s.Integer(7)))
@pytest.mark.parametrize("x", (s.Integer(0), s.Rational(1, 3), s.Integer(1)))
def test_independent_exact_full_domain_majorants(species, mu, ratio, x):
    n = mu * ratio
    j2, j0 = forward.numerator_jets(n, mu, x)
    assert abs(j2) <= 6 * mu / n
    assert 0 < j0 <= 34 * mu / n
    w2 = spectral.weight(species, 2, x)
    w0 = spectral.weight(species, 0, x)
    assert w2 >= 0 and w0 >= 0
    value = forward.nonlocal_integrand(species, n, mu, x)
    assert abs(value) <= mu * (w2 / 16 + 17 * w0 / 1152) / n


@pytest.mark.parametrize("species", ("scalar", "vector"))
@pytest.mark.parametrize("beta", (s.Rational(1, 3), s.Rational(3, 5)))
@pytest.mark.parametrize("angle", (-s.S.One, s.S.Zero, s.Rational(2, 3), s.S.One))
def test_independent_all_angle_gaussian_physical_cut(species, beta, angle):
    n, mu, kappa = map(s.Integer, (2, 1, 3))
    energy = 4 * n / (1 - beta * beta)
    transfer = -(energy - 4 * mu) * (1 - angle) / 2
    crossed = 4 * mu - energy - transfer
    tt, trace = gaussian.tensor_shapes(energy, transfer, crossed, mu)
    actual = (
        spectral.weight(species, 2, beta) * tt / (32 * s.pi * beta)
        + spectral.weight(species, 0, beta) * trace / (4608 * s.pi * beta)
    ) / kappa**2
    if species == "scalar":
        y = s.Symbol("pair_angle", real=True)
        tree = spectral.scalar_tree(energy, y, mu, n, kappa)
        a0 = s.integrate(tree, (y, -1, 1)) / 2
        a2 = 5 * s.integrate(tree * s.legendre(2, y), (y, -1, 1)) / 2
        expected = beta * (a0 * a0 + a2 * a2 * s.legendre(2, angle) / 5) / (32 * s.pi)
    else:
        expected = gaussian.proca.whole_cut(energy, angle, mu, n, kappa)
    assert s.simplify(actual - expected) == 0


@pytest.mark.parametrize("species", ("scalar", "vector"))
def test_explicit_upper_feynman_boundary_and_zero_remainder(species):
    mass = s.Integer(2)
    assert tuple(
        value.doit() for value in gaussian.effective_curvature(species, mass, 0)
    ) == (0, 0)
    values = gaussian.feynman_boundary_curvature(species, mass, 9)
    for value in values:
        assert isinstance(value, s.Limit)
        assert value.args[3] == s.Symbol("+")
        assert value.args[0].has(s.I)
        assert not value.has(s.Float)


@pytest.mark.parametrize("bad", (True, False, 1.0, s.Float(1), None, "1", s.I, s.oo))
def test_bound_domains_reject_inexact_or_unspecified_masses(bad):
    with pytest.raises((TypeError, ValueError)):
        forward.absolute_nonlocal_bound("scalar", bad, 1)
    with pytest.raises((TypeError, ValueError)):
        forward.absolute_nonlocal_bound("vector", 2, bad)


@pytest.mark.parametrize("bad", (-1, 0, s.Rational(1, 2)))
def test_bound_domain_rejects_light_loop_below_external_mass(bad):
    with pytest.raises((TypeError, ValueError)):
        forward.absolute_nonlocal_bound("scalar", bad, 1)


@pytest.mark.parametrize("bad", ("Phi", "graviton", "", 1, True, None))
def test_unsupported_gaussian_species_rejected(bad):
    with pytest.raises((TypeError, ValueError)):
        gaussian.effective_curvature(bad, 2, 1)
    with pytest.raises((TypeError, ValueError)):
        forward.nonlocal_integrand(bad, 2, 1)


def test_exact_light_integral_and_positive_entire_polynomial():
    x = s.Symbol("independent_light_velocity", real=True)
    value = forward.nonlocal_integrand("scalar", 1, 1, x)
    primitive = s.integrals.rationaltools.ratint(value, x)
    exact = s.simplify(primitive.subs(x, 1) - primitive.subs(x, 0))
    assert exact == s.Rational(319, 4800) - 23 * s.pi / 1280
    P = 128 + 120 * (1 - x * x) + 4 * (1 - x * x) ** 2 + 3 * (1 - x * x) ** 4
    assert s.factor(value - x * x * P / (1920 * (1 + x * x) ** 3)) == 0


def test_actual_heavy_interval_uses_whole_bound_not_its_series():
    k = source.KAPPA
    lo, hi = forward.fixed_HP_interval()
    bound = forward.absolute_nonlocal_bound("scalar", source.HEAVY_MASS2, 1, k)
    bound += forward.absolute_nonlocal_bound("vector", source.VECTOR_MASS2, 1, k)
    center = -(s.log(source.HEAVY_MASS2) + 2) / (320 * s.pi**2 * k * k)
    assert s.simplify((lo + hi) / 2 - center) == 0
    assert s.simplify((hi - lo) / 2 - bound) == 0
    assert bound * s.pi**2 * k * k < s.Rational(1, 10**7)
    assert all(
        v > 0 for v in forward.data()["whole_exact_positive_bound_margins"].values()
    )


def test_unmatched_light_constants_are_explicit_and_not_set_to_heavy_scheme():
    data = forward.data()
    coeff = data["whole_light_Gaussian_with_unmatched_local"]
    names = {str(x) for x in coeff.free_symbols}
    assert "unmatched_Phi_R_squared" in names
    assert "unmatched_Phi_Weyl_squared" in names
    assert "not the full-source" in data["actual_fixed_HP_verdict"]


def test_integer_math_apis_preserve_exact_rational_arithmetic():
    values = (
        *local.projection_matrix(4),
        *local.nullspace_columns(4),
        *gaussian.tensor_shapes(6, -1, -1, 1),
        gaussian.local_amplitude(1, 2, 6, -1, 1, 3),
        *gaussian.effective_curvature("scalar", 2, 1),
        *forward.numerator_jets(2, 1, 0),
        forward.nonlocal_integrand("scalar", 2, 1, 0),
        forward.nonlocal_integrand("vector", 2, 1, 1),
        forward.absolute_nonlocal_bound("scalar", 2, 1, 3),
    )
    assert all(isinstance(v, s.Basic) and not v.has(s.Float) for v in values)


def test_parent_cache_warmth_does_not_mutate_the_new_packet():
    before = dict(audit.residuals())
    audit.previous.packets()
    assert audit.residuals() == before
    assert source.data()["checks"] is not source.previous.data()["checks"]
    base.certify_residuals(before)


def test_whole_original_frontier_matching_and_qualifications_retained():
    assert len(audit.frontier()) == 9
    assert len(audit.matching()) == 145
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.qualifications()) == 6
    assert audit.validate_scope(audit.frontier(), audit.matching())
    assert any(row["status"].startswith("REJECTED_") for row in audit.matching())
    assert "Regge" in audit.observable()["not_established"]
    assert "P8" in audit.observable()["not_established"]
