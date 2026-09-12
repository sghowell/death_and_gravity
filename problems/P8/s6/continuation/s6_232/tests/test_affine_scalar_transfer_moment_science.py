"""Independent transfer-cut convolution, positive measures and rational matching controls."""

import pytest
import sympy as s
from p8_vacuum_affine_scalar_transfer_moment import atomic, audit, moment
from p8_vacuum_affine_scalar_transfer_moment import coefficients as c


@pytest.mark.parametrize(
    "name,value", list(audit.residuals().items()), ids=list(audit.residuals())
)
def test_all_exact_residuals(name, value):
    assert value == 0, name


@pytest.mark.parametrize(
    "name,value", list(audit.gates().items()), ids=list(audit.gates())
)
def test_all_proof_gates(name, value):
    assert value is True, name


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[row[0] for row in audit.bad_cases()]
)
def test_every_unsupported_input_rejected(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


def test_unchanged_original_frontier_and_non_completion_control():
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 88
    assert audit.frontier() == audit.previous.frontier()
    assert all(row["status"] != "COMPLETE" for row in audit.frontier())
    assert "NOT_FULL_LOOP_BOUND" in audit.ITEM["status"]
    assert audit.require_tolerances(3, 1, s.Rational(1, 2)) == (3, 1, s.Rational(1, 2))


def test_independent_original_crossing_polynomial():
    v, t = s.symbols("v t")
    ss, uu = 2 - t / 2 + v, 2 - t / 2 - v
    A = (
        2 * c.lam * ((ss - 2) ** 2 + (t - 2) ** 2 + (uu - 2) ** 2)
        + 3 * c.gam * ss * t * uu
        - 8 * c.gam
    )
    assert s.expand(A).coeff(v, 2) == 4 * c.lam - 3 * c.gam * t
    assert s.diff(A, v, 4) == 0
    assert s.diff(A, v, 2, t).subs({v: 0, t: 0}) / 2 == -3 * c.gam
    # At this first jet fixed-s and fixed-v agree by crossing-evenness.
    S = s.Symbol("s")
    fixed_s = A.subs(v, S + t / 2 - 2)
    assert s.expand(s.diff(fixed_s, S, 2, t).subs({S: 2, t: 0}) / 2 + 3 * c.gam) == 0


def test_independent_full_two_angle_elastic_cut_convolution():
    y, z = s.symbols("intermediate_cosine scattering_cosine", real=True)
    a, b, beta = s.symbols("a b beta", real=True)
    # Azimuth average of (n_out.n)^2, with n_in.n=y.
    out_square = z * z * y * y + (1 - z * z) * (1 - y * y) / 2
    product = (a + b * y * y) * (a + b * out_square)
    literal = beta / (64 * s.pi) * s.integrate(s.expand(product), (y, -1, 1))
    wanted = beta / (32 * s.pi) * ((a + b / 3) ** 2 + 4 * b * b * s.legendre(2, z) / 45)
    assert s.expand(literal - wanted) == 0
    mapped = literal.subs(
        {
            a: c.original.a,
            b: c.original.b,
            beta: c.original.beta,
            z: 1 + 2 * c.transfer / (c.S - 4),
        },
        simultaneous=True,
    )
    packet = c.data()
    assert (
        s.simplify(mapped - packet["complete_first_elastic_nonforward_coefficient"])
        == 0
    )
    assert (
        s.simplify(
            s.diff(mapped, c.transfer).subs(c.transfer, 0)
            - packet["complete_first_elastic_transfer_slope"]
        )
        == 0
    )
    # Squaring a pointwise tree is not the nonforward cut convolution.
    assert s.expand(literal - beta * (a + b * z * z) ** 2 / (32 * s.pi)) != 0


@pytest.mark.parametrize("ell", range(11))
def test_independent_rodrigues_positive_endpoint_derivative(ell):
    x = s.Symbol("x")
    P = s.diff((x * x - 1) ** ell, x, ell) / (2**ell * s.factorial(ell))
    value = s.diff(P, x).subs(x, 1)
    assert value == s.Rational(ell * (ell + 1), 2) and value >= 0


def test_crossing_half_shift_and_inverse_fourth_sign():
    mu, t = s.symbols("mu t", positive=True)
    r0, r1 = s.symbols("rho rho_t")
    shifted = 2 * (r0 + t * r1) / (s.pi * (mu + t / 2 - 2) ** 3)
    wanted = 2 / s.pi * (r1 / (mu - 2) ** 3 - s.Rational(3, 2) * r0 / (mu - 2) ** 4)
    assert s.cancel(s.diff(shifted, t).subs(t, 0) - wanted) == 0
    unshifted = 2 * (r0 + t * r1) / (s.pi * (mu - 2) ** 3)
    assert (
        s.cancel(
            s.diff(unshifted, t).subs(t, 0) - wanted - 3 * r0 / (s.pi * (mu - 2) ** 4)
        )
        == 0
    )


@pytest.mark.parametrize("ells", [(0, 0, 0, 0), (0, 2, 4, 6), (2, 0, 2, 0)])
def test_independent_positive_discrete_measure_split(ells):
    mus = (s.Integer(5), s.Integer(10), s.Integer(25), s.Integer(60))
    weights = (s.Rational(1, 3), s.Rational(2, 5), s.Rational(5, 7), s.Rational(7, 11))
    K = s.Integer(16)
    low0 = high0 = low4 = high4 = positive = 0
    for mu, w, ell in zip(mus, weights, ells):
        B = 2 * w / (mu - 2) ** 3
        J = 2 * w / (mu - 2) ** 4
        positive += 2 * w * ell * (ell + 1) / ((mu - 4) * (mu - 2) ** 3)
        if mu <= K:
            low0 += B
            low4 += J
        else:
            high0 += B
            high4 += J
    b0 = low0 + high0
    b1 = positive - s.Rational(3, 2) * (low4 + high4)
    full = b1 + 3 * b0 / (2 * (K - 2)) + 3 * low4 / 2
    decomposition = (
        positive
        + 3 * low0 / (2 * (K - 2))
        + s.Rational(3, 2) * (high0 / (K - 2) - high4)
    )
    assert s.cancel(full - decomposition) == 0
    assert full >= 0 and high4 <= high0 / (K - 2)


def test_erasing_low_cut_breaks_the_required_inequality():
    mu, K = s.Integer(5), s.Integer(100)
    b0 = 2 / (mu - 2) ** 3
    J4 = 2 / (mu - 2) ** 4
    b1 = -3 / (mu - 2) ** 4
    assert b1 + 3 * b0 / (2 * (K - 2)) < 0
    assert b1 + 3 * b0 / (2 * (K - 2)) + 3 * J4 / 2 > 0


def test_independent_first_elastic_integral_upper_and_positive_dropped_terms():
    S, K = s.symbols("S K", positive=True)
    l, g = s.symbols("lambda gamma", positive=True)
    literal = (
        3
        / s.pi**2
        * s.integrate(
            16 * l * l + s.Rational(9, 16) * g * g * S * S + 64 * g * g / S**4,
            (S, 4, K),
        )
    )
    upper = (48 * l * l * K + s.Rational(9, 16) * g * g * K**3 + g * g) / s.pi**2
    difference = (192 * l * l + 36 * g * g + 64 * g * g / K**3) / s.pi**2
    assert s.cancel(upper - literal - difference) == 0
    assert difference.is_positive
    z = s.Symbol("S_minus_four", nonnegative=True)
    assert s.expand((S - 2) - S / 2).subs(S, z + 4) == z / 2


def test_named_mandatory_weight_and_first_elastic_margin_are_exact():
    K = s.Integer(10) ** 198
    required = c.GAMMA - 8 * c.LAMBDA / (K - 2)
    upper = (
        48 * c.LAMBDA * c.LAMBDA * K
        + s.Rational(9, 16) * c.GAMMA * c.GAMMA * K**3
        + c.GAMMA * c.GAMMA
    ) / 9
    assert required > c.GAMMA / 5
    assert 0 < upper < c.GAMMA / 10**200
    assert required / upper > 2 * 10**199
    assert K > 4 and c.LAMBDA > 0 and c.GAMMA > 0


@pytest.mark.parametrize(
    "d0,d1", [(0, 0), (1, s.Rational(1, 2)), (2, s.Rational(1, 4))]
)
def test_general_error_tradeoff_with_positive_full_moment_budget(d0, d1):
    eps = s.Rational(1, 10)
    denominator = c.GAMMA * (2 * (1 - d1) - eps)
    assert denominator > 0
    threshold = 2 + 4 * c.LAMBDA * (1 + d0) / denominator
    # At equality the necessary lower moment equals the selected FULL budget.
    lower = 2 * c.GAMMA * (1 - d1) - 4 * c.LAMBDA * (1 + d0) / (threshold - 2)
    assert s.cancel(lower - eps * c.GAMMA) == 0


def test_independent_positive_atom_two_coefficients_and_extra_v4():
    l, g = s.symbols("lambda gamma", positive=True)
    gap = 2 * l / g
    mass2 = 2 + gap
    coupling2 = g * gap**4
    v, t = s.symbols("v t")
    amplitude = coupling2 * (
        1 / (mass2 - (2 - t / 2 + v)) + 1 / (mass2 - t) + 1 / (mass2 - (2 - t / 2 - v))
    )
    b20 = s.diff(amplitude, v, 2).subs({v: 0, t: 0}) / 2
    b21 = s.diff(amplitude, v, 2, t).subs({v: 0, t: 0}) / 2
    b40 = s.diff(amplitude, v, 4).subs({v: 0, t: 0}) / 24
    assert s.cancel(b20 - 4 * l) == 0
    assert s.cancel(b21 + 3 * g) == 0
    assert s.cancel(b40 - g * g / l) == 0
    assert s.cancel(2 * coupling2 / gap**4 - 2 * g) == 0
    assert atomic.COUPLING2 == s.Rational(1, 2**26)
    assert 4 < atomic.MASS2 < moment.SPLIT


def test_rational_atom_is_real_tree_not_exact_elastic_unitarity():
    # A finite proxy of the separately named rational family makes the scope
    # issue explicit. It is not an approximation of the original couplings.
    energy2 = s.Integer(8)
    pole_mass2 = s.Integer(20)
    g2 = s.Rational(1, 100)
    x = s.Symbol("x", real=True)
    t, u = -(energy2 - 4) * (1 - x) / 2, -(energy2 - 4) * (1 + x) / 2
    A = g2 * (1 / (pole_mass2 - energy2) + 1 / (pole_mass2 - t) + 1 / (pole_mass2 - u))
    assert A.subs(x, 0) > 0
    # Positive nonzero real angular amplitude has real positive t0 but no
    # absorptive part, whereas its exact elastic optical RHS is positive.
    assert s.im(A) == 0
    assert A.subs(x, s.Rational(1, 2)) ** 2 > 0
    assert (
        "lacks the light elastic cut required by exact unitarity"
        in atomic.data()["scope"]
    )


def test_full_tree_size_and_no_automatic_transfer_premise():
    K = moment.SPLIT
    upper = 4 * c.LAMBDA * K * K + 3 * c.GAMMA * K**3 / 4 + 8 * c.GAMMA
    assert upper < s.Rational(1, 10**202)
    assert "differentiable in t" in moment.data()["explicit_unproved_physical_premises"]
    assert "not a UV cutoff" in moment.data()["full_low_moment_definition"]


@pytest.mark.parametrize("masses", [(5, 7), (6, 9), (8, 21)])
def test_independent_full_moment_Gram_lower_for_negative_transfer(masses):
    mu0, mu1 = map(s.Integer, masses)
    w0, w1 = s.Rational(2, 3), s.Rational(5, 7)
    b0 = 2 * w0 / (mu0 - 2) ** 3 + 2 * w1 / (mu1 - 2) ** 3
    J = 2 * w0 / (mu0 - 2) ** 4 + 2 * w1 / (mu1 - 2) ** 4
    b4 = 2 * w0 / (mu0 - 2) ** 5 + 2 * w1 / (mu1 - 2) ** 5
    b1 = -s.Rational(3, 2) * J
    assert b0 * b4 - J * J > 0
    assert b1 < 0
    assert b4 > 4 * b1 * b1 / (9 * b0) > 0


def test_independent_b40_kernel_and_named_physical_lower():
    mu, v = s.symbols("mu v")
    crossing_pair = 1 / (mu - 2 - v) + 1 / (mu - 2 + v) - 2 / (mu - 2)
    coefficient = s.diff(crossing_pair, v, 4).subs(v, 0) / (24 * s.pi)
    assert s.cancel(coefficient - 2 / (s.pi * (mu - 2) ** 5)) == 0
    lower = c.GAMMA**2 / (8 * c.LAMBDA)
    assert lower > 0
    assert s.cancel(4 * (-3 * c.GAMMA / 2) ** 2 / (9 * 8 * c.LAMBDA) - lower) == 0
    assert atomic.data()["diagnostic_nonzero_v4_coefficient"] == 8 * lower
    assert s.diff(c.data()["complete_original_crossing_tree"], c.v, 4) == 0


def test_positive_atom_Gram_saturation_not_exact_unitarity_claim():
    b0 = 4 * c.LAMBDA
    b4 = atomic.data()["diagnostic_nonzero_v4_coefficient"]
    J = atomic.data()["diagnostic_full_atom_low_moment"]
    assert b0 * b4 - J * J == 0
    assert (
        atomic.data()["gates"]["moment_bound_attainment_not_exact_S_matrix_optimality"]
        is True
    )
