"""Independent labelled vertices, identical normalization and exact matching margins."""

from itertools import permutations

import pytest
import sympy as s
from p8_vacuum_affine_scalar_matching_scale import amplitude as a
from p8_vacuum_affine_scalar_matching_scale import audit, dispersion, unitarity


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


def test_unchanged_frontier_and_conditional_not_cutoff():
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 87
    assert audit.frontier() == audit.previous.frontier()
    assert all(row["status"] != "COMPLETE" for row in audit.frontier())
    assert "CONDITIONAL_MATCHING_DISJUNCTION" in audit.ITEM["status"]
    assert "NOT_PHYSICAL_CUTOFF" in audit.ITEM["status"]
    assert audit.require_tolerances(6, s.Rational(1, 2), 1) == (6, s.Rational(1, 2), 1)


def test_independent_twenty_four_labelled_original_vertices():
    S, x = a.S, a.x
    T, U = -(S - 4) * (1 - x) / 2, -(S - 4) * (1 + x) / 2
    av, bv, cv = (S - 2) / 2, (T - 2) / 2, (U - 2) / 2
    dots = s.Matrix(
        [[1, av, bv, cv], [av, 1, cv, bv], [bv, cv, 1, av], [cv, bv, av, 1]]
    )
    lo = 0
    high = 0
    for i, j, k, l in permutations(range(4)):
        lo += a.lam * dots[i, j] * dots[k, l]
        high += (
            2
            * a.gam
            * (
                -dots[i, j] * dots[j, k] * dots[k, l]
                + dots[i, i] * dots[j, k] * dots[k, l]
            )
        )
    actual = s.expand(lo + high - s.factorial(4) * a.gam / 3)
    assert s.expand(actual - a.TREE) == 0
    assert s.expand(actual - (a.TREE + 8 * a.gam)) != 0


@pytest.mark.parametrize("ell", range(9))
def test_independent_rodrigues_projection(ell):
    x = a.x
    P = s.diff((x * x - 1) ** ell, x, ell) / (2**ell * s.factorial(ell))
    projected = s.expand(s.integrate(a.TREE * P, (x, -1, 1)))
    target = 2 * a.Q0 if ell == 0 else 4 * a.b / 15 if ell == 2 else 0
    assert s.expand(projected - target) == 0


def test_independent_exact_three_node_angular_norm():
    z = s.sqrt(s.Rational(3, 5))
    quadrature = (
        s.Rational(5, 9) * (a.TREE.subs(a.x, z) ** 2 + a.TREE.subs(a.x, -z) ** 2)
        + s.Rational(8, 9) * a.TREE.subs(a.x, 0) ** 2
    )
    integral = s.integrate(s.expand(a.TREE * a.TREE), (a.x, -1, 1))
    assert s.expand(quadrature - integral) == 0
    assert s.expand(integral - 2 * a.Q0 * a.Q0) == s.expand(8 * a.b * a.b / 45)


def test_independent_phase_space_and_identical_bubble_factor():
    beta, c = s.symbols("beta c", positive=True)
    # Full dOmega integration, one final-state factorial, and 2 Im A.
    imA = (
        s.Rational(1, 2) * s.Rational(1, 2) * beta / (32 * s.pi**2) * (4 * s.pi) * c * c
    )
    assert imA == beta * c * c / (32 * s.pi)
    t = beta * c / (32 * s.pi)
    assert s.simplify(beta * imA / (32 * s.pi) - t * t) == 0
    wrong = beta * c / (16 * s.pi)
    assert s.simplify(beta * imA / (16 * s.pi) - wrong * wrong) == -wrong * wrong / 2
    assert s.simplify(2 * imA - beta / (32 * s.pi**2) * (4 * s.pi) * c * c) != 0


@pytest.mark.parametrize(
    "real,imag",
    [
        (s.Rational(1, 2), s.Rational(1, 2)),
        (0, 0),
        (0, 1),
        (s.Rational(3, 10), s.Rational(2, 5)),
    ],
)
def test_direct_unitarity_disk(real, imag):
    partial = real + s.I * imag
    S = 1 + 2 * s.I * partial
    assert s.Abs(S) ** 2 <= 1
    assert imag >= real * real + imag * imag
    assert abs(real) <= s.Rational(1, 2)


def test_inelastic_channel_and_nonzero_tree_controls():
    matrix = s.Matrix([[3, 4], [-4, 3]]) / 5
    assert matrix.T * matrix == s.eye(2)
    partial = (matrix[0, 0] - 1) / (2 * s.I)
    assert s.im(partial) > s.Abs(partial) ** 2
    tree = s.Rational(1, 1000)
    assert s.Abs(1 + 2 * s.I * tree) ** 2 > 1
    assert tree < s.Rational(1, 2)


def test_independent_massive_endpoint_and_central_polynomials():
    z = s.Symbol("nonnegative", nonnegative=True)
    plus = s.Poly(s.expand((a.a + a.b).subs(a.S, z + 4)), z)
    minus = s.Poly(s.expand((a.a - a.b).subs(a.S, z + 4)), z)
    for p in (plus, minus):
        assert p.coeff_monomial(1) == 24 * a.lam - 8 * a.gam
        assert all(c.is_positive for c in p.all_coeffs()[:-1])
    central = s.expand(s.Rational(9, 16) * a.S * (a.S - 4) ** 2 - 8 - a.S**3 / 4)
    assert all(c > 0 for c in s.Poly(central.subs(a.S, z + 16), z).all_coeffs())
    assert all(
        v > 0 for v in unitarity.data()["all_positive_rational_margins"].values()
    )


@pytest.mark.parametrize(
    "energy,relation",
    [(10**133, "below"), (2 * 10**133, "above"), (10**134, "correction")],
)
def test_independent_exact_scale_brackets(energy, relation):
    S = s.Integer(energy) ** 2
    Q = (
        a.LAMBDA * (10 * S * S - 32 * S + 40) / 3
        + a.GAMMA * S * (S - 4) ** 2 / 2
        - 8 * a.GAMMA
    )
    if relation == "below":
        assert Q / (32 * 3) < s.Rational(1, 2)
    else:
        assert 1 - 4 / S > s.Rational(999, 1000) ** 2
        bound = s.Rational(999, 1000) * Q / (32 * s.Rational(22, 7))
        assert bound > (3 if relation == "above" else 50000)
    assert a.LAMBDA > 0 and a.GAMMA > 0


def test_independent_forward_half_second_derivative():
    v = s.Symbol("v")
    S, T, U = 2 + v, 0, 2 - v
    forward = (
        2 * a.lam * ((S - 2) ** 2 + (T - 2) ** 2 + (U - 2) ** 2)
        + 3 * a.gam * S * T * U
        - 8 * a.gam
    )
    assert s.diff(forward, v, 2) / 2 == 4 * a.lam
    assert s.expand(forward - (4 * a.lam * v * v + 8 * a.lam - 8 * a.gam)) == 0
    spectral = s.Symbol("spectral", positive=True)
    pair = 1 / (spectral - 2 - v) + 1 / (spectral - 2 + v) - 2 / (spectral - 2)
    assert (
        s.cancel(
            s.diff(pair, v, 2).subs(v, 0) / (2 * s.pi)
            - 2 / (s.pi * (spectral - 2) ** 3)
        )
        == 0
    )


def test_independent_exact_annulus_optical_lower_and_tolerances():
    K, eta = s.symbols("K eta", positive=True)
    sqnorm = a.GAMMA**2 * a.S**6 / 16
    absorptive = s.Rational(3, 4) * (1 - eta) ** 2 * sqnorm / (64 * s.pi)
    bound = s.integrate(2 * absorptive / (s.pi * a.S**3), (a.S, K / 2, K))
    assert (
        s.simplify(bound - 45 * (1 - eta) ** 2 * a.GAMMA**2 * K**4 / (131072 * s.pi**2))
        == 0
    )
    chosen = bound.subs({K: dispersion.AGREEMENT_ENERGY**2, eta: s.Rational(1, 2)})
    assert s.simplify(chosen - 90 * a.LAMBDA / s.pi**2) == 0
    assert 9 * a.LAMBDA > 8 * a.LAMBDA
    # Halving the cut or confusing b2 with B'' would lose this named contradiction.
    assert s.Rational(9, 2) * a.LAMBDA < 8 * a.LAMBDA


@pytest.mark.parametrize(
    "eta", [s.S.Zero, s.Rational(1, 4), s.Rational(1, 2), s.Rational(3, 4)]
)
def test_full_error_triangle_sharpness_and_complex_error(eta):
    tree = 1 + 3 * a.x * a.x
    norm = s.integrate(tree * tree, (a.x, -1, 1))
    closest = (1 - eta) * tree
    assert s.integrate((closest - tree) ** 2, (a.x, -1, 1)) == eta * eta * norm
    assert s.integrate(closest * closest, (a.x, -1, 1)) == (1 - eta) ** 2 * norm
    complex_amplitude = (1 + s.I * eta) * tree
    squared = s.expand(complex_amplitude * s.conjugate(complex_amplitude))
    assert s.integrate(squared, (a.x, -1, 1)) >= (1 - eta) ** 2 * norm


def test_small_tree_at_dispersion_scale_and_tensor_units_only():
    E = dispersion.AGREEMENT_ENERGY
    bound = 4 * a.LAMBDA * E**4 + 3 * a.GAMMA * E**6 / 4 + 8 * a.GAMMA
    assert bound < s.Rational(1, 10**46)
    assert E < unitarity.ENERGY_LO < unitarity.ENERGY_HI < 10**398
    assert (
        "No uniform quantum decoupling" in dispersion.data()["tensor_pole_comparison"]
    )
