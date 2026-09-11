"""Independent full-function germs, analytic bounds and function-class controls."""

from fractions import Fraction
from functools import cache

import mpmath as mp
import pytest
import sympy as s
from p8_exceptional_vacuum import family
from p8_offshell_vacuum import band
from p8_offshell_vacuum import jets as parent_jets
from p8_vacuum_full_analytic_target_match import (
    audit,
    bounds,
    calibration,
    holomorphic,
    jets,
    transport,
)
from p8_vacuum_two_loop_elastic_cut import audit as previous


@pytest.mark.parametrize("name", list(audit.residuals()))
def test_named_exact_identity(name):
    assert audit.residuals()[name] == 0


@pytest.mark.parametrize("name", list(audit.gates()))
def test_each_proof_gate(name):
    assert audit.gates()[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_rejected_input(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@cache
def tree_evaluator():
    return s.lambdify(
        (family.u, family.X), family.data()["original_retuned_tree_scalar"], "mpmath"
    )


def number(value):
    value = s.Rational(value)
    return mp.mpf(int(value.p)) / int(value.q)


def full_density(n, kappa, lam, t, base):
    """Independent evaluation of the full rational/exponential target."""
    p, y, l3, l4, l5 = map(number, base)
    kappa, lam = number(kappa), number(lam)
    u = t * p / mp.sqrt(kappa)
    x = t * t * y / kappa
    denominator = x**n + (1 - x) ** n
    step = x**n / denominator
    step_prime = n * x ** (n - 1) * (1 - x) ** (n - 1) / denominator**2
    bump = n * x * x * mp.exp(-n * x * x)
    bump_prime = 2 * n * x * mp.exp(-n * x * x) * (1 - n * x * x)
    B = step + (1 - step) * bump
    Bprime = step_prime * (1 - bump) + (1 - step) * bump_prime
    R = 1 + B * (x - 1) / (1 + u * u) ** 3
    RX = (Bprime * (x - 1) + B) / (1 + u * u) ** 3
    A3 = RX / x if x else -2 * n / (1 + u * u) ** 3
    A4 = -A3 - mp.mpf(7) * RX * RX / (4 * R)
    A5 = RX * RX / (R * x) if x else mp.mpf(0)
    tree = tree_evaluator()(u, x)
    F0 = -mp.mpf(28) / 25
    Fv = (x - u * u) / 2 + (lam * kappa - n * F0) * x * x + (-mp.mpf(n) / 3 - F0) * u**4
    # The identity S=1-B and expm1 prevent subtraction of two unit
    # numbers before a tiny calibrated field correction is resolved.
    defect = -mp.expm1(-(u**4)) + mp.exp(-(u**4)) * B
    F = Fv + defect * (tree - Fv)
    return kappa * F + t**4 * (A3 * l3 + A4 * l4) / kappa + t**6 * A5 * l5 / kappa**2


@pytest.mark.parametrize(
    "n,kappa,coupling", [(8, 11, s.Rational(2, 7)), (12, 17, s.Rational(3, 19))]
)
def test_full_function_Cauchy_extraction_of_sixth_and_eighth_germs(n, kappa, coupling):
    base = (
        s.Rational(2, 3),
        s.Rational(-3, 5),
        s.Rational(5, 7),
        s.Rational(2, 9),
        s.Rational(4, 11),
    )
    substitutions = dict(zip(jets.VARIABLES, base))
    substitutions.update({jets.N: n, jets.KAPPA: kappa, jets.LAM: coupling})
    with mp.workdps(150):
        radius = mp.mpf(1) / 128
        count = 64
        angles = [2 * mp.pi * j / count for j in range(count)]
        values = [
            full_density(n, kappa, coupling, radius * mp.exp(1j * a), base)
            for a in angles
        ]
        for degree in range(9):
            observed = mp.fsum(
                v * mp.exp(-1j * degree * a) for v, a in zip(values, angles)
            ) / (count * radius**degree)
            expected = (
                number(
                    jets.data()["canonical_density_field_degrees"]
                    .get(degree, 0)
                    .subs(substitutions)
                )
                if degree in (2, 4, 6, 8)
                else mp.mpf(0)
            )
            assert abs(observed - expected) < mp.mpf("1e-100")
        # Independent extraction, not a test of only pre-truncated germs.
        assert abs(values[0]) > 0


@pytest.mark.parametrize("t", [1, -1, s.I, s.I / 2])
def test_actual_full_calibrated_function_and_all_higher_tail(t):
    base = (
        s.Rational(2, 3),
        s.Rational(-3, 5),
        s.Rational(5, 7),
        s.Rational(2, 9),
        s.Rational(4, 11),
    )
    n, kappa, coupling = 1024, 10**800, s.Rational(1, 10**600)
    substitutions = dict(zip(jets.VARIABLES, base))
    substitutions.update({jets.N: n, jets.KAPPA: kappa, jets.LAM: coupling})
    with mp.workdps(2450):
        z = mp.mpc(number(s.re(t)), number(s.im(t)))
        actual = full_density(n, kappa, coupling, z, base)
        terms = {
            degree: number(piece.subs(substitutions)) * z**degree
            for degree, piece in jets.data()["canonical_density_field_degrees"].items()
        }
        correction = actual - terms[2] - terms[4]
        pointwise = number(
            calibration.data()[
                "all_target_field_degrees_above_four_pointwise_W_squared_majorant"
            ]
        )
        assert abs(correction) <= pointwise
        assert abs(correction) > mp.mpf("1e-1600") * abs(z) ** 6
        tail = number(
            calibration.data()["all_even_field_degrees_at_least_ten_Cauchy_majorant"]
        )
        assert abs(correction - terms[6] - terms[8]) <= tail * abs(z) ** 10
        assert abs(actual - full_density(n, kappa, coupling, -z, base)) < mp.mpf(
            "1e-2400"
        )


@pytest.mark.parametrize("seed", range(10))
def test_literal_Lorentz_contraction_cube_caps(seed):
    eta = s.diag(1, -1, -1, -1)
    grad = s.Matrix([s.Rational(((seed + 3 * j) % 11) - 5, 5) for j in range(4)])
    Hessian = s.Matrix(
        4,
        4,
        lambda i, j: s.Rational(((seed + 2 * min(i, j) + 3 * max(i, j)) % 9) - 4, 4),
    )
    raised = eta * grad
    Y = (grad.T * raised)[0]
    Z = (raised.T * Hessian * raised)[0]
    box = s.trace(eta * Hessian)
    L3 = box * Z
    L4 = (raised.T * Hessian * eta * Hessian * raised)[0]
    assert abs(Y) <= 4
    assert abs(L3) <= 64
    assert abs(L4) <= 64
    assert abs(Z * Z) <= 256


@pytest.mark.parametrize("angle", range(8))
def test_full_calibrated_density_on_declared_complex_amplitude_circle(angle):
    base = (
        s.Rational(2, 3),
        s.Rational(-3, 5),
        s.Rational(5, 7),
        s.Rational(2, 9),
        s.Rational(4, 11),
    )
    with mp.workdps(150):
        t = mp.mpf(10) ** 300 * mp.exp(2j * mp.pi * angle / 8)
        observed = full_density(1024, 10**800, s.Rational(1, 10**600), t, base)
        assert abs(observed) < mp.mpf(10) ** 811
        assert abs(observed) > 0


def test_independent_literal_quartic_target():
    j = parent_jets.data()
    gamma = s.Symbol("quartic_gamma", real=True)
    coupling = s.Symbol("quartic_lambda", real=True)
    observed = jets.data()["canonical_density_field_degrees"][4].subs(
        jets.N, jets.KAPPA * gamma
    )
    observed = observed.subs(
        {
            jets.PHI: j["phi"],
            jets.Y: j["X"],
            jets.L3: j["L3"],
            jets.L4: j["L4"],
            jets.LAM: coupling,
        },
        simultaneous=True,
    )
    assert s.expand(observed - j["target_quartic"]) == 0


@pytest.mark.parametrize("n", [6, 8, 12])
def test_full_rational_switch_division_is_removable(n):
    x = s.Symbol("x")
    T = x**n / (x**n + (1 - x) ** n)
    assert s.limit(T / x, x, 0) == 0
    assert s.limit(s.diff(T, x) / x, x, 0) == 0
    assert s.limit(s.diff(n * x * x * s.exp(-n * x * x), x) / x, x, 0) == 2 * n


@pytest.mark.parametrize("angle", range(8))
def test_rational_step_and_bump_complex_caps_in_independent_larger_domain(angle):
    with mp.workdps(100):
        n = 1024
        x = mp.mpf(1) / 4 * mp.exp(2j * mp.pi * angle / 8)
        r = x / (1 - x)
        step = r**n / (1 + r**n)
        prime = n * r ** (n - 1) / (1 - x) ** 2 / (1 + r**n) ** 2
        assert abs(step) <= 2 * mp.mpf(3) ** (-n)
        assert abs(step / x) <= 8 * mp.mpf(3) ** (-n)
        assert abs(prime / x) <= 100 * n * mp.mpf(3) ** (-n)
        # The bump needs its own, much smaller, declared argument domain.
        x = mp.sqrt(mp.mpf(1) / (4 * n)) * mp.exp(2j * mp.pi * angle / 8)
        bump = n * x * x * mp.exp(-n * x * x)
        prime = 2 * n * x * mp.exp(-n * x * x) * (1 - n * x * x)
        assert abs(bump) <= 2 * n * abs(x) ** 2
        assert abs(prime / x) <= 8 * n


def test_retuned_tree_bound_from_literal_coefficient_sum():
    u, x = family.u, family.X
    num, den = s.fraction(s.factor(family.data()["original_retuned_tree_scalar"]))
    assert s.expand(den - 200 * (1 + u * u) ** 12) == 0
    cap = sum(abs(c) for c in s.Poly(num, u, x).coeffs()) / (
        200 * s.Rational(15, 16) ** 12
    )
    assert cap == s.Rational(354778132704900677632, 3243658447265625)
    assert cap < 10**8


@pytest.mark.parametrize("radius,last", [(2, 0), (2, 8), (3, 6), (10, 12)])
def test_Cauchy_even_tail_against_exact_geometric_function(radius, last):
    # f(t)=1/(1-a t^2), |a| rho^2=1/2 gives circle bound M=2.
    a = s.Rational(1, 2 * radius**2)
    true_tail = a ** (last // 2 + 1) / (1 - a)
    assert true_tail > 0
    assert true_tail <= bounds.even_tail(2, radius, last)
    explicit = sum(a**j for j in range(last // 2 + 1, 100))
    assert explicit < true_tail


@pytest.mark.parametrize("W", [0, s.Rational(1, 10), s.Rational(1, 2), 1])
def test_homogeneous_remainder_can_be_integrated_without_volume_constant(W):
    e6, e8, e10 = s.Rational(2, 7), s.Rational(3, 11), s.Rational(5, 13)
    observed = e6 * W**6 + e8 * W**8 + e10 * W**10
    assert observed <= (e6 + e8 + e10) * W * W
    jets_values = [W * s.Rational(j, 20) for j in range(21)]
    assert max(v * v for v in jets_values) <= sum(v * v for v in jets_values)
    assert max(jets_values) == W
    if W == 0:
        assert observed == 0


def test_majorant_exact_rational_inputs_and_mixed_signs():
    p = jets.PHI**6 - 2 * jets.PHI**2 * jets.Y**2 + 3 * jets.L3 * jets.Y
    expected = 1 + 2 * 4**2 + 3 * 64 * 4
    assert bounds.weighted_majorant(p) == expected
    assert bounds.even_tail(Fraction(2, 3), s.Integer(2)) == s.Rational(2, 3) / (
        2**10 * s.Rational(3, 4)
    )
    assert bounds.compose(Fraction(1, 3), s.Rational(2, 7)) == s.Rational(13, 21)


def test_error_ownership_and_actual_common_class():
    d = calibration.data()
    E6 = d["actual_sixth_density_coefficient_majorant"]
    E8 = d["actual_eighth_density_coefficient_majorant"]
    E10 = d["all_even_field_degrees_at_least_ten_Cauchy_majorant"]
    assert E6 == s.Rational(2858465025, 10**1603)
    assert d[
        "full_target_action_minus_quartic_target_common_class_coefficient"
    ] == 21 * (E6 + E8 + E10)
    t = transport.data()
    assert (
        t["previous_quartic_target_action_error_coefficient"]
        == band.data()["combined_tree_action_error_coefficient"]
    )
    assert t["full_analytic_target_stationary_action_match_coefficient"] < s.Rational(
        1, 10**800
    )
    assert t["rank_regular_target_X_absolute_upper_on_real_common_class"] < s.Rational(
        1, 4096
    )
    assert t["mapped_field_Fourier_radius"] == 3
    assert t["mapped_quadratic_source_Fourier_radius"] == 6
    assert t["centered_box_multiplier_radius"] == 38


def test_field_amplitude_circle_is_not_the_clock_or_momentum_cutoff():
    d = holomorphic.data()
    assert d["dimensionless_field_modulus_upper"] == s.Rational(1, 10**100)
    assert d["dimensionless_gradient_invariant_modulus_upper"] == s.Rational(4, 10**200)
    assert d["fixed_kappa"] != 1
    assert (
        d["complex_field_amplitude_radius"]
        != band.data()["centered_box_spectral_radius"]
    )
    assert "not a momentum cutoff" in d["scope"]
    assert "not_inferred" in transport.data()


def test_frontier_is_only_append_and_original_P8_open():
    assert audit.frontier() == previous.frontier()
    assert audit.matching()[:-1] == previous.matching()
    assert len(audit.matching()) == len(previous.matching()) + 1
    assert "NOT_QUANTUM_TARGET" in audit.matching()[-1]["status"]
    assert audit.controls()["original_P8_not_closed"] is True
