"""Root-authored physical-source, endpoint-domain and Fourier-cost audit."""
from functools import cache

import pytest
import sympy as sp
from p8_preparation_cost import construction, spectral


@cache
def source_pairing_residual():
    # Generic actual canonical derivative chart, independently differentiated.
    a, b, e, omega, omega_prime, jl, jh, sigma = sp.symbols(
        "A B E omega omega_prime jL jH sigma", real=True)
    x = sp.Matrix(sp.symbols("l v Q w", real=True))
    y = sp.Matrix(sp.symbols("L V R W", real=True))
    matrix = sp.Matrix([[0, 1, 0, 0], [-a, 0, -e+2*omega_prime, 2*omega],
                        [0, 0, 0, 1], [-e, -2*omega, -b, 0]])
    forcing = sp.Matrix([0, jl*sigma, 0, jh*sigma])
    pairing = (x[0]*(y[1]-2*omega*y[2])-(x[1]-2*omega*x[2])*y[0]
               +x[2]*y[3]-x[3]*y[2])
    rate = (sum(sp.diff(pairing, x[i])*(matrix*x+forcing)[i] for i in range(4))
            +sum(sp.diff(pairing, y[i])*(matrix*y)[i] for i in range(4))
            +sp.diff(pairing, omega)*omega_prime)
    fs, fr, w2, volume = sp.symbols("fs fr w2 volume", positive=True)
    actual = rate.subs({jl: volume/(2*fs), jh: -volume*w2/(2*fr)})
    g_mode = y[0]/fs-w2*y[2]/fr
    return sp.expand(actual+volume*sigma*g_mode/2)


def test_full_canonical_pairing_has_the_actual_physical_source_weight():
    assert source_pairing_residual() == 0


@pytest.mark.parametrize("momentum", [1, sp.Rational(5, 2), 4])
@pytest.mark.parametrize("endpoint", [(1, 0, 0, 0), (0, 1, 0, 0),
                                    (0, 0, 1, 0), (0, 0, 0, 1)])
def test_source_and_first_derivative_vanish_at_both_physical_endpoints(momentum, endpoint):
    d = construction.derive()
    x, u, ell = d["x"], d["u"], d["ell"]
    polynomial = construction.hermite_polynomial(endpoint, momentum)
    # Independent evaluation of the source jets as Leibniz sums, including
    # all coefficient derivatives and the physical u versus x scale.
    for side, location in ((0, d["left"]), (1, d["right"])):
        f_jets = [sp.diff(polynomial, x, j).subs(x, side)/ell**j for j in range(6)]
        for source_order in (0, 1):
            value = sum(sp.binomial(source_order, k)
                        *sp.diff(coefficient, u, source_order-k).subs(
                            {u: location, d["K"]: momentum})*f_jets[j+k]
                        for j, coefficient in enumerate(d["source_coefficients"])
                        for k in range(source_order+1))
            assert sp.cancel(value) == 0


@pytest.mark.parametrize("degree", [2, 3, 4, 6])
def test_smooth_domain_can_have_unattained_minimum_cost(degree):
    # For C sigma=int_0^1 sigma, target1, the unique L2 minimizer is1.
    # These exact H0^2 approximants have correct loading and cost strictly
    # above1. The constant optimizer fails the boundary value traces.
    x = sp.Symbol("x", real=True)
    bump = (1-x**degree)**2*(1-(1-x)**degree)**2
    area = sp.integrate(bump, (x, 0, 1))
    source = bump/area
    assert sp.integrate(source, (x, 0, 1)) == 1
    for endpoint in (0, 1):
        assert source.subs(x, endpoint) == 0
        assert sp.diff(source, x).subs(x, endpoint) == 0
    assert sp.integrate(source**2, (x, 0, 1)) > 1
    assert sp.Integer(1).subs(x, 0) != 0


def test_zero_frequency_is_not_a_strict_positive_lowpass_operator():
    result = spectral.lowpass_bounds(0)
    assert result["lowpass_norm_upper"] == 0
    assert result["tail_energy_fraction_lower"] == 1
    assert spectral.sinc_remainder(0)["operator_error_upper"] == 0


def test_fixed_physical_source_units_retain_the_clock_jacobian():
    result = spectral.source_units(3, 2, 100)
    assert result["physical_frequency_edge"] == 50
    assert result["physical_L2_squared_multiplier"] == sp.Rational(81, 8)
    assert result["physical_Fourier_energy_multiplier"] == sp.Rational(81, 8)
    assert result["physical_L1_multiplier"] == sp.Rational(9, 2)


def test_correct_target_frame_changes_the_finite_delta_comparison():
    delta = sp.Rational(1, 10**17)
    limiting = spectral.finite_delta_error(delta, 1, 10**6, sp.Rational(1, 10000))
    prepared = spectral.finite_delta_error(delta, 1, 10**6, sp.Rational(1, 10000), target="prepared")
    assert limiting-prepared == 200*delta
    assert limiting == 42*sp.Rational(1, 10000)+delta*(8600+12600000*10**6)


def test_inconsistent_or_non_strict_budgets_are_not_certified():
    with pytest.raises(ValueError):
        spectral.minimum_cost_envelope(100, 1, 1, 2, 1)
    with pytest.raises(ValueError):
        spectral.minimum_cost_envelope(100, 1, 1, 2, 10, 2)


def test_calibrated_tail_floor_is_not_a_computed_optimizer():
    cal = spectral.calibration()
    assert cal["unit_target_envelope"]["tail_cost_infimum_lower"] == sp.Rational(722, 3)
    assert cal["even_tail_cost_lower"] == 540000
    assert cal["unit_target_envelope"]["strict_squared_budget_margin"] == sp.Rational(10**12, 3)
    assert not cal["unit_target_envelope"]["optimizer_numerically_computed"]
    assert not cal["unit_target_envelope"]["H0_2_attainment_claimed"]
    assert not cal["sinc"]["actual_band_moments_computed"]
