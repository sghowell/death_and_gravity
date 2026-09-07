"""Exact source-cost envelopes, not a numerically solved band optimizer.

The actual loading map is proved elsewhere in this child. This module
replays Fourier/operator arithmetic conditional on its verified Gramian
and source-construction constants; passing numbers does not verify a map.
All frequencies use physical g proper time u=T/tau, not a canonical source
weight or an instantaneous heavy-mode eigenvalue.
"""

from fractions import Fraction
from functools import cache

import sympy as sp

LEFT = -sp.Rational(1, 50)
RIGHT = -sp.Rational(1, 100)
LENGTH = RIGHT-LEFT
OMEGA_MAX = sp.Integer(100)
DELTA_MAX = sp.Rational(1, 10**9)
LOADING_NORM_UPPER = sp.Rational(4, 5)


def rational(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, sp.Rational, Fraction)):
        raise TypeError(f"{name} must be an exact integer or rational")
    if isinstance(value, Fraction):
        return sp.Rational(value.numerator, value.denominator)
    return sp.Rational(value)


def nonnegative(value, name):
    value = rational(value, name)
    if value < 0:
        raise ValueError(f"{name} must be nonnegative")
    return value


def positive(value, name):
    value = rational(value, name)
    if value <= 0:
        raise ValueError(f"{name} must be positive")
    return value


def frequency(omega):
    omega = nonnegative(omega, "Omega")
    if omega > OMEGA_MAX:
        raise ValueError("The calibrated frequency domain is 0<=Omega<=100")
    return omega


def order_number(order):
    order = nonnegative(order, "Taylor order")
    if order.q != 1:
        raise ValueError("Taylor order must be an integer")
    return int(order)


def lowpass_bounds(omega=OMEGA_MAX):
    """Bounds for the actual zero-extended source's Fourier energy.

    The rational bounds use pi>3. At Omega=0 the lowpass operator is zero;
    for positive Omega the displayed rational norm ceiling is strict.
    """
    omega = frequency(omega)
    q = omega*LENGTH/3
    return {"omega": omega, "length": LENGTH,
            "exact_trace": omega*LENGTH/sp.pi,
            "lowpass_norm_upper": q, "tail_coercivity_lower": 1-q,
            "band_energy_fraction_upper": q,
            "tail_energy_fraction_lower": 1-q}


def sinc_remainder(omega=OMEGA_MAX, order=4, loading_norm_upper=LOADING_NORM_UPPER):
    """Continuous finite-rank Fourier error, not computed source moments."""
    omega, order = frequency(omega), order_number(order)
    loading_norm_upper = nonnegative(loading_norm_upper, "loading norm upper bound")
    error = (omega*LENGTH)**(2*order+3)/(3*sp.factorial(2*order+3))
    return {"order": order, "operator_error_upper": error,
            "band_gramian_error_upper": loading_norm_upper**2*error,
            "moment_degrees_required": 2*order,
            "actual_band_moments_computed": False}


def moment_weights(omega=OMEGA_MAX, order=4):
    """Weights of m_(2n-j) m_j^T, m_j=int_J (u-midpoint)^j h(u) du.

    Actual moment vectors are not computed here. All coefficients include
    the physical Fourier factor 1/pi; one must not count it a second time.
    """
    omega, order = frequency(omega), order_number(order)
    return tuple((2*n-j, j,
                  omega/sp.pi*(-1)**(n+j)*omega**(2*n)*sp.binomial(2*n, j)
                  /sp.factorial(2*n+1))
                 for n in range(order+1) for j in range(2*n+1))


def minimum_cost_envelope(omega, target_norm, gramian_lower, gramian_upper,
                          source_budget, additional_norm_lower=0):
    """Certified envelope for an infimum, with a strict optimizer budget.

    Require g_- I<=C C*<=g_+ I and any independently proved additional
    minimum-source-norm lower bound. The sufficient strict-budget test
    guarantees that the unconstrained L2 spectral optimizer is feasible
    and that its value equals the H0^2 infimum. Its actual value/source is
    not evaluated, and attainment in H0^2 is not asserted.
    """
    band = lowpass_bounds(omega)
    target = nonnegative(target_norm, "target norm")
    lower, upper = positive(gramian_lower, "Gramian lower"), positive(gramian_upper, "Gramian upper")
    budget = nonnegative(source_budget, "source L2 budget")
    extra = nonnegative(additional_norm_lower, "additional source-norm lower")
    if lower > upper or (target == 0 and extra > 0):
        raise ValueError("Inconsistent Gramian or zero-target norm bounds")
    norm_floor = max(target**2/upper, extra**2)
    norm_ceiling = target**2/lower
    if norm_floor > norm_ceiling:
        raise ValueError("The supplied minimum-norm bounds are inconsistent")
    optimizer_norm_ceiling = norm_ceiling/band["tail_coercivity_lower"]
    slack = budget**2-optimizer_norm_ceiling
    if target > 0 and slack <= 0:
        raise ValueError("The sufficient strict spectral-optimizer budget is not verified")
    return {"source_minimum_norm_squared_lower": norm_floor,
            "source_minimum_norm_squared_upper": norm_ceiling,
            "spectral_optimizer_norm_squared_upper": optimizer_norm_ceiling,
            "strict_squared_budget_margin": slack,
            "tail_cost_infimum_lower": band["tail_coercivity_lower"]*norm_floor,
            "tail_cost_infimum_upper": norm_ceiling,
            "optimizer_numerically_computed": False,
            "H0_2_attainment_claimed": False}


def constructive_bounds(omega, target_norm, source_l2_upper, source_second_derivative_l2_upper):
    """An actual H0^2 right-inverse control, not the spectral optimizer."""
    omega = frequency(omega)
    target = nonnegative(target_norm, "target norm")
    l2 = positive(source_l2_upper, "right-inverse L2 norm")*target
    d2 = positive(source_second_derivative_l2_upper, "right-inverse second-derivative norm")*target
    tail = l2**2 if omega == 0 else min(l2**2, d2**2/omega**4)
    return {"source_L2_upper": l2, "source_L1_upper": l2/10,
            "source_second_derivative_L2_upper": d2,
            "source_sup_Sobolev_upper": d2/3000,
            "tail_energy_upper": tail,
            "claims_mostly_low_band": False}


def finite_delta_error(delta, target_norm, source_l2, loading_error=0,
                       target="limit", momentum_squared=1):
    """Fixed-source output bound in the actual P_delta endpoint frame.

    target='limit' compares to Lv; 'prepared' compares to the actual
    analytic prepared H_delta^+ v. No physical map or past state is reset.
    """
    delta = positive(delta, "delta")
    kk = rational(momentum_squared, "K")
    if delta > DELTA_MAX or not 1 <= kk <= 4:
        raise ValueError("Require 0<delta<=10^-9 and 1<=K<=4")
    if target not in ("limit", "prepared"):
        raise ValueError("Target must be 'limit' or 'prepared'")
    norm = nonnegative(target_norm, "target norm")
    source = nonnegative(source_l2, "source L2")
    error = nonnegative(loading_error, "normalized loading error")
    light = 8600 if target == "limit" else 8400
    return 42*error+delta*(light*norm+12600000*source)


def source_units(planck_mass, time_scale, omega=0):
    """Pi(T)=M^2/tau^2 sigma(T/tau); Fourier transform has no prefactor."""
    mass, tau = positive(planck_mass, "M"), positive(time_scale, "tau")
    omega = frequency(omega)
    return {"physical_frequency_edge": omega/tau,
            "physical_L2_squared_multiplier": mass**4/tau**3,
            "physical_Fourier_energy_multiplier": mass**4/tau**3,
            "physical_L1_multiplier": mass**2/tau,
            "physical_sup_multiplier": mass**2/tau**2}


@cache
def calibration():
    """Arithmetic bridge; companion modules certify the primitive inputs."""
    lower, upper = sp.Rational(9, 4*10**12), sp.Rational(16, 25)
    envelope = minimum_cost_envelope(100, 1, lower, upper, 10**6, 19)
    hermite = constructive_bounds(100, 1, 2*10**8, 4*10**14)
    error = finite_delta_error(sp.Rational(1, 10**17), 1, 10**6)
    return {"gramian_lower_input": lower, "gramian_upper_input": upper,
            "regular_light_norm_lower_input": sp.Integer(19),
            "regular_even_norm_lower_input": sp.Integer(900),
            "band": lowpass_bounds(), "sinc": sinc_remainder(),
            "unit_target_budget": sp.Integer(10**6), "unit_target_envelope": envelope,
            "even_tail_cost_lower": sp.Rational(2, 3)*900**2,
            "Hermite_control": hermite,
            "delta_example": sp.Rational(1, 10**17), "limit_target_error": error,
            "actual_prepared_error": finite_delta_error(sp.Rational(1, 10**17), 1, 10**6,
                                                         target="prepared")}


@cache
def checks():
    x, y = sp.symbols("x y", real=True)
    omega = sp.Integer(100)
    weights = sum(coefficient*x**left*y**right for left, right, coefficient in moment_weights(omega))
    direct = omega/sp.pi*sum((-1)**n*omega**(2*n)*(x-y)**(2*n)/sp.factorial(2*n+1)
                             for n in range(5))
    r1, r2, lam, target = sp.symbols("r1 r2 lambda target", positive=True)
    denominator = r1+r2+2*lam
    sigma = sp.Matrix([(r2+lam)*target/denominator, (r1+lam)*target/denominator])
    sigma_prime = sigma.diff(lam)
    a, b, c, q, budget = sp.symbols("a b c q budget", positive=True)
    n2 = a*a+b*b
    tail_root, band_root = (b*c-a*q)/n2, (a*c+b*q)/n2
    residuals = {
        "finite_rank_moment_kernel": sp.expand(weights-direct),
        "physical_source_L1_from_L2": sp.sqrt(LENGTH)-sp.Rational(1, 10),
        "two_coordinate_KKT_constraint": sp.factor(sum(sigma)-target),
        "KKT_budget_monotonicity": sp.factor(sp.diff(sigma.dot(sigma), lam)
                                              +2*(sigma_prime.T*sp.diag(r1+lam, r2+lam)*sigma_prime)[0]),
        "dual_moment_constraint": sp.factor(a*band_root+b*tail_root-c),
        "dual_moment_sphere": sp.factor((band_root**2+tail_root**2-budget**2)
                                         .subs(q*q, n2*budget**2-c*c)),
        "dual_moment_positive_root_threshold": sp.factor(b*b*c*c-a*a*(n2*budget**2-c*c)
                                                          -n2*(c*c-a*a*budget**2)),
        "exact_prepared_target_coefficient": sp.Integer(42)*200-8400,
    }
    cal = calibration()
    margins = {
        "positive_band_coercivity": cal["band"]["tail_coercivity_lower"],
        "strict_unit_optimizer_budget": cal["unit_target_envelope"]["strict_squared_budget_margin"],
        "positive_unit_cost_floor": cal["unit_target_envelope"]["tail_cost_infimum_lower"],
        "finite_sinc_operator_error_below_10_minus_8": sp.Rational(1, 10**8)-cal["sinc"]["operator_error_upper"],
        "finite_delta_example_below_one_thousandth": sp.Rational(1, 1000)-cal["limit_target_error"],
        "distinct_actual_prepared_target_error": cal["limit_target_error"]-cal["actual_prepared_error"],
    }
    if any(value != 0 for value in residuals.values()):
        raise ValueError("An exact spectral-cost identity failed")
    if any(value.is_positive is not True for value in margins.values()):
        raise ValueError("A continuous spectral-cost margin failed")
    return {"residuals": residuals, "margins": margins}
