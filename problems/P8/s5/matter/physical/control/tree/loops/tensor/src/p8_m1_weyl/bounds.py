"""Exact sufficient, centre-data finite-window tensor perturbation bounds.

These estimates concern the chosen order-reduced two-data representative,
not arbitrary four-data solutions of the resummed fourth-order equation.
"""

from functools import cache

import sympy as sp
from p8_m1_loops.background import exact_finite_parameter

from . import reduction as r

LOWER_K = sp.Integer(10)**11
UPPER_K = sp.Integer(10)**12
K_MAX = 8*UPPER_K
Q_MAX = K_MAX**2
HALF_WINDOW = sp.Rational(1, 100)
EPS_MAX = sp.Rational(10)**(-30)
PERTURBATION_NORM = sp.Integer(10)**16
EVOLUTION_ERROR = sp.Integer(10)**15
BORN_REMAINDER = sp.Integer(10)**30
ORIGINAL_RESIDUAL = sp.Integer(10)**38
G_NORM = sp.Rational(4, 3)
MTAU_EXAMPLE = sp.Integer(10)**324


def exponential_upper(argument, order=12):
    """Positive Taylor polynomial and geometric remainder, exact rational."""
    argument = sp.Rational(argument)
    if argument < 0 or argument >= order+2:
        raise ValueError("Exponential enclosure requires 0<=argument<order+2")
    return (sum(argument**j/sp.factorial(j) for j in range(order+1))
            + argument**(order+1)/sp.factorial(order+1)/(1-argument/(order+2)))


@cache
def compact_extrema():
    """Exhaustive exact critical-point/endpoint checks on compact intervals."""
    v = sp.Symbol("v", real=True)
    specs = {
        "Hdot": (4-8*v, v, 0, 1, 4),
        "B": (-56*r.X+96*r.X**3, r.X, -1, 1, 40),
        "friction_F": (-448*r.X+768*r.X**3, r.X, -1, 1, 320),
        "physical_speed_shift": (-64+128*v, v, 0, 1, 64),
        "canonical_G_non_q": (2688*v-4608*v**2, v, 0, 1, 1920),
    }
    result = {}
    for name, (polynomial, variable, lower, upper, target) in specs.items():
        derivative = sp.diff(polynomial, variable)
        roots = sp.solve(derivative, variable)
        if len(roots) != sp.degree(derivative, variable):
            raise ValueError("Expected the complete distinct stationary root list")
        points = [sp.Integer(lower), sp.Integer(upper)]
        for root in roots:
            if sp.simplify(derivative.subs(variable, root)) != 0:
                raise ValueError("A stationary root failed exact substitution")
            if root.is_real is not True:
                raise ValueError("This polynomial-extrema proof expects real roots")
            if lower <= root <= upper:
                points.append(root)
        values = [sp.simplify(polynomial.subs(variable, point)) for point in points]
        margins = [sp.simplify(target-abs(value)) for value in values]
        if any(margin.is_nonnegative is not True for margin in margins) or not any(margin == 0 for margin in margins):
            raise ValueError("A sharp compact coefficient bound failed")
        result[name] = {"polynomial": str(polynomial), "variable": str(variable),
                        "interval": [str(lower), str(upper)], "absolute_sup": str(target),
                        "derivative": str(derivative), "stationary_roots": list(map(str, roots)),
                        "candidate_points": list(map(str, points)), "candidate_values": list(map(str, values)),
                        "nonnegative_margins": list(map(str, margins))}
    return result


@cache
def constant_checks():
    norm_bound = 2560+4096*UPPER_K+sp.Rational(122880)/LOWER_K
    full_duration = 2*HALF_WINDOW
    return {
        "scale_factor_coarse_upper": 2-1/(1-4*HALF_WINDOW/(1-HALF_WINDOW)),
        "ell_ratio_coarse_lower": 1-HALF_WINDOW-sp.Rational(1, 2),
        "baseline_positive_W": LOWER_K**2/16-60,
        "baseline_derivative_energy_rate": LOWER_K**2/16-132,
        "baseline_energy_ratio": sp.Rational(11, 7)-exponential_upper(sp.Rational(4, 11)),
        "baseline_energy_norm": G_NORM**2-sp.Rational(11, 7),
        "perturbation_norm_including_k_enhancement": PERTURBATION_NORM-norm_bound,
        "low_energy_higher_derivative_parameter": sp.Rational(1, 100)-EPS_MAX*Q_MAX,
        "positive_physical_two_derivative_speed": sp.Rational(1, 2)-256*EPS_MAX,
        "positive_y_action_kinetic_coefficient": sp.Rational(1, 2)-1024*EPS_MAX,
        "extra_energy_exponent": sp.Rational(1, 100)-PERTURBATION_NORM*EPS_MAX/25,
        "total_energy_ratio": sp.Rational(5, 3)-1/(1-sp.Rational(4, 11)-sp.Rational(1, 100)),
        "perturbed_energy_norm": G_NORM**2-sp.Rational(5, 3),
        "Duhamel_evolution_error": EVOLUTION_ERROR-G_NORM**2*PERTURBATION_NORM*full_duration,
        "first_Born_remainder": BORN_REMAINDER-G_NORM**3*(PERTURBATION_NORM*full_duration)**2/2,
        "metric_to_old_canonical_phase_norm": 4-2-(1+48/LOWER_K)**2,
        "original_equation_residual": ORIGINAL_RESIDUAL-8*sp.Integer(20000000000)*Q_MAX,
    }


def validate_constants():
    compact_extrema()
    values = constant_checks()
    if any(value.is_positive is not True for value in values.values()):
        raise ValueError("A sufficient tensor-window inequality failed")
    return values


def estimates(c_c_abs, mtau):
    """Input is an independently supplied bound, never inferred from running."""
    c_c_abs, mtau = map(exact_finite_parameter, (c_c_abs, mtau))
    if c_c_abs.is_nonnegative is not True or mtau.is_positive is not True:
        raise ValueError("Need a nonnegative coefficient bound and positive M*tau")
    epsilon = c_c_abs/mtau**2
    if (EPS_MAX-epsilon).is_nonnegative is not True:
        raise ValueError("This sufficient tensor gate requires |cC|/(M*tau)^2<=10^-30")
    return {"epsilon_upper": epsilon,
            "relative_energy_norm_evolution_error": EVOLUTION_ERROR*epsilon,
            "relative_first_Born_remainder": BORN_REMAINDER*epsilon**2,
            "relative_original_equation_residual": ORIGINAL_RESIDUAL*epsilon**2,
            "higher_derivative_spatial_parameter": epsilon*Q_MAX}


def residual_majorants():
    """Polynomial l1 bounds using |ell0^(n+1) H^(n)|<=2^(n+3)n!.

    Records bound the gamma and velocity coefficients at each extra beta
    power.  q remains a positive variable until the final q<=Q_MAX step.
    """
    jet_bounds = {jet: sp.Integer(2)**(n+3)*sp.factorial(n) for n, jet in enumerate(r.H_JETS)}
    result = {}
    for power, record in r.residual_coefficients().items():
        result[power] = {}
        for component, coefficient in record.items():
            polynomial = sp.Poly(sp.expand(coefficient), *r.H_JETS, r.Q)
            majorant = 0
            for exponents, value in polynomial.terms():
                monomial = abs(value)*r.Q**exponents[-1]
                for jet, exponent in zip(r.H_JETS, exponents[:-1]):
                    monomial *= jet_bounds[jet]**exponent
                majorant += monomial
            result[power][component] = sp.expand(majorant)
    return result


def residual_majorant_checks():
    record = residual_majorants()
    target = {
        0: {"gamma": 102400*r.Q, "velocity": 1024000+10240*r.Q},
        1: {"gamma": 27525120*r.Q+65536*r.Q**2, "velocity": 275251200+1310720*r.Q},
        2: {"gamma": 1677721600*r.Q, "velocity": 16777216000},
    }
    identities = {f"coefficient_{power}_{component}": sp.expand(value-target[power][component])
                  for power, entry in record.items() for component, value in entry.items()}
    # With q>=1, |epsilon|<=1, |epsilon|q<=1, divide gamma coefficients
    # by sqrt(q), and bound the full phase forcing by the indicated C*q.
    summed = 1136640+304152576+18454937600
    identities["residual_phase_majorant_arithmetic"] = sp.Integer(summed)-18760226816
    if any(value != 0 for value in identities.values()):
        raise ValueError("The exact fourth-order residual majorant changed")
    if summed >= 20000000000:
        raise ValueError("Residual phase majorant exceeds its conservative bound")
    return identities


def report():
    checks = validate_constants()
    residual_majorant_checks()
    example = estimates(1, MTAU_EXAMPLE)
    return {"domain": {"centre_data": "same physical gamma(t0), gamma_dot(t0), hence same old Y(t0), P(t0)",
                        "window": "abs(t-t0)<=ell0/100; ell0=tau*sqrt(1+u0^2); a(t0)=1",
                        "centre_momenta": "10^11<=k_com*ell0<=4*10^12, both TT polarizations",
                        "coarse_window_geometry": "1/2<=a<=2 and 1/2<=ell/ell0<=2",
                        "matching_input": "constant cC; epsilon=abs(cC)/(M*tau)^2<=10^-30",
                        "not_initial_data": "not equal redefined y data or a newly chosen corrected quantum vacuum"},
            "exact_sufficient_margins": {key: str(value) for key, value in checks.items()},
            "exact_compact_coefficient_extrema": compact_extrema(),
            "perturbation_energy_norm_generator_bound": "ell0*(abs(f)+abs(DeltaW)/sqrt(W0))<=10^16*abs(cC)/(M*ell0)^2",
            "retained_frequency_enhancement": str(2560+4096*UPPER_K+sp.Rational(122880)/LOWER_K),
            "relative_energy_ratio_about_centre": ["3/5", "5/3"],
            "relative_same_data_evolution_error": "10^15*epsilon",
            "relative_first_Born_remainder": "10^30*epsilon^2; Taylor remainder of the chosen second-order representative only",
            "original_fourth_order_residual_bound": "ell0^2*abs((M/2)*a^(3/2)*Efull[gamma_red])<=10^38*epsilon^2*N_t(Y,P)",
            "residual_coefficient_majorants": {str(power): {key: str(value) for key, value in entry.items()}
                                                for power, entry in residual_majorants().items()},
            "conditional_example": {"assumption": "independently supplied abs(cC)<=1 and M*tau=10^324",
                                     **{key: str(value) for key, value in example.items()}},
            "limits": ["no actual four-data fourth-order approximation/branch existence theorem",
                       "no unknown matching coefficient or omitted operator error is bounded",
                       "no scalar sector, corrected quantum state, complete loop or front-velocity verdict"]}
