"""Exact conservative hard-band two-scalar propagation and phase-map bounds.

The high powers of q are retained. This is a sufficient bound on one chosen
first-order reduced Hamiltonian, not a bound on arbitrary omitted operators.
"""

from functools import cache

import sympy as sp
from p8_m1_control import bounds as free_bounds
from p8_m1_loops.background import exact_finite_parameter
from p8_m1_tree import majorant

from . import reduction as r

LOWER_K = sp.Integer(10)**11
UPPER_K = sp.Integer(10)**12
Q_MAX = sp.Integer(10)**27
HALF_WINDOW = sp.Rational(1, 100)
CS = sp.Rational(42828424431, 2000000)
ROW_Q = sp.Integer(432)
ROW_P = sp.Integer(172)
NORMALIZED_ROW_Y = sp.Integer(4)*10**8
NORMALIZED_ROW_P = sp.Integer(8000)
DUAL_ENERGY_ROW_SQUARED = sp.Integer(4)*10**17
GENERATOR = sp.Integer(10)**114
PHASE_MAP = sp.Integer(10)**84
EPS_MAX = sp.Rational(10)**(-116)
G_NORM = sp.Rational(4, 3)
EVOLUTION_ERROR = sp.Integer(10)**115
BORN_REMAINDER = sp.Integer(10)**229
PULLED_PHASE_ERROR = sp.Integer(10)**116
MTAU_EXAMPLE = sp.Integer(10)**324


def exact_row(chart):
    """CD row with w=-l*Lambda exposed before taking absolute values."""
    d, lam, j, q, h, theta, l = r.DELTA, r.LAMBDA, r.J, r.Q, r.H, r.THETA, r.L
    if chart == "unitary":
        return sp.Matrix([-2-2*d*lam*q/j, 3*d*l*theta/j, d*theta/j+h/q, d*l*lam/j])
    if chart == "gamma":
        return sp.Matrix([-2*d*theta*q/j-2*h, 3*d*l*theta/j, -1/q-d*lam/j, d*l*lam/j])
    raise ValueError("Use unitary or gamma")


def row_checks():
    return {f"{chart}_{index}": r.cd(value)
            for chart in ("unitary", "gamma")
            for index, value in enumerate(exact_row(chart)-r.chart_data(chart)["lensing_row"])}


def row_majorants():
    """Positive Laurent monomial bound with q>=1 and Jfixed>=1/40.

    It sums absolute row coefficients, hence also bounds their Euclidean
    norm. Only the q powers are enlarged: q^-1<=1<=q.
    """
    j_inverse = sp.Symbol("J_inverse", positive=True)
    values = {r.H: 8, r.THETA: 8, r.LAMBDA: 1, r.DELTA: sp.Rational(1, 2),
              r.L: sp.Rational(1, 5), j_inverse: 40}
    output = {}
    for chart in ("unitary", "gamma"):
        row = exact_row(chart).subs(r.J, 1/j_inverse)
        output[chart] = {}
        for label, indices, q_degree in (("coordinate_row_over_q", (0, 1), 1), ("momentum_row", (2, 3), 0)):
            bound = sp.Integer(0)
            for index in indices:
                for term in sp.Add.make_args(sp.expand(row[index])):
                    power = term.as_powers_dict().get(r.Q, sp.Integer(0))
                    if power.is_Integer is not True or power > q_degree:
                        raise ValueError("Unexpected momentum power in the lensing row")
                    coefficient = sp.expand(term/r.Q**power)
                    polynomial = sp.Poly(coefficient, *values, domain=sp.QQ)
                    for powers, number in polynomial.terms():
                        bound += abs(number)*sp.prod(values[variable]**degree for variable, degree in zip(values, powers))
            output[chart][label] = bound
    return output


@cache
def lineage_bounds():
    previous = majorant.domain_checks()
    if sp.Rational(previous["common_boundary_bound_over_q"]) != CS:
        raise ValueError("Pinned exact canonical boundary majorant changed")
    free = free_bounds.build_bounds()
    if sp.Integer(free["q_threshold"]) != 10**20:
        raise ValueError("Pinned free high-q threshold changed")
    return {"common_boundary_bound_over_q": str(CS), "q_threshold": free["q_threshold"],
            "free_energy_ratio": ["7/11", "11/7"],
            "source": "replayed S5.7 and S5.8, including full coupled W, Omega and S"}


def constant_checks():
    rows = row_majorants()
    # ||rP||<=46*172; ||rY||<20*432*q+2*q*CS*46*172.
    actual_y = 20*ROW_Q+2*CS*46*ROW_P
    actual_p = 46*ROW_P
    phase_coefficient_unitary = sp.Rational(4, 3)*92*(1+CS)*10**9
    phase_coefficient_gamma = sp.Rational(8, 3)*20*10**9
    checks = {
        "scalar_window_high_q": LOWER_K**2/16-10**20,
        "fixed_q_at_least_one": LOWER_K**2/4-1,
        "local_and_fixed_q_upper": Q_MAX-(16*UPPER_K)**2,
        "T_fixed_bound": 46**2-4*512,
        "T_fixed_inverse_bound": 20**2-4*100,
        "unitary_coordinate_row": ROW_Q-rows["unitary"]["coordinate_row_over_q"],
        "gamma_coordinate_row": ROW_Q-rows["gamma"]["coordinate_row_over_q"],
        "unitary_momentum_row": ROW_P-rows["unitary"]["momentum_row"],
        "gamma_momentum_row": ROW_P-rows["gamma"]["momentum_row"],
        "normalized_lensing_Y_row": NORMALIZED_ROW_Y-actual_y,
        "normalized_lensing_P_row": NORMALIZED_ROW_P-actual_p,
        "dual_free_energy_row_squared": DUAL_ENERGY_ROW_SQUARED-2*NORMALIZED_ROW_Y**2-NORMALIZED_ROW_P**2,
        "dual_free_energy_row_sqrt": sp.Integer(10)**18-DUAL_ENERGY_ROW_SQUARED,
        "rank_one_energy_generator_coefficient": sp.Integer(3)*10**18-sp.Rational(16, 3)*DUAL_ENERGY_ROW_SQUARED,
        "q_power_seven_halves": sp.Integer(10)**190-Q_MAX**7,
        "finite_band_generator": GENERATOR-sp.Integer(3)*10**113,
        "unitary_phase_map_coefficient": sp.Integer(3)*10**15-phase_coefficient_unitary,
        "gamma_phase_map_coefficient": sp.Integer(3)*10**15-phase_coefficient_gamma,
        "q_power_five_halves": sp.Integer(10)**136-Q_MAX**5,
        "finite_band_phase_map": PHASE_MAP-sp.Integer(3)*10**83,
        "free_energy_norm": G_NORM**2-sp.Rational(11, 7),
        "interaction_exponent_covers_full_window": 1-sp.Rational(11, 7)*2*HALF_WINDOW,
        "small_interaction_exponent": sp.Rational(1, 100)-GENERATOR*EPS_MAX,
        "phase_map_neumann_smallness": sp.Rational(1, 2)-PHASE_MAP*EPS_MAX,
        "reduced_momentum_Hessian_positive": sp.Rational(1, 2)-sp.Rational(8, 3)*EPS_MAX*Q_MAX**2*NORMALIZED_ROW_P**2,
        "low_energy_spatial_parameter": sp.Rational(1, 100)-EPS_MAX*Q_MAX,
        "same_reduced_data_evolution_error": EVOLUTION_ERROR-G_NORM*GENERATOR/(1-sp.Rational(1, 100)),
        "first_Born_remainder": BORN_REMAINDER-G_NORM*GENERATOR**2/(2*(1-sp.Rational(1, 100))),
        "same_pulled_old_phase_error": PULLED_PHASE_ERROR-6*GENERATOR-8*PHASE_MAP,
    }
    return {key: sp.sympify(value) for key, value in checks.items()}


def validate_constants():
    if any(sp.simplify(value) != 0 for value in row_checks().values()):
        raise ValueError("An exact CD lensing row identity failed")
    values = constant_checks()
    if any(value.is_nonnegative is not True for value in values.values()):
        raise ValueError("A scalar-band sufficient inequality failed")
    # The inverse-T bound is non-strict at the coarse endpoints; all uses
    # below accept <=20. Exact old alpha bounds are in fact strict there.
    return values


def estimates(c_c_abs, mtau):
    """Conditional matching input only; Floats and unproved finite values fail."""
    c_c_abs, mtau = map(exact_finite_parameter, (c_c_abs, mtau))
    if c_c_abs.is_nonnegative is not True or mtau.is_positive is not True:
        raise ValueError("Need nonnegative |cC| and positive M*tau")
    epsilon = c_c_abs/mtau**2
    if (EPS_MAX-epsilon).is_nonnegative is not True:
        raise ValueError("Scalar sufficient gate requires |cC|/(M*tau)^2<=10^-116")
    return {"epsilon_upper": epsilon, "interaction_exponent_upper": GENERATOR*epsilon,
            "relative_same_reduced_data_evolution_error": EVOLUTION_ERROR*epsilon,
            "relative_first_Born_remainder": BORN_REMAINDER*epsilon**2,
            "branch_phase_map_norm": PHASE_MAP*epsilon,
            "chosen_branch_phase_inverse_norm": 1/(1-PHASE_MAP*epsilon),
            "relative_same_pulled_old_phase_error": PULLED_PHASE_ERROR*epsilon,
            "spatial_higher_derivative_parameter": Q_MAX*epsilon}


def controls():
    return {"tensor_epsilon_gate_is_not_this_scalar_gate": GENERATOR*sp.Rational(10)**(-30)-sp.Rational(1, 100),
            "omitting_momentum_enhancement_changes_generator": Q_MAX**7-1,
            "zero_matching_coefficient_has_zero_bound": estimates(0, 1)["relative_same_reduced_data_evolution_error"]}


def report():
    checks = validate_constants()
    return {"domain": {"window": "abs(t-t0)<=ell0/100; ell0=tau*sqrt(1+u0^2); a(t0)=1",
                        "centre_momenta": "10^11<=k_com*ell0<=4*10^12, both coupled scalar channels",
                        "one_chart_per_window": "gamma if abs(x0)<=9/50, unitary otherwise; covers all finite centres",
                        "coarse_geometry": "1/2<=a<=2 and 1/2<=ell/ell0<=2; qfixed,qbar<=10^27; qbar>=10^20",
                        "matching_gate": "independently supplied constant |cC|/(M*tau)^2<=10^-116"},
            "pinned_free_and_boundary_bounds": lineage_bounds(),
            "row_majorants": {chart: {key: str(value) for key, value in row.items()} for chart, row in row_majorants().items()},
            "exact_sufficient_margins": {key: str(value) for key, value in checks.items()},
            "fixed_unit_normalized_lensing_rows": {"Y_row": "4*10^8*qfixed", "P_row": "8000", "dual_energy_row_squared": "4*10^17*qfixed"},
            "perturbation_generator_bound": "||A1||_E <= 3*10^18*qfixed^(7/2) < 10^114; full connection retained in A0",
            "branch_phase_map_bound": "||Fmap||_E <= 3*10^15*qfixed^(5/2) < 10^84",
            "reduced_momentum_Hessian": "eigenvalues >1/2 on this finite band under the sufficient coefficient gate",
            "same_reduced_data_error": "10^115*epsilon, epsilon=|cC|/(M*tau)^2",
            "first_Born_remainder": "10^229*epsilon^2, of the specified reduced Hamiltonian only",
            "same_chosen_pulled_old_phase_error": "10^116*epsilon, for X_old=(I+epsilon0*Fmap)X_red treated as a chosen pullback",
            "conditional_example": {"assumption": "independently supplied |cC|<=1 and M*tau=10^324",
                                     **{key: str(value) for key, value in estimates(1, MTAU_EXAMPLE).items()}},
            "limits": ["the chosen finite-order pullback is only formally matched to the original action through first order",
                       "equal reduced data, equal pulled old phase data and equal physical metric/lapse/shift data are distinct",
                       "no exact higher-derivative branch approximation or unknown O(epsilon^2) action bound",
                       "no global fixed-comoving or unbounded momentum estimate; high powers of q are not removed",
                       "no finite matching value or complete loop/nonlinear/quantum causality conclusion"]}
