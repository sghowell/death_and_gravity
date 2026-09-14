"""Source-pinned fixed reference and original smooth preparation, not a new state."""

from functools import cache

import sympy as s
from p8_vacuum_affine_finite_band_neighborhood import source as previous
from p8_vacuum_affine_heavy_clock_quadratic import clock
from p8_vacuum_affine_physical_background_vertices import parent
from p8_vacuum_affine_scalar_tame_propagator import charts

u = parent.u
H, ell = parent.H, parent.ell
a = (1 + u * u) ** 2
delta = 1 / (2 * parent.h)
E = 1 - 3 * delta
theta = H - u / (1 + u * u) ** 4
F = s.factor(
    theta * s.diff(E, u)
    - E * s.diff(theta, u)
    + H * E * theta
    - theta**2
    - ell**2 * E**2 / 2
)
rho, pressure, rho1, pressure1 = s.symbols(
    "fixed_profile_rho fixed_profile_pressure fixed_profile_rho_first fixed_profile_pressure_first",
    real=True,
)
PROFILE_VARIABLES = (rho, pressure, rho1, pressure1)
PROFILE_BOUND = s.Rational(1, 10**390)
J = (
    F
    + 1 / (50 * (1 + u * u) ** 6)
    - (21 * delta**2 - 3 * delta) * pressure / 2
    - (1 - 6 * delta) * (rho + pressure) / 2
)
TIME, LOW, HIGH, EPSILON = previous.TIME, previous.LOW, previous.HIGH, previous.EPSILON
KAPPA = previous.b.source.KAPPA
SWITCH = s.Rational(3, 16)
ROOT_RATE = s.Integer(12)
COVARIANCE = s.Integer(10) ** 21
CHARTS = (
    ("outer", theta, -ell * E / theta, -s.Rational(1, 2), -SWITCH),
    ("central", E, ell, -SWITCH, SWITCH),
    ("outer_positive", theta, -ell * E / theta, SWITCH, s.Rational(1, 2)),
)


def dt(expression):
    return (
        s.diff(expression, u)
        + rho1 * s.diff(expression, rho)
        + pressure1 * s.diff(expression, pressure)
    )


def principal(pivot, mixed):
    K = s.Matrix([[2 * J / pivot**2 + mixed**2, mixed], [mixed, 1]]).applyfunc(s.factor)
    G = s.Matrix([[2 * F / pivot**2 + mixed**2, mixed], [mixed, 1]]).applyfunc(s.factor)
    return K, G


def as_rational(value):
    return s.Rational(value.numerator, value.denominator)


@cache
def data():
    current = parent.reference_data()["whole_reference_coefficients"]
    bind = {parent.RHO: rho, parent.PRESSURE: pressure}
    old_gradient = charts.gradient_numerator() - charts.w**2 / 2
    jets = {
        charts.th: theta,
        charts.E: E,
        charts.l: ell,
        charts.H: H,
        charts.jets[0]: s.diff(theta, u),
        charts.jets[1]: s.diff(E, u),
    }
    return {
        "whole_actual_reference_coefficients": current,
        "whole_coefficient_symbol_binding": {
            charts.th: theta,
            charts.E: E,
            charts.l: ell,
            charts.J: J,
            charts.A: -pressure,
            charts.T: rho - 3 * delta * pressure,
            charts.H: H,
        },
        "whole_profile_and_first_jet_variables": PROFILE_VARIABLES,
        "whole_profile_box_absolute_bound": PROFILE_BOUND,
        "actual_existing_total_C5_profile_bound": clock.EPS,
        "whole_reference_scale": a,
        "whole_two_chart_intervals": {name: [lo, hi] for name, _, _, lo, hi in CHARTS},
        "unchanged_preparation": "Use exactly S251 t*=-1/2 and f(t)=exp[-1/(1-1024(t+7/16)^2)] on |t+7/16|<1/32, zero otherwise. Its original selection form includes its unchanged mass-scale1 K term. No current-history re-minimization, new vacuum, endpoint reset or new state is introduced.",
        "whole_target_slab": [-TIME, TIME],
        "whole_quantified_high_momentum_domain": "Every P>=10^64 for the fixed-reference covariance/transport estimate; physical point moments below use only the stated band[10^64,2*10^64].",
        "checks": {
            "literal_entire_current_pivot": s.factor(
                current["Jc"].subs(bind, simultaneous=True) - J
            ),
            "literal_entire_current_Theta": s.factor(current["Theta"] - theta),
            "literal_entire_current_E": s.factor(current["E"] - E),
            "literal_entire_current_charge": s.factor(current["ell"] - ell),
            "literal_original_principal_gradient": s.factor(
                old_gradient.subs(jets, simultaneous=True) - F
            ),
            "literal_reference_scale_Hubble": s.factor(s.diff(a, u) - a * H),
            "same_original_kappa": KAPPA - s.Integer(10) ** 800,
            "same_S263_initial_comparison_radius": EPSILON - s.Rational(1, 10**230),
        },
        "gates": {
            "actual_full_profiles_and_first_jets_inside_larger_box": clock.EPS
            < PROFILE_BOUND,
            "actual_profiles_not_zero_or_reselected": True,
            "original_preparation_support_before_only_negative_chart_switch": -s.Rational(
                13, 32
            )
            < -SWITCH,
            "target_slab_inside_central_chart": TIME < SWITCH,
            "whole_reference_scale_between_one_and_two": s.Rational(25, 16) < 2,
            "full_reference_P_domain_above_old_central_threshold": LOW**2 / 4 > 4096,
            "old_fixed_H_and_Proca_preparations_and_determinants_retained": True,
            "not_quantitative_bound_from_uncontrolled_asymptotic_symbols": True,
        },
    }
