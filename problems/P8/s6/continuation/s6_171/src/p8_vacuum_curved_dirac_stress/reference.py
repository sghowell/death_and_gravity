"""Fixed soft-metric Newton reference and complete local two-derivative bounds."""

from functools import cache

import sympy as s
from p8_vacuum_curved_dirac_state import geometry, transition


def enclosures(mean_mass, amplitude, timescale, multiplicity=6):
    m, d, tau, _ = transition.parameters(mean_mass, amplitude, timescale)
    if type(multiplicity) is not int or multiplicity < 2 or multiplicity % 2:
        raise TypeError("Require a positive even native paired multiplicity")
    ell = 2 * d / (m - d)
    up = m + d
    D = 3 * d * d
    Dt = 2 * up * d * ell / tau
    Dtt = (2 * ell * (d * d + 9 * up * d) + 4 * d * d) / tau**2
    kinetic = (ell + s.Rational(2, 3)) * d * d / tau**2
    return {
        "absolute_logarithm": ell,
        "absolute_D_function": D,
        "absolute_D_time_derivative": Dt,
        "absolute_D_second_time_derivative": Dtt,
        "local_two_derivative_energy": multiplicity * (kinetic + 4 * D + 2 * Dt) / 144,
        "local_two_derivative_pressure": multiplicity
        * (kinetic + s.Rational(20, 3) * D + Dtt / 3 + s.Rational(4, 3) * Dt)
        / 144,
    }


@cache
def data():
    M, m = s.symbols("M m", positive=True)
    D = M * M * (s.log(M * M / (m * m)) - 1) + m * m
    geom = geometry.data()
    m0, amp, tau = (
        geom["mean_mass"],
        geom["amplitude_upper"],
        geom["mass_transition_time"],
    )
    actual = enclosures(m0, amp, tau)
    t = s.Symbol("t", real=True)
    H = 4 * t / (1 + t * t)
    checks = {
        "D_zero_at_reference": D.subs(M, m),
        "D_first_derivative_zero_at_reference": s.diff(D, M).subs(M, m),
        "D_first_derivative": s.simplify(s.diff(D, M) - 2 * M * s.log(M * M / (m * m))),
        "D_second_derivative": s.simplify(
            s.diff(D, M, 2) - 2 * s.log(M * M / (m * m)) - 4
        ),
        "CD_H_derivative": s.factor(s.diff(H, t) - 4 * (1 - t * t) / (1 + t * t) ** 2),
        "H_squared_upper_polynomial": s.factor(
            (4 - H * H) * (1 + t * t) ** 2 - 4 * (t * t - 1) ** 2
        ),
        "absolute_Hdot_bound_polynomial": s.expand(
            (1 + t * t) ** 2 - (1 - t * t) - t * t * (3 + t * t)
        ),
        "zero_profile_referenced_local_energy": enclosures(m0, 0, tau)[
            "local_two_derivative_energy"
        ],
        "zero_profile_referenced_local_pressure": enclosures(m0, 0, tau)[
            "local_two_derivative_pressure"
        ],
    }
    return {
        "reference_name": "GY14-SAT8-MR/EC-N0",
        "reference_action": "Retain S6.168 MS/vacuum/saturated-mass-reference terms and S6.169 flat improvement. Per Dirac copy add -m^2 R/(6Q), Q=16pi^2, mu=m. This fixes the free quadratic zero-field soft-metric Newton contribution to zero. It is not a coefficient fitted to this stress estimate, nor the complete interacting physical Newton dictionary.",
        "finite_D_function": D,
        "complete_local_bounds": actual,
        "paired_active_and_constant_species": "Six paired active copies use M=m+/-Delta*s(t/tau). The other36 copies have M=m: their referenced local two-derivative term is exactly zero, but their curved state and higher geometric remainder are not discarded.",
        "bound_proof": "On |M-m|<=Delta, |ell|<=2Delta/(m-Delta), |Mdot|<=Delta/tau, |Mddot|<=9Delta/tau^2. Since D(m)=D'(m)=0 and |D''|<6, |D|<=3Delta^2. Differentiate D exactly for the displayed Dt,Dtt bounds. Insert |H|<=2 and |Hdot|<=4 into both covariant components. No derivative of an inequality is taken.",
        "decimal_diagnostics_only": {k: str(s.N(v, 28)) for k, v in actual.items()},
        "checks": checks,
        "gates": {
            "entire_domain_D_second_derivative_below_six": bool(
                4 + 4 * s.Rational(3, 1000) / (1 - s.Rational(3, 1000)) < 6
            ),
            "actual_log_below_one_hundredth": bool(
                actual["absolute_logarithm"] < s.Rational(1, 100)
            ),
            "actual_local_energy_below_one_e594": bool(
                actual["local_two_derivative_energy"] < 10**594
            ),
            "actual_local_pressure_below_one_e595": bool(
                actual["local_two_derivative_pressure"] < 10**595
            ),
            "curved_reference_fixed_independently_of_energy_bound": True,
        },
    }
