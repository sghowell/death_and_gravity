"""Sharper out-particle energy and a same-operator in/out T00 difference."""

from functools import cache

import sympy as s
from p8_exceptional_vacuum import analytic
from p8_vacuum_flat_dirac_hadamard import symbols
from p8_vacuum_flat_dirac_production import calibration as original

from . import tube


def parameters(mean_mass, amplitude, timescale):
    m, d, tau = map(symbols.rational, (mean_mass, amplitude, timescale))
    floor = s.Rational(99, 100) * m
    if m <= 0 or d < 0 or tau <= 0:
        raise ValueError("Require m,tau>0 and Delta>=0")
    if d / floor >= s.Rational(1, 100) or tau * floor < tube.SCALE:
        raise ValueError("Require Delta/(.99m)<1/100 and .99m*tau>=2^20")
    return m, d, tau, floor


def amplitude_upper(momentum, mean_mass, amplitude, timescale):
    p = symbols.rational(momentum)
    if p < 0:
        raise ValueError("Require nonnegative momentum magnitude")
    _m, d, tau, floor = parameters(mean_mass, amplitude, timescale)
    return tube.TRANSITION_CONSTANT * p * d / (tau**4 * (p * p + floor * floor) ** 3)


def enclosures(mean_mass, amplitude, timescale, multiplicity=6):
    _m, d, tau, floor = parameters(mean_mass, amplitude, timescale)
    if type(multiplicity) is not int or multiplicity < 1:
        raise TypeError("Require a positive native color/flavor multiplicity")
    C = tube.TRANSITION_CONSTANT
    return {
        "mass_floor": floor,
        "uniform_exact_transition_amplitude_upper": C * d / (2 * tau**4 * floor**5),
        "complete_free_out_particle_energy_upper": s.Rational(
            32 * multiplicity, 315 * 9
        )
        * C**2
        * d
        * d
        / (tau**8 * floor**6),
        "uniform_in_out_local_energy_density_difference_upper": s.Rational(
            8 * multiplicity, 3 * 9
        )
        * C
        * d
        / (tau**4 * floor),
    }


@cache
def data():
    x = s.Symbol("x", positive=True)
    out_integral = s.integrate(x**4 / (1 + x * x) ** s.Rational(11, 2), (x, 0, s.oo))
    difference_integral = s.integrate(
        x**3 / (1 + x * x) ** s.Rational(5, 2), (x, 0, s.oo)
    )
    p = original.data()
    m, d, tau = (
        p["quadratic_subsystem_mean_mass"],
        p["profile_amplitude_rational_upper"],
        p["same_profile_time_scale"],
    )
    values = enclosures(m, d, tau)
    old = p["complete_quadratic_transition_and_out_energy_enclosure"]
    C = s.Integer(tube.TRANSITION_CONSTANT)
    N, Delta, time, floor = s.symbols("N Delta tau m0", positive=True)
    out = 32 * N * C * C * Delta * Delta / (315 * s.pi**2 * time**8 * floor**6)
    diff = 8 * N * C * Delta / (3 * s.pi**2 * time**4 * floor)
    checks = {
        "out_energy_radial_integral": out_integral - s.Rational(8, 315),
        "state_difference_radial_integral": difference_integral - s.Rational(2, 3),
        "all_out_energy_factors": s.simplify(
            out
            - 4
            * N
            / s.pi**2
            * C
            * C
            * Delta
            * Delta
            / time**8
            * out_integral
            / floor**6
        ),
        "all_energy_difference_factors": s.simplify(
            diff - 4 * N / s.pi**2 * C * Delta / time**4 * difference_integral / floor
        ),
        "same_clock_transition_scale": tau - s.Rational(1, 10**100),
        "same_mass_profile_amplitude_upper": d - 3 * 10**197,
        "same_mass_reference": m - 10**200,
        "same_protected_mass_floor": values["mass_floor"] - old["strict_mass_lower"],
        "zero_momentum_exact_transition": amplitude_upper(0, m, d, tau),
        "zero_mass_variation_exact_transition": amplitude_upper(1, m, 0, tau),
    }
    return {
        "same_quadratic_operator_mean_mass": m,
        "same_mass_profile_amplitude_rational_upper": d,
        "same_profile_time_scale": tau,
        "active_color_flavor_multiplicity": 6,
        "actual_complete_exact_frame_enclosures": values,
        "new_out_energy_over_named_kappa_scale": values[
            "complete_free_out_particle_energy_upper"
        ]
        / analytic.KAPPA,
        "uniform_local_state_difference_over_named_kappa_scale": values[
            "uniform_in_out_local_energy_density_difference_upper"
        ]
        / analytic.KAPPA,
        "actual_bound_decimal_diagnostics": {
            k: str(s.N(v, 25)) for k, v in values.items()
        },
        "out_energy_proof": "rho_out=4N/(2pi^2) integral p^2 omega_out |beta|^2 dp. Use omega_out<=2E, |beta|<=C pDelta/(tau^4 E^6), and integral p^4/E^11=8/(315m0^6). These are all exact quadratic transition orders.",
        "local_energy_difference_proof": "Use the SAME free operator and the S6.166 Hadamard in/out states. The symmetric Dirac T00 is the one-particle Hamiltonian bilinear. The covariance difference has trace norm4|beta| per color/flavor, invariant under time evolution. Since ||H(t)||<=2E, |rho_in(t)-rho_out_state(t)|<=4N/pi^2 integral p^2 E|beta| dp<=8NC Delta/(3pi^2 tau^4 m0) for every real t.",
        "renormalization_scope": "Both expectation values use one identical local state-independent subtraction and finite counterterm prescription for this same external mass. Those terms cancel in the state difference. No absolute local renormalized energy, curved/interacting state or backreaction bound is inferred. The named kappa comparison is not a relative error against the zero bounce density.",
        "external_source_energy_exchange": "For the homogeneous prescribed mass, the state-difference identity is d_t Delta rho=M'(t) Delta<bar psi psi>; it is not an isolated conserved fermion energy. The background source remains external in this calculation.",
        "cutoff_scope": "The estimate is for the explicitly specified quadratic operator on its full momentum space, not a proof that the physical SAT8 interacting EFT has a cutoff above those momenta.",
        "checks": checks,
        "gates": {
            "actual_adiabatic_gap_exceeds_fixed_threshold": bool(
                tau * values["mass_floor"] > tube.SCALE
            ),
            "actual_uniform_transition_below_one_e_minus_390": bool(
                values["uniform_exact_transition_amplitude_upper"]
                < s.Rational(1, 10**390)
            ),
            "actual_out_particle_energy_below_one_e_20": bool(
                values["complete_free_out_particle_energy_upper"] < 10**20
            ),
            "actual_local_state_energy_difference_below_one_e_411": bool(
                values["uniform_in_out_local_energy_density_difference_upper"] < 10**411
            ),
            "out_energy_reference_ratio_below_one_e_minus_780": bool(
                values["complete_free_out_particle_energy_upper"] / analytic.KAPPA
                < s.Rational(1, 10**780)
            ),
            "state_difference_reference_ratio_below_one_e_minus_389": bool(
                values["uniform_in_out_local_energy_density_difference_upper"]
                / analytic.KAPPA
                < s.Rational(1, 10**389)
            ),
            "sharper_out_energy_allowance_than_frozen_parent": bool(
                values["complete_free_out_particle_energy_upper"]
                < old["complete_free_out_particle_energy_density_upper"]
            ),
            "positive_real_time_frequency_below_twice_floor": bool(
                m + d < 2 * values["mass_floor"]
            ),
        },
    }
