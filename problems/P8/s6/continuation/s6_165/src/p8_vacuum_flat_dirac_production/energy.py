"""Finite free out-particle energy density, with explicit degeneracies."""

from functools import cache

import sympy as s

from . import transition


def enclosure(mass, amplitude, timescale, multiplicity=6):
    m, d, tau = map(transition.exact, (mass, amplitude, timescale))
    if min(m, tau) <= 0 or d < 0:
        raise ValueError("Require positive mass/time and nonnegative amplitude")
    if type(multiplicity) is not int or multiplicity < 1:
        raise TypeError("Require a positive native active color/flavor multiplicity")
    floor = s.Rational(99, 100) * m
    if d / floor >= s.Rational(1, 100):
        raise ValueError("Require amplitude/lower_mass<1/100")
    first = s.Rational(8 * multiplicity, 9) * 5 * d * d / (tau**4 * floor**2)
    remainder = s.Rational(8 * multiplicity, 9) * d**6 / (225 * floor**2)
    beta_uniform = 5 * d / (2 * tau * tau * floor**3) + d**3 / (40 * floor**3)
    return {
        "strict_mass_lower": floor,
        "max_integrated_mixing_upper": d / (2 * floor),
        "first_transition_part_of_energy_upper": first,
        "higher_transition_part_of_energy_upper": remainder,
        "complete_free_out_particle_energy_density_upper": first + remainder,
        "uniform_exact_transition_amplitude_upper": beta_uniform,
    }


@cache
def data():
    p, m, d, tau, N = s.symbols("p m0 Delta tau N", positive=True)
    E = s.sqrt(p * p + m * m)
    A = 5 * p * d / (tau * tau * E**4)
    B = p**3 * d**3 / (5 * E**6)
    x = s.Symbol("x", positive=True)
    Ia = s.integrate(x**4 / (1 + x * x) ** s.Rational(7, 2), (x, 0, s.oo))
    Ib = s.integrate(x**8 / (1 + x * x) ** s.Rational(11, 2), (x, 0, s.oo))
    rho = 8 * N / s.pi**2 * (5 * d * d / (tau**4 * m * m) + d**6 / (225 * m * m))
    return {
        "transition_integral_upper": A,
        "all_higher_transition_upper": B,
        "complete_free_out_energy_upper": rho,
        "mode_measure": "rho_out=4N/(2pi^2) integral_0^infinity p^2 omega_out |beta_p|^2 dp. N=6 active color/flavor components, two spin/helicity states and particles plus antiparticles. The twelve inert flavors have constant flat masses and beta=0.",
        "energy_argument": "omega_out<=2E, E=sqrt(p^2+m0^2), and (A+B)^2<=2A^2+2B^2. The radial integrals are integral p^4/E^7 dp=1/(5m0^2) and integral p^8/E^11 dp=1/(9m0^2). Both upper terms decay as p^-3 in the transition amplitude, making the out-particle energy density finite.",
        "scope": "Energy of the free out-particle Hamiltonian normal ordered against its free out-vacuum, per unit spatial volume. It is not a finite total energy, a bound on transient renormalized local stress/coherence, an interacting gauge observable or curved-space backreaction.",
        "checks": {
            "first_radial_integral": Ia - s.Rational(1, 5),
            "higher_radial_integral": Ib - s.Rational(1, 9),
            "complete_radial_energy_composition": s.simplify(
                rho
                - 8
                * N
                / s.pi**2
                * (25 * d * d * Ia / (tau**4 * m * m) + d**6 * Ib / (25 * m * m))
            ),
            "spherical_measure_with_Dirac_degeneracy": 4
            * N
            * 4
            * s.pi
            / (2 * s.pi) ** 3
            - 2 * N / s.pi**2,
            "two_norm_inequality_positive_difference": s.expand(
                2 * A * A + 2 * B * B - (A + B) ** 2 - (A - B) ** 2
            ),
            "out_frequency_and_two_norm_factors": 2 * N * 2 * 2 - 8 * N,
            "first_amplitude_large_momentum_degree": 1 - 4 + 3,
            "higher_amplitude_large_momentum_degree": 3 - 6 + 3,
        },
    }
