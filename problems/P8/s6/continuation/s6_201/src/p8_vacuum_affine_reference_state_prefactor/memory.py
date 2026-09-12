"""Complete overlapping-time prefactor memory correction and tail."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_metric_noise import reference
from p8_vacuum_affine_retarded_state_remainder import memory as original

from . import prefactor


@cache
def constants():
    C, B, m = reference.PAIR_REF, prefactor.BETA, reference.MASS
    J10 = original.constants()["J10_pi_lower_bound"]
    return {
        "full_J10_bound": J10,
        "full_memory_H1_coefficient": 72 * C**2 * B**2 * J10,
        "low_external_K7_numerator": 72 * C**2 * B**2 * 100,
        "low_external_K_numerator": 72 * C**2 * B**2 * 100 / m**6,
        "high_external_K_numerator": s.Rational(36, 100) * C**2 * B**2 * J10,
    }


@cache
def data():
    nu, mu, L, J = s.symbols("nu mu L J10", positive=True)
    C, B = reference.PAIR_REF, prefactor.BETA
    c = constants()
    K, r = s.symbols("K r", positive=True)
    checks = {
        "complete_internal_weight_symmetry": s.expand(
            nu * mu * (nu**-12 + mu**-12) - (mu / nu) * nu**-10 - (nu / mu) * mu**-10
        ),
        "nine_pair_two_leg_weight_and_internal_integral": 9
        * 2
        * C**2
        * B**2
        * 2
        * L
        * J
        - 36 * C**2 * B**2 * L * J,
        "unit_time_and_spatial_H1_factor": 2 * 36 * C**2 * B**2 * J
        - 72 * C**2 * B**2 * J,
        "both_leg_half_band_tenth_radial_tail": s.integrate(r**-8, (r, K / 2, s.oo))
        - 128 / (7 * K**7),
    }
    return {
        "definition": "R_alpha,memory=J_memory[alpha_initial W8]-J_memory[unit W8], using the same complete nine-pair finite-regulator formula and the unchanged theta(t-s).",
        "pointwise": "The exact weight increment multiplies the full unit-reference product, whose modulus is<=Cref^2 nu mu ||Dhat||||Gammahat||. There is no derivative of alpha,beta or replacement of the retarded triangle.",
        "all_internal": "For P=k+l and L=1+|P|/(Amax m), mu<=L nu and nu<=L mu. The full integral of nu mu(nu^-12+mu^-12) is<=2L J10. With nine pairs the bound is36 Cref^2 B^2 L J10.",
        "full_memory": "Unit-time Cauchy and L<2(1+|P|^2) give |R_alpha,memory|<1e12 M[D]M[Gamma], M[f]^2=||f||L2^2+||grad_x f||L2^2.",
        "low_external_tail": "For |P|<=K/2 the removed both-leg union forces both internal momenta above K/2. J10_tail(K/2)<100/K^7 gives the recorded K^-7 numerator, bounded by the recorded1/K numerator for K>=m.",
        "high_external_tail": "For |P|>K/2 use the full internal integral and L/(1+|P|^2)<1/(100K). This retains high-low pairs at every external momentum.",
        "constants": c,
        "checks": checks,
        "gates": {
            "strict_memory_display": c["full_memory_H1_coefficient"] < 10**12,
            "low_external_tail_display": c["low_external_K_numerator"] < 3 * 10**16,
            "high_external_tail_display": c["high_external_K_numerator"] < 10**10,
            "original_massive_J10_bound": 0
            < c["full_J10_bound"]
            < 1 / reference.MASS**7,
        },
    }
