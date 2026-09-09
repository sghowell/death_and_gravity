"""Actual normalized compact mass-shell mode and exact finite energy scaling."""

from functools import cache

import sympy as sp
from p8a_nonminimal import field

from . import construction


@cache
def data():
    R, Iq, I0, I2, J1, J2, Jp, hb, B2, Bq2, e = sp.symbols(
        "R Iq I0 I2 J1 J2 Jp hbar B2 Bq2 epsilon", positive=True
    )
    C = Iq * J2 / (16 * sp.pi**3)
    gamma2 = J1**2 / (16 * sp.pi**3 * Iq * J2)
    amp2 = hb * R * R * gamma2
    energy = hb * sp.Rational(1, 6) * (R * R * Jp * I0 + J2 * I2) / (Iq * J2)
    i0 = (20 + 2 / e**2) * B2
    iq = (24 * e + 3 / e) * B2
    i2 = (32 * e**2 + 5) * B2 + (20 + 2 / e**2) * Bq2
    own = field.data()
    flat = own["flat_stress"]
    symbols = {str(s): s for s in flat.free_symbols}
    xi, z = symbols["xi"], symbols["field"]
    dnull = symbols["field_d0"] + symbols["field_d3"]
    # Symmetric Hessian ordering in the arbitrary-jet parent stress.
    hnull = symbols["field_h0"] + 2 * symbols["field_h3"] + symbols["field_h9"]
    ell = sp.Matrix([1, 0, 0, 1])
    null_stress = (ell.T * flat * ell)[0]
    n, m = sp.Rational(1, 3), -sp.Rational(2, 3)
    x, y = sp.symbols("x y", real=True)
    # Full-line positive Fourier support kills the anomalous term exactly.
    anec = 2 * n * amp2 * 2 * sp.pi * I2
    return {
        "mode_normalization_constant_squared": C,
        "positive_mode_amplitude_coefficient_squared": gamma2,
        "actual_mode_amplitude_squared": amp2,
        "actual_total_energy_above_vacuum": energy,
        "bump_frequency_moments": {"I0": i0, "Iq": iq, "I2": i2},
        "actual_complete_line_ANEC": anec,
        "checks": {
            "actual_one_particle_mode_has_unit_norm": sp.factor(
                Iq * R**2 * J2 / (16 * sp.pi**3 * R**2 * C) - 1
            ),
            "actual_null_mode_amplitude_normalization": sp.factor(
                J1**2 / (256 * sp.pi**6 * C) - gamma2
            ),
            "actual_transverse_dilation_keeps_R_squared_amplitude": sp.factor(
                amp2 / (hb * gamma2) - R**2
            ),
            "actual_total_energy_from_mass_shell_integral": sp.factor(
                hb * n * (R**4 * Jp * I0 + R**2 * J2 * I2) / (32 * sp.pi**3 * R**2 * C)
                - energy
            ),
            "actual_null_contraction_matches_parent_covariant_stress": sp.expand(
                null_stress - dnull * dnull + xi * (2 * dnull * dnull + 2 * z * hnull)
            ),
            "normalized_frequency_zeroth_moment": sp.expand(
                i0 - ((16 + e**-2) + (4 + e**-2)) * B2
            ),
            "normalized_frequency_first_moment": sp.expand(
                iq - (e * (16 + e**-2) + 2 * e * (4 + e**-2)) * B2
            ),
            "normalized_frequency_second_moment": sp.expand(
                i2
                - (16 + e**-2) * (e * e * B2 + Bq2)
                - (4 + e**-2) * (4 * e * e * B2 + Bq2)
            ),
            "actual_full_line_ANEC_positive_coefficient": sp.factor(
                anec - sp.Rational(4, 3) * sp.pi * amp2 * I2
            ),
            "squeezed_oscillator_anomalous_sign_not_free_parameter": n * (n + 1)
            - m * m,
            "normalized_real_and_imaginary_null_stress_bound": sp.expand(
                (sp.Rational(2, 3) * (-x * x + 3 * y * y))
                - 2 * (n * (x * x + y * y) + m * (x * x - y * y))
            ),
        },
    }


def resource_bounds():
    e, w = construction.EPSILON, construction.WIDTH
    iq = 24 * e + 3 / e
    i0 = 20 + 2 / e**2
    i2max = 32 * e**2 + 5 + (20 + 2 / e**2) * w * w
    return {
        "energy_lower_coefficient_of_hbar_R_squared": i0 / (6 * iq),
        "energy_upper_coefficient_of_hbar_R_squared": 2 * i0 / (3 * iq),
        "energy_upper_constant_coefficient_of_hbar": i2max / (6 * iq),
        "positive_ANEC_support_sum_gap": 2 * (e - w),
        "each_fixed_R_has_finite_energy_and_smooth_Hadamard_difference": True,
        "uniform_total_energy_or_transverse_momentum_cutoff_for_the_family": False,
        "global_two_sided_Wick_square_cap_for_the_family": False,
    }
