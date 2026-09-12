"""Uniform joint complex time and external-momentum frequency domain."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_metric_noise import reference
from p8_vector_state.comparison import AMAX

DELTA = s.Rational(1, 10**6)
OUTER = reference.OUTER
INNER = reference.RADIUS
MASS = reference.MASS


@cache
def constants():
    ea = 6 * OUTER / (1 - 6 * OUTER)
    ep = 2 * AMAX * DELTA + DELTA**2
    ew = ea + (1 + ea) * ep
    return {
        "outer_time_radius": OUTER,
        "inner_readout_time_radius": INNER,
        "spatial_radius_over_nu": DELTA,
        "relative_time_scale_defect": ea,
        "additional_momentum_frequency_defect": ep,
        "complete_relative_squared_frequency_defect": ew,
        "complex_z_modulus_upper": (1 + ew) / (1 - ew),
        "physical_complex_spatial_momentum_over_nu": s.Rational(101, 100)
        * (AMAX + DELTA),
        "Wprime_over_nu": 3 / INNER,
        "canonical_momentum_amplitude_over_sqrt_nu": 3 + (3 / INNER + 5) / MASS,
    }


@cache
def data():
    d = constants()
    rows = reference.data()["actual_W8_complex_coefficient_majorants"]
    delta, k, nu = s.symbols("delta k nu", positive=True)
    checks = {
        "full_squared_momentum_defect": s.expand(
            (k + delta * nu) ** 2 - k * k - 2 * k * delta * nu - delta**2 * nu**2
        ),
        "both_time_and_spatial_relative_defects": s.expand(
            (1 + d["relative_time_scale_defect"])
            * (1 + d["additional_momentum_frequency_defect"])
            - 1
            - d["complete_relative_squared_frequency_defect"]
        ),
        "actual_outer_and_inner_time_radii": OUTER - 2 * INNER,
        "same_fixed_mass": MASS - 1000,
    }
    return {
        "domain": "Fix real |k|>=m, nu=sqrt(m^2+|k|^2/Amax^2), l=-k+p, ||p||complex_Euclidean<=delta nu with delta1e-6. Combine with |u-t|<=1e-4 about every real t in[-1/2,1/2]. Nu is tied to the original real k, not a complex norm.",
        "frequency": "The time scale defect is ea=6r/(1-6r). Since |k|<=Amax nu, the extra squared-momentum defect is<=ep=2Amax delta+delta^2. Relative to the real-time omega(k)^2, the full squared-frequency defect is<=ea+(1+ea)ep<1e-3 for both legs.",
        "branches": "Both squared frequencies lie in the right half-plane. Their positive-real square-root continuations have Re omega>=.99nu and |omega|<2nu. Also |z|=|1-m^2/omega^2|<1.01, and |l/a|Euclidean<2nu.",
        "W8": "The complete four W8 coefficient polynomials are the same source-pinned rational functions as S186. The unchanged |u|<.51, |1+u^2|>.7, |z|<1.01 and |omega|>=.99nu bounds keep both full relative W defects below1e-3. Hence Re W>.5nu and |W|<3nu jointly.",
        "time_derivative": "On the inner time disc1e-4/2, Cauchy from the outer disc gives |Wprime|<=60000nu. With |dT|,|dL|<5, |f_amp|<=nu^-1/2 and |p_amp|<64sqrt(nu), for either phase sign. No initial alpha or beta appears.",
        "constants": d,
        "checks": checks,
        "gates": {
            "complete_joint_frequency_defect_below_one_thousandth": d[
                "complete_relative_squared_frequency_defect"
            ]
            < s.Rational(1, 1000),
            "joint_z_modulus_below_one_point_zero_one": d["complex_z_modulus_upper"]
            < s.Rational(101, 100),
            "joint_omega_lower": s.Rational(99, 100) ** 2
            < 1 - d["complete_relative_squared_frequency_defect"],
            "joint_omega_upper": AMAX**2
            * (1 + d["complete_relative_squared_frequency_defect"])
            < 4,
            "inverse_scale_modulus_below_one_point_zero_one": 1
            + d["relative_time_scale_defect"]
            < s.Rational(101, 100) ** 2,
            "full_complex_spatial_momentum_below_two_nu": d[
                "physical_complex_spatial_momentum_over_nu"
            ]
            < 2,
            "unchanged_W8_complex_polynomial_bounds": all(
                row["ratio_defect"] < s.Rational(1, 1000) for row in rows.values()
            ),
            "canonical_momentum_below_sixty_four": d[
                "canonical_momentum_amplitude_over_sqrt_nu"
            ]
            < 64,
            "high_band_nu_below_twice_real_momentum": 1 + 1 / AMAX**2 < 4,
        },
    }
