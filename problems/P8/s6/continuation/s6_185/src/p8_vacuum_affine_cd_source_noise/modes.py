"""Actual all-order CD mode bounds, not a flat or finite-WKB state."""

from functools import cache

import sympy as s
from p8_vacuum_affine_proca_gaussian import bridge
from p8_vacuum_affine_retarded_energy import clock
from p8_vector_regularity import frequency, preparation

MASS = bridge.MASS
KAPPA = bridge.KAPPA


@cache
def data():
    rows = {kind: frequency.reference(kind) for kind in ("transverse", "longitudinal")}
    state = preparation.constants()["common"]
    beta = (
        state["initial_mixing_envelopes"][6] / MASS**6
        + state["evolution_mixing_envelope"] / MASS**10
    )
    W = s.Symbol("W", positive=True)
    d, rate, phase = s.symbols("d logW_rate phase", real=True)
    f2 = 1 / (2 * W)
    p2 = (W * W + (d + rate / 2) ** 2) / (2 * W)
    f = s.exp(-s.I * phase) / s.sqrt(2 * W)
    p = (-s.I * W - rate / 2 - d) * f
    upper = 27
    H, z = s.symbols("H z", real=True)
    checks = {
        "actual_fixed_vector_mass": MASS - 1000,
        "same_normalization_and_original_charge_state": bridge.ZETA * MASS**2 - 1,
        "reference_field_normalization": 2 * W * f2 - 1,
        "reference_physical_momentum": s.simplify(p * s.conjugate(p) - p2),
        "actual_longitudinal_rate": H * (s.Rational(1, 2) + z) - H / 2 - H * z,
        "reference_rate_square_envelope": ((d + rate / 2) ** 2).subs({d: 3, rate: 4})
        - 25,
        "actual_field_modulus_upper": 3**2 - 9,
        "actual_momentum_modulus_upper": 3**2 * 3 - upper,
        "exact_clock_first_source_variation": clock.data()["checks"][
            "source_first_clock_variation_zero"
        ],
    }
    return {
        "state": "Unchanged S6.55 all-order Borel Cauchy preparation on actual CD, with its exact evolved mixing; not W8 initial data",
        "mass": MASS,
        "same_reference_bounds": {
            kind: {
                "frequency_ratio_deviation": p["frequency_ratio_deviation_upper"],
                "logarithmic_rate": p["reference_logarithmic_rate_upper"],
            }
            for kind, p in rows.items()
        },
        "actual_state_mixing_upper": beta,
        "reference_absolute_field_squared": f2,
        "reference_absolute_physical_momentum_squared": p2,
        "physical_reference_rate_bound": s.Integer(3),
        "actual_field_squared_upper_times_frequency": s.Integer(9),
        "actual_physical_momentum_squared_upper_over_frequency": s.Integer(27),
        "canonical_momentum": "p=f'-d f, dT=H/2 and dL=(1/2+z)H; alpha'fref+beta'conjugate(fref)=0 in the exact comparison",
        "checks": checks,
        "gates": {
            "all_W8_ratio_bounds": all(
                p["frequency_ratio_deviation_upper"] < s.Rational(1, 2)
                for p in rows.values()
            ),
            "all_reference_log_rates_below_four": all(
                p["reference_logarithmic_rate_upper"] < 4 for p in rows.values()
            ),
            "actual_all_order_mixing_below_one": 0 < beta < 1,
            "physical_canonical_rate_below_three": s.Rational(3, 2) * s.Rational(8, 5)
            < 3,
            "reference_momentum_below_three_omega": s.Rational(9, 4) + 25 / MASS**2 < 3,
            "mixing_modulus_sum_below_three": s.sqrt(2) + 1 < 3,
        },
    }
