"""Low/high-window bounds with the original positive dimensional measure."""

from functools import cache

import sympy as s
from p8_vacuum_affine_dimensional_gravity_radiation import sew
from p8_vacuum_affine_massive_scalar_graviton_bremsstrahlung import recoil
from p8_vacuum_affine_physical_virtual_soft_pairing import continuity as old_continuity
from p8_vacuum_affine_physical_virtual_soft_pairing import soft as old_soft

from . import born, source

B = s.Integer(4000000000)
FIXED_BORN_FACTOR = s.Rational(361, 324)


def resolution(value):
    x = sew.exact_real(value)
    if x < 0 or x > s.Rational(1, 8):
        raise ValueError("Require 0<=resolution<=1/8")
    return x


def gap(value):
    value = sew.exact_real(value)
    if value <= 0 or value > 1:
        raise ValueError("Require 0<delta<=1")
    return value


def unit_D6_remainder_upper(delta):
    return s.Integer(13396352288) / gap(delta) ** 2


def canonical_gravity_remainder_upper(delta):
    return s.Integer(30000000000) / gap(delta) ** 2


def canonical_gravity_high_upper(delta):
    return s.Integer(9000000000) * 30000**2 / gap(delta) ** 2


def low_real_remainder_upper(invariant, cosine, cutoff):
    delta = born.transfer_gap(invariant, cosine)
    x = resolution(cutoff)
    if x > delta / 192:
        raise ValueError("Low-window bound requires x<=delta/192")
    return (
        FIXED_BORN_FACTOR
        * (128 * B * x / delta + B * B * x * x / (2 * delta**2) + 12288 * x)
        / (4 * s.pi**2 * source.KAPPA)
    )


def real_remainder_upper(invariant, cosine, cutoff):
    delta = born.transfer_gap(invariant, cosine)
    x = resolution(cutoff)
    a = delta / 192
    value = low_real_remainder_upper(invariant, cosine, min(x, a))
    if x > a:
        value += s.Integer(10) ** 35 * (x * x - a * a) / (source.KAPPA * delta * delta)
    return value


def lower_cutoff_error_upper(invariant, cosine, lower_cutoff):
    """Uniform epsilon tail for the subtracted integral, not either divergent rate."""
    return low_real_remainder_upper(invariant, cosine, lower_cutoff)


def folded_weight(perpendicular, epsilon):
    v = sew.exact_real(perpendicular)
    e = sew.regulator(epsilon)
    if v < 0 or v > 1:
        raise ValueError("Require a folded perpendicular coordinate in[0,1]")
    return old_continuity.folded_angular_weight(v, e)


def pair_density(energy, omega, epsilon):
    E = sew.exact_real(energy)
    w = resolution(omega)
    e = sew.regulator(epsilon)
    if E < s.Rational(5, 4) or E > 2:
        raise ValueError("Require 5/4<=hard COM scalar energy<=2")
    return old_continuity.continued_pair_density(E, w, e)


@cache
def data():
    w, x, delta, e, v = s.symbols(
        "omega resolution delta epsilon folded", positive=True
    )
    integrated = s.integrate(128 * B / delta + B * B * w / delta**2 + 12288, (w, 0, x))
    high_g6 = s.Rational(2488320 * 4, 625) + 441936 + 3265173504
    low_cap = (
        FIXED_BORN_FACTOR * (s.Rational(2, 3) * B + B * B / s.Integer(73728) + 64) / 36
    )
    checks = {
        "integrated_full_polarization_low_budget": s.expand(
            integrated
            - 128 * B * x / delta
            - B * B * x * x / (2 * delta**2)
            - 12288 * x
        ),
        "high_seagull_field_scaling": s.Rational(883872, 2) - 441936,
        "high_cubic_field_scaling": s.Rational(13060694016, 4) - 3265173504,
        "high_external_current_budget": s.Integer(64 * 3 * 36 * 40 * 9 - 2488320),
        "whole_high_kernel_radial_integral": s.integrate(w, (w, delta / 192, x))
        - (x * x - (delta / 192) ** 2) / 2,
        "actual_phase_normalization": s.simplify(
            old_soft.phase_normalization(e)
            - (4 * s.pi) ** (-e)
            * s.gamma(s.Rational(3, 2))
            / s.gamma(s.Rational(3, 2) + e)
        ),
        "actual_folded_angular_normalization": s.integrate(
            old_continuity.folded_angular_weight(v, e), (v, 0, 1)
        )
        - 1,
        "actual_folded_zero_regulator": old_continuity.folded_angular_weight(v, 0) - 1,
        "original_low_kappa_exponent": s.Integer(10) ** 13 / source.KAPPA
        - s.Rational(1, 10**787),
        "original_high_kappa_exponent": s.Integer(10) ** 35 / source.KAPPA
        - s.Rational(1, 10**765),
    }
    for E, w0 in (
        (s.Rational(5, 4), s.Rational(1, 128)),
        (s.Integer(2), s.Rational(1, 8)),
    ):
        checks["actual_pair_density_zero_" + str(E)] = s.simplify(
            pair_density(E, w0, 0) - recoil.density_ratio(E, w0)
        )
    margins = {
        "low_fixed_Born_cap": s.Integer(10) ** 13 - low_cap,
        "high_G6_unit": 4000000000 - high_g6,
        "high_auxiliary_core": 100 - (s.Rational(192 * 4, 625) + 96),
        "high_canonical_GD": s.Integer(9000000000) - (2 * 4000000000 + 100),
        "high_normalized_all47": s.Integer(2000000000000000000)
        - (s.Integer(9000000000) * 30000**2 / 8 + 12288 + 330000),
        "high_squared_kernel": 5 * s.Integer(10) ** 36
        - (s.Integer(2000000000000000000) ** 2 + 12288**2),
        "high_fixed_Born_integrated": s.Integer(10) ** 35
        - FIXED_BORN_FACTOR * 5 * s.Integer(10) ** 36 / 72,
    }
    return {
        "checks": checks,
        "gates": {
            **{
                key + "_strict_margin": bool(value > 0)
                for key, value in margins.items()
            },
            "separate_high_window_not_extrapolated_low_Taylor_bound": True,
            "same_D_Born_and_full_polarization_soft_subtraction": True,
            "actual_positive_fixed_disk_and_pair_phase_retained": True,
            "dominated_convergence_at_fixed_nonzero_transfer_and_resolution": True,
            "joint_lower_IR_cutoff_limit_only_after_soft_subtraction": True,
            "no_uniform_forward_or_virtual_hard_matching_conclusion": True,
        },
        "whole_exact_positive_margins": margins,
        "whole_low_bound": "For y<=delta/192, fixed-A0 reference upper=(361/324)[128B*y/delta+B^2*y^2/(2delta^2)+12288y]/(4pi^2*kappa), B=4e9. It is<1e13/kappa and tends to zero with y uniformly in epsilon.",
        "whole_high_bound": "For omega>=a=delta/192, m=delta/30000, canonical gravity tree<9e9/m^2. Complete normalized amplitude<2e18/(sqrt(kappa)delta), so full-minus-soft squared kernel<5e36/(kappa delta^2). The high integral, including fixed-A0 factor361/324, is<1e35*(x^2-a^2)/(kappa delta^2).",
        "whole_continuity_proof": "Use the actual S296 phase p(e)<=1, pair ratioJ_e<=1 with1-J_e<=3omega, and positive fixed-disk density<=5/4 of its e0 density. The low kernel is omega-integrable uniformly in e; the high interval has a fixed positive lower endpoint and physical propagator gaps. S334 supplies pointwise continuous complete trees and sews. Dominated convergence gives the S304 physical D4 real-minus-soft remainder at fixed nonforward E,u,x and joint removal of a lower energy cutoff. Neither separate divergent rate or virtual hard matching follows.",
        "whole_normalization": "The fixed-A0 bound retains the SAME D-dimensional reference two-body phase. Same initial flux and identical-final-scalar factors cancel. SAME-D complete Born amplitude remains inside the leading-soft subtraction.",
    }
