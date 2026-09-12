"""A separately named canonical polynomial model with a bounded classical potential."""

from functools import cache

import sympy as s
from p8_vacuum_affine_scalar_matching_scale import amplitude as original
from p8_vacuum_affine_scalar_transfer_moment import atomic as atom

NAME = "V2S-T1"
LAMBDA, GAMMA = original.LAMBDA, original.GAMMA
GAP = 2 * LAMBDA / GAMMA
MASS2 = GAP + 2
G2 = GAMMA * GAP**4
G = s.Rational(1, 8192)
CONTACT = -G2 * (3 / GAP - 2 / GAP**2)
D, g = s.symbols("positive_heavy_gap positive_trilinear", positive=True)
phi, H = s.symbols("light_phi heavy_H", real=True)
C = -g * g * (3 / D - 2 / D**2)
V = phi * phi / 2 + (D + 2) * H * H / 2 - g * H * phi * phi / 2 - C * phi**4 / 24
VALLEY_QUARTIC = g * g * (D - 1) / (6 * D * D * (D + 2))


@cache
def data():
    vp, vH = s.symbols("velocity_phi velocity_H", real=True)
    kinetic = (vp * vp + vH * vH) / 2
    square = (
        phi * phi / 2
        + (D + 2) * (H - g * phi * phi / (2 * (D + 2))) ** 2 / 2
        + VALLEY_QUARTIC * phi**4
    )
    coercive = (
        phi * phi / 3
        + (3 * D + 5) * (H - 3 * g * phi * phi / (2 * (3 * D + 5))) ** 2 / 6
        + g * g * (9 * D - 10) * phi**4 / (24 * D * D * (3 * D + 5))
    )
    checks = {
        "full_potential_positive_square_identity": s.cancel(V - square),
        "explicit_global_quadratic_coercivity_identity": s.cancel(
            V - (phi * phi + H * H) / 6 - coercive
        ),
        "canonical_light_velocity_entry": s.diff(kinetic, vp, 2) - 1,
        "canonical_heavy_velocity_entry": s.diff(kinetic, vH, 2) - 1,
        "canonical_velocity_cross_zero": s.diff(kinetic, vp, vH),
        "light_vacuum_mass_squared": s.diff(V, phi, 2).subs({phi: 0, H: 0}) - 1,
        "heavy_vacuum_mass_squared": s.diff(V, H, 2).subs({phi: 0, H: 0}) - (D + 2),
        "vacuum_mass_cross_zero": s.diff(V, phi, H).subs({phi: 0, H: 0}),
        "literal_cubic_vertex_factor": s.diff(g * H * phi * phi / 2, H, phi, phi) - g,
        "literal_contact_vertex_factor": s.diff(C * phi**4 / 24, phi, 4) - C,
        "positive_square_remaining_quartic": s.cancel(
            -C / 24 - g * g / (8 * (D + 2)) - VALLEY_QUARTIC
        ),
        "actual_fixed_coupling_squared": G2 - s.Rational(1, 2**26),
        "actual_positive_coupling_root": G * G - G2,
        "unchanged_positive_atom_mass": MASS2 - atom.MASS2,
        "unchanged_positive_atom_coupling": G2 - atom.COUPLING2,
    }
    return {
        "separately_named_model": NAME,
        "Minkowski_classical_Lagrangian": "L=(partial phi)^2/2-phi^2/2+(partial H)^2/2-M_H^2 H^2/2+(g/2)H phi^2+(C/24)phi^4, signature+---. No original affine action or coefficient is changed.",
        "canonical_kinetic_matrix": s.hessian(kinetic, (vp, vH)),
        "full_classical_potential": V,
        "positive_square_potential": square,
        "explicit_coercive_square_remainder": coercive,
        "actual_global_coercivity": "At the actual D>2, V>=(phi^2+H^2)/6 globally, with the displayed remainder a sum of nonnegative terms. This controls the classical flat potential, not a quantum or cosmological stability problem.",
        "generic_parameter_domain": "D>1,g>0; the actual parameters obey this strict domain.",
        "actual_parameters": {
            "target_lambda": LAMBDA,
            "target_gamma": GAMMA,
            "D": GAP,
            "M_H_squared": MASS2,
            "g": G,
            "g_squared": G2,
            "C": CONTACT,
        },
        "classical_vacuum_result": "The potential is coercive and nonnegative with unique global minimum phi=H=0. The canonical kinetic matrix is positive identity and the vacuum masses are1 and sqrt(D+2)>2. This is a classical flat two-field statement, not a quantum vacuum/matching theorem or a cosmological solution.",
        "model_scope": "This separate local polynomial model has only the named light scalar and heavy scalar. It is not the original Proca/DHOST affine parent, does not match its full independent functions, and is not asserted to be an exact quantum UV completion. Classical power counting, a bounded potential and tree positivity do not close original V/G/B.",
        "checks": checks,
        "gates": {
            "actual_gap_greater_than_one": GAP > 1,
            "actual_coercive_square_all_coefficients_positive": GAP > 2,
            "actual_light_and_heavy_vacuum_masses_positive": MASS2 > 4,
            "actual_contact_is_negative_and_tiny": 0
            < -CONTACT
            < s.Rational(1, 10**204),
            "actual_dimensionless_exchange_strength_tiny": 0
            < G2 / MASS2
            < s.Rational(1, 10**204),
            "actual_remaining_quartic_strictly_positive": VALLEY_QUARTIC.subs(
                {D: GAP, g: G}
            )
            > 0,
        },
    }
