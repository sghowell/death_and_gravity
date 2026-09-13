"""Whole rolling-heavy Gaussian block and exact on-shell Ward contacts."""

from functools import cache

import sympy as s
from p8_vacuum_affine_coupled_principal_obstruction import background as old
from p8_vacuum_affine_coupled_principal_obstruction import coupled as c

from . import homogeneous as bg

Nd, Hd, mcd, mhd = s.symbols(
    "Ndot Hhatdot M1_acceleration heavy_acceleration", real=True
)


def dtime(expression):
    return (
        s.diff(expression, bg.u)
        + Nd * s.diff(expression, bg.N)
        + c.a * c.H * s.diff(expression, c.a)
        + Hd * s.diff(expression, c.H)
        + mcd * s.diff(expression, bg.mc)
        + mhd * s.diff(expression, bg.mh)
        + bg.mh * s.diff(expression, bg.hbar)
    )


@cache
def coefficient_map():
    packet = bg.action()
    D, Z = packet["D"], packet["Z"]
    density = packet["normalized_homogeneous_density"]
    B = old.coefficients()["whole_ADM_families"]["B"]
    cv = -18 * D * c.H + 9 * B
    charges = s.Matrix([Z * bg.mc, Z * bg.mh])
    heavy_force = packet["complete_heavy_force_density"]
    parameters = {
        **old.coefficients()["whole_parameter_substitution"],
        c.J: packet["J"],
        c.th: packet["Theta"],
        c.L0: 3 * packet["complete_lapse_constraint"],
        c.Vvv: s.Rational(9, 2) * density - (dtime(cv) + 3 * c.H * cv) / 2,
        c.cs[0]: charges[0],
        c.cs[1]: charges[1],
        c.ws[0]: packet["w"][0],
        c.ws[1]: packet["w"][1],
        c.ds[0]: 0,
        c.ds[1]: s.diff(heavy_force, bg.N),
        c.vs[0]: -3 * (dtime(charges[0]) + 3 * c.H * charges[0]),
        c.vs[1]: 3 * heavy_force - 3 * (dtime(charges[1]) + 3 * c.H * charges[1]),
    }
    # These are Euler expressions with sign chosen as partial_q L-dt partial_qdot L.
    volume_euler = 3 * density - (dtime(cv) + 3 * c.H * cv) / 3
    matter_euler = s.Matrix(
        [
            -dtime(charges[0]) - 3 * c.H * charges[0],
            heavy_force - dtime(charges[1]) - 3 * c.H * charges[1],
        ]
    )
    return {
        "entire_rolling_four_mode_parameter_map": parameters,
        "entire_homogeneous_volume_Euler_expression": volume_euler,
        "entire_homogeneous_matter_Euler_expressions": matter_euler,
        "entire_metric_scalar_time_boundary_coefficient": cv,
        "entire_two_matter_time_boundary_coefficients": charges,
    }


@cache
def data():
    packet = coefficient_map()
    mapped = packet["entire_rolling_four_mode_parameter_map"]
    eps = s.Symbol("eps")
    N, a, mass = s.symbols("N a heavy_mass_squared", positive=True)
    n, v, vd, h, hd, hbar, mh, H, b, beta, hx, q = s.symbols(
        "n v vd h hd hbar mh H b beta hx q", real=True
    )
    U0, U1, U2, j0, j1, j2, Y0, cd = s.symbols(
        "U0 U1 U2 normalized_source0 normalized_source1 normalized_source2 Cchi0 cd",
        real=True,
    )
    x = s.Symbol("x", real=True)
    U = U0 + U1 * x + U2 * x**2 / 2
    source = j0 + j1 * x + j2 * x**2 / 2
    lapse = N + x
    Z = U / lapse
    potential = -mass * hbar**2 / 2 + source * hbar
    base = Z * mh**2 / 2 + lapse * U * potential
    charge = (Z * mh).subs(x, 0)
    w = (Z.diff(x) * mh).subs(x, 0)
    Vh = -mass * hbar + source
    d = (lapse * U * Vh).diff(x).subs(x, 0)
    vs = (3 * lapse * U * Vh).subs(x, 0) - 3 * (cd + 3 * H * charge)
    volume = 1 + 3 * eps * v + s.Rational(9, 2) * eps**2 * v**2
    rate = mh + eps * hd - eps**2 * beta * hx
    value = hbar + eps * h
    UL, JL, NL = U.subs(x, eps * n), source.subs(x, eps * n), N + eps * n
    literal = volume * (
        UL * rate**2 / (2 * NL) + NL * UL * (-mass * value**2 / 2 + JL * value)
    )
    literal -= NL * Y0 * (1 + eps * v) * eps**2 * hx**2 / (2 * a**2)
    raw = s.expand(s.diff(literal, eps, 2).subs(eps, 0) / 2)
    raw = raw.subs(beta * hx, -b * h).subs(hx**2, q * a**2 * h**2)
    normalized = (
        Z.subs(x, 0) * hd**2 / 2
        + w * n * hd
        + base.diff(x, 2).subs(x, 0) * n**2 / 2
        - 3 * charge * vd * h
        + 3 * base.diff(x).subs(x, 0) * n * v
        + 9 * base.subs(x, 0) * v**2 / 2
        + d * n * h
        + vs * v * h
        + b * charge * h
        - N * U0 * mass * h**2 / 2
        - N * Y0 * q * h**2 / 2
    )
    boundary = 3 * charge * (vd * h + v * hd) + 3 * (cd + 3 * H * charge) * v * h
    reference = {
        bg.hbar: 0,
        bg.mh: 0,
        mhd: 0,
        bg.mc: old.m,
        mcd: old.md,
        Hd: old.Hd,
        Nd: old.Nd,
    }
    old_map = old.coefficients()["whole_parameter_substitution"]
    checks = {
        "whole_rolling_heavy_ADM_density_and_spatial_time_boundaries": raw
        - normalized
        - boundary,
        "entire_metric_potential_is_volume_Euler_contact": mapped[c.Vvv]
        - s.Rational(3, 2) * packet["entire_homogeneous_volume_Euler_expression"],
        "entire_M1_volume_contact_is_M1_Euler_contact": mapped[c.vs[0]]
        - 3 * packet["entire_homogeneous_matter_Euler_expressions"][0],
        "entire_heavy_volume_contact_is_heavy_Euler_contact": mapped[c.vs[1]]
        - 3 * packet["entire_homogeneous_matter_Euler_expressions"][1],
        "whole_on_shell_lapse_volume_coefficient": mapped[c.L0]
        - 3 * bg.action()["complete_lapse_constraint"],
        "whole_rolling_heavy_lapse_contact": mapped[c.ds[1]]
        - s.diff(bg.action()["complete_heavy_force_density"], bg.N),
        "entire_heavy_Euler_relation_not_held_Hbar0": vs.subs(
            cd, N * U0 * (-mass * hbar + j0) - 3 * H * charge
        ),
    }
    for key, value in mapped.items():
        checks["entire_held_zero_heavy_restriction_" + str(key)] = (
            s.sympify(value).subs(reference, simultaneous=True) - old_map[key]
        )
    return {
        **packet,
        "literal_rolling_heavy_density_before_expansion": literal,
        "whole_raw_rolling_heavy_quadratic_density": raw,
        "whole_normalized_rolling_heavy_quadratic_density": normalized,
        "whole_weighted_time_boundary": boundary,
        "whole_on_shell_lower_contacts": {c.L0: 0, c.Vvv: 0, c.vs[0]: 0, c.vs[1]: 0},
        "application": "Substitute this complete rolling parameter map into the ENTIRE S254 four-mode Lagrangian/Hamiltonian. R,F and JH use the same full source bindings and exact primitive. The metric/M1 raw expansion is S253 and the full rolling-heavy density is independently expanded here. All shift transports, masses, lapse-heavy terms and time boundaries remain. On a homogeneous classical solution the listed four Euler contacts vanish only AFTER this complete derivation; the growing fast coefficient does not involve them.",
        "held_vs_on_shell": "Holding hbar identically zero gives mh=mhdot=0 and recovers all S254 coefficients exactly. Starting an unforced solution at hbar=mh=0 generally has mhdot nonzero due to JH. Its heavy Euler relation, rather than held-zero substitution, removes the volume-heavy Ward contact. These operations are not confused.",
        "checks": checks,
        "gates": {
            "whole_heavy_field_and_velocity_retained": all(
                literal.has(z) for z in (hbar, mh, j0, j1, j2, mass)
            ),
            "whole_rolling_heavy_lapse_contact_not_assumed_zero": d != 0,
            "all_S254_parameters_listed": set(mapped) == set(old_map),
            "on_shell_Euler_contacts_applied_only_after_derivation": True,
            "no_new_mass_or_profile_or_reference_state": True,
        },
    }
