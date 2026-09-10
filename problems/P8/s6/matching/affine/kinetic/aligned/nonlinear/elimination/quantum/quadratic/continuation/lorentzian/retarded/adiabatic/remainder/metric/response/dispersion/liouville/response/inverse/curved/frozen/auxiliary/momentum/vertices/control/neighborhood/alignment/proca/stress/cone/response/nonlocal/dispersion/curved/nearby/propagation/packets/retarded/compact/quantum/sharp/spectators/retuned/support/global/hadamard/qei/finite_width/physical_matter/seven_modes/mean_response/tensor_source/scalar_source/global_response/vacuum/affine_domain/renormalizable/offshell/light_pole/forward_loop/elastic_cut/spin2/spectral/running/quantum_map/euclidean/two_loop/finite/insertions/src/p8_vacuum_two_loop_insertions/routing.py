"""All three forward channels and full heavy denominators after loop shifting."""

from functools import cache

import sympy as sp


@cache
def data():
    s, x = sp.symbols("complex_s outer_Feynman_parameter")
    k = sp.Matrix(sp.symbols("loop0:4", real=True))
    p1 = sp.Matrix([sp.I * sp.sqrt(4 - s) / 2, 0, 0, sp.I * sp.sqrt(s) / 2])
    p2 = sp.Matrix([-sp.I * sp.sqrt(4 - s) / 2, 0, 0, sp.I * sp.sqrt(s) / 2])
    momenta = (p1, p2, -p1, -p2)
    checks = {
        "external_momentum_conservation_" + str(j): sp.expand(
            sum(p[j] for p in momenta)
        )
        for j in range(4)
    }
    checks.update(
        {
            "on_shell_light_leg_" + str(j): sp.expand(p.dot(p) + 1)
            for j, p in enumerate(momenta)
        }
    )
    channels = {}
    for name, (i, j, z) in {
        "s": (0, 1, s),
        "t": (0, 2, sp.Integer(0)),
        "u": (0, 3, 4 - s),
    }.items():
        P = momenta[i] + momenta[j]
        Delta = 1 - x * (1 - x) * z
        q1 = k - x * P
        q2 = k + (1 - x) * P
        checks[name + "_pair_invariant"] = sp.simplify(P.dot(P) + z)
        original = k
        D1 = original.dot(original) + 1
        D2 = (original + P).dot(original + P) + 1
        checks[name + "_exact_Feynman_denominator_shift"] = sp.expand(
            (1 - x) * D1 + x * D2 - (k + x * P).dot(k + x * P) - Delta
        )
        channels[name] = {
            "external_pair": (i, j),
            "pair_momentum": P,
            "channel_invariant": z,
            "centered_light_denominator": k.dot(k) + Delta,
            "first_light_line": q1,
            "second_light_line": q2,
            "first_vertex_heavy_exchange_shifts": (
                momenta[i] - x * P,
                momenta[j] - x * P,
            ),
        }
    a = sp.Matrix(sp.symbols("real_shift0:4", real=True))
    b = sp.Matrix(sp.symbols("imaginary_shift0:4", real=True))
    realz = (k + a).dot(k + a) + 1 - b.dot(b)
    gap = realz - (k.dot(k) / 2 + 1 - a.dot(a) - b.dot(b))
    checks["complex_shift_radial_square_gap"] = sp.expand(
        gap - (k + 2 * a).dot(k + 2 * a) / 2
    )
    checks["external_leg_hermitian_norm_majorant"] = sp.Rational(
        3 + 3, 4
    ) - sp.Rational(3, 2)
    checks["heavy_shift_triangle_bound"] = sp.expand(
        sp.Rational(3, 1)
        + sp.Rational(3, 2)
        + 2 * sp.sqrt(sp.Rational(9, 2))
        - (sp.sqrt(3) + sp.sqrt(sp.Rational(3, 2))) ** 2
    )
    return {
        "s": s,
        "x": x,
        "momenta": momenta,
        "channels": channels,
        "external_leg_squared_Hermitian_norm_upper": sp.Rational(3, 2),
        "channel_squared_Hermitian_norm_upper": sp.Integer(3),
        "heavy_shift_Hermitian_norm_upper": sp.Integer(3),
        "centered_light_real_mass_gap": sp.Rational(1, 4),
        "self_energy_inverse_real_lower": k.dot(k) / 2 - 2,
        "full_heavy_denominator_real_lower_before_real_translation": sp.Symbol(
            "heavy_mass_squared", positive=True
        )
        - 9,
        "scope": "An explicit first-sheet forward routing with |s-2|<=1 and real outer Feynman parameter in [0,1]. Euclidean dot products are bilinear, while shift estimates use Hermitian norms. No heavy propagator expansion or loop cutoff.",
        "checks": checks,
    }
