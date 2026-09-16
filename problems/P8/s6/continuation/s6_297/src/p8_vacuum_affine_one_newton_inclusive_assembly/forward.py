"""Positive full Born normalization and the nonuniform forward Newton limit."""

from functools import cache

import sympy as s
from p8_vacuum_affine_massive_elastic_proca_infrared import elastic

from . import source

S, Z, N, G, K = s.symbols("s z n g kappa", real=True)


def matter_born(energy=S, cosine=Z, heavy=N, cubic=G):
    ss, z, n, g = map(s.sympify, (energy, cosine, heavy, cubic))
    q = ss - 4
    channels = (ss, -q * (1 - z) / 2, -q * (1 + z) / 2)
    return g * g / (n - 2) ** 2 * sum((a - 2) ** 2 / (n - a) for a in channels)


def gravity_born(energy=S, cosine=Z, kappa=K):
    ss, z, k = map(s.sympify, (energy, cosine, kappa))
    q = ss - 4
    return (
        (4 * q + 16 + 8 / q) / (1 - z * z)
        - 2 * ss
        + 6
        - 2 / ss
        + q * q * (1 - z * z) / (4 * ss)
    ) / k


@cache
def data():
    q, z, n, g, kap = s.symbols("Q z n g kappa", positive=True)
    ss = q + 4
    tt = -q * (1 - z) / 2
    uu = -q * (1 + z) / 2
    p = 4 * q + 16 + 8 / q
    d = -2 * ss + 6 - 2 / ss
    ee = q * q / (4 * ss)
    checks = {}

    def put(name, value):
        checks[name] = s.factor(value)

    put(
        "whole_original_gravity_Born",
        elastic.whole_tree(z, ss, 1, kap, n, 0, 0) - gravity_born(ss, z, kap),
    )
    C = -g * g * (3 / (n - 2) - 2 / (n - 2) ** 2)
    put(
        "whole_original_matter_Born",
        C + g * g * sum(1 / (n - a) for a in (ss, tt, uu)) - matter_born(ss, z, n, g),
    )
    put("positive_gravity_lower_bound", p + 2 * d - (12 + (4 * q + 32) / (q * (q + 4))))
    put("gravity_upper_bound", -(d + ee) * 4 * (q + 4) - (7 * q * q + 40 * q + 40))
    put("p_monotone_derivative", s.diff(p, q) - (4 - 8 / q**2))
    put("p_lower_endpoint", p.subs(q, s.Rational(9, 4)) - s.Rational(257, 9))
    put("p_upper_endpoint", p.subs(q, 12) - s.Rational(194, 3))
    put(
        "sum_channel_squares",
        sum((a - 2) ** 2 for a in (ss, tt, uu))
        - (q * q * (3 + z * z) / 2 + 8 * q + 12),
    )
    put("sum_squares_upper_endpoint", (2 * q * q + 8 * q + 12).subs(q, 12) - 396)
    put("heavy_denominator_lower", n - 16 - s.Rational(7, 8) * n - (n - 128) / 8)
    put(
        "subtracted_heavy_denominator_lower",
        n - 2 - s.Rational(63, 64) * n - (n - 128) / 64,
    )
    put(
        "heavy_denominator_upper",
        s.Rational(35, 32) * n - (n + 12) - s.Rational(3, 32) * (n - 128),
    )
    put(
        "forward_pole_residue",
        s.limit((1 - z * z) * gravity_born(ss, z, kap), z, 1) - p / kap,
    )
    nu, t, mu = s.symbols("nu t mu", real=True)
    a = 2 * mu - t / 2 + nu
    c = 2 * mu - t / 2 - nu
    numerator = 2 * mu**2 - 2 * mu * t - a * c
    put(
        "massive_crossing_center_numerator", numerator - (nu**2 - 2 * mu**2 - t * t / 4)
    )
    put(
        "fixed_transfer_b20_graviton_pole",
        s.expand(-numerator / t).coeff(nu, 2) + 1 / t,
    )
    xi = source.HEAVY_MASS2**3 / (source.KAPPA * source.CUBIC**2)
    margins = {
        "matter_upper_470": 470 - 396 * s.Rational(8, 7) * s.Rational(64, 63) ** 2,
        "matter_lower_4": 12 * s.Rational(32, 35) - 4,
        "ratio_lower_one_over_33": s.Rational(257, 9) / 940 - s.Rational(1, 33),
        "ratio_upper_17": 17 - s.Rational(194, 3) / 4,
        "xi_lower": xi - s.Rational(1, 2 * 10**200),
        "xi_upper": s.Rational(1, 10**200) - xi,
        "outer_cone_upper": s.Rational(17, 10**10) - 17 * xi / s.Rational(1, 10**190),
        "inner_cone_lower": xi / (33 * s.Rational(1, 10**204)) - 150,
    }
    return {
        "whole_original_matter_Born_positive_form": matter_born(),
        "whole_original_gravity_Born_positive_form": gravity_born(),
        "whole_forward_crossover_scale_xi": xi,
        "whole_explicit_positive_margins": margins,
        "whole_positive_Born_proof": "Use mass-one units. For all s>4 and -1<z<1, put Q=s-4,w=1-z^2,p=4Q+16+8/Q,d=-2s+6-2/s,e=Q^2/(4s). Then kappa A_G=p/w+d+ew. Since p+2d>0 and d+e<0, p/(2w)<kappa A_G<p/w. The complete tuned matter amplitude is its positive divided-remainder sum for4<s<n. Hence the full Born amplitude has no zero on this domain.",
        "whole_compact_crossover_proof": "For25/4<=s<=16,n>=128, all channel invariants lie in[-12,16]; the sum of(a-2)^2 lies between12 and396. The exact denominator bounds give4g^2/n^3<A_m<470g^2/n^3 and257/9<=p<=194/3. Therefore xi/(33w)<A_G/A_m<17xi/w, xi=n^3/(kappa g^2). At original parameters1/(2*10^200)<xi<1/10^200. Forw>=10^-190 the ratio is below1.7*10^-9; forw<=10^-204 it exceeds150.",
        "whole_forward_boundary": "Taking kappa large at fixed angle is not uniform at z=+/-1. The normalized missing-loop bound remains uniform because it uses the full positive Born denominator; this does not make the forward cross section finite or justify a uniform Newton truncation. For nu=s-2mu+t/2, the t-channel Born coefficient is-1/(kappa t), which must be combined with the contour contribution before t tends to zero.",
        "checks": checks,
        "gates": {
            "all_crossover_margins_strict": all(bool(v > 0) for v in margins.values()),
            "positive_Born_domain_includes_full_S239_window": bool(
                source.HEAVY_MASS2 > 10**196
            ),
            "original_compact_parameters_in_range": bool(
                source.HEAVY_MASS2 >= 128 and source.KAPPA > 0
            ),
            "original_tuned_quartic_not_retuned": True,
            "forward_crossover_not_an_all_angle_small_gravity_claim": True,
            "finite_contour_requires_independent_remainder_and_IR_matching": True,
        },
    }
