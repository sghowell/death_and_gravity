"""Variable complex disks and an arbitrary finite number of exact rotations."""

from functools import cache

import sympy as s

RADIUS = s.Rational(1, 1024)
ORDER = 20


def order(value):
    if type(value) is not int or value < 1:
        raise TypeError("Require a positive native rotation order")
    return value


def factor(count):
    return 8192 * order(count)


def radius(count, index):
    n = order(count)
    if type(index) is not int or not 0 <= index <= n:
        raise TypeError("Require a native index from zero through N")
    return RADIUS - s.Rational(index, 2048 * n)


def floor(count, index):
    radius(count, index)
    return s.Rational(9, 10) * (1 - s.Rational(1, 4 * count)) ** index


@cache
def data():
    n = ORDER
    r = 2 * RADIUS
    geom = 2 * RADIUS + RADIUS**2
    mass = s.Rational(3, 1000)
    checks = {
        "actual_twenty_frame_factor": factor(n) - 163840,
        "all_order_initial_radius": radius(n, 0) - s.Rational(1, 1024),
        "all_order_final_radius": radius(n, n) - s.Rational(1, 2048),
        "Cauchy_step_factor_slack": factor(n) - 4 * 2048 * n,
        "geometry_relative_square_defect": (s.Rational(256, 255)) ** 4
        - 1
        - (256**4 - 255**4) / s.Integer(255**4),
    }
    for j in range(n):
        checks[f"radius_loss_{j}"] = (
            radius(n, j) - radius(n, j + 1) - s.Rational(1, 2048 * n)
        )
        checks[f"frequency_floor_recursion_{j}"] = floor(n, j + 1) - floor(n, j) * (
            1 - s.Rational(1, 4 * n)
        )
    return {
        "order": n,
        "Cauchy_factor": factor(n),
        "nested_relative_radii": {j: radius(n, j) for j in range(n + 1)},
        "complex_frequency_floors": {j: floor(n, j) for j in range(n + 1)},
        "local_lengths": "L=sqrt(t^2+tau^2), D=1+t^2, q=p/D^2, E=sqrt(q^2+m0^2), m0=.99m; 0<tau<=1, Delta/m<=.003.",
        "complex_profile_bounds": "On |z-t|<=L/1024 the physical SAT8 branch has |s(z/tau)|<2 and |s'(z/tau)|<2w(t/tau). The momentum ratio is near one and |omega(z)|>=.9E.",
        "initial_connection_majorant": "G0=16q/E^2*(Delta/tau*w(t/tau)+mL/D). It satisfies G0/E<=32/(LE). The expansion-induced term is retained even if Delta=0.",
        "exact_induction": "For any positive integer N, K_N=8192N, and LE>=2K_N, e_(j+1)=sqrt(e_j^2+g_j^2), g_(j+1)=(-1)^j partial_t atan(g_j/e_j)/2. On disk j, |g_j|<=G0*(K_N/(LE))^j and |e_j|>=.9*(1-1/(4N))^j E>.5E. Every connection and all exact transition orders are retained.",
        "all_order_high_momentum_condition": "min_t LE>=min(tau*p/4,(m0^3*p/4)^(1/4)). For each fixed N it suffices that p>=max(8K_N/tau,4(2K_N)^4/m0^3). This does NOT require the fixed actual mass-timescale to exceed unbounded N.",
        "checks": checks,
        "gates": {
            "metric_relative_defect_below_one_over_256": bool(
                geom < s.Rational(1, 256)
            ),
            "momentum_square_relative_defect_below_one_over_16": bool(
                s.Rational(256, 255) ** 4 - 1 < s.Rational(1, 16)
            ),
            "profile_eighth_power_defect_below_one_over_32": bool(
                8 * r * (1 + r) ** 7 < s.Rational(1, 32)
            ),
            "profile_complex_modulus_below_two": bool(
                (1 + r) ** 8 * s.Rational(32, 31) < 2**8
            ),
            "profile_derivative_weight_below_two": bool(s.Rational(32, 31) ** 9 < 2**8),
            "mass_square_relative_error_below_one_over_24": bool(
                (4 * mass + 4 * mass**2 + 1 - s.Rational(99, 100) ** 2)
                / s.Rational(99, 100) ** 2
                < s.Rational(1, 24)
            ),
            "total_frequency_defect_below_one_ninth": bool(
                s.Rational(1, 16) + s.Rational(1, 24) < s.Rational(1, 9)
            ),
            "frequency_modulus_above_nine_tenths": bool(
                1 - s.Rational(1, 9) > s.Rational(9, 10) ** 2
            ),
            "all_twenty_one_frequency_floors_above_one_half": all(
                floor(n, j) > s.Rational(1, 2) for j in range(n + 1)
            ),
            "all_order_Bernoulli_floor_above_one_half": bool(
                s.Rational(9, 10) * s.Rational(3, 4) > s.Rational(1, 2)
            ),
            "complex_ratio_cap_allows_sqrt_induction": bool(
                s.Rational(1, 256 * n) ** 2 < s.Rational(1, 4 * n)
            ),
            "initial_majorant_dimensionless_slack": bool(
                16 * (2 * mass + 1) / s.Rational(99, 100) < 32
            ),
        },
    }
