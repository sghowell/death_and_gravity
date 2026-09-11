"""A C2 partition and an explicit one-band Fourier-to-L1 inequality."""

from functools import cache

import sympy as s

y, t, h = s.symbols("y t h", positive=True)


@cache
def data():
    chi = 1 - 10 * t**3 + 15 * t**4 - 6 * t**5
    R = s.symbols("R", nonnegative=True)
    u = s.symbols("u", real=True)
    G0, G2 = s.symbols("G0 G2", positive=True)
    checks = {
        "chi_left": chi.subs(t, 0) - 1,
        "chi_right": chi.subs(t, 1),
        "chi_left_first": s.diff(chi, t).subs(t, 0),
        "chi_right_first": s.diff(chi, t).subs(t, 1),
        "chi_left_second": s.diff(chi, t, 2).subs(t, 0),
        "chi_right_second": s.diff(chi, t, 2).subs(t, 1),
        "chi_monotone": s.factor(s.diff(chi, t) + 30 * t * t * (1 - t) ** 2),
        "chi_first_max": s.Rational(30, 16) - s.Rational(15, 8),
        "chi_second_u_form": s.expand(
            s.diff(chi, t, 2).subs(t, (1 - u) / 2) + 15 * u * (1 - u * u)
        ),
        "chi_second_square_completion": s.expand(
            s.Rational(4, 27)
            - R * (1 - R) ** 2
            - (R - s.Rational(1, 3)) ** 2 * (s.Rational(4, 3) - R)
        ),
        "partition_first_upper": 2 * s.Rational(15, 8) - s.Rational(15, 4),
        "partition_second_upper": s.Integer(4) * 6 - 24,
        "band_second_product_upper": 1 + 2 * s.Rational(15, 4) + 24 - s.Rational(65, 2),
        "Fourier_min_split_integral": s.factor(
            G0 * s.sqrt(G2 / G0) + G2 / s.sqrt(G2 / G0) - 2 * s.sqrt(G0 * G2)
        ),
        "sine_band_norm_square": 16 * s.Rational(3, 2) ** 2 * s.Rational(65, 2) - 1170,
        "small_band_geometric": s.factor(1 / (1 - 1 / s.sqrt(2)) - (2 + s.sqrt(2))),
        "tail_inverse_square_integral": s.integrate(1 / t**2, (t, 1, s.oo)) - 1,
        "finite_partition_telescope": s.Integer(0),
    }
    # This is a symbolic finite telescoping check for arbitrary function values.
    symbols = s.symbols("c0:9")
    checks["finite_partition_telescope"] = s.expand(
        sum(symbols[j] - symbols[j - 1] for j in range(1, len(symbols)))
        - symbols[-1]
        + symbols[0]
    )
    return {
        "chi_transition": chi,
        "definition": "chi(y)=1 for y<=1; chi(y)=1-10(y-1)^3+15(y-1)^4-6(y-1)^5 for 1<y<2; chi(y)=0 for y>=2",
        "partition": "h_j=2^j; phi_j(x)=chi(x/h_j)-chi(2x/h_j); sum over every integer j is1",
        "band_support": "[h/2,2h], length3h/2; the two transition derivatives have disjoint interiors",
        "partition_derivative_bounds": (s.Integer(1), s.Rational(15, 4), s.Integer(24)),
        "hypothesis": "On [h/2,2h], max(|a|,h|a'|,h^2|a''|)<=A_j",
        "g_L1_upper_over_h_A": s.Rational(3, 2),
        "g_second_L1_upper_times_h_over_A": s.Rational(195, 4),
        "band_kernel": "K_j(t)=-2 Im[e^(2it) integral_0^infinity a(x)phi_j(x)e^(itx) dx]",
        "band_kernel_L1_upper_over_A": s.Integer(35),
        "checks": checks,
        "gates": {
            "chi_second_strict_upper": s.Rational(100, 3) < 36,
            "band_constant_rounded_outward": 1170 < 35**2,
            "geometric_sum_below_four": 2 + s.sqrt(2) < 4,
            "sqrt_two_times_four_below_six": 4 * s.sqrt(2) < 6,
        },
    }
