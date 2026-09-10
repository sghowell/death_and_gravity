"""The local-quartic outer bubble after its one overall constant subtraction."""

from functools import cache

import sympy as sp
from p8_vacuum_global_one_loop_insertions import halfplane

N = sp.Integer(6)


def enclosure(mass, yukawa_upper, quartic, Q_lower):
    m, Y, L, Q = map(halfplane.exact, (mass, yukawa_upper, quartic, Q_lower))
    if m < 2 or Y < 0 or Q <= 0:
        raise ValueError("Require m>=2, Y>=0, finite real L and Q>0")
    T = 4 * m * m
    C = 2 * N * Y / Q
    E = 8 * C * L * L / (3 * Q * T * T)
    return {
        "fermion_mass": m,
        "Y_upper": Y,
        "quartic": L,
        "Q_lower": Q,
        "spectral_threshold": T,
        "spectral_factor_upper": C,
        "channel_subtracted_linear_envelope_coefficient": 4 * C * L * L / (Q * T),
        "forward_second_coefficient_upper": E,
        "strictly_positive_family_coefficient": bool(Y > 0 and L != 0),
        "scope": "Only the six local-quartic bubble contributions with one fully subtracted fermion self-energy insertion, not a complete new-model two-loop amplitude.",
    }


@cache
def data():
    x, u, z, v = sp.symbols("x u z v")
    m, Y, L, Q = sp.symbols("m Y L Q", positive=True)
    A = x * (1 - x)
    Delta = x * u + 1 - x - A * z
    Delta2 = Delta.subs(z, 2)
    Bsub = -sp.log(Delta / Delta2)
    C = 2 * N * Y / Q
    T = 4 * m * m
    E = 8 * C * L * L / (3 * Q * T * T)
    checks = {
        "bubble_value_reference_zero": Bsub.subs(z, 2),
        "bubble_first_derivative": sp.factor(sp.diff(Bsub, z) - A / Delta),
        "bubble_second_derivative": sp.factor(sp.diff(Bsub, z, 2) - A * A / Delta**2),
        "halfplane_parameter_gap_margin": sp.expand(
            Delta - x * (u - 4) - (1 - x) - A * (4 - z) - 4 * x * x
        ),
        "center_parameter_gap": sp.expand(
            Delta2 - ((1 - x) ** 2 + x * (u - 1) + x * x)
        ),
        "first_derivative_parameter_weight": sp.integrate(1 - x, (x, 0, 1))
        - sp.Rational(1, 2),
        "second_derivative_parameter_weight": sp.integrate((1 - x) ** 2, (x, 0, 1))
        - sp.Rational(1, 3),
        "crossed_pair_cancels_linear_term": sp.diff(
            Bsub.subs(z, 2 + v) + Bsub.subs(z, 2 - v), v
        ).subs(v, 0),
        "crossed_pair_second_coefficient": sp.factor(
            sp.diff(Bsub.subs(z, 2 + v) + Bsub.subs(z, 2 - v), v, 2).subs(v, 0) / 2
            - sp.diff(Bsub, z, 2).subs(z, 2)
        ),
        "zero_transfer_has_no_forward_derivative": sp.diff(Bsub.subs(z, 0), v, 2),
        "second_coefficient_bound_normalization": sp.factor(
            E - N * Y * L * L / (3 * Q * Q * m**4)
        ),
        "once_subtracted_integral_finite_majorant": sp.factor(
            (L * L / Q) * (4 * C) * (1 / T) - 4 * C * L * L / (Q * T)
        ),
        "second_derivative_integral_finite_majorant": sp.factor(
            (L * L / Q) * (16 * C / 3) * (1 / (2 * T * T)) - E
        ),
    }
    return {
        "symbols": {"x": x, "u": u, "z": z, "L": L, "Q": Q},
        "mixed_mass_parameter_denominator": Delta,
        "once_subtracted_bubble_integrand": Bsub,
        "forward_second_coefficient_integrand": A * A / Delta2**2,
        "analytic_domain": "Re(z)<=4, u>=4mF^2>=16; Delta has strictly positive real part. This covers |z-2|<=2 and t=0.",
        "derivative_majorants": "|B'(z;u)|<=1/[2(u-4)], |B''(z;u)|<=1/[3(u-4)^2].",
        "renormalized_channel": "A_channel(z)-A_channel(2)=(L^2/Q) integral_T^infinity w(u)[B(z;u)-B(2;u)] du.",
        "family_forward_coefficient": "b2_family=(L^2/Q) integral_T^infinity w(u) integral_0^1 A^2/[x u+1-x-2A]^2 dx du >0 for Y>0 and L!=0.",
        "uniform_forward_second_coefficient_upper": E,
        "scope": "Finite once-subtracted outer bubble and its strictly positive forward coefficient for the specified insertion family. Local constant conversion to the full reference scheme has zero b2 here, not in every other graph family.",
        "checks": checks,
    }
