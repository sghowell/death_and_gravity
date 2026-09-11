"""Complete curved-mode transition bounds, including all geometric species."""

from functools import cache

import sympy as s
from p8_vacuum_flat_dirac_hadamard.symbols import rational

from . import geometry, tube


def parameters(mean_mass, amplitude, timescale, count=20):
    n = tube.order(count)
    m, d, tau = map(rational, (mean_mass, amplitude, timescale))
    m0 = s.Rational(99, 100) * m
    if m <= 0 or not 0 <= d <= s.Rational(3, 1000) * m or not 0 < tau <= 1:
        raise ValueError("Require m>0, 0<=Delta/m<=.003 and 0<tau<=1")
    if tau * m0 < 2 * tube.factor(n):
        raise ValueError("The specified finite-order uniform gap condition fails")
    return m, d, tau, m0


def bounds(momentum, mean_mass, amplitude, timescale):
    p = rational(momentum)
    if p < 0:
        raise ValueError("Require nonnegative comoving momentum")
    m, d, tau, m0 = parameters(mean_mass, amplitude, timescale)
    k = s.Integer(tube.factor(20)) ** 20
    compact = 32 * k * p * (d + m * tau) / (tau**20 * (p * p / 16 + m0 * m0) ** 11)
    low = 4 * m * k * p / m0**22
    high = 0 if p == 0 else s.Rational(2048, 15) * m * k / (p**5 * m0**16)
    tail = 0 if p == 0 else min(low, high)
    return {
        "compact": compact,
        "tail_low_momentum": low,
        "tail_high_momentum": high,
        "tail": tail,
        "total": compact + tail,
    }


def uniform(mean_mass, amplitude, timescale):
    m, d, tau, m0 = parameters(mean_mass, amplitude, timescale)
    k = s.Integer(tube.factor(20)) ** 20
    return (
        128 * k * (d + m * tau) / (tau**20 * m0**21)
        + s.Rational(2048, 15) * m * k / m0**21
    )


@cache
def data():
    z = s.Symbol("z", positive=True)
    n = s.Symbol("N", integer=True, positive=True)
    radial = s.integrate(z**5 / (1 + z * z) ** 11, (z, 0, s.oo))
    d = geometry.data()
    m, amp, tau = d["mean_mass"], d["amplitude_upper"], d["mass_transition_time"]
    checks = {
        "full_twenty_frame_tail_radial": radial - s.Rational(1, 720),
        "tail_substitution_factor": 96 * 4**5 * radial - s.Rational(2048, 15),
        "twenty_frame_small_momentum_tail_factor": s.Rational(96, 20 + 4) - 4,
        "all_order_tail_UV_exponent": s.simplify(1 - (n + 4) / 4 + n / 4),
        "all_order_tail_mass_exponent": s.simplify(1 + n / 4 - (n + 2) + 3 * n / 4 + 1),
        "compact_twenty_frame_UV_exponent": 1 - 22 + 21,
        "zero_momentum_exact_transition": bounds(0, m, amp, tau)["total"],
        "uniform_majorant_includes_constant_mass_geometry": uniform(m, 0, tau)
        - 128 * tube.factor(20) ** 20 * m / (tau**19 * (s.Rational(99, 100) * m) ** 21)
        - s.Rational(2048, 15)
        * m
        * tube.factor(20) ** 20
        / (s.Rational(99, 100) * m) ** 21,
    }
    actual = uniform(m, amp, tau)
    return {
        "all_order_complete_transition": "For each finite N satisfying the uniform LE condition, |beta_p|<=integral_R G0*(K_N/(LE))^N dt. The same bound controls either half-line transition relative to the exact N-frame instantaneous projector.",
        "compact_bound_all_N": "32K_N^N p(Delta+m*tau)/(tau^N Ec^(N+2)), Ec=sqrt(p^2/16+m0^2), on |t|<=1.",
        "two_tail_bound_all_N": "96mK_N^N(4/p)^(N/4)*integral_0^(p/4) y^(N/4)/(y^2+m0^2)^((N+2)/2)dy. The integral converges at infinity for every N>=1; thus O(p^(-N/4)).",
        "twenty_frame_tail": "min(4mpK20/m0^22,2048mK20/(15p^5m0^16)), with the continuous value zero at p=0 and K20=(163840)^20.",
        "actual_uniform_transition_upper": actual,
        "actual_uniform_transition_diagnostic": str(s.N(actual, 28)),
        "all_order_vs_finite_order": "N20 gives the explicit numerical bound. Arbitrarily high N, with an N-dependent high-momentum threshold, supplies rapid UV comparison for Hadamard regularity. No fixed finite adiabatic order implies Hadamard.",
        "checks": checks,
        "gates": {
            "actual_twenty_frame_uniform_gap": bool(
                tau * s.Rational(99, 100) * m > 2 * tube.factor(20)
            ),
            "actual_uniform_beta_below_one_e_minus_1890": bool(
                actual < s.Rational(1, 10**1890)
            ),
            "constant_mass_geometric_bound_nonzero": bool(uniform(m, 0, tau) > 0),
        },
    }
