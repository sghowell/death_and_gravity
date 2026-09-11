"""Exact instantaneous-projector comparison after four rotating frames."""

from functools import cache

import sympy as s
from p8_vacuum_superadiabatic_state_energy import energy as previous
from p8_vacuum_superadiabatic_state_energy import tube


def enclosures(mean_mass, amplitude, timescale, multiplicity=6):
    m, d, tau, floor = previous.parameters(mean_mass, amplitude, timescale)
    old = previous.enclosures(m, d, tau, multiplicity)
    return {
        "state_to_fourth_frame_energy_upper": old[
            "uniform_in_out_local_energy_density_difference_upper"
        ],
        "fourth_to_first_frame_energy_upper": s.Rational(
            24 * multiplicity * tube.FACTOR**2, 5 * 9
        )
        * d**2
        / (tau**4 * floor**2),
        "first_frame_second_order_Taylor_remainder_upper": s.Rational(
            multiplicity, 672 * 9
        )
        * d**4
        / (tau**4 * floor**4),
    }


@cache
def data():
    a, b, c, d = s.symbols("theta0 theta1 theta2 theta3", real=True)
    rx = lambda t: s.Matrix(
        [[1, 0, 0], [0, s.cos(t), s.sin(t)], [0, -s.sin(t), s.cos(t)]]
    )
    ry = lambda t: s.Matrix(
        [[s.cos(t), 0, s.sin(t)], [0, 1, 0], [-s.sin(t), 0, s.cos(t)]]
    )
    z = (rx(a) * ry(b) * rx(c) * ry(d) * s.Matrix([0, 0, 1]))[2]
    expected = (
        s.cos(a) * s.cos(b) * s.cos(c) * s.cos(d)
        - s.cos(a) * s.sin(b) * s.sin(d)
        - s.sin(a) * s.sin(c) * s.cos(d)
    )
    x = s.Symbol("x", positive=True)
    q, omega, p, rate = s.symbols("q omega p mass_rate", positive=True)
    first = -omega / s.sqrt(1 + q * q)
    radial1 = s.integrate(x**4 / (1 + x * x) ** s.Rational(7, 2), (x, 0, s.oo))
    radial2 = s.integrate(x**6 / (1 + x * x) ** s.Rational(11, 2), (x, 0, s.oo))
    r = s.Rational(tube.FACTOR, tube.SCALE)
    return {
        "four_frame_Bloch_z": z,
        "first_frame_energy_per_helicity": first,
        "physical_energy_operator": "In the physical mass eigenbasis the energy operator is omega sigma3. The time-dependent basis connection g0 sigma2/tau belongs to the evolution generator, not to that physical energy observable. Nor is -e4/tau the local Dirac energy.",
        "subtracted_energy_per_color_flavor": "-2 omega + p^2 Mdot^2/(4 omega^5)",
        "exact_state_comparison": "In the fourth exact frame the half-line unitary Duhamel integral bounds the in- or out-state projector distance from the instantaneous negative projector by integral_halfline |g4| <= 2^42 pDelta/(tau^4 E^6). The comparison projector need not be an exact solution or a Hadamard state. The exact in/out states are unchanged.",
        "real_frame_angle_bound": "Real e_j>=tau omega>=tau E. With A=2pDelta/(tau E^3), r=1024/(tau E)<=1/1024 and w<=1, |theta_j|<=|g_j/e_j|<=A r^j, j=0..3.",
        "Bloch_comparison_bound": "|z4-cos(theta0)| <= |theta0 theta2|+theta1^2+theta2^2+theta3^2 <= A^2(2r^2+r^4+r^6) <3 A^2 r^2. The exact alternating axes remove a would-be theta0 theta1 term.",
        "first_frame_Taylor_bound": "For u=q0^2>=0, Taylor's integral remainder for (1+u)^(-1/2) has magnitude <=3u^2/8. Here q0=p Mdot/(2omega^3), |Mdot|<=Delta/tau. This is a bound on the full first-frame projector, not a truncation of the exact state.",
        "all_momentum_remainder": "For N color/flavor copies and two occupied helicities, use omega<=2E and radial measure1/(2pi^2). The three absolute allowances are 8NC Delta/(3pi^2 tau^4 m0), 24N1024^2 Delta^2/(5pi^2 tau^4 m0^2), and NDelta^4/(672pi^2 tau^4 m0^4). They bound the full exact-state energy after the displayed two UV terms are subtracted, uniformly for every real time.",
        "checks": {
            "exact_alternating_Bloch_z": s.expand(z - expected),
            "no_second_axis_cross_with_first": s.diff(z, a, b).subs(
                {a: 0, b: 0, c: 0, d: 0}
            ),
            "nonzero_first_third_axis_cross": s.diff(z, a, c).subs(
                {a: 0, b: 0, c: 0, d: 0}
            )
            + 1,
            "nonzero_second_fourth_axis_cross": s.diff(z, b, d).subs(
                {a: 0, b: 0, c: 0, d: 0}
            )
            + 1,
            "first_frame_full_energy_series": s.expand(
                s.series(first, q, 0, 6).removeO()
                + omega
                - omega * q * q / 2
                + 3 * omega * q**4 / 8
            ),
            "evolution_connection_not_a_physical_energy_term": s.simplify(
                -omega * s.sqrt(1 + q * q) - first + omega * q * q / s.sqrt(1 + q * q)
            ),
            "second_order_two_helicity_subtraction": s.simplify(
                (omega * q * q).subs(q, p * rate / (2 * omega**3))
                - p * p * rate * rate / (4 * omega**5)
            ),
            "fourth_first_radial_integral": radial1 - s.Rational(1, 5),
            "first_frame_radial_integral": radial2 - s.Rational(2, 63),
            "fourth_first_full_factor": 2 * 2 * 3 * 4 * s.Rational(1, 2) * radial1
            - s.Rational(24, 5),
            "first_frame_full_factor": 2
            * 2
            * s.Rational(3, 8)
            * s.Rational(1, 16)
            * s.Rational(1, 2)
            * radial2
            - s.Rational(1, 672),
            "state_remainder_mass_dimension": 1 + 4 - 1 - 4,
            "frame_remainder_mass_dimension": 2 + 4 - 2 - 4,
            "Taylor_remainder_mass_dimension": 4 + 4 - 4 - 4,
        },
        "gates": {
            "real_angle_series_ratio_at_most_one_over_1024": bool(
                r <= s.Rational(1, 1024)
            ),
            "Bloch_sum_strictly_below_three": bool(2 + r * r + r**4 < 3),
            "all_remainder_radial_integrals_finite": bool(radial1 > 0 and radial2 > 0),
        },
    }
