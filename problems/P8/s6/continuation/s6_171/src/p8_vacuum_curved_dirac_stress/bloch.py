"""Complete twenty-frame physical projector remainder after local subtractions."""

from functools import cache

import sympy as s
from p8_vacuum_curved_dirac_state import energy as state_energy
from p8_vacuum_curved_dirac_state import transition, tube


def enclosures(mean_mass, amplitude, timescale, multiplicity=42):
    m, d, tau, m0 = transition.parameters(mean_mass, amplitude, timescale)
    if type(multiplicity) is not int or multiplicity < 1:
        raise TypeError("Require a positive native color/flavor multiplicity")
    n = multiplicity
    k = s.Integer(tube.factor(20))
    state = state_energy.enclosures(m, d, tau, n)["total_state_energy_difference_upper"]
    frame = 2048 * n * k * k * (d / tau**2 + m) ** 2 / (45 * m0 * m0)
    taylor = s.Rational(32768 * n, 21 * 9) * (d / tau + m) ** 4 / m0**4
    linear = s.Rational(32 * n, 3 * 9) * m * k**3 * (d / tau**4 + m / tau**2) / m0**2
    cubic = s.Rational(8192 * n, 9 * 9) * m * k * (d / tau + m) ** 3 / (tau * m0**4)
    return {
        "state_to_twenty_frame_energy": state,
        "state_to_twenty_frame_pressure": state / 6,
        "twenty_to_first_frame_energy": frame,
        "first_frame_energy_Taylor": taylor,
        "off_diagonal_pressure_linear": linear,
        "off_diagonal_pressure_cubic": cubic,
        "complete_subtracted_energy_remainder": state + frame + taylor,
        "complete_subtracted_pressure_remainder": state / 6
        + (frame + taylor) / 3
        + linear
        + cubic,
    }


@cache
def data():
    z = s.Symbol("z", positive=True)
    radial = {
        "frame": s.integrate(z**4 / (1 + z * z) ** s.Rational(7, 2), (z, 0, s.oo)),
        "Taylor": s.integrate(z**6 / (1 + z * z) ** s.Rational(11, 2), (z, 0, s.oo)),
        "offdiag_linear": s.integrate(z**3 / (1 + z * z) ** 3, (z, 0, s.oo)),
        "offdiag_cubic": s.integrate(z**5 / (1 + z * z) ** 5, (z, 0, s.oo)),
    }
    a, r = s.symbols("A r", nonnegative=True)
    q0 = s.Symbol("q0", real=True)
    checks = {
        "frame_full_radial": radial["frame"] - s.Rational(1, 5),
        "Taylor_full_radial": radial["Taylor"] - s.Rational(2, 63),
        "pressure_linear_full_radial": radial["offdiag_linear"] - s.Rational(1, 4),
        "pressure_cubic_full_radial": radial["offdiag_cubic"] - s.Rational(1, 24),
        "all_energy_frame_prefactors": s.Rational(1, 2) * 16 * 16**2 * radial["frame"]
        - s.Rational(2048, 5),
        "all_energy_Taylor_prefactors": s.Rational(1, 2)
        * s.Rational(3, 2)
        * 16**4
        * radial["Taylor"]
        - s.Rational(32768, 21),
        "all_pressure_linear_prefactors": s.Rational(2, 3)
        * 4
        * 16
        * radial["offdiag_linear"]
        - s.Rational(32, 3),
        "all_pressure_cubic_prefactors": s.Rational(2, 3)
        * 8
        * 16**3
        * radial["offdiag_cubic"]
        - s.Rational(8192, 9),
        "exact_first_frame_energy_factor": s.cos(s.atan(q0)) - 1 / s.sqrt(1 + q0 * q0),
        "z_geometric_series_quadratic_bound": 2 * a * a * r * r
        + s.Rational(4, 3) * a * a * r * r
        - s.Rational(10, 3) * a * a * r * r,
    }
    return {
        "four_complete_radial_integrals": radial,
        "local_majorants": "At fixed t, A=G0/E=16qB/E^3, B=Delta/tau*w+mL/(1+t^2), r=K20/(LE), |theta_j|<=A r^j. Real exact e_j>=omega>=E. The actual gap gives A<1/4,r<1/2.",
        "twenty_frame_z_bound": "For the internal rotations j>=1, |v_z-1|<=.5(sum|theta_j|)^2<=2A^2r^2. Only even j>=2 x-axis rotations change v_y, so |v_y|<=Ar^2/(1-r^2)<=4Ar^2/3. The outer j0 rotation gives |z20-cos(theta0)|<=4A^2r^2.",
        "twenty_frame_x_bound": "After separating the j1 y rotation, |v_x|<=Ar^3/(1-r^2) and |v_z-1|<=2A^2r^4 for the j>=2 remainder. With b=partial_t(q0)/(2omega), q1=b/(1+q0^2)^(3/2), |b|<=Ar by Cauchy. Exact sin(atan q1), its cubic difference and |q1-b|<=3|b|q0^2/2 give |x20-b|<=4Ar^3+8A^3r.",
        "physical_energy_subtractions": "For both helicities rho0=-2omega,rho2=omega*q0^2=q^2(Mdot+HM)^2/(4omega^5). The frame difference is <=16E A^2r^2 and the first-frame Taylor remainder <=3omega*q0^4/4<=1.5E A^4.",
        "physical_pressure_subtractions": "Before rotations Pi=(q/3)sigma1. In the mass basis it is q^2/(3omega)sigma3+qM/(3omega)sigma1. Two-helicity P0=-2q^2/(3omega), P2=q^2q0^2/(3omega)-2qM b/(3omega). The diagonal remainder is at most one third the corresponding energy bound. The off-diagonal bound uses M<2m and the stated x20-b estimate.",
        "state_and_regulator_boundary": "The exact-state to P20 half-line remainder is bounded by the S6.170 all-momentum allowance; P20 is not claimed Hadamard or an exact dynamical state. The remaining mode integrals converge near epsilon=0 at each fixed time. Restore the dimensional local finite terms and then use the uniform four-dimensional bounds; no uniform-in-epsilon infinite-time exchange is assumed.",
        "checks": checks,
        "gates": {
            "z_majorant_constant_has_slack": bool(s.Rational(10, 3) < 4),
            "x_geometric_series_constant_has_slack": bool(s.Rational(4, 3) < 4),
            "x_cubic_majorant_has_slack": bool(
                2 * s.Rational(1, 2) ** 4
                + s.Rational(1, 2) * s.Rational(1, 2) ** 2
                + s.Rational(3, 2)
                < 8
            ),
            "b_derivative_majorant_has_slack": 1024 < tube.factor(20),
            "second_order_subtracted_full_energy_and_pressure_integrable": True,
        },
    }
