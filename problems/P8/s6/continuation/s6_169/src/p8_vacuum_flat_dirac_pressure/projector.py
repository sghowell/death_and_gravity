"""The physical pressure observable and the full fourth-frame x component."""

from functools import cache

import sympy as s
from p8_vacuum_superadiabatic_state_energy import tube


@cache
def data():
    a, b, c, d = s.symbols("theta0 theta1 theta2 theta3", real=True)
    rx = lambda t: s.Matrix(
        [[1, 0, 0], [0, s.cos(t), s.sin(t)], [0, -s.sin(t), s.cos(t)]]
    )
    ry = lambda t: s.Matrix(
        [[s.cos(t), 0, s.sin(t)], [0, 1, 0], [-s.sin(t), 0, s.cos(t)]]
    )
    nx = (rx(a) * ry(b) * rx(c) * ry(d) * s.Matrix([0, 0, 1]))[0]
    p, M, v, w = s.symbols("p M Mdot Mddot", positive=True)
    omega = s.sqrt(p * p + M * M)
    q0 = p * v / (2 * omega**3)
    D = lambda expr: s.diff(expr, M) * v + s.diff(expr, v) * w
    leading = D(q0) / (2 * omega)
    explicit = p / (4 * omega**4) * (w - 3 * M * v * v / omega**2)
    q1 = leading / (1 + q0 * q0) ** s.Rational(3, 2)
    pressure0 = -2 * p * p / (3 * omega)
    pressure2 = p * p * q0 * q0 / (3 * omega) - 2 * p * M * leading / (3 * omega)
    expected = -p * p * M * w / (6 * omega**5) + p * p * (p * p + 6 * M * M) * v * v / (
        12 * omega**7
    )
    theta = s.Symbol("theta", real=True)
    sigma1 = s.Matrix([[0, 1], [1, 0]])
    sigma2 = s.Matrix([[0, -s.I], [s.I, 0]])
    sigma3 = s.diag(1, -1)
    U = s.cos(theta / 2) * s.eye(2) - s.I * s.sin(theta / 2) * sigma2
    mass_rotated = s.simplify(U.H * sigma1 * U)
    r = s.Rational(tube.FACTOR, tube.SCALE)
    Amax = 2 * s.Rational(1, 100) / tube.SCALE
    return {
        "physical_pressure_per_helicity_in_mass_basis": "Pi=p^2/(3omega) sigma3+pM/(3omega) sigma1. It is not the rotating-frame evolution generator or a frequency eigenvalue.",
        "fourth_frame_Bloch_x": nx,
        "leading_second_derivative_angle": leading,
        "exact_first_angle_ratio": q1,
        "two_helicity_zeroth_pressure_subtraction": pressure0,
        "two_helicity_second_pressure_subtraction": pressure2,
        "exact_angle_dictionary": "q0=pMdot/(2omega^3), b=q0dot/(2omega). The exact q1=b/(1+q0^2)^(3/2), theta1=atan(q1). This retains the derivative connection and the first frequency correction.",
        "uniform_x_remainder_proof": "With A=2pDelta/(tau E^3), r=1024/(tau E)<=1/1024, |theta_j|<=A r^j, |q1|<=Ar and |q0|<=A<1/4. Use |nx-sin(theta1)|<=|theta3|+|theta1|(theta2^2+theta3^2)/2, |sin(theta1)-theta1|<=|theta1|^3/6, |atan(q1)-q1|<=|q1|^3/3, and |b-q1|<=2|q1|q0^2. Thus |nx-b|<=A r^3+3A^3 r.",
        "diagonal_and_state_remainders": "The physical pressure operator has norm p/3<=omega/3. Its exact-state versus fourth-frame allowance, and its z-component frame/Taylor allowance, are each at most one third of the corresponding frozen S6.168 energy allowance.",
        "off_diagonal_all_momentum_bound": "Use |pM/omega|<=M<2m, not a constant multiple of E. Both helicities and radial measure give the x-remainder allowances Nm1024^3Delta/(3pi^2tau^4m0^2) and 2Nm1024Delta^3/(3pi^2tau^4m0^4). They follow from integral p^3/E^6=1/(4m0^2), integral p^5/E^10=1/(24m0^4).",
        "checks": {
            "exact_fourth_frame_x": s.expand(
                nx - s.cos(b) * s.sin(d) - s.sin(b) * s.cos(c) * s.cos(d)
            ),
            "first_x_rotation_does_not_change_Bloch_x": s.diff(nx, a),
            "physical_mass_rotation_of_pressure_matrix": sum(
                s.simplify(v) ** 2
                for v in (mass_rotated - s.cos(theta) * sigma1 - s.sin(theta) * sigma3)
            ),
            "leading_derivative_angle_exact": s.simplify(leading - explicit),
            "exact_q1_frequency_and_connection": s.simplify(
                D(s.atan(q0)) / (2 * omega * s.sqrt(1 + q0 * q0)) - q1
            ),
            "complete_second_pressure_subtraction": s.simplify(pressure2 - expected),
            "negative_instantaneous_vacuum_pressure": s.simplify(
                pressure0 + 2 * omega / 3 - 2 * M * M / (3 * omega)
            ),
            "pressure_mode_derivative_dimensions": 2 + 1 + 3 - 5 - 1,
        },
        "gates": {
            "whole_domain_angle_cap_below_one_quarter": bool(Amax < s.Rational(1, 4)),
            "real_power_derivative_below_two": bool(
                s.Rational(9, 4) * (1 + s.Rational(1, 16)) < 4
            ),
            "entire_x_remainder_coefficient_below_three": bool(
                2 + (r * r + r**4 + r**6) / 2 < 3
            ),
            "fixed_ratio_below_one_over_1024": bool(r <= s.Rational(1, 1024)),
        },
    }
