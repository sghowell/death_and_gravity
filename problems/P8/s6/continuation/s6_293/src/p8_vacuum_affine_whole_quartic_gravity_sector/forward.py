"""Complete crossed finite coefficient and explicit original-domain regulator bound."""

from functools import cache

import sympy as s

from . import proper, source

MU, K, EP = source.MU, source.K, source.EP
E2 = s.Symbol("formal_resolution_squared", positive=True)
X = proper.X


def center_integrands(epsilon=EP, mass=MU):
    epsilon, mass = map(s.sympify, (epsilon, mass))
    y = X * (1 - X)
    A = mass * (1 - 2 * y)
    T = (1 - 2 * epsilon * y) * A ** (epsilon - 1) - mass**2 * (epsilon - 1) * (
        epsilon - 2
    ) * y * y * A ** (epsilon - 3) / (1 + epsilon)
    E = -2 * epsilon * y * y * A ** (epsilon - 1) + 2 * mass * (
        2 + epsilon
    ) * epsilon * (epsilon - 1) * y**3 * A ** (epsilon - 2) / (1 + epsilon)
    return T, E


def center_jets(epsilon=EP, mass=MU):
    return tuple(
        s.Integral(value, (X, 0, 1)) for value in center_integrands(epsilon, mass)
    )


def soft_b20(mass=MU):
    mass = s.sympify(mass)
    return 3 * s.pi / (8 * mass)


def proper_first_jet(mass=MU):
    mass = s.sympify(mass)
    return (
        3 * s.pi * s.log(2 * mass) / 8 - 3 * s.Catalan / 2 + s.Rational(5, 4) - s.pi / 4
    ) / mass


def endpoint_first_jet(mass=MU):
    mass = s.sympify(mass)
    return -(1 - s.pi / 4) / mass


def raw_b20(epsilon=EP, mass=MU, quartic=proper.C, kappa=K, scale_squared=proper.NU2):
    epsilon, mass, quartic, kappa, scale_squared = map(
        s.sympify, (epsilon, mass, quartic, kappa, scale_squared)
    )
    T, E = center_jets(epsilon, mass)
    return (
        quartic
        * s.gamma(-epsilon)
        * (4 * s.pi * scale_squared) ** (-epsilon)
        * (-2 * T + E)
        / (16 * s.pi**2 * kappa)
    )


def formal_soft_divided_b20(
    epsilon=EP,
    mass=MU,
    quartic=proper.C,
    kappa=K,
    scale_squared=proper.NU2,
    resolution_squared=E2,
):
    epsilon, mass, quartic, kappa, scale_squared, resolution_squared = map(
        s.sympify, (epsilon, mass, quartic, kappa, scale_squared, resolution_squared)
    )
    return raw_b20(epsilon, mass, quartic, kappa, scale_squared) - quartic * soft_b20(
        mass
    ) * (resolution_squared / scale_squared) ** epsilon / (
        8 * s.pi**2 * kappa * epsilon
    )


def finite_b20(mass=MU, quartic=proper.C, kappa=K, resolution_squared=E2):
    mass, quartic, kappa, resolution_squared = map(
        s.sympify, (mass, quartic, kappa, resolution_squared)
    )
    return (
        quartic
        * (
            3
            * s.pi
            * (s.EulerGamma + s.log(mass / (2 * s.pi * resolution_squared)) - 1)
            / 8
            - 3 * s.Catalan / 2
            + s.Rational(7, 4)
        )
        / (8 * s.pi**2 * kappa * mass)
    )


def require_resolution(resolution_squared):
    value = source.require_mass(resolution_squared)
    if not s.Rational(1, 4) <= value <= 1:
        raise ValueError("Require original mu=nu^2=1 and1/4<=E^2<=1")
    return value


def require_quartic(value):
    if isinstance(value, bool) or not isinstance(value, (int, s.Rational)):
        raise TypeError("Require an exact real rational quartic")
    return s.Rational(value)


def finite_magnitude_bound(quartic, kappa, resolution_squared):
    quartic = require_quartic(quartic)
    kappa = source.require_mass(kappa)
    require_resolution(resolution_squared)
    return abs(quartic) / (4 * kappa)


def regulator_remainder(epsilon, quartic, kappa, resolution_squared):
    epsilon = source.require_mass(epsilon)
    if epsilon > s.Rational(1, 8):
        raise ValueError("Require0<epsilon<=1/8 in the original massive domain")
    quartic = require_quartic(quartic)
    kappa = source.require_mass(kappa)
    require_resolution(resolution_squared)
    return 8 * epsilon * abs(quartic) / kappa


@cache
def data():
    a, mu, e, x, z = s.symbols(
        "channel mass_squared epsilon parameter root_coordinate", positive=True
    )
    st, tt = s.symbols("s t", real=True)
    uu = 4 * mu - st - tt
    y = x * (1 - x)
    A = mu - a * y
    V = (a - 2 * mu) ** 2 - 2 * mu**2 / (1 + e)
    AA = s.Symbol("positive_denominator", positive=True)
    pair = V * AA ** (e - 1) / 2 + (a - 2 * mu) * AA**e
    endpoint = (a + 2 * mu / (1 + e)) * y * AA**e
    differentiate = lambda f: s.diff(f, a) - y * s.diff(f, AA)
    checks = {}
    T2 = differentiate(differentiate(pair)).subs(a, 2 * mu)
    center = mu * (1 - 2 * y)
    T2base = (1 - 2 * e * y) * AA ** (e - 1) - mu**2 * (e - 1) * (e - 2) / (
        1 + e
    ) * y * y * AA ** (e - 3)
    T2expected = T2base.subs(AA, center)
    checks["entire_D_twice_differentiated_proper_moment"] = s.factor(
        (T2 - T2base) / AA**e
    )
    E2 = differentiate(differentiate(endpoint)).subs(a, 2 * mu)
    E2base = -2 * e * y * y * AA ** (e - 1) + 2 * mu * (2 + e) / (1 + e) * e * (
        e - 1
    ) * y**3 * AA ** (e - 2)
    E2expected = E2base.subs(AA, center)
    checks["entire_D_twice_differentiated_endpoint"] = s.factor((E2 - E2base) / AA**e)
    T20 = 1 / center - 2 * mu**2 * y * y / center**3
    T21 = (
        s.log(center) / center
        - 2 * y / center
        + mu**2 * y * y * (5 - 2 * s.log(center)) / center**3
    )
    checks["whole_proper_zeroth_forward_jet"] = s.factor(T2expected.subs(e, 0) - T20)
    checks["whole_proper_first_forward_jet"] = s.factor(
        s.diff(T2expected, e).subs(e, 0) - T21
    )
    checks["whole_endpoint_zero_forward_jet"] = E2expected.subs(e, 0)
    E21 = -2 * y * y / center - 4 * mu * y**3 / center**2
    checks["whole_endpoint_first_forward_jet"] = s.factor(
        s.diff(E2expected, e).subs(e, 0) - E21
    )
    checks["complete_endpoint_UV_is_constant_on_shell"] = s.expand(
        sum((q + 2 * mu) / 6 for q in (st, tt, uu)) - 5 * mu / 3
    )
    Mst, Mtt, Muu = s.symbols("M_s M_t M_u")
    Tzero = (
        sum(
            2 * ((q - 2 * mu) ** 2 - 2 * mu**2) * m + (q - 2 * mu)
            for q, m in ((st, Mst), (tt, Mtt), (uu, Muu))
        )
        + mu
    )
    soft = (
        2
        * sum(
            ((q - 2 * mu) ** 2 - 2 * mu**2) * m
            for q, m in ((st, Mst), (tt, Mtt), (uu, Muu))
        )
        - mu
    )
    checks["whole_six_pair_proper_pole_is_massive_four_leg_soft"] = s.expand(
        Tzero - soft
    )
    checks["whole_quartic_pair_contact_LSZ_UV"] = 4 * mu - 12 * mu + 8 * mu
    D = s.Symbol("dimension", positive=True)
    trace = 2 * mu + (D - 2) * a / 2
    qp = -a * trace + (D - 1) * a * trace / (D - 2)
    checks["complete_D_two_endpoint_projector"] = s.factor(
        2 * qp / a - a - 4 * mu / (D - 2)
    )
    B = s.Symbol("entire_B0")
    Delta = mu - a * y
    A0 = 2 * Delta * B / (D - 2)
    loop2 = A0 + Delta * B
    eta = (2 / D - 1) * loop2 + (mu + a * y) * B
    checks["complete_bubble_transverse_eta"] = s.factor(eta - 2 * a * y * B)
    checks["literal_full_quartic_bubble_phase_half"] = (
        s.I * (-s.I) * s.I**2 * s.I / 2 + s.I / 2
    )
    fB = (1 + 6 * z * z + z**4) / (1 + z * z) ** 3
    checks["whole_soft_forward_rational_kernel"] = s.factor(
        mu * T20.subs(x, (1 + z) / 2) - fB
    )
    primitive = 3 * s.atan(z) / 2 - z * (1 - z * z) / (2 * (1 + z * z) ** 2)
    checks["whole_soft_forward_primitive"] = s.factor(s.diff(primitive, z) - fB)
    checks["whole_soft_forward_integral"] = (
        primitive.subs(z, 1) - primitive.subs(z, 0) - 3 * s.pi / 8
    )
    rat = s.Rational(5, 2) * (1 - z * z) ** 2 / (1 + z * z) ** 3 - (1 - z * z) / (
        1 + z * z
    )
    checks["whole_proper_log_and_rational_decomposition"] = s.expand_log(
        s.factor(
            mu * T21.subs(x, (1 + z) / 2) - fB * s.log(mu * (1 + z * z) / 2) - rat
        ),
        force=True,
    )
    checks["whole_rational_remainder_integral"] = s.integrate(rat, (z, 0, 1)) - (
        1 - 3 * s.pi / 16
    )
    Dkernel = (1 - z * z) ** 2 / (2 * (1 + z * z) ** 2)
    checks["whole_endpoint_forward_rational_kernel"] = s.factor(
        -mu * E21.subs(x, (1 + z) / 2) - Dkernel
    )
    checks["whole_endpoint_forward_integral"] = s.integrate(Dkernel, (z, 0, 1)) - (
        1 - s.pi / 4
    )
    theta = s.Symbol("theta", real=True)
    checks["whole_log_angular_substitution"] = s.trigsimp(
        fB.subs(z, s.tan(theta)) / s.cos(theta) ** 2
        - (s.Rational(3, 2) - s.cos(4 * theta) / 2)
    )
    ip = s.sin(2 * theta) / 4 - theta / 4 - s.sin(4 * theta) / 16
    checks["whole_log_integration_by_parts_primitive"] = s.trigsimp(
        s.diff(ip, theta) - s.sin(4 * theta) * s.tan(theta) / 4
    )
    checks["whole_log_integration_by_parts_value"] = (
        ip.subs(theta, s.pi / 4) - ip.subs(theta, 0) - (s.Rational(1, 4) - s.pi / 16)
    )
    Icos = -s.pi * s.log(2) / 4 + s.Catalan / 2
    Ilog = -3 * s.pi * s.log(2) / 8 - 3 * Icos + s.Rational(1, 4) - s.pi / 16
    T21closed = (
        3 * s.pi * s.log(2 * mu) / 8 - 3 * s.Catalan / 2 + s.Rational(5, 4) - s.pi / 4
    )
    checks["whole_complete_proper_finite_forward_constant"] = s.expand_log(
        s.expand(3 * s.pi * s.log(mu) / 8 + Ilog + 1 - 3 * s.pi / 16 - T21closed),
        force=True,
    )
    nu2, E2scale = s.symbols("loop_scale_squared resolution_squared", positive=True)
    total = (
        T21closed
        + (s.EulerGamma - s.log(4 * s.pi * nu2)) * 3 * s.pi / 8
        + (1 - s.pi / 4) / 2
        - 3 * s.pi * s.log(E2scale / nu2) / 8
    )
    closed = (
        3 * s.pi * (s.EulerGamma + s.log(mu / (2 * s.pi * E2scale)) - 1) / 8
        - 3 * s.Catalan / 2
        + s.Rational(7, 4)
    )
    checks["whole_known_soft_divided_quartic_forward_coefficient"] = s.expand_log(
        s.expand(total - closed), force=True
    )
    checks["whole_known_soft_divided_loop_scale_cancels"] = s.diff(total, nu2)
    checks["whole_proper_log_and_rational_decomposition"] = s.simplify(
        checks["whole_proper_log_and_rational_decomposition"]
    )
    e = s.Symbol("epsilon", nonnegative=True)
    x = s.Symbol("unit_parameter", real=True)
    a = s.Symbol("positive_power", positive=True)
    y = x * (1 - x)
    A = 1 - 2 * y
    g = (e - 1) * (e - 2) / (1 + e)
    h = (e + 2) * (e - 1) / (1 + e)
    bound_checks = {
        "whole_g_rational_division": s.factor(g - (e - 4 + 6 / (1 + e))),
        "whole_two_minus_g_margin": s.factor(2 - g - e * (5 - e) / (1 + e)),
        "whole_g_first_derivative": s.factor(s.diff(g, e) - 1 + 6 / (1 + e) ** 2),
        "whole_g_second_derivative": s.factor(s.diff(g, e, 2) - 12 / (1 + e) ** 3),
        "whole_h_rational_division": s.factor(h - e + 2 / (1 + e)),
        "whole_h_plus_two_margin": s.factor(h + 2 - e * (e + 3) / (1 + e)),
        "whole_h_first_derivative": s.factor(s.diff(h, e) - 1 - 2 / (1 + e) ** 2),
        "whole_h_second_derivative": s.factor(s.diff(h, e, 2) + 4 / (1 + e) ** 3),
        "whole_positive_center_gap": s.expand(
            A - s.Rational(1, 2) - 2 * (x - s.Rational(1, 2)) ** 2
        ),
        "whole_quartic_original_contact_majorant": s.S.Zero,
    }
    n = s.Symbol("heavy_mass_squared", positive=True)
    bound_checks["whole_quartic_original_contact_majorant"] = s.factor(
        4 / n - (3 / (n - 2) - 2 / (n - 2) ** 2) - (n - 4) ** 2 / (n * (n - 2) ** 2)
    )
    for j in range(3):
        bound_checks["whole_low_Gamma_log_moment_" + str(j)] = s.diff(1 / a, a, j) * (
            -1
        ) ** j - s.factorial(j) / a ** (j + 1)
    gamma = [
        s.factorial(j) / (s.Rational(7, 8)) ** (j + 1) + s.factorial(j)
        for j in range(3)
    ]
    margins = {
        "gamma_zero_below3": 3 - gamma[0],
        "gamma_first_below3": 3 - gamma[1],
        "gamma_second_below5": 5 - gamma[2],
        "exp3_above16_by_finite_positive_sum": sum(
            s.Rational(3) ** j / s.factorial(j) for j in range(5)
        )
        - 16,
        "exp2_above4_by_finite_positive_sum": sum(
            s.Rational(2) ** j / s.factorial(j) for j in range(3)
        )
        - 4,
        "exp3_above8": sum(s.Rational(3) ** j / s.factorial(j) for j in range(4)) - 8,
        "proper_integrand_base_below2": 2 - (1 + s.Rational(1, 16) + s.Rational(1, 2)),
        "proper_first_base_below2": 2 - (s.Rational(1, 2) + s.Rational(5, 4)),
        "endpoint_first_base_below_half": s.Rational(1, 2)
        - (s.Rational(1, 4) + s.Rational(3, 128)),
        "endpoint_second_base_below_half": s.Rational(1, 2)
        - (s.Rational(3, 8) + s.Rational(1, 32)),
        "endpoint_first_derivative_below2": 2
        - 2 * (s.Rational(1, 32) + s.Rational(1, 2)),
        "endpoint_second_derivative_below4": 4
        - 2 * (s.Rational(1, 32) + 1 + s.Rational(1, 2)),
        "combined_uniform_second_below1024": 1024
        - (50 * 9 + 2 * 12 * 18 + 3 * 40 + 12),
        "finite_bracket_below8": 8 - (6 + s.Rational(7, 4)),
    }
    actualg = s.Rational(1, 8192)
    actualn = s.Rational(10**200, 512) + 2
    actualk = s.Integer(10) ** 800
    actualC = -(actualg**2) * (3 / (actualn - 2) - 2 / (actualn - 2) ** 2)
    margins["actual_contact_bound"] = 4 * actualg**2 / actualn - abs(actualC)
    margins["actual_known_b20_below_1e_minus1005"] = s.Rational(
        1, 10**1005
    ) - actualg**2 / (actualk * actualn)
    checks.update(
        {"explicit_bound_" + name: value for name, value in bound_checks.items()}
    )
    variable = s.Symbol("crossed_forward_variable", real=True)
    checks["whole_constant_quartic_counterterm_zero_b20"] = s.diff(
        s.Symbol("unmatched_constant_quartic"), variable, 2
    )
    checks["whole_constant_curvature_anchor_zero_b20"] = s.diff(
        (4 * source.MU + variable) + (4 * source.MU - variable) + 2 * source.MU,
        variable,
        2,
    )
    return {
        "whole_D_second_channel_integrands": center_integrands(),
        "whole_D_raw_forward_coefficient": raw_b20(),
        "whole_formal_S278_soft_divided_coefficient": formal_soft_divided_b20(),
        "whole_exact_first_coefficients": (
            soft_b20(),
            proper_first_jet(),
            endpoint_first_jet(),
        ),
        "whole_complete_known_finite_quartic_coefficient": finite_b20(),
        "whole_actual_finite_magnitude_majorant": finite_magnitude_bound(
            source.CONTACT, source.KAPPA, s.S.One
        ),
        "whole_actual_uniform_regulator_majorant": 8
        * source.EP
        * abs(source.CONTACT)
        / source.KAPPA,
        "whole_exact_positive_bound_margins": margins,
        "whole_Catalan_integral_proof": "The center has A=mu(1-2x(1-x))>=mu/2. Reflection and x=(1+z)/2 give the complete soft kernel(1+6z^2+z^4)/(1+z^2)^3 with integral3pi/8. With z=tan(theta), its logarithmic moment reduces to int_0^(pi/4)(3/2-cos4theta/2)ln(sec(theta)^2/2). The defining Catalan integral gives intlncos=-pi ln2/4+Catalan/2; integration by parts gives intcos4theta lncos=1/4-pi/16. Together with the entire rational remainder this fixes T2'(0). The endpoint derivative is-(1-pi/4)/mu, including its full D coefficient.",
        "whole_fixed_soft_prescription_and_scale_proof": "Subtract only the linear C coefficient of the already inherited S278 analytic factor logW_E=Bsoft*(E/nu)^(2e)/(8pi^2 kappa e). This fixes the displayed formal finite expression; nu cancels and E remains explicit. It is not a newly chosen finite counterterm or a proof of physical factorization/unitarity. Changing E obeys the inherited resolution law and can change this subset's sign.",
        "whole_explicit_uniform_remainder_proof": "For original mu=nu^2=1,1/4<=E^2<=1 and0<e<=1/8, let R_e=Gamma(1-e)(4pi)^(-e)[2T2(e)-E2(e)]-2B2*(E^2)^e. R0=0. On the complete parameter interval A>=1/2, the zeroth/first/second e derivative bounds are T2:(4,8,18), E2:(1,2,4). Splitting the Gamma integral at1 gives|Gamma^(j)(1-e)|<j!/(7/8)^(j+1)+j!, hence bounds(3,3,5). Since ln4pi<3, its scale-weighted bounds are(3,12,50). Therefore|R''|<50*9+2*12*18+3*40+12=1014<1024. Taylor's theorem bounds the complete finite-regulator error by32e*absC/(pi^2 kappa)<8e*absC/kappa. This is an evaluated uniform majorant, not an unspecified derivative supremum or physical detector error.",
        "whole_actual_finite_bound_proof": "For the stated original resolution interval, -3<ln[1/(2piE^2)]<0,0<EulerGamma<1 and0<Catalan<1. The absolute finite bracket is<8, soabsb20<absC/(4kappa). The exact original contact obeysabsC<4g^2/n, with differenceg^2(n-4)^2/[n(n-2)^2]. Its unchanged numerical parameters giveg^2/(kappa*n)<10^-1005. No E is chosen to obtain a sign.",
        "whole_remaining_matching_and_physical_scope": "Finite constant C and constant R Phi^2 matching havezero crossed b20 at this insertion but remain in the full amplitude. Higher-derivative EFT matching, all other coupling sectors, full physical Coulomb/inclusive/dressed IR, fixed-transfer Regge/all-loop bounds and original V/G/B/P8 remain open. This known formal coefficient is not an isolated positivity verdict.",
        "checks": checks,
        "gates": {
            "whole_D_center_jets_before_epsilon_expansion": True,
            "entire_Catalan_and_endpoint_finite_integrals": True,
            "existing_S278_soft_prescription_not_new_finite_choice": True,
            "explicit_uniform_regulator_bound_not_unevaluated_supremum": all(
                value > 0 for value in margins.values()
            ),
            "actual_finite_bound_without_resolution_sign_choice": True,
            "constant_anchors_zero_b20_not_full_matching_closure": True,
        },
    }
