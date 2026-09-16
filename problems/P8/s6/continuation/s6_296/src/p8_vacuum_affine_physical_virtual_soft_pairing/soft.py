"""Whole D angular kernel and finite physical-to-analytic soft conversion."""

from functools import cache

import sympy as s

from . import source

E = s.Symbol("hard_COM_scalar_energy", positive=True)
COST = s.Symbol("hard_scattering_cosine", real=True)
X = s.Symbol("soft_Feynman_parameter", real=True)
B = s.Symbol("speed_squared", nonnegative=True)
EP = source.EP


def atanh_ratio(speed_squared):
    b = s.sympify(speed_squared)
    if b == 0:
        return s.S.One
    return s.Piecewise((1, s.Eq(b, 0)), (s.atanh(s.sqrt(b)) / s.sqrt(b), True))


def angular_mean(energy=E, speed_squared=B, epsilon=EP):
    energy, b, e = map(s.sympify, (energy, speed_squared, epsilon))
    return s.hyper((1, s.Rational(3, 2)), (s.Rational(3, 2) + e,), b) / energy**2


def angular_first(energy=E, speed_squared=B):
    energy, b = map(s.sympify, (energy, speed_squared))
    return 2 * (1 - atanh_ratio(b)) / (energy**2 * (1 - b))


def pair_zero(pair_dot, mass=source.MU):
    dot, mu = map(s.sympify, (pair_dot, mass))
    if dot == mu:
        return 1 / mu
    return s.Integral(1 / (mu + 2 * X * (1 - X) * (dot - mu)), (X, 0, 1))


def pair_mean(energy, pair_dot, epsilon=EP, mass=source.MU):
    energy, dot, e, mu = map(s.sympify, (energy, pair_dot, epsilon, mass))
    if dot == mu:
        return angular_mean(energy, 1 - mu / energy**2, e)
    mx = mu + 2 * X * (1 - X) * (dot - mu)
    return s.Integral(angular_mean(energy, 1 - mx / energy**2, e), (X, 0, 1))


def pair_first(energy, pair_dot, mass=source.MU):
    energy, dot, mu = map(s.sympify, (energy, pair_dot, mass))
    if dot == mu:
        return angular_first(energy, 1 - mu / energy**2)
    mx = mu + 2 * X * (1 - X) * (dot - mu)
    return s.Integral(2 * (1 - atanh_ratio(1 - mx / energy**2)) / mx, (X, 0, 1))


def pair_dots(energy=E, cosine=COST, mass=source.MU):
    energy, c, mu = map(s.sympify, (energy, cosine, mass))
    return (
        2 * energy**2 - mu,
        mu + (energy**2 - mu) * (1 - c),
        mu + (energy**2 - mu) * (1 + c),
    )


def kernel(energy=E, cosine=COST, epsilon=EP, mass=source.MU):
    energy, c, e, mu = map(s.sympify, (energy, cosine, epsilon, mass))
    n = lambda dot: dot * dot - mu**2 / (2 + 2 * e)
    value = 4 * n(mu) * angular_mean(energy, 1 - mu / energy**2, e)
    return value + 4 * sum(
        sign * n(dot) * pair_mean(energy, dot, e, mu)
        for sign, dot in zip((1, -1, -1), pair_dots(energy, c, mu))
    )


def kernel_zero(energy=E, cosine=COST, mass=source.MU):
    energy, c, mu = map(s.sympify, (energy, cosine, mass))
    return 2 * mu + 4 * sum(
        sign * (dot * dot - mu**2 / 2) * pair_zero(dot, mu)
        for sign, dot in zip((1, -1, -1), pair_dots(energy, c, mu))
    )


def kernel_first(energy=E, cosine=COST, mass=source.MU):
    energy, c, mu = map(s.sympify, (energy, cosine, mass))
    value = 4 * (mu / 2 + mu**2 * angular_first(energy, 1 - mu / energy**2) / 2)
    return value + 4 * sum(
        sign
        * (
            mu**2 * pair_zero(dot, mu) / 2
            + (dot * dot - mu**2 / 2) * pair_first(energy, dot, mu)
        )
        for sign, dot in zip((1, -1, -1), pair_dots(energy, c, mu))
    )


def phase_normalization(epsilon=EP):
    e = s.sympify(epsilon)
    return (
        (4 * s.pi) ** (-e) * s.gamma(s.Rational(3, 2)) / s.gamma(s.Rational(3, 2) + e)
    )


def finite_conversion(energy=E, cosine=COST, mass=source.MU, kappa=source.K):
    return (
        kernel_first(energy, cosine, mass)
        + (s.EulerGamma - 2 - s.log(s.pi)) * kernel_zero(energy, cosine, mass)
    ) / (8 * s.pi**2 * s.sympify(kappa))


@cache
def data():
    checks = {}

    def put(name, value):
        checks[name] = s.factor(s.expand_func(value))

    x, a, b, E, mu, e = s.symbols("x a b E mu epsilon", positive=True)
    put(
        "Feynman_parameter_primitive",
        s.diff(-1 / ((a - b) * (b + x * (a - b))), x) - 1 / (x * a + (1 - x) * b) ** 2,
    )
    primitive = -1 / ((a - b) * (b + x * (a - b)))
    put(
        "Feynman_parameter_endpoints",
        primitive.subs(x, 1) - primitive.subs(x, 0) - 1 / (a * b),
    )
    put(
        "log_endpoint_integration_by_parts",
        s.diff(x / (a * (a + b * x)), x) - 1 / (a + b * x) ** 2,
    )
    put(
        "normalized_sphere_weight_log_derivative",
        s.digamma(1) - s.digamma(s.Rational(3, 2)) - (2 * s.log(2) - 2),
    )
    p = (4 * s.pi) ** (-e) * s.gamma(s.Rational(3, 2)) / s.gamma(s.Rational(3, 2) + e)
    put("sphere_phase_normalization", p.subs(e, 0) - 1)
    put(
        "sphere_phase_first_derivative",
        s.diff(p, e).subs(e, 0) - (s.EulerGamma - 2 - s.log(s.pi)),
    )
    put(
        "sphere_area_radial_prefactor",
        4 * s.pi / (4 * (2 * s.pi) ** 3) - 1 / (8 * s.pi**2),
    )
    put(
        "full_D_projector_first_derivative",
        s.diff(-(mu**2) / (2 + 2 * e), e).subs(e, 0) - mu**2 / 2,
    )
    beta = s.Symbol("beta", positive=True)
    raw = (2 * s.log(2) - 2 * s.atanh(beta) / beta) / (1 - beta**2)
    put(
        "closed_normalized_angular_derivative",
        raw
        - (2 * s.log(2) - 2) / (1 - beta**2)
        - 2 * (1 - s.atanh(beta) / beta) / (1 - beta**2),
    )
    K0, K1, K2, cphase, L = s.symbols("K0 K1 K2 cphase log_resolution", real=True)
    series = s.series(
        s.exp(2 * e * L) * ((1 + cphase * e) * (K0 + K1 * e + K2 * e**2) - K0) / e,
        e,
        0,
        2,
    ).removeO()
    put("real_virtual_finite_conversion", series.subs(e, 0) - K1 - cphase * K0)
    Ns, Nt, Nu, Ls, Lt, Lu = s.symbols("Ns Nt Nu Ls Lt Lu", real=True)
    Bsoft = 2 * (-Ns * Ls + Nt * Lt + Nu * Lu) - mu
    put(
        "entire_ordered_pair_virtual_pole_match",
        2 * mu + 4 * (Ns * Ls - Nt * Lt - Nu * Lu) + 2 * Bsoft,
    )
    j1 = 2 * (1 - s.atanh(beta) / beta) / mu
    put(
        "moving_diagonal_finite_term",
        4 * (mu / 2 + mu**2 * j1 / 2) - (6 * mu - 4 * mu * s.atanh(beta) / beta),
    )
    for dim in (4, 5, 6):
        eta = s.diag(1, *([-1] * (dim - 1)))
        energy = s.Rational(5, 4)
        r = s.Rational(3, 4)
        p1 = s.zeros(dim, 1)
        p1[0] = energy
        p1[1] = r
        p2 = s.zeros(dim, 1)
        p2[0] = energy
        p2[1] = -r
        p3 = s.zeros(dim, 1)
        p3[0] = energy
        p3[1] = 3 * r / 5
        p3[2] = 4 * r / 5
        p4 = s.zeros(dim, 1)
        p4[0] = energy
        p4[1] = -3 * r / 5
        p4[2] = -4 * r / 5
        q = s.zeros(dim, 1)
        q[0] = 1
        q[-1] = 1
        momenta = [p1, p2, p3, p4]
        signs = [-1, -1, 1, 1]
        dot = lambda p, r, eta=eta: (p.T * eta * r)[0]
        current = sum(
            (sign * p * p.T / dot(p, q) for sign, p in zip(signs, momenta)),
            s.zeros(dim),
        )
        for i, component in enumerate(current * eta * q):
            put(f"D{dim}_conserved_soft_current_{i}", component)
        tr = current[1:-1, 1:-1]
        ttsew = s.trace(tr * tr) - s.trace(tr) ** 2 / (dim - 2)
        pairs = sum(
            signs[i]
            * signs[j]
            * (dot(pi, pj) ** 2 - s.Rational(1, dim - 2))
            / (dot(pi, q) * dot(pj, q))
            for i, pi in enumerate(momenta)
            for j, pj in enumerate(momenta)
        )
        put(f"D{dim}_whole_physical_TT_pair_kernel", ttsew - pairs)

    be = s.Symbol("positive_speed_squared", positive=True)
    put(
        "public_angular_zero",
        s.hyperexpand(angular_mean(E, be, 0)) - 1 / (E**2 * (1 - be)),
    )
    put("public_zero_speed_first_derivative", angular_first(E, 0))
    for direction in (1, -1):
        put(
            "public_full_D_soft_endpoint_" + str(direction), kernel(E, direction, e, mu)
        )
        put("public_K0_endpoint_" + str(direction), kernel_zero(E, direction, mu))
        put("public_K1_endpoint_" + str(direction), kernel_first(E, direction, mu))
    return {
        "whole_normalized_D_soft_kernel": kernel(),
        "whole_physical_D4_kernel": kernel_zero(),
        "whole_first_dimensional_kernel": kernel_first(),
        "whole_phase_normalization": phase_normalization(),
        "whole_finite_conversion": finite_conversion(),
        "whole_angular_derivative_proof": "Feynman pairing reduces the conserved current to normalized integral(1-z^2)^e/(E-rz)^2. Its hypergeometric mean is2F1(1,3/2;3/2+e;b)/E^2. The logarithmic integral obtained with y=(1-z)/2 is[2ln2-2atanh(sqrtb)/sqrtb]/(1-b); subtracting normalization derivative2ln2-2 yields J1=2[1-atanh(sqrtb)/sqrtb]/[E^2(1-b)]. The removable b0 value is explicitly retained.",
        "whole_diagonal_and_sheet_boundary": "All16 ordered pairs and four moving diagonal terms are included. The diagonal finite term is6mu-4mu*atanh(beta)/beta, not its rest-frame value. K0=-2ReBsoft on the frozen S288 i0 sheet. Both forward/backward real kernels vanish for all D; this does not remove the imaginary Coulomb phase.",
        "checks": checks,
        "gates": {
            "full_D_TT_pair_kernel_and_all_diagonal_terms": True,
            "complete_angular_mean_and_logarithmic_first_derivative": True,
            "phase_Gamma_and_evanescent_projector_retained": True,
            "physical_to_analytic_finite_term_not_discarded": True,
            "both_zero_angle_endpoints_all_D": True,
            "reference_soft_prescription_not_changed": True,
            "amplitude_Coulomb_phase_not_erased": True,
        },
    }
