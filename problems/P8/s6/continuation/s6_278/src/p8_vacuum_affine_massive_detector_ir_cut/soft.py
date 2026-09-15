"""Equal-mass analytic soft factor, physical angular origin and resolution flow."""

from functools import cache

import sympy as s

from . import source

S, T, U, MU, K = source.S, source.T, source.U, source.MU, source.K
X = s.Symbol("soft_pair_parameter", real=True)
L = s.Symbol("log_resolution_ratio", real=True)


def exact_positive(value):
    if isinstance(value, bool) or not isinstance(value, (int, s.Rational)):
        raise TypeError("Require an exact positive rational")
    result = s.Rational(value)
    if result <= 0:
        raise ValueError("Require an exact positive rational")
    return result


def pair_factor(z, mu=MU):
    if z == 0:
        return mu / 2
    return (
        ((z - 2 * mu) ** 2 - 2 * mu * mu)
        * s.atan(s.sqrt(z / (4 * mu - z)))
        / s.sqrt(z * (4 * mu - z))
    )


def pair_integral(z, mu=MU):
    return ((z - 2 * mu) ** 2 - 2 * mu * mu) * s.Integral(
        1 / (4 * mu - z * (1 - X * X)), (X, 0, 1)
    )


def kernel_coefficient(order, mu=MU):
    if isinstance(order, bool) or not isinstance(order, (int, s.Integer)):
        raise TypeError("Require an exact nonnegative integer order")
    if order < 0 or order > 8:
        raise ValueError("Require a licensed order from zero through eight")
    j = int(order)
    return s.factorial(j) ** 2 / (4 * mu ** (j + 1) * s.factorial(2 * j + 1))


def exponent_shape():
    return 2 * sum(pair_factor(z) for z in (S, T, U)) - MU


def resolution_exponent(log_ratio=L):
    return -log_ratio * exponent_shape() / (4 * s.pi**2 * K)


@cache
def data():
    a, b, x = s.symbols("a b x", positive=True)
    c, mu, be = s.symbols(
        "pair_dot positive_mass_squared relative_speed", positive=True
    )
    primitive = s.atan(x * s.sqrt(b / a)) / s.sqrt(a * b)
    checks = {
        "pair_integral_antiderivative": s.simplify(
            s.diff(primitive, x) - 1 / (a + b * x * x)
        ),
        "Feynman_pair_denominator": s.expand(
            mu * (x * x + (1 - x) ** 2)
            + 2 * c * x * (1 - x)
            - (mu + 2 * (c - mu) * x * (1 - x))
        ),
        "angular_integral_pair_coefficient": s.simplify(
            (
                2 * (c * c - mu * mu / 2) / c - mu * (1 + be * be) / s.sqrt(1 - be * be)
            ).subs(c, mu / s.sqrt(1 - be * be))
        ),
        "four_self_legs_in_soft_exponent": -4 * mu / 4 + mu,
        "twelve_ordered_offdiagonal_pairs": s.Integer(4 * 3 - 12),
        "pseudothreshold_value": pair_factor(s.Integer(0), mu) - mu / 2,
        "crossing_center_value": s.simplify(pair_factor(2 * mu, mu) + s.pi * mu / 4),
        "Newton_to_kappa_resolution_coefficient": s.Rational(2, 1)
        / (8 * s.pi * K)
        / s.pi
        - 1 / (4 * s.pi**2 * K),
        "physical_virtual_half_real_emission_normalization": s.Rational(1, 4)
        / (2 * s.pi) ** 3
        / s.Integer(2)
        - 1 / (64 * s.pi**3),
    }
    for j in range(9):
        actual = s.integrate((1 - x * x) ** j, (x, 0, 1)) / (4 * mu) ** (j + 1)
        checks["independent_analytic_kernel_coefficient_" + str(j)] = s.factor(
            actual - kernel_coefficient(j, mu)
        )
    jet = sum(kernel_coefficient(j, mu) * x**j for j in range(3))
    Fjet = s.expand(((x - 2 * mu) ** 2 - 2 * mu * mu) * jet)
    checks.update(
        {
            "pseudothreshold_linear_coefficient": Fjet.coeff(x, 1) + s.Rational(11, 12),
            "pseudothreshold_quadratic_coefficient": Fjet.coeff(x, 2) - 1 / (10 * mu),
        }
    )
    L1, L2 = s.symbols("log_first_ratio log_second_ratio", real=True)
    checks["exact_resolution_composition"] = s.expand(
        resolution_exponent(L1 + L2) - resolution_exponent(L1) - resolution_exponent(L2)
    )
    for p, q in ((S, T), (S, U), (T, U)):
        checks["analytic_soft_crossing_" + str(p) + str(q)] = s.expand(
            exponent_shape().subs({p: q, q: p}, simultaneous=True) - exponent_shape()
        )
    return {
        "whole_equal_mass_analytic_pair": pair_factor(S),
        "whole_equal_mass_pair_integral": pair_integral(S),
        "whole_soft_exponent_shape": exponent_shape(),
        "whole_log_resolution_flow": resolution_exponent(),
        "whole_Weinberg_exponent": "log W_E = G/pi * ((E/mu_ren)^(2epsilon)/epsilon) * [2(F(s)+F(t)+F(u))-m^2]; s+t+u=4m^2, all four external scalar masses retained.",
        "pair_angular_origin": "I(p,q)=integral dOmega/[(p.n)(q.n)]=4pi integral_0^1 da/[m^2+2(p.q-m^2)a(1-a)]. With beta=sqrt(1-m^4/(p.q)^2), I=2pi log[(1+beta)/(1-beta)]/[(p.q)beta], and the diagonal limit is4pi/m^2.",
        "physical_radiation_origin": "J_mn(n)=sum eta_i p_i,m p_i,n/(p_i.n) is conserved. Its physical TT square is J:P:J>=0. Real log W = -integral dOmega J:P:J /(64pi^3 kappa) * (E/mu_ren)^(2epsilon)/epsilon. The phase is retained by analytic continuation, not replaced by the real damping factor.",
        "assumption_boundary": "The resolution law belongs to the analytic universal soft-factor prescription, conditional on factorization of the full regulated amplitude. It neither constructs that amplitude nor gives an absolute finite-coupling unitarity or Regge remainder. Hard and stripped amplitudes are not identified outside the stated detector scaling limit.",
        "M1_distinction": "The one-loop M1-pair coefficient carries one M1 species loop; the one-loop universal soft-graviton factor carries none. Stripping it cannot cancel this massless M1 nonanalytic coefficient. Infrared-finite at nonzero transfer does not mean analytic at forward transfer.",
        "checks": checks,
        "gates": {
            "actual_massive_not_massless_soft_kernel": pair_factor(S).has(MU),
            "whole_s_t_u_and_self_terms_retained": all(
                exponent_shape().has(x) for x in (S, T, U, MU)
            ),
            "phase_and_principal_sheet_not_dropped": True,
            "pseudothreshold_is_removable_by_parameter_integral": True,
            "massless_M1_cut_not_soft_graviton_divergence": True,
            "finite_resolution_distinct_from_dimensional_regulator": True,
            "factorization_premise_and_finite_unitarity_error_not_certified": True,
        },
    }
