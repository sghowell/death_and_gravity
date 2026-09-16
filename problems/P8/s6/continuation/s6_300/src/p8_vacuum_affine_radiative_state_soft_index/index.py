"""Entire finite massless-pair soft index and its angular-energy representation."""

from functools import cache

import sympy as s

from . import recoil, source

SIGNS = (-1, -1, 1, 1)
D = s.Symbol("future_massive_pair_dot", positive=True)
X = source.previous.soft.X
B = s.Symbol("null_direction_gap", nonnegative=True)


def massive_pair(pair_dot, mass=source.MU):
    d, mu = map(s.sympify, (pair_dot, mass))
    if d == mu:
        return mu / 2
    return (d * d - mu * mu / 2) * s.acosh(d / mu) / s.sqrt(d * d - mu * mu)


def massive_pair_integral(pair_dot, mass=source.MU):
    d, mu = map(s.sympify, (pair_dot, mass))
    return (d * d - mu * mu / 2) * s.Integral(
        1 / (mu + 2 * X * (1 - X) * (d - mu)), (X, 0, 1)
    )


def null_pair(gap=B):
    b = s.sympify(gap)
    if b == 0:
        return s.S.Zero
    return b * s.log(2 * b)


def kernel(points, quanta, mass=s.S.One):
    if len(points) != 4:
        raise ValueError("Require the four generated massive scalar momenta")
    mu = s.sympify(mass)
    value = 2 * mu
    for i in range(4):
        for j in range(i + 1, 4):
            d = SIGNS[i] * SIGNS[j] * recoil.dot(points[i], points[j])
            value += 2 * SIGNS[i] * SIGNS[j] * massive_pair(s.factor(d), mu)
    for q in quanta:
        for i, p in enumerate(points):
            a = s.factor(SIGNS[i] * recoil.dot(p, q) / q[0])
            value += 2 * SIGNS[i] * q[0] * a * s.log(2 * a / s.sqrt(mu))
    for i, q in enumerate(quanta):
        for r in quanta[i + 1 :]:
            b = s.factor(recoil.dot(q, r) / (q[0] * r[0]))
            value += 2 * q[0] * r[0] * null_pair(b)
    return value


@cache
def data():
    checks = {}

    def put(name, value):
        checks[name] = s.factor(s.expand_func(s.expand_log(value)))

    d, c = s.symbols("d c", real=True)
    denom = 1 + c * (d - 1)
    put(
        "massive_pair_derivative_numerator",
        2 * d * denom
        - c * (d * d - s.Rational(1, 2))
        - (2 * d + c * (d * d - 2 * d + s.Rational(1, 2))),
    )
    y = s.Symbol("positive_pair_gap_root", positive=True)
    t = s.Symbol("angular_integration_coordinate", real=True)
    primitive = s.atanh(t * y / s.sqrt(2 + y * y)) / (y * s.sqrt(2 + y * y))
    put(
        "massive_Feynman_primitive",
        s.diff(primitive, t) - 1 / (2 + y * y - y * y * t * t),
    )
    put("massive_equal_velocity_value", massive_pair(1, 1) - s.Rational(1, 2))
    put(
        "massive_equal_velocity_integral",
        massive_pair_integral(1, 1).doit() - s.Rational(1, 2),
    )
    z = s.Symbol("large_positive_argument", positive=True)
    put("massless_arcosh_log_limit", s.limit(s.acosh(z) - s.log(2 * z), z, s.oo))
    r = s.Symbol("positive_fictitious_mass_squared", positive=True)
    for value in (s.Rational(1, 3), s.Integer(2)):
        raw = (
            (value * value - r / 2)
            * s.acosh(value / s.sqrt(r))
            / s.sqrt(value * value - r)
        )
        put(
            "massless_pair_finite_limit_" + str(value),
            s.limit(raw - value * s.log(2 * value / s.sqrt(r)), r, 0, dir="+"),
        )
    b12, b13, b23 = s.symbols(
        "radiation_pair_12 radiation_pair_13 radiation_pair_23", real=True
    )
    pair = {(0, 1): b12, (0, 2): b13, (1, 2): b23}
    hard = [-b12 - b13, -b12 - b23, -b13 - b23]
    for i in range(3):
        put(
            "each_massless_collinear_log_coefficient_" + str(i),
            hard[i] + sum(v for key, v in pair.items() if i in key),
        )
    put("common_massless_log_coefficient", sum(hard) + 2 * sum(pair.values()))
    ell = s.symbols("boost_log_energy_0:3", real=True)
    shift = sum(hard[i] * ell[i] for i in range(3)) + sum(
        v * (ell[i] + ell[j]) for (i, j), v in pair.items()
    )
    put("Lorentz_boost_log_energy_cancellation", shift)
    w1, w2, b, a, e = s.symbols("w1 w2 b a epsilon", positive=True)
    put(
        "massive_null_log_rewrite",
        s.expand_log(s.log(2 * w1 * a / (e * w1)), force=True)
        - s.log(2 * a)
        + s.log(e),
    )
    put(
        "null_null_log_rewrite",
        s.expand_log(s.log(2 * w1 * w2 * b / (e * e * w1 * w2)), force=True)
        - s.log(2 * b)
        + 2 * s.log(e),
    )
    put("null_pair_collinear_continuity", s.limit(b * s.log(2 * b), b, 0, dir="+"))
    put("null_pair_exact_collinear_value", null_pair(0))
    w, v, alpha, H = s.symbols("w v alpha kernel", real=True)
    put("collinear_split_linear_term", alpha * w * H + (1 - alpha) * w * H - w * H)
    put(
        "collinear_split_all_other_pairs",
        2 * alpha * w * v * H + 2 * (1 - alpha) * w * v * H - 2 * w * v * H,
    )
    put("collinear_split_self_pair", 2 * alpha * (1 - alpha) * w * w * null_pair(0))
    put("diagonal_massive_contribution", 4 * s.Rational(1, 2) - 2)
    put("pair_double_integral_count", 2 * s.Rational(1, 2) - 1)
    put("N0_original_kernel_domain", s.Integer(len(SIGNS)) - 4)
    # Positivity is a full physical TT statement, not positivity of each
    # signed pair in this collinear-finite representation.
    E0 = s.Rational(5, 4)
    u = s.Matrix([0, s.Rational(4, 5), s.Rational(3, 5)])
    points, _quanta, born = recoil.momenta(E0, [], u)
    put(
        "zero_radiation_exact_Born_state",
        s.Matrix.hstack(*points) - s.Matrix.hstack(*born),
    )
    literal = 2 + 2 * sum(
        SIGNS[i]
        * SIGNS[j]
        * massive_pair_integral(
            s.factor(SIGNS[i] * SIGNS[j] * recoil.dot(points[i], points[j])), 1
        )
        for i in range(4)
        for j in range(i + 1, 4)
    )
    put(
        "S296_complete_elastic_index_normalization",
        (literal - source.previous.soft.kernel_zero(E0, u[2], 1))
        / (4 * s.pi**2 * source.K),
    )
    return {
        "whole_massive_pair": massive_pair(D),
        "whole_massive_pair_integral": massive_pair_integral(D),
        "whole_null_pair_kernel": null_pair(),
        "whole_angular_energy_measure_formula": "For sigma=sum omega_r delta(n_r),R=sigma(S2),Q=(R,integral n dsigma),K=2mu+2sum_i<j eta_i eta_j F_mu(pi.pj)+2sum_i eta_i integral a_i(n)ln[2a_i(n)/sqrt(mu)]dsigma(n)+integral integral b(n,n')ln[2b(n,n')]dsigma(n)dsigma(n'),a_i=pi0-pivec.n,b=1-n.n'. Allmassive pi are future directed;0ln0=0. Atomic self-pairs vanish and the double integral counts each distinct unordered pair twice.",
        "whole_massless_and_boost_proof": "Start from the massive pair kernel and take auxiliary masses m_r=epsilon*omega_r. Its leading massless pair is d ln[2d/(m_i*m_r)]. Each ln m_r coefficient is qr.(sum_i eta_i*pi+sum_s qs)=0. The finite terms are exactly the angular-energy formula, including every radiation-radiation pair. A boost changes each reference omega_r inside the logs; its coefficient is the same zero, establishing frame independence. The four-scalar-only signed current is not conserved after recoil and cannot replace the full current.",
        "whole_continuum_boundary": "The mixed kernels are continuous by the uniform massive Doppler gap, and b ln(2b) extends continuously at b0 on the compact double sphere. The smooth recoil map and weak convergence of finite positive angular energy measures thus extend K continuously at fixed bounded total energy. Positivity follows from discrete physical TT sums and their limit. This is a kinematic coefficient extension, not a quantum measure or all-amplitude convergence theorem.",
        "checks": checks,
        "gates": {
            "all_massive_null_and_null_null_ordered_pairs": True,
            "each_collinear_mass_log_cancels_by_conservation": True,
            "full_conserved_TT_current_not_massive_subset": True,
            "boost_reference_energy_logs_cancel": True,
            "arbitrary_collinear_split_invariance": True,
            "finite_energy_measure_continuity_without_N_bound": True,
            "no_full_quantum_state_or_rate_inferred": True,
        },
    }
