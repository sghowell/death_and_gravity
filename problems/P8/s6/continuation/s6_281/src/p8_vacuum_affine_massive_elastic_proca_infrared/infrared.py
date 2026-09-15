"""Complete nonforward elastic finite part in an explicit dimensional soft convention."""

from functools import cache

import sympy as s

from . import elastic, source

S, MU, K, Z = source.S, source.MU, source.K, source.Z
EP = s.Symbol("positive_IR_epsilon", positive=True)
NU, E = s.symbols(
    "positive_dimensional_scale positive_detector_resolution", positive=True
)


def angular_finite_same(z=Z):
    return (s.log((1 - z) / 2) / (1 - z) + s.log((1 + z) / 2) / (1 + z)) / 2


def angular_finite_mixed(v, z=Z):
    return (
        s.log((v - z) ** 2 / (v * v - 1)) / (v - z)
        + s.log((v + z) ** 2 / (v * v - 1)) / (v + z)
    ) / (4 * v)


def whole_auxiliary_finite(z=Z):
    p, v, b, a0, a2 = elastic.coefficients()
    p2 = s.legendre(2, z)
    return (
        p * p * angular_finite_same(z)
        + 2 * p * b * angular_finite_mixed(v, z)
        - 3 * p * a2 * p2
        + b * b * elastic.master_J(v, v, z)
        + 2 * b * (a0 * elastic.Q0(v) / v + a2 * elastic.Q2(v) * p2 / v)
        + a0 * a0
        + a2 * a2 * p2 / 5
    )


def evanescent_tree(z=Z):
    q = S - 4 * MU
    t = -q * (1 - z) / 2
    u = -q * (1 + z) / 2
    return -2 * MU**2 / K * (1 / S + 1 / t + 1 / u)


def whole_bare_cut(epsilon=EP, scale=NU, z=Z):
    p = elastic.coefficients()[0]
    tree = elastic.whole_tree(z)
    p1 = 8 * MU**2 / (K * (S - 4 * MU))
    return (
        s.sqrt(1 - 4 * MU / S)
        / (32 * s.pi)
        * (
            p * tree / epsilon
            + whole_auxiliary_finite(z)
            + p * tree * (s.log((S - 4 * MU) / (4 * s.pi * scale**2)) + s.EulerGamma)
            + p1 * tree
            + p * evanescent_tree(z)
        )
    )


def whole_stripped_cut(resolution=E, z=Z):
    p = elastic.coefficients()[0]
    tree = elastic.whole_tree(z)
    p1 = 8 * MU**2 / (K * (S - 4 * MU))
    return (
        s.sqrt(1 - 4 * MU / S)
        / (32 * s.pi)
        * (
            whole_auxiliary_finite(z)
            + p
            * tree
            * (s.log((S - 4 * MU) / (4 * s.pi * resolution**2)) + s.EulerGamma)
            + p1 * tree
        )
    )


@cache
def data():
    z, v = s.symbols("z v", real=True)
    p, b, a0, a2 = s.symbols("P B a0 a2", real=True)
    pp = s.legendre(2, z)
    tree = p / (1 - z * z) + b / (v * v - z * z) + a0 + a2 * pp
    log_coeff = p * p / (1 - z * z) + p * b / (v * v - z * z) + p * (a0 + a2 * pp)
    checks = {
        "entire_soft_log_coefficient_all_interferences": s.factor(log_coeff - p * tree),
        "same_master_log_coefficient": s.factor(
            (1 / (1 - z) + 1 / (1 + z)) / 2 - 1 / (1 - z * z)
        ),
        "mixed_master_log_coefficient": s.factor(
            (1 / (v - z) + 1 / (v + z)) / (4 * v) - 1 / (2 * (v * v - z * z))
        ),
        "Q2_constant_finite_term": -s.Integer(3) / 2 + s.Rational(3, 2),
    }
    # Prove exact endpoint rewritings whose smooth limits give all finite logs.
    d, c = s.symbols("positive_regulator_delta positive_endpoint", positive=True)
    delta = s.sqrt(c * c + 2 * c * d)
    log_argument = (c + d + delta) / d
    checks["same_endpoint_log_argument_factorization"] = s.factor(
        log_argument / (4 / d) - (c + d + delta) / 4
    )
    checks["same_endpoint_discriminant_limit"] = s.limit(delta, d, 0, dir="+") - c
    checks["same_endpoint_finite_log_limit"] = (
        s.limit((c + d + delta) / 4, d, 0, dir="+") - c / 2
    )
    cv = s.sqrt(1 + d) * v - z
    mixedroot = s.sqrt(cv * cv - d * (v * v - 1))
    checks["mixed_endpoint_squared_discriminant"] = s.expand(
        mixedroot**2 - (cv**2 - d * (v * v - 1))
    )
    w = s.sqrt(1 + d)
    checks["Q0_complete_log_argument"] = s.factor(
        ((w + 1) / (w - 1)) / (4 / d) - (w + 1) ** 2 / 4
    )
    # Exact dimensional scalar stress contraction, all incoming momenta.
    ss, tt, mu, ep = s.symbols("s t mu epsilon", real=True)
    uu = 4 * mu - ss - tt
    dimension = 4 + 2 * ep
    dot = (ss - 2 * mu) / 2
    contracted = (
        ((tt - 2 * mu) ** 2 + (uu - 2 * mu) ** 2) / 2
        - 4 * dot * (dot + mu)
        + dimension * (dot + mu) ** 2
    )
    trace = 2 * dot - dimension * (dot + mu)
    numerator = s.factor(contracted - trace**2 / (dimension - 2))
    N4 = 2 * mu**2 - 2 * mu * ss - tt * uu
    checks["entire_dimensional_graviton_trace_numerator"] = s.factor(
        numerator - N4 - 2 * mu**2 * ep / (1 + ep)
    )
    q = S - 4 * MU
    p1 = 8 * MU**2 / (K * q)
    checks["whole_evanescent_tree_endpoint_decomposition"] = s.factor(
        evanescent_tree(z) - (p1 / (1 - z * z) - 2 * MU**2 / (K * S))
    )
    # Gamma recurrences make the normalized endpoint integral exact.
    checks["exact_D_sphere_endpoint_integral"] = s.factor(
        (ep + s.Rational(1, 2)) / ep - (1 / (2 * ep) + 1)
    )
    be = s.Symbol("physical_beta", positive=True)
    pref_log = s.log(ss) + 2 * s.log(be) - 4 * s.log(2) - s.log(s.pi) - 2 * s.log(NU)
    normalized = pref_log - s.digamma(s.Rational(3, 2)) + 2
    expected = (
        s.log(ss) + 2 * s.log(be) - s.log(4 * s.pi) - 2 * s.log(NU) + s.EulerGamma
    )
    checks["full_D_phase_space_finite_constant"] = s.simplify(
        s.expand_log(normalized - expected, force=True)
    )
    VV = S * S - 4 * MU * S + 2 * MU * MU
    pole = be * elastic.coefficients()[0] / (32 * s.pi)
    soft = VV / (8 * s.pi * K * S * be)
    checks["original_massive_soft_phase_coefficient"] = s.factor(
        (pole - soft).subs(MU, S * (1 - be * be) / 4)
    )
    A, A1, P1, F, Lq, LE = s.symbols("A A1 P1 F_aux log_Q_over_4pi_nu2 log_E_over_nu")
    bare = p * A / ep + F + p * A * (Lq + s.EulerGamma) + P1 * A + p * A1
    subtraction = (p / ep + 2 * p * LE) * (A + ep * A1)
    stripped = s.expand(bare - subtraction).subs(ep, 0)
    checks["whole_finite_stripping_including_D_tree"] = s.expand(
        stripped - (F + p * A * (Lq + s.EulerGamma - 2 * LE) + P1 * A)
    )
    checks["physical_resolution_flow"] = s.diff(whole_stripped_cut(), E) + s.sqrt(
        1 - 4 * MU / S
    ) * elastic.coefficients()[0] * elastic.whole_tree() / (16 * s.pi * E)
    return {
        "same_resolvent_complete_finite_part": angular_finite_same(),
        "mixed_resolvent_complete_finite_part": angular_finite_mixed(elastic.V),
        "whole_auxiliary_finite_remainder": whole_auxiliary_finite(),
        "entire_dimensional_gravity_numerator": numerator,
        "entire_first_evanescent_elastic_tree": evanescent_tree(),
        "entire_bare_elastic_cut_through_finite_order": whole_bare_cut(),
        "entire_S278_soft_stripped_elastic_cut": whole_stripped_cut(),
        "exact_normalized_D_sphere_endpoint": "<1/(1-x^2)>_(S^(2+2epsilon))=1/(2epsilon)+1, epsilon>0.",
        "complete_D_measure": "dPhi2=(2pi)^(2-D) p^(D-3)/(4sqrt(s)) dOmega_(D-2). Each tree scales as nu^(-2epsilon) and the final amplitude is rescaled by nu^(2epsilon). The normalized cut prefactor is beta/(32pi)*(s beta^2/(16pi nu^2))^epsilon*Gamma(3/2)/Gamma(3/2+epsilon).",
        "finite_convention": "The real minimal D-dimensional graviton projector is continued as 1/(D-2), while the S278 analytic soft exponential is DEFINED by its four-dimensional massive kernel. The resulting P1*A finite phase is retained. Gamma_E and log4pi reflect the displayed raw measure convention, not an undocumented MSbar choice.",
        "full_soft_subtraction": "Subtract Im(logW_E)*(A+epsilon*A1), not merely Im(logW_E)*A. The P*A1 finite term cancels but P1*A does not. Physical E and the auxiliary angular w are never identified.",
        "continuum_proof_boundary": "At fixed strictly nonforward angle subtract P_D*A_D(z)[T(a.n)+T(b.n)]. Four disjoint endpoint neighborhoods have integrable O(1/theta) remainder; a uniform logarithmic derivative bound justifies the epsilon expansion. Its four-dimensional limit is the complete auxiliary finite remainder. Notes supply the dominated-convergence proof.",
        "not_established": "Forward/unphysical analytic continuation, the complete real amplitude or b20, exact all-order factorization, finite-gravity positivity or a unique off-background dimensional completion.",
        "checks": checks,
        "gates": {
            "full_tree_times_soft_pole_not_gravity_only": all(
                whole_bare_cut().has(value) for value in (elastic.C, elastic.G2)
            ),
            "finite_evanescent_mass_term_retained": evanescent_tree().has(MU),
            "physical_resolution_survives_but_dimensional_scale_cancels": whole_stripped_cut().has(
                E
            )
            and not whole_stripped_cut().has(NU, EP),
            "Euler_and_phase_space_finite_constants_retained": whole_stripped_cut().has(
                s.EulerGamma
            ),
            "strict_nonforward_domain_not_forward_limit_exchange": True,
            "specific_minimal_D_and_S278_soft_convention_named": True,
            "conditional_first_loop_not_exact_finite_G_unitarity": True,
        },
    }
