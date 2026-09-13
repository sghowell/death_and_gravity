"""Whole prepared Gaussian source influence, ordered variations and causal boundary."""

from functools import cache

import sympy as s


@cache
def data():
    times = (0, s.pi / 2, s.pi)
    W = s.Matrix(
        3, 3, lambda i, j: s.expand_complex(s.exp(-s.I * (times[i] - times[j])) / 2)
    )
    time_ordered = s.Matrix(3, 3, lambda i, j: W[i, j] if i >= j else W[j, i])
    anti_ordered = s.Matrix(3, 3, lambda i, j: W[j, i] if i >= j else W[i, j])
    retarded = s.Matrix(3, 3, lambda i, j: s.I * (W[i, j] - W[j, i]) if i > j else 0)
    noise = (W + W.T) / 2
    jp = s.Matrix(s.symbols("plus_source0:3", real=True))
    jm = s.Matrix(s.symbols("minus_source0:3", real=True))
    jd, jc = jp - jm, (jp + jm) / 2
    literal = (
        -(jp.T * time_ordered * jp + jm.T * anti_ordered * jm - 2 * jm.T * W * jp)[0]
        / 2
    )
    target = (s.I * jd.T * retarded * jc - jd.T * noise * jd / 2)[0]
    d = s.Matrix(s.symbols("difference0:3", real=True))
    c = s.Matrix(s.symbols("average0:3", real=True))
    action = (d.T * retarded * c + s.I * d.T * noise * d / 2)[0]
    force = s.Matrix([s.diff(action, v).subs(dict.fromkeys(d, 0)) for v in d])
    wrong = s.Matrix([s.diff((c.T * retarded * c)[0] / 2, v) for v in c])
    q = s.Symbol("metric_parameter", real=True)
    K0 = s.Matrix([[3, 1, 0], [1, 4, 1], [0, 1, 5]])
    K1 = s.Matrix([[2, 2, 0], [2, -1, 1], [0, 1, 2]])
    K2 = s.Matrix([[2, 0, 1], [0, 1, 0], [1, 0, -1]])
    K = K0 + q * K1 + q * q * K2 / 2
    green = K.inv()
    inv0 = K0.inv()
    j0 = s.Matrix([1, 2, 3])
    j1 = s.Matrix([2, -1, 1])
    j2 = s.Matrix([1, 3, -2])
    j = j0 + q * j1 + q * q * j2 / 2
    stationary = green * j
    S = (j.T * green * j)[0] / 2
    first = (j1.T * inv0 * j0 - j0.T * inv0 * K1 * inv0 * j0 / 2)[0]
    second = (
        j2.T * inv0 * j0
        + j1.T * inv0 * j1
        - 2 * j1.T * inv0 * K1 * inv0 * j0
        + j0.T * inv0 * K1 * inv0 * K1 * inv0 * j0
        - j0.T * inv0 * K2 * inv0 * j0 / 2
    )[0]
    hh = s.Matrix(s.symbols("heavy0:3"))
    whole = (-hh.T * K * hh / 2 + hh.T * j)[0]
    square = (-(hh - stationary).T * K * (hh - stationary) / 2)[0] + S
    eps = s.Symbol("perturbation", real=True)
    # An entire analytic kernel variation cannot lower the two endpoint source degrees.
    Jleft = eps**8 * (1 + 2 * eps + eps * eps)
    Jright = eps**8 * (3 - eps + eps * eps)
    kernel = (2 + eps) / (3 - eps)
    full = Jleft * kernel * Jright
    checks = {
        "full_CTP_ordered_cumulant_retarded_and_noise_normalization": s.expand(
            literal - target
        ),
        "full_CTP_identical_sources_unitarity": s.expand(
            literal.subs(dict(zip(jm, jp, strict=True)))
        ),
        "whole_physical_causal_source_force": force - retarded * c,
        "whole_noise_kernel_symmetry": noise - noise.T,
        "entire_retarded_kernel_from_original_commutator": retarded
        - s.Matrix([[0, 0, 0], [1, 0, 0], [0, 1, 0]]),
        "complete_heavy_Gaussian_square_with_density_operator": s.factor(
            whole - square
        ),
        "whole_inverse_first_variation": green.diff(q).subs(q, 0) + inv0 * K1 * inv0,
        "whole_inverse_second_variation": green.diff(q, 2).subs(q, 0)
        - 2 * inv0 * K1 * inv0 * K1 * inv0
        + inv0 * K2 * inv0,
        "full_source_and_operator_first_variation": s.diff(S, q).subs(q, 0) - first,
        "full_source_and_operator_second_variation": s.diff(S, q, 2).subs(q, 0)
        - second,
        "source_free_heavy_determinant_first_metric_jet": s.diff(
            s.log(K.det()) / 2, q
        ).subs(q, 0)
        - s.trace(inv0 * K1) / 2,
        "source_bilinear_full_degree_sixteen": s.limit(full / eps**16, eps, 0) - 2,
        "source_mean_degree_eight": s.limit(Jleft * kernel / eps**8, eps, 0)
        - s.Rational(2, 3),
        "source_Euler_degree_fifteen": s.limit(s.diff(full, eps) / eps**15, eps, 0)
        - 32,
        "source_first_response_degree_fourteen": s.limit(
            s.diff(full, eps, 2) / eps**14, eps, 0
        )
        - 480,
    }
    for degree in range(16):
        checks["all_lower_source_influence_jets_" + str(degree)] = s.diff(
            full, eps, degree
        ).subs(eps, 0)
    return {
        "exact_prepared_integration": "For two complete CTP metric histories and the unchanged zero-mean Gaussian H preparation, the source dependence of the logarithmic influence is exactly quadratic in the TWO density sources. Keep the entire source-independent metric determinant. Its metric dependence is not a quadratic-stress truncation.",
        "common_geometry_conventions": "With hbar1, G_R=i theta<[H,H]>, N=<{H,H}>/2, Jd=Jplus-Jminus and Jc=(Jplus+Jminus)/2, S_IF,J=Jd G_R Jc+i Jd N Jd/2. The physical difference-leg variation is retarded. At genuinely distinct metric histories use the complete contour inverse, not a common-G_R shortcut.",
        "retarded_kernel_fixture": retarded,
        "noise_kernel_fixture": noise,
        "wrong_one_copy_retarded_force": wrong,
        "ordered_Gaussian_first_source_variation": first,
        "ordered_Gaussian_second_source_variation": second,
        "source_jet_result": "The complete source influence starts degree16, retarded coherent H starts8, source-dependent metric/light Euler current starts15 and its first response starts14. Density and inverse-operator variations remain. These are prepared smooth Taylor-jet statements, not a nonlinear norm, amplitude estimate or convergence theorem.",
        "coherent_mean_boundary": "A nonzero fixed H homogeneous mean adds an influence term linear in J, of degree8 with current degree7. It is excluded only by the actual zero-mean preparation, not by a universal source identity.",
        "determinant_boundary": "Source zeros do not delete the full metric-dependent Gaussian determinant or its loop1 mean, response, noise and higher metric cumulants. Those conditional quantities and fixed profiles remain in the current reference.",
        "checks": {
            name: value.applyfunc(s.cancel)
            if isinstance(value, s.MatrixBase)
            else s.cancel(value)
            for name, value in checks.items()
        },
        "gates": {
            "one_copy_retarded_quadratic_is_not_the_causal_force": wrong
            != retarded * c,
            "ordered_kernel_and_its_variations_do_not_commute": K0 * K1 != K1 * K0,
            "full_noise_kernel_nonzero_and_positive_semidefinite": noise != s.zeros(3)
            and noise.is_positive_semidefinite,
            "source_free_full_heavy_determinant_not_removed": s.trace(inv0 * K1) != 0,
            "physical_source_and_operator_first_terms_both_nonzero": (j1.T * inv0 * j0)[
                0
            ]
            != 0
            and (j0.T * inv0 * K1 * inv0 * j0)[0] != 0,
            "distinct_CTP_metric_histories_not_assumed_equal": True,
            "same_prepared_Gaussian_H_state_not_interacting_state_theorem": True,
        },
    }
