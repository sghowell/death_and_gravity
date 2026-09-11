"""Exact flat operators, common-regulator Gaussian identity and source contact."""

from functools import cache

import sympy as s


@cache
def data():
    z = s.symbols("zeta", positive=True)
    p = s.Matrix(s.symbols("p0:4", real=True))
    eta, eye = s.diag(1, -1, -1, -1), s.eye(4)
    p2 = (p.T * eta * p)[0]
    outer = p * p.T * eta
    K, D = -p2 * eye + outer, -outer
    O, N, L = eye + z * K, eye + z * D, (1 - z * p2) * eye
    inverse = (eye - z * outer) / (1 - z * p2)
    pe2 = (p.T * p)[0]
    Ke = pe2 * eye - p * p.T
    Oe = eye + z * Ke
    inverse_e = (eye + z * p * p.T) / (1 + z * pe2)
    W, S = s.Matrix(s.symbols("W0:4")), s.Matrix(s.symbols("S0:4"))
    shift = inverse * S
    original = (
        (W.T * eta * O * W)[0] - 2 * (W.T * eta * S)[0] + (S.T * eta * S)[0]
    ) / 2
    centered = (
        ((W - shift).T * eta * O * (W - shift))[0]
        + (S.T * eta * (eye - inverse) * S)[0]
    ) / 2
    checks = {
        "K_gradient_zero": K * p,
        "divergence_K_zero": p.T * eta * K,
        "mixed_KD_zero": K * D,
        "mixed_DK_zero": D * K,
        "ON_normally_hyperbolic": O * N - L,
        "NO_normally_hyperbolic": N * O - L,
        "Lorentz_left_inverse": O * inverse - eye,
        "Lorentz_right_inverse": inverse * O - eye,
        "Lorentz_pairing_self_adjoint": eta * O - O.T * eta,
        "sourced_divergence_identity": p.T * eta * O - p.T * eta,
        "closed_source_identity": inverse * p - p,
        "retarded_response_algebra": inverse - eye + z * K / (1 - z * p2),
        "Euclidean_inverse": Oe * inverse_e - eye,
        "Euclidean_gradient_unchanged": inverse_e * p - p,
        "source_contact_complete_square": s.cancel(original - centered),
    }
    return {
        "Lorentz_operator": O,
        "Lorentz_inverse_off_pole": inverse,
        "Euclidean_operator": Oe,
        "Euclidean_inverse": inverse_e,
        "source_action_kernel": "kappa/2 <S,(I-O_inverse)S>; the local +S^2/2 contact is retained",
        "fixed_metric_vector_determinant": "At a specified common finite regulator the W translation has unit Jacobian. The vector quadratic determinant depends on g,zeta,kappa but not u; it is not a renormalized stress tensor or a bound on light/gravity loops.",
        "real_time_boundary": "A retarded inverse solves the sourced equation. A single-branch quadratic action requires its specified variational/in-out boundary inverse; varying S G_ret S symmetrizes G_ret and does not give a purely retarded equation.",
        "covariant_identity": "K=div d, D=grad div, KD=DK=0, O=I+zeta K, N=I+zeta D, L=ON=NO. L/zeta is normally hyperbolic on one-forms after raising/lowering by the physical metric.",
        "checks": {
            k: v.applyfunc(s.cancel) if isinstance(v, s.MatrixBase) else s.cancel(v)
            for k, v in checks.items()
        },
    }


@cache
def variation():
    a, b, c, d = s.symbols("a b c d", real=True)
    x = s.Matrix(s.symbols("s0:2", real=True))
    G = s.Matrix([[a, b], [c, d]])
    action = (x.T * G * x)[0] / 2
    derivative = s.Matrix([s.diff(action, v) for v in x])
    eps, source2, source4, kernel0 = s.symbols("epsilon source2 source4 kernel0")
    clock = s.expand((eps**2 * source2) ** 2 * kernel0 / 2)
    vacuum = s.expand((eps**4 * source4) ** 2 * kernel0 / 2)
    return {
        "clock_first_possible_source_functional_degree": 4,
        "vacuum_first_possible_source_functional_degree": 8,
        "single_eight_scalar_vertex_first_possible_light_loop_orders": {
            "four_point": 2,
            "two_point": 3,
            "vacuum": 4,
        },
        "retarded_single_branch_error": derivative - G * x,
        "checks": {
            "single_branch_kernel_is_symmetrized": (
                derivative - (G + G.T) * x / 2
            ).applyfunc(s.expand),
            "clock_degree_four": clock.coeff(eps, 4) - source2**2 * kernel0 / 2,
            "clock_degrees_below_four": sum(clock.coeff(eps, i) for i in range(4)),
            "vacuum_degree_eight": vacuum.coeff(eps, 8) - source4**2 * kernel0 / 2,
            "vacuum_degrees_below_eight": sum(vacuum.coeff(eps, i) for i in range(8)),
            "four_point_single_vertex_two_loops": (8 - 4) // 2 - 2,
            "two_point_single_vertex_three_loops": (8 - 2) // 2 - 3,
            "vacuum_single_vertex_four_loops": 8 // 2 - 4,
        },
    }
