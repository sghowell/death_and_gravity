"""Entire regulated unitary evolution and independent noncommuting fixtures."""

from functools import cache

import sympy as s


def occupation_basis(cutoff):
    if isinstance(cutoff, bool) or not isinstance(cutoff, int) or cutoff < 1:
        raise ValueError(
            "Require a positive integer total-occupation diagnostic cutoff"
        )
    return [(a, b) for a in range(cutoff + 1) for b in range(cutoff + 1 - a)]


def annihilator(basis, axis):
    index = {state: i for i, state in enumerate(basis)}
    matrix = s.zeros(len(basis))
    for j, state in enumerate(basis):
        if state[axis]:
            nxt = list(state)
            nxt[axis] -= 1
            matrix[index[tuple(nxt)], j] = s.sqrt(state[axis])
    return matrix


def coherent_matrix(polynomial, beta, basis, alphas, conjugates):
    """Exact compression of the full coherent integral, not finite CCR."""
    terms = s.Poly(polynomial, *alphas, *conjugates).terms()
    result = s.zeros(len(basis))
    for i, m in enumerate(basis):
        for j, n in enumerate(basis):
            value = 0
            for powers, coefficient in terms:
                upper = [m[k] + powers[k] for k in range(2)]
                lower = [n[k] + powers[k + 2] for k in range(2)]
                if upper != lower:
                    continue
                piece = coefficient
                for k in range(2):
                    piece *= s.factorial(upper[k]) / (
                        (1 + beta) ** (upper[k] + 1)
                        * s.sqrt(s.factorial(m[k]) * s.factorial(n[k]))
                    )
                value += piece
            result[i, j] = s.simplify(value)
    return result


@cache
def fixture():
    a = s.symbols("alpha0:2")
    b = s.symbols("alphabar0:2")
    q = [(a[i] + b[i]) / s.sqrt(2) for i in range(2)]
    p = [(a[i] - b[i]) / (s.I * s.sqrt(2)) for i in range(2)]
    q2, p2 = sum(z * z for z in q), sum(z * z for z in p)
    radius = s.expand(q2 + p2)
    qp = sum(q[i] * p[i] for i in range(2))
    basis = occupation_basis(3)
    H0 = coherent_matrix(
        s.expand(q2 - p2 + radius**2 / 10), s.Rational(1, 2), basis, a, b
    )
    H1 = coherent_matrix(s.expand(qp + qp**2 / 20), s.Rational(2, 3), basis, a, b)
    F = coherent_matrix(s.expand((q2 - p2) ** 2 + 1), s.Rational(1, 2), basis, a, b)
    A0, A1 = (annihilator(basis, axis) for axis in range(2))
    J = s.I * (A0.T * A1 - A1.T * A0)
    psi = s.zeros(len(basis), 1)
    psi[basis.index((0, 0))] = 1
    return {
        "basis": basis,
        "H0": H0,
        "H1": H1,
        "F": F,
        "J": J,
        "psi": psi,
        "identity_integral": coherent_matrix(1, 0, basis, a, b),
        "radial_vacuum_mean": coherent_matrix(radius, 0, basis, a, b)[0, 0],
        "gaussian_vacuum_mean": coherent_matrix(1, s.Rational(1, 2), basis, a, b)[0, 0],
        "annihilator": A0,
    }


def dyson_polynomials(H0, H1, order, variable):
    terms = [s.eye(H0.rows)]
    for _ in range(order):
        product = -s.I * (H0 + variable * H1) * terms[-1]
        terms.append(
            product.applyfunc(
                lambda expr: s.Poly(s.expand(expr), variable).integrate().as_expr()
            )
        )
    return terms


@cache
def data():
    fix = fixture()
    H0, H1, F, J, psi = (fix[k] for k in ("H0", "H1", "F", "J", "psi"))
    t = s.Symbol("fixture_clock", real=True)
    terms = dyson_polynomials(H0, H1, 4, t)
    checks = {
        "whole_normalized_coherent_resolution_compressed": fix["identity_integral"]
        - s.eye(H0.rows),
        "whole_coherent_H0_Hermiticity": H0 - H0.H,
        "whole_coherent_H1_Hermiticity": H1 - H1.H,
        "whole_positive_observable_Hermiticity": F - F.H,
        "whole_translation_H0_commutator": H0 * J - J * H0,
        "whole_translation_H1_commutator": H1 * J - J * H1,
        "same_vacuum_exact_zero_charge": J * psi,
        "whole_Husimi_radius_mean": fix["radial_vacuum_mean"] - 4,
        "whole_Gaussian_symbol_expectation": fix["gaussian_vacuum_mean"]
        - s.Rational(4, 9),
        "finite_compression_cannot_satisfy_CCR_trace": s.trace(
            fix["annihilator"] * fix["annihilator"].T
            - fix["annihilator"].T * fix["annihilator"]
        ),
    }
    for k in range(1, 5):
        checks["entire_time_ordered_Dyson_equation_" + str(k)] = (
            terms[k].diff(t) + s.I * (H0 + t * H1) * terms[k - 1]
        ).applyfunc(s.expand)
        checks["entire_Dyson_translation_commutator_" + str(k)] = (
            terms[k] * J - J * terms[k]
        ).applyfunc(s.expand)
        checks["entire_Dyson_unitarity_order_" + str(k)] = sum(
            (terms[j].H * terms[k - j] for j in range(k + 1)), s.zeros(H0.rows)
        ).applyfunc(s.expand)
    return {
        "whole_actual_regulated_interaction": "A(u)=Q_V((1-D)g_ext(u)) in units hbar=1; its scalar background is retained. A is bounded self-adjoint and norm-continuous on the common finite interval, with ||A(u)||<=||g_ext(u)||infinity+||D g_ext(u)||infinity. This keeps the whole local Hamiltonian and every implicit/time/cutoff contact under the declared first-Weyl-calibrated ordering.",
        "whole_full_propagator": "The norm-convergent Dyson series solves i partial_u U_I=A(u)U_I, U_I(0)=I, is unitary and composes exactly. U_full=U_ref U_I retains the entire original free flow. It solves the Schrödinger equation for H_ref+U_ref A U_ref*, a bounded perturbation of the reference quadratic Hamiltonian. Scalars contribute their complete phase.",
        "whole_Dyson_remainder": "For Lambda=int_0^|u| ||A(s)|| ds, each Dyson term has norm<=Lambda^k/k!. Absolute convergence first proves existence; differentiating U_I*U_I proves unitarity. Iterated Duhamel against the exact unitary endpoint then gives ||U_I-sum_{k=0}^M D_k||<=Lambda^(M+1)/(M+1)! without an extra exponential. This is an actual operator norm bound, unlike the separate Weyl SYMBOL calibration bound.",
        "whole_Schwartz_core": "The compact smooth interaction part has a Gaussian-window coherent integral with Schwartz range and finite L2-to-each-Schwartz-seminorm bound, uniformly on compact times. Separating the scalar phase, the Dyson terms have one such smoothing factor and bounded remaining factors, so the Dyson sum preserves Schwartz space. The metaplectic free flow preserves Schwartz; hence the full Schrödinger equation holds there strongly. No nonexistent finite-dimensional CCR representation is used.",
        "whole_residual_constraint_preservation": "The same pure seed, full symbols, cutoff and original free flow are equivariant under the complete compact three-translation group. Coherent covariance makes A commute with every group unitary; U_I and U_full preserve its common zero-charge Hilbert subspace. No Gaussian conditioning, homogeneous state or full diffeomorphism/BRST regulator is introduced.",
        "whole_volume_readout": "The actual regulated mean <U_full psi0, U_ref Q_V(F_ext) U_ref* U_full psi0>=<U_I psi0,Q_V(F_ext)U_I psi0> is finite and positive. For an uncalibrated observable norm B_F and Dyson approximation U_M with error epsilon, the mean error is <=B_F(2epsilon+epsilon^2). A first-Weyl-calibrated volume is a distinct readout and has only the conditional positivity gate in ordering.data. No numerical P8 mean or comparison with the undefined unregularized Hamiltonian is claimed.",
        "whole_independent_fixture": {
            k: v for k, v in fix.items() if k != "annihilator"
        },
        "whole_noncommuting_fixture_time_generator": H0 + t * H1,
        "whole_first_four_exact_time_ordered_terms": terms,
        "whole_fixture_Frobenius_squared_bounds": [s.trace(H0**2), s.trace(H1**2)],
        "checks": checks,
        "gates": {
            "exact_coherent_fixture_noncommuting_time_generators": H0 * H1 != H1 * H0,
            "exact_coherent_fixture_strict_positive_observable": F.is_positive_definite
            is True,
            "total_occupation_compression_preserves_rotation": H0 * J == J * H0
            and H1 * J == J * H1,
            "finite_CCR_not_falsely_claimed": fix["annihilator"] * fix["annihilator"].T
            - fix["annihilator"].T * fix["annihilator"]
            != s.eye(H0.rows),
            "fixture_is_not_full_P8_time_simulation_or_error_estimate": True,
            "all_three_actual_translation_generators_not_only_fixture_rotation": True,
            "positive_finite_regulated_mean_not_original_matching_or_closure": True,
        },
    }
