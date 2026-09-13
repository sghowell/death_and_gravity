"""Whole fixed-state Weyl-quadratic mean, response and noise at a finite regulator."""

from functools import cache

import sympy as s
from p8_vacuum_affine_coupled_gaussian_state import gaussian


def omega(dimension):
    if type(dimension) is not int or dimension not in (2, 4, 8):
        raise ValueError(
            "The diagnostic interface supports complete2,4,8-dimensional phases"
        )
    half = dimension // 2
    out = s.zeros(dimension)
    out[:half, half:] = s.eye(half)
    out[half:, :half] = -s.eye(half)
    return out


def exact_matrix(value, *, dimension=None, real=True, symmetric=False):
    M = s.Matrix(value)
    if M.rows != M.cols or M.rows not in (2, 4, 8):
        raise ValueError("Require a complete exact finite even phase matrix")
    if dimension is not None and M.shape != (dimension, dimension):
        raise ValueError("Inconsistent full phase dimensions")
    if any(x.has(s.Float) or x.free_symbols for x in M):
        raise ValueError("The finite diagnostic interface requires fixed exact entries")
    if real and any(x.is_real is not True for x in M):
        raise ValueError("Require provably real entries")
    if symmetric and M != M.T:
        raise ValueError("Require a symmetric Weyl-quadratic coefficient")
    return M


def pure_covariance(value):
    V = exact_matrix(value, symmetric=True)
    O = omega(V.rows)
    if not all(
        s.factor(V[:j, :j].det()).is_positive is True for j in range(1, V.rows + 1)
    ):
        raise ValueError("Require strictly positive Cauchy covariance")
    if (V * O * V - O / 4).applyfunc(s.cancel) != s.zeros(V.rows):
        raise ValueError("Require the full exact pure-state CCR covariance")
    return V


def symplectic(value, dimension):
    S = exact_matrix(value, dimension=dimension)
    O = omega(dimension)
    if (S * O * S.T - O).applyfunc(s.cancel) != s.zeros(dimension):
        raise ValueError("Require the complete canonical symplectic propagator")
    return S


def kernels(A, B, left, right, initial, contact=None):
    """No time step, cutoff removal, or physical-background map is implicit."""
    V = pure_covariance(initial)
    dimension = V.rows
    A = exact_matrix(A, dimension=dimension, symmetric=True)
    B = exact_matrix(B, dimension=dimension, symmetric=True)
    L = symplectic(left, dimension)
    R = symplectic(right, dimension)
    C = (
        s.zeros(dimension)
        if contact is None
        else exact_matrix(contact, dimension=dimension, symmetric=True)
    )
    O = omega(dimension)
    W = L * (V + s.I * O / 2) * R.T
    VA, VB = L * V * L.T, R * V * R.T
    pulled_A, pulled_B = L.T * A * L, R.T * B * R
    connected = s.expand(s.trace(A * W * B * W.T) / 2)
    susceptibility = (
        s.trace((pulled_A * O * pulled_B - pulled_B * O * pulled_A) * V) / 2
    )
    return {
        "mean_A": s.cancel(s.trace(A * VA) / 2),
        "mean_B": s.cancel(s.trace(B * VB) / 2),
        "full_unequal_time_two_point": W,
        "ordered_connected_Wick": s.simplify(connected),
        "symmetric_noise": s.simplify(s.re(connected)),
        "observable_retarded_susceptibility_before_step": s.cancel(susceptibility),
        "effective_action_retarded_part_before_step": s.cancel(-susceptibility),
        "effective_action_seagull": s.cancel(-s.trace(C * VA) / 2),
    }


def coefficient_fixture(case, dimension=4):
    if type(case) is not int or not 0 <= case < 12:
        raise ValueError("Unknown exact full quadratic coefficient fixture")
    omega(dimension)
    m = s.Integer(case + 1)
    return s.Matrix(
        dimension,
        dimension,
        lambda i, j: (
            (-1) ** (i + j) * (m + i + j + 1) / s.Integer(7 + i * j + i + j)
            + (m + i if i == j else 0)
        ),
    )


@cache
def data():
    checks = {}
    O = omega(4)
    for case in range(4):
        prepare = gaussian.symplectic_fixture(case)
        V = prepare * prepare.T / 2
        left = gaussian.symplectic_fixture(case + 2)
        right = gaussian.symplectic_fixture(case + 5)
        A, B, C = [coefficient_fixture(case + j) for j in range(3)]
        result = kernels(A, B, left, right, V, C)
        W = result["full_unequal_time_two_point"]
        left_variation = left * right.inv() * O * B * right
        V_variation = left_variation * V * left.T + left * V * left_variation.T
        independent = s.trace(A * V_variation) / 2
        reverse = kernels(B, A, right, left, V)
        checks[f"entire_Wick_commutator_matches_quadratic_Lie_algebra_{case}"] = (
            2 * s.im(result["ordered_connected_Wick"])
            - result["observable_retarded_susceptibility_before_step"]
        )
        checks[f"whole_Duhamel_kick_response_matches_commutator_{case}"] = (
            independent - result["observable_retarded_susceptibility_before_step"]
        )
        checks[f"ordered_connected_reverse_is_conjugate_{case}"] = reverse[
            "ordered_connected_Wick"
        ] - s.conjugate(result["ordered_connected_Wick"])
        checks[f"noise_symmetric_under_full_time_vertex_swap_{case}"] = (
            reverse["symmetric_noise"] - result["symmetric_noise"]
        )
        checks[f"full_unequal_time_CCR_{case}"] = (
            W - W.conjugate() - s.I * left * O * right.T
        )
        checks[f"full_local_second_vertex_retained_{case}"] = (
            result["effective_action_seagull"] + s.trace(C * left * V * left.T) / 2
        )
        U = gaussian.symplectic_fixture(case + 7)
        transformed = kernels(
            U.inv().T * A * U.inv(),
            U.inv().T * B * U.inv(),
            U * left * U.inv(),
            U * right * U.inv(),
            U * V * U.T,
            U.inv().T * C * U.inv(),
        )
        for name in (
            "mean_A",
            "ordered_connected_Wick",
            "symmetric_noise",
            "observable_retarded_susceptibility_before_step",
            "effective_action_seagull",
        ):
            checks[f"entire_common_canonical_basis_invariance_{name}_{case}"] = (
                transformed[name] - result[name]
            )
    vacuum = s.eye(4) / 2
    stationary = kernels(s.eye(4), s.eye(4), s.eye(4), s.eye(4), vacuum)
    W0 = stationary["full_unequal_time_two_point"]
    checks["vacuum_total_number_energy_has_zero_connected_variance"] = stationary[
        "symmetric_noise"
    ]
    checks["ground_Wick_pair_uses_transpose_not_adjoint"] = W0 * W0.T
    tensor_O = omega(2)
    for case in range(3):
        x = s.Rational(case + 1, 7)
        T = s.Matrix([[1, x], [x / 3, 1 + x * x / 3]])
        A, B = coefficient_fixture(case, 2), coefficient_fixture(case + 2, 2)
        result = kernels(A, B, T, s.eye(2), s.eye(2) / 2)
        checks[f"tensor_full_Wick_commutator_{case}"] = (
            2 * s.im(result["ordered_connected_Wick"])
            - result["observable_retarded_susceptibility_before_step"]
        )
        checks[f"tensor_full_symplectic_transport_{case}"] = (
            T * tensor_O * T.T - tensor_O
        )
    return {
        "whole_fixed_state_W": "W(t,s)=S(t,t*)[V0+iOmega/2]S(s,t*)^T, with the entire S251 scalar covariance or each complete tensor covariance; all initial data are fixed under perturbation.",
        "entire_Weyl_quadratic_mean": "<O_A(t)>=Tr[A(t)V(t,t)]/2, O_A=z^T A z/2 in Weyl order. The effective-action current for a Hamiltonian coefficient is its negative.",
        "ordered_connected_quadratic_Wick": "<O_A(t)O_B(s)>c=Tr[A W(t,s) B W(t,s)^T]/2. The last factor is an ordinary transpose, not a conjugate transpose.",
        "noise": "N_AB(t,s)=Re <O_A(t)O_B(s)>c. Positivity follows for arbitrary finite real test combinations from the squared norm of the centered Hermitian quadratic operator on the fixed Gaussian state.",
        "retarded_observable_response": "For t>s, chi_AB=-i<[O_A(t),O_B(s)]>=Tr[(A_* Omega B_*-B_* Omega A_*)V0]/2, A_*=S(t,t*)^T A S(t,t*), B_*=S(s,t*)^T B S(s,t*). Multiply by theta(t-s); no advanced response is included.",
        "effective_action_second_variation": "Gamma_,AB=-<O_HAB>delta(t-s)-theta(t-s)chi_HA,HB. The entire second Hamiltonian vertex H_AB is the seagull. Initial-state derivatives are zero in the fixed canonical representation; a time-dependent chart instead transports the same density and its boundary phase.",
        "finite_regulator_scope": "Exact Gaussian operator/closed-time-path identities for the whole physical quadratic reference and specified coefficient probes at a common finite regulator. They do not by themselves define full parent metric/light vertices, a covariant gauge measure, local finite counterterms, regulator removal or a renormalized physical loop mean.",
        "checks": {
            name: value.applyfunc(s.cancel)
            if isinstance(value, s.MatrixBase)
            else s.cancel(value)
            for name, value in checks.items()
        },
        "gates": {
            "anomalous_mixing_retained_in_fixed_reference": gaussian.symplectic_fixture(
                0
            )[:2, 2:]
            != s.zeros(2),
            "transpose_not_adjoint_negative_control": s.trace(W0 * W0.conjugate().T) / 2
            != 0,
            "causal_step_and_seagull_are_both_required": True,
            "state_not_reminimized_on_perturbed_history": True,
            "coefficient_current_not_automatically_covariant_metric_stress": True,
            "finite_cutoff_does_not_bound_or_renormalize_full_loop": True,
        },
    }
