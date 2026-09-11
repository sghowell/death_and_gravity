"""Exact two-momentum Gaussian response before the ultraviolet limit."""

from functools import cache

import sympy as s


def pair_source(H, Ck, Cq, J):
    return J * H * Cq + Ck * H * J.T


def pair_block(H):
    zero = s.zeros(H.rows)
    return zero.row_join(H).col_join(H.T.row_join(zero))


@cache
def data():
    J = s.Matrix([[0, 1], [-1, 0]])

    def symmetric(prefix):
        a, b, c = s.symbols(prefix + "0:3", real=True)
        return s.Matrix([[a, b], [b, c]])

    Ck, Cq = symmetric("Ck"), symmetric("Cq")
    H = s.Matrix(2, 2, s.symbols("H0:4", real=True))
    D = s.Matrix(2, 2, s.symbols("D0:4", real=True))
    Jall = s.diag(J, J)
    C = s.diag(Ck, Cq)
    M, N = pair_block(H), pair_block(D)
    full = Jall * M * C + C * M * Jall.T
    x, y, z, w = s.symbols("x y z w", real=True)
    Uk = s.Matrix([[1, x], [y, 1 + x * y]])
    Uq = s.Matrix([[1, z], [w, 1 + z * w]])
    U = s.diag(Uk, Uq)
    transported = U.T * N * U
    comm = transported * Jall * M - M * Jall * transported
    checks = {
        "full_covariance_off_diagonal_pair_source": full[:2, 2:]
        - pair_source(H, Ck, Cq, J),
        "reversed_pair_source_transpose": full[2:, :2] - pair_source(H, Ck, Cq, J).T,
        "full_two_mode_symplectic_transport": U * Jall * U.T - Jall,
        "unequal_propagators_full_Kubo_trace_identity": s.trace(N * U * full * U.T)
        - s.trace(comm * C),
        "pair_current_counts_both_reverse_blocks": s.trace(N * U * full * U.T)
        - 2 * s.trace(D.T * Uk * pair_source(H, Ck, Cq, J) * Uq.T),
    }
    return {
        "Fourier_covariance": "C(k,q)=<sym Z(k) Z(q)^dagger>, with unperturbed C(k,q)=delta(k-q)C_k. Complex amplitudes are the complexification of the real-field canonical algebra; reversed pairs obey the corresponding adjoint relation.",
        "exact_pair_tangent": "delta C_kq'=A_k delta C_kq+delta C_kq A_q^t+J M_Gamma(k,q)C_q+C_k M_Gamma(k,q)J^t, delta C_kq(t0)=0. The initial state is unchanged because the prescribed external perturbation vanishes in the common preparation neighborhood.",
        "retarded_solution": "Integrate U_k(t,s)[J M_Gamma(k,q;s)C_q(s)+C_k(s)M_Gamma(k,q;s)J^t]U_q(t,s)^t from t0 to t. The propagators and physical momenta are generally different.",
        "full_current": "delta J_D(t)=-1/2 integral tr[M_D(q,k;t)delta C_kq(t)]dkdq-1/2 integral tr[M_DGamma(k,k;t)C_k(t)]dk, with the indicated Fourier measures and external-momentum pairing. The second term is a local-in-time metric contact.",
        "Kubo": "On every common finite real-mode regulator the complete expression equals -<H_DGamma>+i integral <[H_D(t),H_Gamma(s)]>ds. The factor and sign follow from varying +H_Gamma and reading out -H_D, not from replacing the commutator by the noise anticommutator.",
        "initial_state_boundary": "External-metric vector response only. No factorized or vanishing quantum metric initial data or switched-on interaction is presumed.",
        "checks": checks,
        "gates": {
            "two_different_momenta_and_propagators": Uk != Uq,
            "contact_retained_separately_from_covariance": True,
            "full_real_regulator_before_continuum": True,
        },
    }
