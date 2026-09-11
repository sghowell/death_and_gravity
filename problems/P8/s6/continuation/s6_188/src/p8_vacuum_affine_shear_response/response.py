"""Exact covariance tangent, Kubo sign and metric-variation contact."""

from functools import cache

import sympy as s

from .hamiltonian import K0, V0


@cache
def data():
    J = s.Matrix([[0, 1], [-1, 0]])

    def sym(prefix):
        a, b, c = s.symbols(prefix + "0:3", real=True)
        return s.Matrix([[a, b], [b, c]])

    A, B, C = sym("readout"), sym("vertex"), sym("covariance")
    tangent = J * B * C - C * B * J
    comm = A * J * B - B * J * A
    W = C + s.I * J / 2
    wick = lambda M, N: s.trace(M * W * N * W.T) / 2
    x, y = s.symbols("x y", real=True)
    S = s.Matrix([[1, x], [y, 1 + x * y]])
    J6 = s.diag(J, J, J)
    Jblock = s.zeros(6)
    Jblock[:3, 3:] = s.eye(3)
    Jblock[3:, :3] = -s.eye(3)
    M0 = s.diag(V0, K0)
    generator = Jblock * M0
    checks = {
        "full_three_mode_generator_symplectic": s.simplify(
            generator * Jblock + Jblock * generator.T
        ),
        "full_three_mode_energy_flow_cancels": s.simplify(
            generator.T * M0 + M0 * generator
        ),
        "covariance_tangent_is_symmetric": s.simplify(tangent - tangent.T),
        "quadratic_commutator_matrix_symmetric": comm - comm.T,
        "covariance_tangent_matches_commutator_trace": s.expand(
            s.trace(A * tangent) - s.trace(comm * C)
        ),
        "Wick_commutator_sign_and_factor": s.expand(
            wick(A, B) - wick(B, A) - s.I * s.trace(comm * C) / 2
        ),
        "independent_symplectic_transport": s.simplify(S * J * S.T - J),
        "transported_CCR_unchanged": s.simplify(S * (W - W.T) * S.T - s.I * J),
        "canonical_phase_dimension": J6.rows - 6,
        "actual_initial_energy_covariance_trace_bound": 3 * (27 + 9) - 108,
    }
    return {
        "canonical_covariance": "C=(1/2)<ZZ^T+(ZZ^T)^T>, C'=A C+C A^T, A=J M; W=<ZZ^T>=C+iJ/2",
        "actual_initial_state": "Use the unchanged S6.55 all-order Cauchy state. For each real Fourier oscillator sector, C0 is built from A_r=f_r/g_r, pi_r=g_r p_r, gT^2=a, gL^2=a m^2/omega^2; p=f'-d f. No W8 reset is made.",
        "initial_energy_trace_upper": "tr[M0(t0)^(1/2) C0 M0(t0)^(1/2)] <=108 nu, nu=sqrt(m^2+|k|^2/a0^2), by the actual-state |f|^2<9/omega, |p|^2<27 omega bounds.",
        "exact_covariance_tangent": "delta C'=A0 delta C+delta C A0^T+eta[J M_D C0-C0 M_D J], with delta C(t0)=0 because the prescribed external metric perturbation vanishes near the common Cauchy surface.",
        "full_current_tangent": "delta J_G(t)=-1/2 tr[M_G(t) delta C(t)]-eta(t)/2 tr[M_GD(t) C0(t)]. The second term is the nonzero metric contact.",
        "Kubo_identity": "Equivalently delta J_G=-<H_GD>eta+i integral_(t0)^t <[H_G(t),H_D(s)]>eta(s)ds, before continuum renormalization.",
        "Wick_identity": "For quadratic H_M=Z^T M Z/2, <H_M H_N>_connected=tr(M W N W^T)/2. The two orderings give the same retarded sign and normalization as the exact covariance tangent.",
        "initial_state_boundary": "This is the response to a prescribed external homogeneous metric history in the common-Cauchy vector state. It does not postulate a factorized or zero quantum metric initial state.",
        "checks": checks,
    }
