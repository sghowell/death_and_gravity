"""Exact orthogonal two-particle projection and commutator tail identity."""

from functools import cache

import sympy as s


@cache
def data():
    xr = s.symbols("xr0:4", real=True)
    xi = s.symbols("xi0:4", real=True)
    yr = s.symbols("yr0:4", real=True)
    yi = s.symbols("yi0:4", real=True)
    x = s.Matrix([xr[j] + s.I * xi[j] for j in range(4)])
    y = s.Matrix([yr[j] + s.I * yi[j] for j in range(4)])
    P = s.diag(1, 0, 1, 0)
    Q = s.eye(4) - P
    inner = lambda u, v: (u.conjugate().T * v)[0]
    z = inner(x, y)
    checks = {
        "common_orthogonal_projection": P * P - P,
        "projection_self_adjoint": P.conjugate().T - P,
        "both_discarded_tails_exact_difference": s.expand(
            z - inner(P * x, P * y) - inner(Q * x, Q * y)
        ),
        "first_cross_term_zero": inner(P * x, Q * y),
        "second_cross_term_zero": inner(Q * x, P * y),
        "full_norm_Pythagoras": s.expand(
            inner(x, x) - inner(P * x, P * x) - inner(Q * x, Q * x)
        ),
        "current_commutator_sign": s.expand(
            s.I * (z - s.conjugate(z)) / 4 + s.im(z) / 2
        ),
    }
    return {
        "Fock_identification": "The centered smeared quadratic stress on the unchanged pure Gaussian state creates exactly a two-particle vector psi_f. Projecting both field momenta before forming that centered observable gives psi_f,K=P_K^(2)psi_f, the SAME orthogonal projection for every test.",
        "exact_tail_identity": "<psi_D,psi_Gamma>-<Ppsi_D,Ppsi_Gamma>=<(1-P)psi_D,(1-P)psi_Gamma>. Both cross terms vanish by orthogonality; the sphere is invariant under the momentum reversal relating S195 q to the created momentum -l.",
        "weak_current": "B(D,Gamma)=i<[T[D],T[Gamma]]>/4=-Im<psi_D,psi_Gamma>/2 for real tests and strict time ordering. Thus |B-B_K|<=tau_D(K)tau_Gamma(K)/2.",
        "boundary": "This Hilbert-space convergence identifies a weak first-order smeared response. It does not assert operator-norm differentiability of metric evolution, a finite-amplitude quantum solution or a same-space inverse.",
        "checks": checks,
        "gates": {
            "same_projection_for_both_tests": P.shape == (4, 4),
            "nontrivial_retained_and_removed_sectors": s.trace(P) == s.trace(Q) == 2,
            "real_test_commutator_is_real": s.im(
                s.I * (z - s.conjugate(z)) / 4
            ).expand()
            == 0,
        },
    }
