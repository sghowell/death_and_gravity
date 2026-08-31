"""Theorem C, S4-1: exact 2nd-order reduction of the 4D linearised system (route (b)).

Everything is exact ``fmpq_mpoly`` arithmetic in the context (t, A, N, W, V, Ap, Np,
Wp, Vp, k), in the style of ``linsys.linear_constraint_propagation``: background
derivatives by Cramer (du_i = S*Delta~ * u_i'), perturbation derivatives from the 4D
linearised system  P(u) p' = [DQ - Psi - kappa P_s] p  (dp_j = S^2 Delta~^2 p_j'),
A eliminated at the end by the background constraint A = Anum/S (``elim``).

Identities established (all verified by exact cancellation, see the tests):
  (G)  the gauge vector  g = (A', N' + kappa N, W', V')  solves the 4D linearised
       system *identically* (zero polynomial residual, no elimination needed) and
       satisfies the linearised momentum constraint modulo the background constraint;
  (R)  with the g-annihilating pair (N-clock gauge-invariant fluid perturbations)
           chi := Delta~ [(N' + kappa N) V_p - V' N_p]   (prop. to the v_rel- and
                                                          V-perturbation, see note),
           eta := Delta~ [(N' + kappa N) W_p - W' N_p],
       the closure identities
           S^2 Delta~^2 dn e . chi' = Achi chi + Bchi eta + Cchi c(p),
           S^2 Delta~^2 dn e . eta' = Aeta chi + Beta eta + Ceta c(p),
       hold modulo the background equations, where dn = (kappa - A) S, e = Delta~
       (N' + kappa N), and c(p) is the linearised momentum constraint functional
       (c = 0 on the constraint surface Sigma).  Hence on Sigma the pair (chi, eta)
       satisfies a closed 2x2 first-order system, and chi satisfies the scalar
       2nd-order equation obtained by eliminating eta (``scalar form'' of the note
       notes/s4-reduction.md: (w chi')' + a1 chi' + a0 chi = 0).
The quotient interpretation: g spans the kernel of (chi, eta) on Sigma, so the pair
is the flow induced on Sigma / span(g); at kappa = kbar the admissible solution IS g
and chi == 0.  Run ``python -m p4.validated.tc_reduce`` for the identity checks.
"""
from __future__ import annotations

from flint import fmpq, fmpq_mpoly_ctx

from .linsys import LinSystem


class Reduction:
    """Exact polynomial data of the reduction; all members are fmpq_mpoly."""

    def __init__(self, L=None):
        L = L or LinSystem()
        sys = L.sys
        ctx = fmpq_mpoly_ctx.get(("t", "A", "N", "W", "V", "Ap", "Np", "Wp", "Vp", "k"))
        self.ctx, self.L = ctx, L
        _t, A, N, W, V, _Ap, _Np, _Wp, _Vp, k = ctx.gens()
        self.A, self.N, self.W, self.V, self.k = A, N, W, V, k
        z = ctx.constant(0)
        pr = lambda q: q.project_to_context(ctx)
        self.P = [[pr(sys.P[r][c]) for c in range(4)] for r in range(4)]
        self.Q = [pr(q) for q in sys.Q]
        dP, dQ = sys.dP(), sys.dQ()
        self.dP = [[[pr(dP[r][i][l]) for l in range(4)] for i in range(4)] for r in range(4)]
        self.dQ = [[pr(dQ[r][l]) for l in range(4)] for r in range(4)]
        self.Ps = [[pr(x) for x in row] for row in L.Ps]
        self.dC = [pr(c) for c in L.dC]
        self.S = pr(L.S)
        P, Q = self.P, self.Q
        self.Delta = P[2][2] * P[3][3] - P[2][3] * P[3][2]                       # = 4 S W D
        self.Dpoly = 3 * N**2 * V**2 - N**2 + 4 * N * V - V**2 + 3
        self.SD = self.S * self.Delta
        self.PW = Q[2] * P[3][3] - Q[3] * P[2][3]                                # Delta~ W'
        self.PV = P[2][2] * Q[3] - P[3][2] * Q[2]                                # Delta~ V'
        self.du = [Q[0] * self.Delta, Q[1] * self.SD, self.PW * self.S, self.PV * self.S]
        self.Psi = [[sum((self.du[i] * self.dP[r][i][l] for i in range(4)), z) for l in range(4)]
                    for r in range(4)]
        self.Gcov = [[self.SD * self.dQ[r][l] - self.Psi[r][l] - k * self.SD * self.Ps[r][l]
                      for l in range(4)] for r in range(4)]
        G = self.Gcov
        self.Mdp = [[G[0][l] * self.Delta for l in range(4)],
                    [G[1][l] * self.SD for l in range(4)],
                    [(P[3][3] * G[2][l] - P[2][3] * G[3][l]) * self.S for l in range(4)],
                    [(P[2][2] * G[3][l] - P[3][2] * G[2][l]) * self.S for l in range(4)]]
        self.Anum = self.S + 2 * W * (1 + V**2 * fmpq(1, 3) + N * V * fmpq(4, 3))
        self.ccov = [(k - A) * self.S, -A * self.dC[1], -A * self.dC[2], -A * self.dC[3]]
        self.dn = (k - A) * self.S
        self.e = self.Delta * (Q[1] + k * N)                                      # Delta~ (N'+kN)
        self.lchi = [z, -self.PV, z, self.e]
        self.leta = [z, -self.PW, self.e, z]

    # -- derivations -------------------------------------------------------------------
    def D1(self, f):
        """S Delta~ df/dx along the background flow (f background-only)."""
        return sum((f.derivative(i + 1) * self.du[i] for i in range(4)), self.ctx.constant(0))

    def elim(self, f):
        """(S^degA f)|_{A -> Anum/S}: exact restriction to the background constraint."""
        dA, out = f.degrees()[1], self.ctx.constant(0)
        for exps, coef in f.terms():
            ee = list(map(int, exps))
            eA, ee[1] = ee[1], 0
            out += self.ctx.from_dict({tuple(ee): coef}) * self.Anum**eA * self.S ** (dA - eA)
        return out

    # -- the gauge solution ------------------------------------------------------------
    def gauge_vector(self):
        """ghat = S Delta~ (A', N' + kappa N, W', V')."""
        return [self.du[0], self.du[1] + self.k * self.N * self.SD, self.du[2], self.du[3]]

    def gauge_residual(self):
        """(S Delta~)^3 [P g' - (DQ - Psi - kappa P_s) g]: 4 polynomials, all == 0
        IDENTICALLY (the gauge identity of s2-theorem-b.md section 3.7, now exact)."""
        g, dSD = self.gauge_vector(), self.D1(self.SD)
        z = self.ctx.constant(0)
        return [sum((self.P[r][c] * (self.SD * self.D1(g[c]) - dSD * g[c]) for c in range(4)), z)
                - self.SD * sum((self.Gcov[r][c] * g[c] for c in range(4)), z) for r in range(4)]

    def gauge_constraint_residual(self):
        """c(ghat): must eliminate to 0 (g lies on the linearised constraint surface)."""
        g = self.gauge_vector()
        return self.elim(sum((self.ccov[i] * g[i] for i in range(4)), self.ctx.constant(0)))

    # -- the closure identities ----------------------------------------------------------
    def closure(self, l):
        """Covector of S^2 Delta~^2 (l . p)' along the 4D linearised flow."""
        z = self.ctx.constant(0)
        return [self.SD * self.D1(l[j]) + sum((l[i] * self.Mdp[i][j] for i in range(4)), z)
                for j in range(4)]

    def reduce_pair(self):
        """Multipliers of the closure identities and their exactness residuals.

        Returns dict with, for each of "chi", "eta": (Acoef, Bcoef, Ccoef, resid) s.t.
            dn e L = Acoef lchi + Bcoef leta + Ccoef ccov  + (resid in the Np-slot),
        L = closure covector; elim(resid) == 0 is the exact reduction identity.  Also
        "D2" = S^2 Delta~^2 dn e (the common denominator of the 2x2 system)."""
        out = {"D2": self.SD * self.SD * self.dn * self.e}
        for name, l in (("chi", self.lchi), ("eta", self.leta)):
            Lc = self.closure(l)
            Ac = self.dn * Lc[3] + self.A * self.dC[3] * Lc[0]
            Bc = self.dn * Lc[2] + self.A * self.dC[2] * Lc[0]
            Cc = self.e * Lc[0]
            resid = self.dn * self.e * Lc[1] + Ac * self.PV + Bc * self.PW \
                + Cc * self.A * self.dC[1]
            out[name] = (Ac, Bc, Cc, resid)
        return out

    def identity_tests(self):
        """(gauge_rows_zero, gauge_constraint_zero, chi_closure_zero, eta_closure_zero)."""
        rp = self.reduce_pair()
        return (all(r.is_zero() for r in self.gauge_residual()),
                self.gauge_constraint_residual().is_zero(),
                self.elim(rp["chi"][3]).is_zero(),
                self.elim(rp["eta"][3]).is_zero())

    def apparent_singularities(self):
        """kappa-dependent zero loci of D2 (apparent singularities of the reduced 2x2
        system): dn = 0 at kappa = A(x) (the constraint slaving's pole, real kappa in
        (1, A_max] -- as for the 3D system, s2-theorem-b.md 2.5(c); there (lchi, leta,
        ccov) become dependent on Sigma) and e = 0 at kappa = -F_N(x) (the N-clock
        chart, real kappa in [kbar, 1)).  Both windows lie inside Theorem B's R; off
        the real axis and for |kappa| > A_max only the sonic zero D(x) = 0 remains.
        Returns the two kappa-affine factors ((k - A) S, Q1 + k N)."""
        return self.dn, self.Q[1] + self.k * self.N


if __name__ == "__main__":
    red = Reduction()
    print("gauge identity (identically zero):", all(r.is_zero() for r in red.gauge_residual()))
    print("gauge constraint residual == 0:  ", red.gauge_constraint_residual().is_zero())
    rp = red.reduce_pair()
    print("chi closure exact:", red.elim(rp["chi"][3]).is_zero())
    print("eta closure exact:", red.elim(rp["eta"][3]).is_zero())
