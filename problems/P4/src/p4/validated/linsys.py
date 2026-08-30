"""Linearised CSS system for complex kappa (Theorem B, Stage 1): exact polynomial form
and ball evaluation of its coefficient matrix on a certified background.

Perturbation h = H_ss(x) + eps h_p(x) e^{kappa s} of the KHA system in the cleared
polynomial form of ``systems.sonic_system`` (P(u) u_x + P_s(u) u_s = Q(u),
u = (A, N, W, V); P_s = s-derivative coefficients of the fluid rows, rows 3-4 of
``p4.css`` times 3WS resp. WS).  Linearised:
    P(u) p' + Psi(u, u') p + kappa P_s(u) p = DQ(u) p,   Psi_{rl} = sum_i u_i' dP_{ri}/du_l,
and the linearised momentum constraint (C~ = 0 on the background, dC~/dA = S)
    (kappa - A) S A_p = A (C~_N N_p + C~_W W_p + C~_V V_p)
eliminates A_p (S1: kappa A_p = dC.y_p with C = A C~/S).  Rows 2-4 contain no A_p';
multiplied by (kappa - A) S they form the linear 3D system for q = (N_p, W_p, V_p):
    Pk(x; kappa) q' = Gk(x; kappa) q,
    Pk = (kappa - A) S P_{rq},      Gk = (kappa - A) S (DQ - Psi - kappa P_s)_{rq} + A DQ_{rA} C~_q,
polynomial in (u, u', kappa) with kappa-degree <= 2:
    Pk = Pk0(u) + kappa Pk1(u),   Gk = sum_j kappa^j [Ga_j(u) + sum_i Gb_{j,i}(u) u_i'].
This is the system S1 integrates (``p4.perturb.linearized_rhs``); row 1 of the 4D
linearisation is implied by the Bianchi identity (S1: 1e-11; ``linsonic`` checks it as a
ball identity order by order).  Indices: 4D u = (A, N, W, V) = 0..3; 3D q = 0..2 <-> 4D 1..3.
"""
from __future__ import annotations

from flint import acb, acb_mat, arb, fmpq

from .arbseries import Series, to_arb
from .polysys import _PolyEvaluator
from .shootsys import eval_box
from .systems import sonic_constraint_poly, sonic_system


def to_acb(x):
    if isinstance(x, acb):
        return x
    return acb(x) if isinstance(x, arb) else acb(to_arb(x))


def kappa_box(centre, width):
    """acb box  centre + [-w, w] + i[-w, w]  (centre: decimal string / arb / acb)."""
    return to_acb(centre) + acb(arb(0, width), arb(0, width))


def abs_up(z):
    """Exact arb upper bound of |z| for acb or arb."""
    return z.abs_upper() if isinstance(z, acb) else to_arb(z).abs_upper()


def plain_from_scaled(x, n, w, v):
    """(A, N, W, V) balls from centre-scaled (n, w, v) = (N e^x, W e^{-2x}, V e^{-x}); A from the
    momentum constraint A - 1 = 2WT/S, T = 1 + V^2/3 + (4/3)NV (as in ``css.Am1_constraint``)."""
    ex = to_arb(x).exp()
    N, W, V = n / ex, w * ex * ex, v * ex
    T = 1 + V * V / 3 + N * V * arb(fmpq(4, 3))
    return [1 + 2 * W * T / (1 - V * V), N, W, V]


class LinSystem:
    """Exact polynomial data (fmpq_mpoly in (t, A, N, W, V)) of the 3D linearised system."""

    def __init__(self):
        sys = sonic_system()
        self.sys, ctx = sys, sys.ctx
        t, A, N, W, V = ctx.gens()
        z = ctx.constant(0)
        S = 1 - V**2
        self.S, self.Apoly = S, A
        self.Ps = [[z] * 4, [z] * 4, [z, z, 3 * S, 4 * V * W], [z, z, 4 * V * S, 4 * (1 + V**2) * W]]
        self.C = sonic_constraint_poly(ctx)
        self.dC = [self.C.derivative(k + 1) for k in range(4)]
        dP, dQ = sys.dP(), sys.dQ()
        R, Qi = range(3), range(3)
        self.Pk = [[[-A * S * sys.P[r + 1][q + 1] for q in Qi] for r in R],
                   [[S * sys.P[r + 1][q + 1] for q in Qi] for r in R]]
        self.Ga = [[[-A * S * dQ[r + 1][q + 1] + A * dQ[r + 1][0] * self.dC[q + 1] for q in Qi] for r in R],
                   [[S * dQ[r + 1][q + 1] + A * S * self.Ps[r + 1][q + 1] for q in Qi] for r in R],
                   [[-S * self.Ps[r + 1][q + 1] for q in Qi] for r in R]]
        self.Gb = [[[[A * S * dP[r + 1][i][q + 1] for i in range(4)] for q in Qi] for r in R],
                   [[[-S * dP[r + 1][i][q + 1] for i in range(4)] for q in Qi] for r in R]]

    # -- point evaluation on background balls ------------------------------------------
    def background_derivative(self, u):
        """u' = (A', N', W', V') at a regular point of the 4D system (Cramer on the fluid block)."""
        args = [arb(0)] + list(u)
        P, Q = self.sys.P, self.sys.Q
        ev = lambda p: eval_box(p, args)
        a, b, c, d = ev(P[2][2]), ev(P[2][3]), ev(P[3][2]), ev(P[3][3])
        q3, q4 = ev(Q[2]), ev(Q[3])
        det = a * d - b * c
        return [ev(Q[0]) / ev(P[0][0]), ev(Q[1]) / ev(P[1][1]), (q3 * d - b * q4) / det, (a * q4 - c * q3) / det]

    def matrices(self, u, kappa, du=None):
        """(Pk, Gk) as 3x3 acb_mat at a point: background balls u (4-list), kappa (acb ball/box)."""
        du = self.background_derivative(u) if du is None else du
        args = [arb(0)] + list(u)
        k = to_acb(kappa)
        Pk = [[k * eval_box(self.Pk[1][r][q], args) + eval_box(self.Pk[0][r][q], args) for q in range(3)]
              for r in range(3)]
        Gk = []
        for r in range(3):
            row = []
            for q in range(3):
                acc = eval_box(self.Ga[0][r][q], args) + k * eval_box(self.Ga[1][r][q], args)
                acc += k * k * eval_box(self.Ga[2][r][q], args)
                for i in range(4):
                    if not self.Gb[0][r][q][i].is_zero():
                        acc += du[i] * eval_box(self.Gb[0][r][q][i], args)
                    if not self.Gb[1][r][q][i].is_zero():
                        acc += k * du[i] * eval_box(self.Gb[1][r][q][i], args)
                row.append(acc)
            Gk.append(row)
        return acb_mat(Pk), acb_mat(Gk)

    def coefficient_matrix(self, u, kappa, du=None):
        """A(x; kappa) with q' = A q  (regular point: (kappa - A) S Delta~ != 0 required)."""
        Pk, Gk = self.matrices(u, kappa, du)
        return Pk.solve(Gk)

    def A_p(self, u, kappa, q):
        """A_p from the linearised momentum constraint, q = (N_p, W_p, V_p) acb."""
        args = [arb(0)] + list(u)
        num = sum((eval_box(self.dC[l + 1], args) * q[l] for l in range(3)), acb(0))
        return u[0] * num / ((to_acb(kappa) - u[0]) * eval_box(self.S, args))

    # -- series evaluation along a background series -----------------------------------
    def series_parts(self, u, cap=None):
        """Real Series matrices (Pser[j], Gser[j]) along the background series u (4 Series):
        Pk = sum_j kappa^j Pser[j] (j = 0, 1), Gk = sum_j kappa^j Gser[j] (j = 0, 1, 2)."""
        u = [ui.with_cap(cap) for ui in u]
        ev = _PolyEvaluator([Series.var(cap)] + u)
        th = [ui.deriv() for ui in u]
        Pser = [[[ev(self.Pk[j][r][q]) for q in range(3)] for r in range(3)] for j in range(2)]
        Gser = []
        for j in range(3):
            Mj = []
            for r in range(3):
                row = []
                for q in range(3):
                    acc = ev(self.Ga[j][r][q])
                    for i in (range(4) if j < 2 else ()):
                        if not self.Gb[j][r][q][i].is_zero():
                            acc = acc + ev(self.Gb[j][r][q][i]) * th[i]
                    row.append(acc)
                Mj.append(row)
            Gser.append(Mj)
        return Pser, Gser


def background_from_state(st, point=False):
    """(A, N, W, V) balls at st.x from a ``tmint.State`` (the A3 tube: point set or interval set)."""
    u = st.u_point() if point else st.u_interval()
    return plain_from_scaled(st.x, u[0], u[1], u[2])
