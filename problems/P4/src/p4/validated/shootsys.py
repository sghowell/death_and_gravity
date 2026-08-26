"""The shooting system for A3: constraint-reduced CSS system in the centre-scaled
variables, made autonomous and polynomial by adjoining t = e^x.

Unknowns u = (n, w, v, T) with (n, w, v) = (N e^x, W e^{-2x}, V e^{-x}) (S1's
scaled variables, all O(1) on the whole shooting range) and T = e^x, independent
variable x, theta = d/dx:
    P(T, n, w, v) . (n', w', v')  = Q(T, n, w, v)      (``systems.centre_system``, whose
                                                       t d/dt is exactly d/dx)
    T'                             = T.
Everything here is exact-polynomial (fmpq_mpoly re-mapped into the 5-generator
context (x, n, w, v, T)); Taylor coefficients of solutions at a regular point,
the coefficients of the fundamental (variational) matrix, and interval
enclosures of f = P^{-1} Q and Df over boxes are computed in ball arithmetic.

Why this system: S1 showed the 4D plain system is off-constraint unstable toward
the centre, so the constraint-reduced 3D system must be used; in the plain
variables N ~ e^{-x} and W ~ e^{2x} become badly scaled by x = -4 (W ~ 1e-3 in
denominators), whereas (n, w, v) stay O(1) and tend to the regular fixed point.
"""
from __future__ import annotations

from flint import arb, arb_mat

from . import recursion
from .arbseries import Series, abs_upper, to_arb
from .polysys import PolySystem, make_ctx
from .systems import centre_system


def remap(poly, ctx2, gmap):
    """Re-express ``poly`` in context ``ctx2``; ``gmap[k]`` = index of old generator k."""
    dd = {}
    for exps, c in poly.terms():
        e2 = [0] * ctx2.nvars()
        for k, e in enumerate(exps):
            if e:
                e2[gmap[k]] += int(e)
        dd[tuple(e2)] = c
    return ctx2.from_dict(dd)


def shoot_system():
    """4D autonomous polynomial system P~ u' = Q~, u = (n, w, v, T), theta = d/dx."""
    cs = centre_system()
    ctx = make_ctx(("x", "n", "w", "v", "T"))
    gmap = {0: 4, 1: 1, 2: 2, 3: 3}                      # old t -> T, (n, w, v) keep positions
    z = ctx.constant(0)
    T = ctx.gens()[4]
    P = [[z] * 4 for _ in range(4)]
    Q = [z] * 4
    for r in range(3):
        for i in range(3):
            P[r][i] = remap(cs.P[r][i], ctx, gmap)
        Q[r] = remap(cs.Q[r], ctx, gmap)
    P[3][3] = 1 + z
    Q[3] = T
    return PolySystem(ctx, P, Q, euler=False)


def regular_level_equations(d=4):
    """E_n = F_{n-1}: the ordinary Taylor recursion at a regular point (M_n = n P_0)."""
    return recursion.LevelEquations([[(arb(1), r, -1)] for r in range(d)])


# ----------------------------------------------------------------------------
# Taylor coefficients of solutions and of the fundamental matrix
# ----------------------------------------------------------------------------
def _mat(rows):
    return arb_mat(rows)


def _tolists(M):
    return [[M[i, j] for j in range(M.ncols())] for i in range(M.nrows())]


def taylor_coefficients(sys, u0, p, eqs=None, blocks=None):
    """Taylor coefficients u_0..u_p (list over i of d-lists of arb) of the solution
    through ``u0`` (arb balls; an interval u0 gives enclosures for every point in it).

    With ``blocks`` (block lower-triangular P) each order is solved by block forward
    substitution:  n P_0 u_n = -F_{n-1}(u_{<n})  (interval LU on the full matrix is
    far too pessimistic on wide boxes)."""
    d = sys.d
    if blocks is None:
        eqs = regular_level_equations(d) if eqs is None else eqs
        coefs, _ = recursion.solve_recursion(sys, eqs, [list(u0)], p)
        return coefs
    coefs = [list(u0)]
    P0 = [[eval_box(sys.P[r][l], [arb(0)] + list(u0)) for l in range(d)] for r in range(d)]
    for n in range(1, p + 1):
        u = recursion.series_from_coefs(coefs, d, cap=n, extra_zero=False)
        F = sys.residual(u, cap=n, t=Series.var(n))
        rhs = [[-F[r][n - 1] / n] for r in range(d)]
        sol = block_solve(P0, rhs, blocks)
        coefs.append([sol[r, 0] for r in range(d)])
    return coefs


def variational_coefficients(sys, coefs, p, blocks=None):
    """Taylor coefficients Y_0..Y_p (arb_mat) of the fundamental matrix Y' = Df(u(t)) Y,
    Y(0) = I, along the solution series ``coefs`` (which must have length >= p+1).

    From P(u) u' = Q(u):  P Y' = (DQ - Psi) Y =: G Y, Psi_{rl} = sum_i (u_i') dP_{ri}/du_l, so
        P_0 i Y_i = sum_{j<i} G_j Y_{i-1-j} - sum_{1<=j<i} P_j (i-j) Y_{i-j}."""
    d = sys.d
    cap = p + 1
    u = recursion.series_from_coefs(coefs[:cap], d, cap=cap, extra_zero=False)
    t = Series.var(cap)
    Pser = sys.P_series(u, cap=cap, t=t)
    dQ = sys.dQ_series(u, cap=cap, t=t)
    Psi = sys.psi_series(u, cap=cap, t=t)
    Pk = [_mat([[Pser[r][l][k] for l in range(d)] for r in range(d)]) for k in range(cap)]
    Gk = [_mat([[dQ[r][l][k] - Psi[r][l][k] for l in range(d)] for r in range(d)]) for k in range(cap)]
    P0l = _tolists(Pk[0])
    Y = [arb_mat(d, d, [arb(1) if i == j else arb(0) for i in range(d) for j in range(d)])]
    for i in range(1, p + 1):
        rhs = Gk[0] * Y[i - 1]
        for j in range(1, i):
            rhs = rhs + Gk[j] * Y[i - 1 - j] - Pk[j] * Y[i - j] * (i - j)
        Yi = Pk[0].solve(rhs) if blocks is None else block_solve(P0l, _tolists(rhs), blocks)
        Y.append(Yi / i)
    return Y


# ----------------------------------------------------------------------------
# interval enclosures of the right-hand side and its Jacobian over a box
# ----------------------------------------------------------------------------
def _pow(x, e):
    """x**e by repeated multiplication (python-flint's arb ** int is NaN on balls
    containing 0 with nonzero radius, e.g. [+/- 1e-20] ** 1)."""
    out = None
    base = x
    while e:
        if e & 1:
            out = base if out is None else out * base
        e >>= 1
        if e:
            base = base * base
    return arb(1) if out is None else out


def eval_box(poly, args):
    """Rigorous ball evaluation of an fmpq_mpoly at arb ``args`` (t first); box-safe."""
    tot = arb(0)
    for exps, coef in poly.terms():
        term = to_arb(coef)
        for k, e in enumerate(map(int, exps)):
            if e:
                term *= _pow(args[k], e)
        tot += term
    return tot


def block_solve(P, B, blocks):
    """Solve P X = B for a block lower-triangular P (lists of arb; ``blocks`` = list of
    (start, stop)) by block forward substitution -- far tighter than interval LU on the
    full matrix.  Returns X as an arb_mat."""
    ncols = len(B[0])
    X = [[None] * ncols for _ in range(len(P))]
    for (a, b) in blocks:
        rhs = [[B[r][c] for c in range(ncols)] for r in range(a, b)]
        for (a2, b2) in blocks:
            if b2 > a:
                break
            for r in range(a, b):
                for c in range(ncols):
                    for l in range(a2, b2):
                        rhs[r - a][c] -= P[r][l] * X[l][c]
        sub = _mat([[P[r][l] for l in range(a, b)] for r in range(a, b)]).solve(_mat(rhs))
        for r in range(a, b):
            for c in range(ncols):
                X[r][c] = sub[r - a, c]
    return _mat(X)


def rhs_enclosure(sys, box, blocks=None):
    """f(box) = P(box)^{-1} Q(box) as a list of arb (rigorous over the box)."""
    args = [arb(0)] + list(box)
    d = sys.d
    blocks = [(0, d)] if blocks is None else blocks
    P = [[eval_box(sys.P[r][l], args) for l in range(d)] for r in range(d)]
    Q = [[eval_box(sys.Q[r], args)] for r in range(d)]
    f = block_solve(P, Q, blocks)
    return [f[r, 0] for r in range(d)], P


def jacobian_enclosure(sys, box, f=None, P=None, blocks=None):
    """Df(box) = P^{-1} (DQ - Psi) as an arb_mat, Psi_{rl} = sum_i f_i dP_{ri}/du_l."""
    args = [arb(0)] + list(box)
    d = sys.d
    blocks = [(0, d)] if blocks is None else blocks
    if f is None:
        f, P = rhs_enclosure(sys, box, blocks)
    dP, dQ = sys.dP(), sys.dQ()
    G = [[arb(0)] * d for _ in range(d)]
    for r in range(d):
        for l in range(d):
            acc = eval_box(dQ[r][l], args)
            for i in range(d):
                if not dP[r][i][l].is_zero():
                    acc -= f[i] * eval_box(dP[r][i][l], args)
            G[r][l] = acc
    return block_solve(P, G, blocks)


def horner_vec(coefs, h):
    """sum_i u_i h^i for a list of d-lists ``coefs``; returns a list of arb."""
    d = len(coefs[0])
    out = []
    for i in range(d):
        acc = arb(0)
        for c in reversed(coefs):
            acc = acc * h + c[i]
        out.append(acc)
    return out


def horner_mat(mats, h):
    acc = mats[-1]
    for M in reversed(mats[:-1]):
        acc = acc * h + M
    return acc


class Hessian:
    """Second derivative of f = P^{-1} Q (w.r.t. u) over a box, from the exact polynomials:
    differentiating P f = Q twice,
        P D2f[i,j] = D2Q[i,j] - sum_l D2P_{.l}[i,j] f_l - sum_l (DP_{.l}[i] Df_{lj} + DP_{.l}[j] Df_{li}).
    ``norm(box, blocks, Sc)`` returns max_r sum_{i,j} |D2f_r[i,j]| Sc[r]/(Sc[i] Sc[j]), so that
    ||Df(z) - Df(z')||_w <= norm * ||z - z'||_w for z, z' in the (convex) box (weighted inf-norms)."""

    def __init__(self, sys):
        self.sys = sys
        d = sys.d
        dP, dQ = sys.dP(), sys.dQ()
        self.d2P = [[[[dP[r][l][i].derivative(j + 1) for j in range(d)] for i in range(d)]
                     for l in range(d)] for r in range(d)]
        self.d2Q = [[[dQ[r][i].derivative(j + 1) for j in range(d)] for i in range(d)] for r in range(d)]
        self.dP, self.dQ = dP, dQ

    def norm(self, box, blocks, Sc):
        sys, d = self.sys, self.sys.d
        args = [arb(0)] + list(box)
        f, P = rhs_enclosure(sys, box, blocks)
        Df = jacobian_enclosure(sys, box, f, P, blocks)
        dPv = [[[eval_box(self.dP[r][l][i], args) if not self.dP[r][l][i].is_zero() else None
                 for i in range(d)] for l in range(d)] for r in range(d)]
        best = arb(0)
        rows = [arb(0)] * d
        for i in range(d):
            for j in range(i, d):
                v = []
                for r in range(d):
                    acc = eval_box(self.d2Q[r][i][j], args) if not self.d2Q[r][i][j].is_zero() else arb(0)
                    for l in range(d):
                        if not self.d2P[r][l][i][j].is_zero():
                            acc -= eval_box(self.d2P[r][l][i][j], args) * f[l]
                        if dPv[r][l][i] is not None:
                            acc -= dPv[r][l][i] * Df[l, j]
                        if dPv[r][l][j] is not None:
                            acc -= dPv[r][l][j] * Df[l, i]
                    v.append([acc])
                x = block_solve(P, v, blocks)
                mult = 1 if i == j else 2
                for r in range(d):
                    rows[r] += abs_upper(x[r, 0]) * (mult * Sc[r] / (Sc[i] * Sc[j]))
        for r in range(d):
            best = best.max(rows[r])
        return best
