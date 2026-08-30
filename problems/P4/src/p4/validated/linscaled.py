"""Theorem B, Stage 2: the linearised systems in the centre-scaled variables, as exact
kappa-graded polynomial matrices in z = (n, w, v, T) and z' = (n', w', v') (fmpq_mpoly).

Scaled perturbation  q~ = Lambda^{-1} q,  Lambda = diag(T^{-1}, T^2, T)  (3D: q = (N_p, W_p, V_p)),
resp. Lambda_4 = diag(T^2, T^{-1}, T^2, T) for the 4D system p = (A_p, N_p, W_p, V_p) [A^_p = A_p e^{-2x}].
With u = (A, N, W, V) = ((S + 2 w T^2 Tt)/S, n/T, w T^2, v T), S = 1 - v^2 T^2, Tt = 1 + v^2T^2/3 + 4nv/3,
and u' expressed through (z, z') (A' by differentiating the constraint), the Stage-1 system
    sum_j kappa^j Pk_j(u) q' = sum_j kappa^j [Ga_j(u) + sum_i Gb_{j,i}(u) u_i'] q
becomes, after multiplying every row by the common nonvanishing factor S^{dA+2} T^{dN+2} (dA, dN the
maximal A- and N-degrees) and stripping the common T-power of each row (exact division),
    sum_j kappa^j P_j(z) q~' = sum_j kappa^j G_j(z, z') q~,     G_j := [...] Lambda - P_j Lambda Gamma,
Gamma = diag(-1, 2, 1) (resp. (2, -1, 2, 1)).  Same solutions (T > 0, S > 0 on the shooting range);
the coefficients are polynomials in the 8-generator context (x, n, w, v, T, dn, dw, dv) (x unused).
The kappa-derivative system (q~, dq~/dkappa) is assembled block-wise (``augmented``).
"""
from __future__ import annotations

from flint import acb_mat, arb, fmpq

from .arbseries import Series, abs_upper
from .linsys import LinSystem
from .polysys import _PolyEvaluator, make_ctx
from .shootsys import eval_box

CTX = make_ctx(("x", "n", "w", "v", "T", "dn", "dw", "dv"))
_X, _n, _w, _v, _T, _dn, _dw, _dv = CTX.gens()
_S = 1 - _v**2 * _T**2
_Tt = 1 + _v**2 * _T**2 * fmpq(1, 3) + _n * _v * fmpq(4, 3)
_Anum = _S + 2 * _w * _T**2 * _Tt                                    # A = Anum / S
_Ttp = _v * (_dv + _v) * _T**2 * fmpq(2, 3) + (_dn * _v + _n * _dv) * fmpq(4, 3)
_Sp = -2 * _v * (_dv + _v) * _T**2
_numAp = 2 * _T**2 * ((_dw + 2 * _w) * _Tt * _S + _w * _Ttp * _S - _w * _Tt * _Sp)   # S^2 A'
# S^2 u_i' * T^{shift_i} as polynomials (index order A, N, W, V)
_DU = [(_numAp, 0), (_S**2 * (_dn - _n), -1), (_S**2 * (_dw + 2 * _w) * _T**2, 0), (_S**2 * (_dv + _v) * _T, 0)]


def _prod(a, exps):
    """prod a[k-1]^{e_k} over the non-t generators (arb)."""
    out = arb(1)
    for k, e in enumerate(map(int, exps)):
        if e and k:
            out *= a[k - 1] ** e
    return out


def _substitute(p, dA, dN):
    """S^{dA} T^{dN} p(u(z)) for a Stage-1 polynomial p in (t, A, N, W, V)."""
    out = CTX.constant(0)
    for exps, c in p.terms():
        et, eA, eN, eW, eV = map(int, exps)
        assert et == 0
        term = _Anum**eA * _S**(dA - eA) * _n**eN * _T**(dN - eN) * (_w * _T**2)**eW * (_v * _T)**eV
        out += term * c
    return out


def _tval(p):
    return min((int(e[4]) for e, _ in p.terms()), default=None)


def _strip_rows(P, G, d):
    for r in range(d):
        vals = [_tval(M[r][c]) for M in P + G for c in range(d)]
        v = min(x for x in vals if x is not None)
        for M in P + G:
            for c in range(d):
                M[r][c] = M[r][c] // _T**v
    return P, G


class ScaledLinearSystem:
    """kappa-graded polynomial matrices: P = sum_j kappa^j P[j], G = sum_j kappa^j G[j] (d x d)."""

    def __init__(self, P, G, d, gamma, name="", weights=None):
        self.P, self.G, self.d, self.gamma, self.name = P, G, d, gamma, name
        self.weights = weights or [1.0] * d          # norm weights for the Groenwall bounds
        self.d0 = None                               # base dimension of a kappa-derivative augmentation

    def box_matrices(self, z, dz):
        args = [arb(0)] + list(z) + list(dz)
        f = lambda M: [[eval_box(M[r][c], args) for c in range(self.d)] for r in range(self.d)]
        return [f(M) for M in self.P], [f(M) for M in self.G]

    def series_matrices(self, z, dz, cap=None):
        """Series matrices along z(s), z'(s) (4 + 3 Series); exact products if cap is None."""
        ev = _PolyEvaluator([Series.var(cap)] + [s.with_cap(cap) for s in list(z) + list(dz)])
        f = lambda M: [[ev(M[r][c]) for c in range(self.d)] for r in range(self.d)]
        return [f(M) for M in self.P], [f(M) for M in self.G]

    def majorant_tails(self, amaj, a, K, h):
        """Bound of sum_{k>K} |M_k| h^k for the exact products M(z_K(s)) from the coefficientwise
        majorant p^abs(A(s)), A(s) = sum_i |z_i| s^i (``amaj``: 7 majorant Series; ``a`` = A(h)):
        p^abs(a) - sum_{k<=K} [p^abs(A(s))]_k h^k  (entrywise arb, nonnegative)."""
        ev = _PolyEvaluator([Series.var(K + 1)] + [s.with_cap(K + 1) for s in amaj])
        hp = [arb(1)]
        for _ in range(K):
            hp.append(hp[-1] * h)

        def tail(p):
            acc = ev.vars[0].zero()
            full = arb(0)
            for exps, coef in p.terms():
                term = None
                for k, e in enumerate(exps):
                    if e:
                        term = ev.power(k, e) if term is None else term * ev.power(k, e)
                c = arb(abs(fmpq(coef)))
                acc = acc + (c if term is None else term * c)
                full += c * _prod(a, exps)
            low = sum((abs_upper(acc[k]) * hp[k] for k in range(K + 1)), arb(0))
            return (full - low).nonnegative_part().abs_upper()
        f = lambda M: [[tail(M[r][c]) for c in range(self.d)] for r in range(self.d)]
        return [f(M) for M in self.P], [f(M) for M in self.G]

    def increments(self, a, eps):
        """Entrywise majorants of |M(z, z') - M(z_K, z_K')| for |z_K| <= a, |z - z_K| <= eps (7-lists)."""
        def inc(p):
            tot = arb(0)
            for exps, coef in p.terms():
                base, shifted = arb(1), arb(1)
                for k, e in enumerate(map(int, exps)):
                    if e and k:
                        base *= a[k - 1] ** e
                        shifted *= (a[k - 1] + eps[k - 1]) ** e
                tot += arb(abs(fmpq(coef))) * (shifted - base).nonnegative_part()
            return tot
        f = lambda M: [[inc(M[r][c]) for c in range(self.d)] for r in range(self.d)]
        return [f(M) for M in self.P], [f(M) for M in self.G]


def _build(Pk, Ga, Gb, d, expo, gamma, name):
    """Pk[j][r][c], Ga[j][r][c], Gb[j][r][c][i] Stage-1 polynomials -> ScaledLinearSystem."""
    allp = [p for M in Pk + Ga for row in M for p in row] + [p for M in Gb for row in M for e in row for p in e]
    dA = max(p.degrees()[1] for p in allp)
    dN = max(p.degrees()[2] for p in allp)
    sub = lambda p: _substitute(p, dA, dN)
    P = [[[sub(Pk[j][r][c]) * _S**2 * _T ** expo[c] for c in range(d)] for r in range(d)] for j in range(len(Pk))]
    G = []
    for j in range(len(Ga)):
        M = []
        for r in range(d):
            row = []
            for c in range(d):
                acc = sub(Ga[j][r][c]) * _S**2 * _T ** expo[c]
                if j < len(Gb):
                    for i in range(4):
                        if not Gb[j][r][c][i].is_zero():
                            dpoly, sh = _DU[i]
                            acc += sub(Gb[j][r][c][i]) * dpoly * _T ** (expo[c] + sh)
                if j < len(P):
                    acc -= P[j][r][c] * gamma[c]
                row.append(acc)
            M.append(row)
        G.append(M)
    P, G = _strip_rows(P, G, d)
    return ScaledLinearSystem(P, G, d, gamma, name)


def reduced_system(L=None):
    """3D constraint-reduced system (Stage 1's, with the (kappa - A) factor), scaled."""
    L = L or LinSystem()
    return _build(L.Pk, L.Ga, L.Gb, 3, (1, 4, 3), (-1, 2, 1), "3D")


def full_system(L=None):
    """4D system P p' = (DQ - Psi - kappa P_s) p with A_p kept (no elimination factor), scaled."""
    L = L or LinSystem()
    sys, dP, dQ = L.sys, L.sys.dP(), L.sys.dQ()
    z = sys.ctx.constant(0)
    Pk = [[[sys.P[r][c] for c in range(4)] for r in range(4)]]
    Ga = [[[dQ[r][c] for c in range(4)] for r in range(4)], [[-L.Ps[r][c] for c in range(4)] for r in range(4)]]
    Gb = [[[[-dP[r][i][c] for i in range(4)] for c in range(4)] for r in range(4)]]
    return _build(Pk, Ga, Gb, 4, (4, 1, 4, 3), (2, -1, 2, 1), "4D")


def augmented(S, lam=2.0**-12):
    """System for (q~, dq~/dkappa): P^ = [[P,0],[P_k,P]], G^ = [[G,0],[G_k,G]], P_k = P_1, G_k = G_1 + 2 kappa G_2."""
    d, z = S.d, CTX.constant(0)
    P = S.P + [[[z] * d for _ in range(d)]] * (2 - len(S.P))
    G = S.G + [[[z] * d for _ in range(d)]] * (3 - len(S.G))
    blk = lambda A, B, C, D: [[(A if r < d and c < d else B if r < d else C if c < d else D)[r % d][c % d]
                               for c in range(2 * d)] for r in range(2 * d)]
    Z = [[z] * d for _ in range(d)]
    twoG2 = [[G[2][r][c] * 2 for c in range(d)] for r in range(d)]
    Ph = [blk(P[0], Z, P[1], P[0]), blk(P[1], Z, Z, P[1])]
    Gh = [blk(G[0], Z, G[1], G[0]), blk(G[1], Z, twoG2, G[1]), blk(G[2], Z, Z, G[2])]
    out = ScaledLinearSystem(Ph, Gh, 2 * d, tuple(S.gamma) * 2, S.name + "+dk", [1.0] * d + [lam] * d)
    out.d0 = d
    return out
