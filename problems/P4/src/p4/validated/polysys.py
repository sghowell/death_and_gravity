"""Polynomial first-order systems  P(t,u) . (theta u) = Q(t,u).

``P`` is a d x d matrix and ``Q`` a d-vector of exact multivariate polynomials
(``flint.fmpq_mpoly``) in the variables (t, u_1, ..., u_d); ``theta`` is either
d/dt (``euler=False``; Taylor expansion in t) or t d/dt (``euler=True``;
expansion in t = e^x at the centre).  The residual of a series u(t) is
    F(u) = P(t,u) . theta u - Q(t,u),
and everything the recursion and the tail bound need is a Taylor coefficient
of F, of P(t,u), of the Jacobians DP, DQ (w.r.t. u), or of
    Psi_{r,l} = sum_i (theta u_i) dP_{r,i}/du_l          (d x d series matrix).

Polynomials are evaluated on ``Series`` objects term by term (cached powers),
which is rigorous ball arithmetic; ``abs_eval`` evaluates the polynomial with
absolute-valued coefficients at nonnegative arguments (the Banach-algebra
majorant used by the tail bound).
"""
from __future__ import annotations

from dataclasses import dataclass, field

from flint import arb, fmpq, fmpq_mpoly_ctx

from .arbseries import Series, to_arb


def make_ctx(names):
    """fmpq_mpoly context with generators ``names`` (first one is t)."""
    return fmpq_mpoly_ctx.get(tuple(names))


@dataclass
class PolySystem:
    ctx: object                        # fmpq_mpoly_ctx, gens = (t, u_1..u_d)
    P: list                            # d x d list of fmpq_mpoly
    Q: list                            # d list of fmpq_mpoly
    euler: bool = False
    _dP: list = field(default=None, repr=False)
    _dQ: list = field(default=None, repr=False)

    @property
    def d(self):
        return len(self.Q)

    @property
    def sigma(self):
        """Index shift of theta: (theta u)_j = s_j u_{j+sigma}."""
        return 0 if self.euler else 1

    def gens(self):
        return self.ctx.gens()

    def theta(self, s):
        return s.euler() if self.euler else s.deriv()

    # -- Jacobians (w.r.t. u only; variable index l+1 in the context) --------
    def dP(self):
        if self._dP is None:
            self._dP = [[[self.P[r][i].derivative(l + 1) for l in range(self.d)]
                         for i in range(self.d)] for r in range(self.d)]
        return self._dP

    def dQ(self):
        if self._dQ is None:
            self._dQ = [[self.Q[r].derivative(l + 1) for l in range(self.d)] for r in range(self.d)]
        return self._dQ

    # -- evaluation on series ----------------------------------------------
    def evaluator(self, t, u):
        return _PolyEvaluator([t] + list(u))

    def residual(self, u, cap=None, t=None):
        """Series list F = P(t,u).theta u - Q(t,u) (truncated at ``cap`` if given)."""
        u = [ui.with_cap(cap) for ui in u]
        t = Series.var(cap) if t is None else t.with_cap(cap)
        ev = self.evaluator(t, u)
        th = [self.theta(ui) for ui in u]
        out = []
        for r in range(self.d):
            acc = -ev(self.Q[r])
            for i in range(self.d):
                if not self.P[r][i].is_zero():
                    acc = acc + ev(self.P[r][i]) * th[i]
            out.append(acc)
        return out

    def P_series(self, u, cap=None, t=None):
        """d x d list of Series of P(t, u(t))."""
        u = [ui.with_cap(cap) for ui in u]
        t = Series.var(cap) if t is None else t.with_cap(cap)
        ev = self.evaluator(t, u)
        return [[ev(self.P[r][i]) for i in range(self.d)] for r in range(self.d)]

    def dQ_series(self, u, cap=None, t=None):
        u = [ui.with_cap(cap) for ui in u]
        t = Series.var(cap) if t is None else t.with_cap(cap)
        ev = self.evaluator(t, u)
        dQ = self.dQ()
        return [[ev(dQ[r][l]) for l in range(self.d)] for r in range(self.d)]

    def psi_series(self, u, cap=None, t=None):
        """Psi_{r,l} = sum_i (theta u_i) * dP_{r,i}/du_l as a d x d list of Series."""
        u = [ui.with_cap(cap) for ui in u]
        t = Series.var(cap) if t is None else t.with_cap(cap)
        ev = self.evaluator(t, u)
        th = [self.theta(ui) for ui in u]
        dP = self.dP()
        out = []
        for r in range(self.d):
            row = []
            for l in range(self.d):
                acc = t.zero()
                for i in range(self.d):
                    if not dP[r][i][l].is_zero():
                        acc = acc + ev(dP[r][i][l]) * th[i]
                row.append(acc)
            out.append(row)
        return out


class _PolyEvaluator:
    """Evaluate fmpq_mpoly's on a fixed tuple of Series with cached powers."""

    def __init__(self, vars_):
        self.vars = vars_
        self.pows = [{0: None} for _ in vars_]

    def power(self, k, e):
        cache = self.pows[k]
        if e not in cache:
            cache[e] = self.vars[k] ** e
        return cache[e]

    def __call__(self, poly):
        acc = self.vars[0].zero()
        for exps, coef in poly.terms():
            term = None
            for k, e in enumerate(exps):
                if e:
                    term = self.power(k, e) if term is None else term * self.power(k, e)
            c = to_arb(coef)
            acc = acc + (c if term is None else term * c)
        return acc


# ----------------------------------------------------------------------------
# majorant evaluation
# ----------------------------------------------------------------------------
def abs_eval(poly, args):
    """sum_terms |coef| prod args_k^{e_k}  for nonnegative arb ``args`` (t first)."""
    tot = arb(0)
    for exps, coef in poly.terms():
        term = to_arb(abs(coef))
        for k, e in enumerate(map(int, exps)):
            if e:
                term *= args[k] ** e
        tot += term
    return tot


def abs_eval_increment(poly, args, eps):
    """Upper bound of  poly^abs(args + eps 1_u) - poly^abs(args)  (t-argument not shifted).

    Computed term by term as  |coef| t^{e_t} [prod (r_k+eps)^{e_k} - prod r_k^{e_k}]  so that
    every term is a nonnegative quantity (no cancellation between terms)."""
    tot = arb(0)
    for exps, coef in poly.terms():
        base = arb(1)
        shifted = arb(1)
        for k, e in enumerate(map(int, exps)):
            if e:
                if k == 0:
                    base *= args[0] ** e
                    shifted *= args[0] ** e
                else:
                    base *= args[k] ** e
                    shifted *= (args[k] + eps) ** e
        tot += to_arb(abs(coef)) * (shifted - base).nonnegative_part()
    return tot


def eval_exact(poly, args):
    """Evaluate at exact rational arguments (fmpq / int); returns fmpq."""
    return poly(*[fmpq(a) for a in args])


def eval_arb(poly, args):
    """Evaluate at arb arguments (t first); rigorous ball."""
    tot = arb(0)
    for exps, coef in poly.terms():
        term = to_arb(coef)
        for k, e in enumerate(map(int, exps)):
            if e:
                term *= args[k] ** e
        tot += term
    return tot
