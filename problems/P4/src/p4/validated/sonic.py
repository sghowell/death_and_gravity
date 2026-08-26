"""Validated sonic-point expansion (Theorem A, item A1) in the KHA variables of S1.

Gauge: sonic point at x = 0.  Input: centre c of the sonic velocity V0 and an
optional half-width w (V0 in [c-w, c+w]).  Output: for every order n <= K
    u_n(V0) = sum_{k<=m} u_{n,k} (V0-c)^k  +  R_n,   |R_n| <= rem_n,
a degree-m Taylor model in delta = V0 - c whose polynomial part is computed at
the exact centre (tiny radii) and whose remainder is the Lagrange remainder,
  rem_n = sup_{xi in [c-w,c+w]} |u_n^{(m+1)}(xi)/(m+1)!| w^{m+1},
obtained from a second run of the same bivariate recursion with the interval
as base point (the interval blow-up of that run is harmless: it multiplies
w^{m+1}).  ``balls`` are the resulting enclosures of u_n over the whole interval.

Zeroth order: KHA99 closed forms.  First order: A_1, N_1 from rows 1-2, W_1
affine in V_1 from the rank-one fluid rows, V_1 = root of the exact quadratic
l.F_1 = 0 with V_1 > 0 (EC branch).  Orders >= 2: ``recursion.solve_recursion``
with level equations (F1_{n-1}, F2_{n-1}, F_fluid_{n-1}, l.F_n).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from flint import arb, fmpq

from . import recursion
from .arbseries import Series, precision, to_arb
from .polysys import _PolyEvaluator, eval_arb
from .systems import sonic_constraint_poly, sonic_constraint_propagation, sonic_system


def _lead(x):
    return x[0] if isinstance(x, Series) else x


def _peval(poly, u0, dm):
    if dm == 1:
        return eval_arb(poly, [arb(0)] + list(u0))
    return _PolyEvaluator([Series([0], dm)] + list(u0))(poly)


def sonic_data(V0):
    """Zeroth-order data (A0, N0, W0, V0); V0 an arb or a delta-Series (V0 = base + delta)."""
    s3 = arb(3).sqrt()
    N0 = (V0 * (-1) + s3) / (V0 * (-s3) + 1)
    A0 = (V0 * (2 * s3) - V0**2 * 3 + 7) / ((V0**2 * (-1) + 1) * 4)
    W0 = (V0 * (-2 / s3) - V0**2 + 1) / (V0**2 * (-1) + 1) * arb(fmpq(3, 8))
    return [A0, N0, W0, V0]


def fluid_null_vector(sys, u0, dm=1):
    """l = (P4_W, -P3_W)(u0)/scale (max |l| = 1 at the base point) and the fluid row used."""
    p3w, p4w = _peval(sys.P[2][2], u0, dm), _peval(sys.P[3][2], u0, dm)
    p3v, p4v = _peval(sys.P[2][3], u0, dm), _peval(sys.P[3][3], u0, dm)
    scale = _lead(p3w).abs_upper().max(_lead(p4w).abs_upper())
    ell = [p4w / scale, -(p3w / scale)]
    row = 2 if (abs(_lead(p3w)) + abs(_lead(p3v))).mid() >= (abs(_lead(p4w)) + abs(_lead(p4v))).mid() else 3
    return ell, row


def level_equations(ell, fluid_row):
    return recursion.LevelEquations([[(arb(1), 0, -1)], [(arb(1), 1, -1)], [(arb(1), fluid_row, -1)],
                                     [(ell[0], 2, 0), (ell[1], 3, 0)]])


def first_order(sys, u0, ell, fluid_row, dm=1, branch="EC"):
    """(u1, info): A1, N1 explicit; W1 = alpha + beta V1; V1 root of the exact quadratic."""
    S0 = _peval(sys.P[0][0], u0, dm)
    A1 = _peval(sys.Q[0], u0, dm) / S0
    N1 = _peval(sys.Q[1], u0, dm)
    r = fluid_row
    pw, pv, q0 = (_peval(sys.P[r][2], u0, dm), _peval(sys.P[r][3], u0, dm), _peval(sys.Q[r], u0, dm))
    alpha, beta = q0 / pw, -(pv / pw)

    def q_of(V1):
        u1 = [A1, N1, alpha + beta * V1, V1]
        u = [recursion.make_series([u0[i], u1[i]], 2, dm) for i in range(4)]
        Fr = sys.residual(u, cap=2, t=recursion.var_series(2, dm))
        return Fr[2][1] * ell[0] + Fr[3][1] * ell[1]

    one = arb(1) if dm == 1 else Series([1], dm)
    qm, q0v, qp = q_of(one * (-1)), q_of(one * 0), q_of(one)
    c0 = q0v
    c2 = (qp + qm) / 2 - q0v
    c1 = (qp - qm) / 2
    disc = c1 * c1 - c2 * c0 * 4
    if not _lead(disc) > 0:
        raise ValueError(f"discriminant not provably positive: {_lead(disc)}")
    sq = disc.sqrt()
    roots = sorted([(-c1 - sq) / (c2 * 2), (-c1 + sq) / (c2 * 2)], key=lambda z: _lead(z).mid())
    V1 = roots[1] if branch == "EC" else roots[0]
    if branch == "EC" and not _lead(V1) > 0:
        raise ValueError(f"EC branch root not provably positive: {_lead(V1)}")
    return [A1, N1, alpha + beta * V1, V1], dict(c2=c2, c1=c1, c0=c0, disc=disc, roots=roots)


def _run(sys, V0, K, dm, branch, check):
    """One recursion run; V0 an arb (dm = 1) or delta-Series base + delta (dm >= 2)."""
    u0 = sonic_data(V0)
    ell, row = fluid_null_vector(sys, u0, dm)
    eqs = level_equations(ell, row)
    info = {}
    if check:
        info["det_fluid_block"] = (_peval(sys.P[2][2], u0, dm) * _peval(sys.P[3][3], u0, dm)
                                   - _peval(sys.P[2][3], u0, dm) * _peval(sys.P[3][2], u0, dm))
        info["ell_dot_Q0"] = _peval(sys.Q[2], u0, dm) * ell[0] + _peval(sys.Q[3], u0, dm) * ell[1]
        info["constraint"] = _peval(sonic_constraint_poly(sys.ctx), u0, dm)
        for k, v in info.items():
            if not recursion.contains_zero(v):
                raise ValueError(f"zeroth-order identity {k} violated: {v}")
        if not _lead(ell[1 if row == 2 else 0]) != 0:
            raise ValueError("null vector component not provably nonzero (row equivalence)")
    u1, qinfo = first_order(sys, u0, ell, row, dm, branch)
    info.update(qinfo)
    coefs, (D, E) = recursion.solve_recursion(sys, eqs, [u0, u1], K, dm=dm)
    if check:
        res = recursion.level_residuals(sys, eqs, coefs, n_from=2, dm=dm)
        bad = [n + 2 for n, r in enumerate(res) if not all(recursion.contains_zero(v) for v in r)]
        if bad:
            raise ValueError(f"level residuals do not contain 0 at orders {bad}")
    return dict(coefs=coefs, D=D, E=E, eqs=eqs, ell=ell, row=row, info=info)


@dataclass
class SonicExpansion:
    centre: arb
    width: float
    m: int
    K: int
    sys: object
    eqs: object
    point: dict                     # run at the exact centre (delta-Series coefficients if m >= 0)
    interval: dict = None           # run with the interval as base point (if width > 0)
    balls: list = None              # enclosures of u_n over the interval (list over n of 4 arb)
    rem: list = None                # Lagrange remainder radii (list over n of 4 float)
    D: object = None                # arb_mat over the interval (or at the centre)
    E: object = None
    cert: object = None
    info: dict = field(default_factory=dict)

    @property
    def coefs(self):
        return self.balls

    def floats(self):
        return np.array([[float(c) for c in cn] for cn in self.balls]).T       # (4, K+1)

    def radii(self):
        return np.array([[float(c.rad()) for c in cn] for cn in self.balls]).T

    def series(self):
        return recursion.series_from_coefs(self.balls, 4, extra_zero=False)

    def delta_poly(self, n, i):
        """Polynomial part (degree <= m) of the Taylor model of u_{n,i} in delta = V0 - c."""
        c = self.point["coefs"][n][i]
        return Series(c.coeffs(self.m + 1), self.m + 1) if isinstance(c, Series) else Series([c], 1)

    def eval(self, x, with_tail=True):
        """Balls for (A, N, W, V)(x) over the V0 interval; tail bound added if certified."""
        x = to_arb(x)
        vals = [s(x) for s in self.series()]
        if with_tail:
            if self.cert is None:
                raise ValueError("expansion not certified; call certify() or use with_tail=False")
            tb = self.cert.tail_bound(abs(x).abs_upper())
            vals = [v + arb(0, tb) for v in vals]
        return vals

    def certify(self, nu=None, **kw):
        from .tailbound import certify_tail
        with precision(kw.pop("prec", 256)):
            self.cert = certify_tail(self.sys, self.eqs, self.balls, self.D, self.E, nu=nu, **kw)
        return self.cert


def constraint_exponent(sys, coefs):
    """gamma = Lambda(u_0) / (S_0 Delta~_1): C~(u(x)) = 0 identically iff gamma is no positive integer
    (given C~(u_0) = 0), since S Delta~ C~' = Lambda C~ along the flow and Delta~ has a simple zero."""
    lam, Delta = sonic_constraint_propagation(sys)
    u = recursion.series_from_coefs(coefs[:2], 4, cap=2, extra_zero=False)
    ev = _PolyEvaluator([Series.var(2)] + u)
    D1 = ev(Delta)[1]
    S0 = ev(sys.P[0][0])[0]
    lam0 = ev(lam)[0]
    return lam0 / (S0 * D1), D1


def sonic_expansion(V0, K=50, width=None, m=None, prec=256, branch="EC", check=True):
    """Certified sonic Taylor coefficients for V0 in [c-w, c+w] (c = ``V0``, w = ``width``)."""
    with precision(prec):
        c = to_arb(V0)
        sys = sonic_system()
        if width is None or width == 0:
            pt = _run(sys, c, K, 1, branch, check)
            ex = SonicExpansion(c, 0.0, -1, K, sys, pt["eqs"], pt, balls=pt["coefs"],
                                D=pt["D"][0], E=pt["E"][0], info=dict(pt["info"]))
        else:
            m = 5 if m is None else m
            dm = m + 2
            pt = _run(sys, Series([c, 1], dm), K, dm, branch, check)
            X = c + arb(0, width)
            iv = _run(sys, Series([X, 1], dm), K, dm, branch, check)
            dw = arb(0, width)
            wpow = arb(width) ** (m + 1)
            balls, rem = [], []
            for n in range(K + 1):
                bn, rn = [], []
                for i in range(4):
                    top = iv["coefs"][n][i][m + 1].abs_upper() * wpow
                    poly = Series(pt["coefs"][n][i].coeffs(m + 1), m + 1)
                    bn.append(poly(dw) + arb(0, top))
                    rn.append(float(top))
                balls.append(bn)
                rem.append(rn)
            # level equations with l enclosed over the interval (delta^0 part of the interval run)
            eqs_X = level_equations([iv["ell"][0][0], iv["ell"][1][0]], iv["row"])
            ex = SonicExpansion(c, float(width), m, K, sys, eqs_X, pt, iv, balls, rem,
                                D=iv["D"][0], E=iv["E"][0], info=dict(iv["info"]))
        if check:
            gam, D1 = constraint_exponent(sys, ex.balls)
            ex.info["constraint_exponent"] = gam
            ex.info["Delta1"] = D1
            if not (D1 != 0 and (gam < 1 or not gam.contains_integer())):
                raise ValueError(f"constraint propagation not certified: gamma={gam}, Delta1={D1}")
        return ex
