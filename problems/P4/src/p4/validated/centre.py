"""Validated regular-centre expansion (Theorem A, item A2) in S1's scaled variables.

Variables n = N e^x, w = W e^{-2x}, v = V e^{-x} (A^ = (A-1) e^{-2x} follows from the
momentum constraint), series in t = e^x:  Y(t) = sum_k Y_k t^k.  The cleared
polynomial system is ``systems.centre_system`` (theta = t d/dt = d/dx):
    P(t, Y) . theta Y = Q(t, Y),   Q(0, Y_0) = 0  <=>  n_0 v_0 = -1/2   (fixed-point set),
and at every order k >= 1 the coefficient Y_k solves  (k P_0 - DQ_0) Y_k = (known),
which is uniquely solvable because P_0^{-1} DQ_0 has spectrum {-3, 0, 0}
(no positive integer): Briot-Bouquet, one analytic solution per fixed point.

Two free parameters (n_inf, w_inf) = (n_0, w_0), v_0 = -1/(2 n_0).  The x-translation
symmetry of the CSS system (x -> x + a) acts as
    n_k -> n_k e^{(k-1)a},  w_k -> w_k e^{(k+2)a},  v_k -> v_k e^{(k+1)a},
so the family is one-parametric modulo translations, with invariant mu = n_inf^2 w_inf.
We compute the *normalised* family n_inf = 1, w_inf = mu as a Taylor model in
delta = mu - c (same bivariate recursion as ``sonic``) and rescale with
e^{a} = 1/n_inf for general parameters (``rescale``).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from flint import arb

from . import recursion
from .arbseries import Series, precision, to_arb
from .systems import centre_Ahat_fraction, centre_system


def fixed_point(nhat, what):
    """Y_0 = (n_inf, w_inf, -1/(2 n_inf)); arguments arb or delta-Series."""
    return [nhat, what, (nhat * 2).inv() * (-1) if isinstance(nhat, Series) else -1 / (2 * nhat)]


def level_equations():
    return recursion.LevelEquations([[(arb(1), 0, 0)], [(arb(1), 1, 0)], [(arb(1), 2, 0)]])


def _run(sys, u0, K, dm, check):
    eqs = level_equations()
    coefs, (D, E) = recursion.solve_recursion(sys, eqs, [u0], K, dm=dm)
    info = {}
    if check:
        res = recursion.level_residuals(sys, eqs, coefs, n_from=0, dm=dm)
        bad = [n for n, r in enumerate(res) if not all(recursion.contains_zero(v) for v in r)]
        if bad:
            raise ValueError(f"level residuals do not contain 0 at orders {bad}")
    return dict(coefs=coefs, D=D, E=E, eqs=eqs, info=info)


@dataclass
class CentreExpansion:
    nhat: arb                       # n_inf
    what: arb                       # w_inf (centre of the mu interval when width > 0 and nhat = 1)
    width: float
    m: int
    K: int
    sys: object
    eqs: object
    point: dict
    interval: dict = None
    balls: list = None              # list over k of [n_k, w_k, v_k] enclosures
    rem: list = None
    D: object = None
    E: object = None
    cert: object = None
    info: dict = field(default_factory=dict)

    @property
    def coefs(self):
        return self.balls

    def floats(self):
        return np.array([[float(c) for c in ck] for ck in self.balls]).T       # (3, K+1)

    def radii(self):
        return np.array([[float(c.rad()) for c in ck] for ck in self.balls]).T

    def series(self):
        return recursion.series_from_coefs(self.balls, 3, extra_zero=False)

    def Ahat_series(self):
        """Truncated series of A^ = 2 w T / S (exact recursion to order K; no tail bound)."""
        num, den = centre_Ahat_fraction(self.sys.ctx)
        u = [s.with_cap(self.K + 1) for s in self.series()]
        t = Series.var(self.K + 1)
        ev = self.sys.evaluator(t, u)
        return ev(num) / ev(den)

    def eval(self, x, with_tail=True):
        """Balls for (n, w, v)(x) = Y(e^x); tail bound added if certified and e^x <= nu."""
        t = to_arb(x).exp()
        vals = [s(t) for s in self.series()]
        if with_tail:
            if self.cert is None:
                raise ValueError("expansion not certified; call certify() or use with_tail=False")
            tb = self.cert.tail_bound(t.abs_upper())
            vals = [v + arb(0, tb) for v in vals]
        return vals

    def certify(self, nu=None, **kw):
        from .tailbound import certify_tail
        with precision(kw.pop("prec", 256)):
            self.cert = certify_tail(self.sys, self.eqs, self.balls, self.D, self.E, nu=nu, **kw)
        return self.cert

    def rescale(self, nhat):
        """Coefficient balls for parameters (n_inf, w_inf) = (nhat, what / nhat^2) from the
        normalised family (requires self.nhat = 1): n_k nhat^{1-k}, w_k nhat^{-(k+2)}, v_k nhat^{-(k+1)}."""
        nhat = to_arb(nhat)
        out = []
        for k, (nk, wk, vk) in enumerate(self.balls):
            out.append([nk * nhat ** (1 - k), wk * nhat ** (-(k + 2)), vk * nhat ** (-(k + 1))])
        return out


def centre_expansion(what, nhat=1, K=30, width=None, m=None, prec=256, check=True):
    """Certified centre coefficients.  ``width`` > 0 gives a Taylor model in delta = mu - what
    for the normalised family (nhat must be 1); otherwise a point computation at (nhat, what)."""
    with precision(prec):
        c = to_arb(what)
        nh = to_arb(nhat)
        sys = centre_system()
        if width is None or width == 0:
            pt = _run(sys, fixed_point(nh, c), K, 1, check)
            return CentreExpansion(nh, c, 0.0, -1, K, sys, pt["eqs"], pt, balls=pt["coefs"],
                                   D=pt["D"][0], E=pt["E"][0], info=dict(pt["info"]))
        if not (nh == 1):
            raise ValueError("parametrised centre family is computed for nhat = 1 (use rescale)")
        m = 5 if m is None else m
        dm = m + 2
        one = Series([1], dm)
        pt = _run(sys, fixed_point(one, Series([c, 1], dm)), K, dm, check)
        X = c + arb(0, width)
        iv = _run(sys, fixed_point(one, Series([X, 1], dm)), K, dm, check)
        dw = arb(0, width)
        wpow = arb(width) ** (m + 1)
        balls, rem = [], []
        for k in range(K + 1):
            bk, rk = [], []
            for i in range(3):
                top = iv["coefs"][k][i][m + 1].abs_upper() * wpow
                poly = Series(pt["coefs"][k][i].coeffs(m + 1), m + 1)
                bk.append(poly(dw) + arb(0, top))
                rk.append(float(top))
            balls.append(bk)
            rem.append(rk)
        return CentreExpansion(nh, c, float(width), m, K, sys, iv["eqs"], pt, iv, balls, rem,
                               D=iv["D"][0], E=iv["E"][0], info=dict(iv["info"]))
