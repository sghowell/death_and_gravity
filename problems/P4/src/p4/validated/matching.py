"""A3: the matching map F(V0, a, mu) and its Krawczyk verification.

Unknowns: V0 (sonic velocity), a (x-translation of the regular-centre family),
mu = n_inf^2 w_inf (its scale invariant); the centre family member is
    C(a, mu)(x) = (e^{-a} n~(t'), e^{2a} w~(t'), e^{a} v~(t')),   t' = e^{a + x},
with (n~, w~, v~) the normalised family of ``centre`` (n_inf = 1).  With
Phi(x; V0) = (n, w, v)(x) the sonic-side solution (A1 Taylor model at x0, then
``tmint`` to x_c) the matching map is
    F(V0, a, mu) = Phi(x_c; V0) - C(a, mu)(x_c)  in R^3,
    F' = [ dPhi/dV0 | -dC/da | -dC/dmu ],
    dC/da = (e^{-a}(theta n~ - n~), e^{2a}(theta w~ + 2 w~), e^{a}(theta v~ + v~)),  theta = t' d/dt',
    dC/dmu = (e^{-a}, e^{2a}, e^{a}) * d(n~, w~, v~)/dmu.
A zero of F is a solution of the reduced CSS system that is the analytic EC-branch
sonic solution for x >= x_c and the regular-centre solution for x <= x_c: the
Evans-Coleman solution.  Krawczyk: K(X) = m - Y F(m) + (I - Y F'(X))(X - m) subset int X
proves existence and uniqueness of the zero in the box X.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from flint import arb, arb_mat

from . import centre, lintail, sonic, tmint, variational as va
from .arbseries import Series, abs_upper, precision, to_arb


# ----------------------------------------------------------------------------
# sonic side: Taylor model in delta = V0 - c at x0, with derivative and tails
# ----------------------------------------------------------------------------
def sonic_initial_state(c, w, x0=-0.05, K=40, m=5):
    """TMState at x0 for V0 in [c-w, c+w] from the certified sonic expansion."""
    ex = sonic.sonic_expansion(str(c) if not isinstance(c, arb) else c, K=K, width=w, m=m)
    cert_u = ex.certify()
    if not cert_u.ok or not (cert_u.nu > abs(x0)):
        raise RuntimeError(f"sonic certificate does not cover x0={x0}: {cert_u}")
    sys8 = va.augment(ex.sys)
    by = va.derivative_balls(ex.point["coefs"], ex.interval["coefs"], m, w)
    ell = ex.interval["ell"]
    eqs8 = va.augmented_level_equations(ex.eqs, 4, extra={3: [(ell[0][1], 2, 0), (ell[1][1], 3, 0)]})
    coefs8 = va.augmented_coefs(ex.balls, by)
    bad = va.check_level_residuals(sys8, eqs8, coefs8, 2)
    if bad:
        raise RuntimeError(f"augmented sonic level residuals fail at orders {bad}")
    cert_y = lintail.linear_tail_certificate(sys8, eqs8, coefs8, 4, cert_u)
    x0a = arb(x0)
    # delta-polynomials of (A, N, W, V)(x0): a_k = sum_n u_{n,k} x0^n  (Series in delta, cap m+2)
    dm = m + 2
    polys = []
    for i in range(4):
        acc = Series([arb(0)], dm)
        for n in reversed(range(K + 1)):
            acc = acc * x0a + ex.point["coefs"][n][i]
        polys.append(acc.coeffs(dm))
    rem_top = [sum((abs_upper(ex.interval["coefs"][n][i][m + 1]) * abs(x0a) ** n for n in range(K + 1)), arb(0))
               for i in range(4)]
    scale = [None, x0a.exp(), (-2 * x0a).exp(), (-x0a).exp()]          # N e^x, W e^{-2x}, V e^{-x}
    tail_u, tail_y = cert_u.tail_bound(abs(x0a)), cert_y.tail_bound(abs(x0a))
    a0, a1, a2, Ru, Ry = [], [], [], [], []
    for j, i in enumerate((1, 2, 3)):
        s = scale[i]
        p = [pk * s for pk in polys[i]]
        a0.append(arb(p[0].mid()))
        a1.append(arb(p[1].mid()))
        a2.append(p[2] * 2)                                   # dy/ddelta at delta = 0 (ball)
        # remainders of the degree-1 truncations of u and of y (the delta^1 term of y is the segment)
        ru = p[0].rad() + p[1].rad() * w + sum((abs_upper(p[k]) * w**k for k in range(2, m + 1)), arb(0))
        ru += rem_top[i] * s * w ** (m + 1) + tail_u * s
        ry = p[1].rad() + abs_upper(p[2].rad()) * 2 * w
        ry += sum((abs_upper(p[k]) * k * w ** (k - 1) for k in range(3, m + 1)), arb(0))
        ry += rem_top[i] * s * (m + 1) * w**m + tail_y * s
        Ru.append(float(abs_upper(ru)))
        Ry.append(float(abs_upper(ry)))
    T = x0a.exp()
    a0.append(arb(T.mid()))
    Ru.append(float(T.rad()))
    # point set: tails and rounding only (delta = 0 exactly); interval set: the segment
    # {(a1, 2 a2) delta} plus the transverse remainder box
    Rp = [float(abs_upper(tail_u * scale[i] + polys[i][0].rad() * scale[i])) for i in (1, 2, 3)]
    Rp += [float(T.rad())]
    Rp += [float(abs_upper(tail_y * scale[i] + polys[i][1].rad() * scale[i])) for i in (1, 2, 3)]
    tau = a1 + [arb(0)] + [arb(c.mid()) for c in a2]
    iv = tmint.Remainder.from_segment(tau, w, Ru + Ry)
    st = tmint.State(x0a, a0 + a1, w, tmint.Remainder(7, Rp), iv)
    st.info = dict(cert_u=cert_u, cert_y=cert_y, expansion=ex)
    return st


# ----------------------------------------------------------------------------
# centre side
# ----------------------------------------------------------------------------
@dataclass
class CentreSide:
    ex: object
    cert_u: object
    cert_y: object
    balls_y: list

    @classmethod
    def build(cls, mu_c, w_mu, K=30, m=5):
        ex = centre.centre_expansion(mu_c, nhat=1, K=K, width=w_mu, m=m)
        cert_u = ex.certify()
        if not cert_u.ok:
            raise RuntimeError(f"centre certificate failed: {cert_u}")
        sys6 = va.augment(ex.sys)
        by = va.derivative_balls(ex.point["coefs"], ex.interval["coefs"], m, w_mu)
        eqs6 = va.augmented_level_equations(ex.eqs, 3)
        coefs6 = va.augmented_coefs(ex.balls, by)
        bad = va.check_level_residuals(sys6, eqs6, coefs6, 0)
        if bad:
            raise RuntimeError(f"augmented centre level residuals fail at orders {bad}")
        cert_y = lintail.linear_tail_certificate(sys6, eqs6, coefs6, 3, cert_u)
        return cls(ex, cert_u, cert_y, by)

    def eval(self, a, x, point=False):
        """(C, dC/da, dC/dmu) at x for a (arb ball) and mu in the certified interval
        (or at the centre mu_c if ``point``)."""
        a = to_arb(a)
        tp = (a + to_arb(x)).exp()
        r = tp.abs_upper()
        if not (r < self.cert_u.nu):
            raise ValueError("t' = e^{a+x} not inside the certified centre radius")
        coefs = self.ex.balls if not point else [[c[0] for c in cn] for cn in self.ex.point["coefs"]]
        ycoefs = self.balls_y if not point else va.derivative_point(self.ex.point["coefs"])
        K = self.ex.K
        vals, thet, dmu = [], [], []
        for i in range(3):
            s = arb(0)
            th = arb(0)
            dy = arb(0)
            for k in reversed(range(K + 1)):
                s = s * tp + coefs[k][i]
                th = th * tp + coefs[k][i] * k
                dy = dy * tp + ycoefs[k][i]
            vals.append(s + arb(0, self.cert_u.tail_bound(r)))
            thet.append(th + arb(0, r * self.cert_u.deriv_tail_bound(r)))
            dmu.append(dy + arb(0, self.cert_y.tail_bound(r)))
        f = [(-a).exp(), (2 * a).exp(), a.exp()]
        C = [f[i] * vals[i] for i in range(3)]
        dCda = [f[0] * (thet[0] - vals[0]), f[1] * (thet[1] + 2 * vals[1]), f[2] * (thet[2] + vals[2])]
        dCdmu = [f[i] * dmu[i] for i in range(3)]
        return C, dCda, dCdmu


# ----------------------------------------------------------------------------
# the matching map, Newton refinement and Krawczyk
# ----------------------------------------------------------------------------
def _fmat(rows):
    return np.array([[float(x.mid()) for x in row] for row in rows])


def krawczyk(st, cs, box_a, m):
    """Krawczyk test on X = [c-w, c+w] x box_a x [mu_c - w_mu, mu_c + w_mu]; ``m`` = (c, a_c, mu_c).

    Returns (ok, K, details)."""
    c, a_c, mu_c = m
    w = st.w
    a_ball = to_arb(a_c) + arb(0, box_a)
    up = st.u_point()
    yi = st.y_interval()
    C0, dCda0, dCdmu0 = cs.eval(to_arb(a_c), st.x, point=True)
    Cb, dCda, dCdmu = cs.eval(a_ball, st.x, point=False)
    Fm = [up[i] - C0[i] for i in range(3)]
    FpX = [[yi[i], -dCda[i], -dCdmu[i]] for i in range(3)]
    yp = st.y_point()
    Fpm = [[yp[i], -dCda0[i], -dCdmu0[i]] for i in range(3)]
    Y = np.linalg.inv(_fmat(Fpm))
    Ya = arb_mat([[arb(float(Y[i, j])) for j in range(3)] for i in range(3)])
    FpXa = arb_mat(FpX)
    I = arb_mat(3, 3, [arb(1) if i == j else arb(0) for i in range(3) for j in range(3)])
    M = I - Ya * FpXa
    Xm = arb_mat([[arb(0, w)], [arb(0, box_a)], [arb(0, cs.ex.width)]])
    Kd = -(Ya * arb_mat([[f] for f in Fm])) + M * Xm             # K(X) - m
    halfw = [w, box_a, cs.ex.width]
    ok = all(abs_upper(Kd[i, 0]) < halfw[i] for i in range(3))
    Kbox = [Kd[i, 0] for i in range(3)]
    return ok, Kbox, dict(Fm=Fm, FpX=FpX, Fpm=Fpm, M=M, newton=[-(Ya * arb_mat([[f] for f in Fm]))[i, 0] for i in range(3)])


# ----------------------------------------------------------------------------
# A4: sign certificate for V = v e^x on (-inf, 0]
# ----------------------------------------------------------------------------
def sonic_v_positive(ex, x0, nsonic=10):
    """The certified A1 series (with tail) has V > 0 on [x0, 0]: ball evaluation on nsonic sub-intervals."""
    ok = True
    for j in range(nsonic):
        xb = arb(x0 * (j + 0.5) / nsonic, abs(x0) * 0.5 / nsonic)
        ok = ok and bool(ex.eval(xb)[3] > 0)
    return ok


def mid_v_one_zero(log):
    """From the per-step ranges of v and v' over the certified step tubes (``tmint`` log): exactly
    one step's v-range contains 0, v' is bounded away from 0 there, v > 0 on the steps before
    (larger x) and v < 0 on the steps after.  Returns (ok, zero step (x_end, h) or None)."""
    zero = [l for l in log if l["J"]["v_range"][0] <= 0 <= l["J"]["v_range"][1]]
    if len(zero) != 1:
        return False, None
    z = zero[0]
    ok = (z["J"]["dv_range"][0] > 0 or z["J"]["dv_range"][1] < 0)
    ok = ok and all(l["J"]["v_range"][0] > 0 for l in log if l["x"] > z["x"])
    ok = ok and all(l["J"]["v_range"][1] < 0 for l in log if l["x"] < z["x"])
    return ok, (z["x"], z["h"])


def centre_v_negative(cs, a_c, x_c):
    """v~(t') <= v~_0 + sum_{k>=1}|v~_k| t_c^k + tail < 0 for all t' <= t_c = e^{a + x_c} (A2 series)."""
    tc = (to_arb(a_c) + to_arb(x_c)).exp().abs_upper()
    K = cs.ex.K
    dev = sum((abs_upper(cs.ex.balls[k][2]) * tc**k for k in range(1, K + 1)), arb(0))
    dev += cs.cert_u.tail_bound(tc)
    return bool(cs.ex.balls[0][2] + dev < 0), float(dev)


def sign_certificate_v(st, cs, a_c, nsonic=10):
    """Certify that V = v e^x has exactly one zero on (-inf, 0] along the EC solution:
    (i) sonic side on [x0, 0], (ii) the integrated range [x_c, x0], (iii) the centre side x <= x_c.
    Returns (ok, details)."""
    x0 = float(st.log[0]["x"] + st.log[0]["h"]) if st.log else float(st.x)
    sonic_ok = sonic_v_positive(st.info["expansion"], x0, nsonic)
    mid_ok, zero_step = mid_v_one_zero(st.log)
    centre_ok, dev = centre_v_negative(cs, a_c, st.x)
    return sonic_ok and mid_ok and centre_ok, dict(sonic_ok=sonic_ok, mid_ok=mid_ok, zero_step=zero_step,
                                                   centre_ok=centre_ok, centre_dev=dev)
