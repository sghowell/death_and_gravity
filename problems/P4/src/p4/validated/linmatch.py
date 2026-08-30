"""Theorem B, Stage 2 (b): the matching function E(kappa) as a ball, with dE/dkappa.

Pipeline for one kappa (point or box kappa_c + [-w, w] + i[-w, w]):
 1. sonic side (Stage 1): (A_p, N_p, W_p, V_p)(x0) from the certified series over the box, scaled
    q~ = (N_p T, W_p/T^2, V_p/T) (3D) or p~ = (A_p/T^2, q~) (4D), T = e^{x0};
 2. validated propagation (``linprop``) through the tube region [x_c, x0] and the centre region
    [x_d, x_c] of the background ``Tube`` (3D reduced or 4D full system);
 3. centre condition: with the 4D regular family r_1, r_2 (``lincentre``) at t_d = e^{x_d} and the
    4D vector p~(x_d) (A^_p from the linearised constraint if propagated in 3D),
        E(kappa) := e^{3 x_d} det [ r_1, r_2, p~(x_d) ] restricted to the rows (A^_p, n~_p, v~_p),
    i.e. det[r_1, r_2, p~, e_W] up to sign: zero iff p~(x_d) lies in the regular plane (e_W is never
    in the constraint-surface solution space since dC~/dW != 0 at the centre), analytic in kappa on
    the whole rectangle (4D data), same zeros as S1's E;
    E_fin(kappa) := e^{x_d} A_p(x_d) is S1's matching function evaluated at x_end = x_d (enclosed
    exactly as a ball; its zeros converge to those of E as x_d -> -infinity).
 4. dE/dkappa: the kappa-derivative systems (``linscaled.augmented``) propagate (q~, dq~/dkappa) with
    initial derivative from the delta-polynomials of Stage 1 (Taylor model in kappa; box remainder
    (m+1) sup|q_{n,m+1}| |delta|^m; derivative of the tail by Cauchy's estimate over the box enlarged
    by rho), and the augmented regular family gives (r_i, dr_i/dkappa) with certified tails.
"""
from __future__ import annotations

import time

from flint import acb, acb_mat, arb

from . import lincentre, linprop, linscaled, linsonic, linstep, sonic
from .arbseries import Series, abs_upper, precision, to_arb
from .linsys import LinSystem, abs_up, kappa_box, plain_from_scaled, to_acb
from .shootsys import eval_box


def box_background(c, w, K=41, m=5):
    """Certified A1 expansion for V0 in [c - w, c + w] packaged as a point-like ``SonicExpansion``
    (coefficient balls over the box, fluid null vector ell as balls from its Taylor model in
    delta = V0 - c with the interval-run remainder) -- the form ``linsonic`` consumes, so the
    linearised sonic series then encloses the solutions for every V0 in the box."""
    bgb = sonic.sonic_expansion(c, K=K, width=w, m=m)
    cert = bgb.certify()
    if not cert.ok:
        raise RuntimeError(f"A1 box certificate failed: {cert}")
    dw, wp = arb(0, w), arb(w) ** (m + 1)
    ell = [Series(pt.coeffs(m + 1))(dw) + arb(0, abs_upper(iv[m + 1]) * wp)
           for pt, iv in zip(bgb.point["ell"], bgb.interval["ell"])]
    return sonic.SonicExpansion(bgb.centre, 0.0, -1, bgb.K, bgb.sys, bgb.eqs, dict(ell=ell, row=bgb.point["row"]),
                                None, bgb.balls, bgb.rem, bgb.D, bgb.E, cert, dict(box=bgb))


def det3(cols, rows=(0, 1, 3)):
    return acb_mat([[cols[c][r] for c in range(3)] for r in rows]).det()


class Matcher:
    def __init__(self, tube, ce, bg, x0=-0.05, nu_c=0.06, K_c=50, K_lin=None, K_sonic=40, prec=256):
        self.tube, self.ce, self.bg, self.x0 = tube, ce, bg, x0
        self.nu_c, self.K_c, self.K_lin, self.K_sonic, self.prec = nu_c, K_c, K_lin, K_sonic, prec
        with precision(prec):
            self.L = LinSystem()
            self.S3, self.S4 = linscaled.reduced_system(self.L), linscaled.full_system(self.L)
            self.S3a, self.S4a = linscaled.augmented(self.S3), linscaled.augmented(self.S4)
            last = tube.steps[-1]
            self.x_d = last.x - arb(last.h)

    # -- sonic side -------------------------------------------------------------------------
    def sonic_data(self, kappa_c, width=0.0, m=5, deriv=False, rho=1e-6):
        ex = linsonic.linear_sonic_expansion(self.bg, kappa_c, width=width, m=m, K=self.K_sonic)
        x0 = to_arb(self.x0)
        cert = ex.certify()
        if cert.ok and not (cert.nu > abs(x0)):                  # nu = |x0| up to rounding: nudge
            cert = ex.certify(nu=float(abs(x0).abs_upper()) * 1.01)
        if not cert.ok or not (cert.nu > abs(x0)):
            raise RuntimeError(f"sonic certificate does not cover x0 = {self.x0}: {cert}")
        vals = ex.eval(x0)
        if not deriv:
            return ex, vals, None
        dw = kappa_box(0, width) if width else acb(0)
        dq = [acb(0)] * 3
        for n in reversed(range(self.K_sonic + 1)):
            cn = [sum((ex.point["coefs"][n][k][i, 0] * k * dw ** (k - 1) for k in range(1, m + 1)), acb(0))
                  for i in range(3)]
            if width:
                top = [abs_up(ex.box["coefs"][n][m + 1][i, 0]) * (m + 1) * abs_up(dw) ** m for i in range(3)]
                cn = [c + acb(arb(0, t), arb(0, t)) for c, t in zip(cn, top)]
            dq = [d * x0 + c for d, c in zip(dq, cn)]
        ex2 = linsonic.linear_sonic_expansion(self.bg, kappa_c, width=width + rho, m=m, K=self.K_sonic, check=False)
        c2 = ex2.certify(nu=float(cert.nu))
        if not c2.ok:
            raise RuntimeError(f"sonic certificate over the enlarged box failed: {c2}")
        tb = abs_upper(c2.tail_bound(abs(x0).abs_upper()) / arb(rho))
        dq = [d + acb(arb(0, tb), arb(0, tb)) for d in dq]
        return ex, vals, dq

    # -- constraint at a point --------------------------------------------------------------
    def _Ap(self, u, kappa, q, dq=None):
        args = [arb(0)] + list(u)
        Cq = [eval_box(self.L.dC[l + 1], args) for l in range(3)]
        S = eval_box(self.L.S, args)
        N = sum((Cq[l] * q[l] for l in range(3)), acb(0))
        den = (kappa - u[0]) * S
        Ap = u[0] * N / den
        if dq is None:
            return Ap, None
        dN = sum((Cq[l] * dq[l] for l in range(3)), acb(0))
        return Ap, u[0] * (dN * (kappa - u[0]) - N) / (den * (kappa - u[0]))

    # -- the matching function ---------------------------------------------------------------
    def E(self, kappa_c, width=0.0, deriv=False, system="3D", m=5, verbose=False):
        t0 = time.time()
        with precision(self.prec):
            kappa = kappa_box(kappa_c, width) if width else to_acb(kappa_c)
            ex, vals, dq0 = self.sonic_data(kappa_c, width, m, deriv)
            T0 = to_arb(self.x0).exp()
            Ap, Np, Wp, Vp = vals
            qt0 = [Np * T0, Wp / (T0 * T0), Vp / T0]
            if system == "3D":
                y0 = qt0 + ([dq0[0] * T0, dq0[1] / (T0 * T0), dq0[2] / T0] if deriv else [])
                S = self.S3a if deriv else self.S3
            else:
                dA0 = self._Ap(self.bg.eval(to_arb(self.x0)), kappa, [Np, Wp, Vp], dq0)[1] if deriv else None
                y0 = [Ap / (T0 * T0)] + qt0
                y0 += [dA0 / (T0 * T0), dq0[0] * T0, dq0[1] / (T0 * T0), dq0[2] / T0] if deriv else []
                S = self.S4a if deriv else self.S4
            st = linstep.LohnerSet(y0)
            _, log = linprop.propagate(self.tube, S, kappa, [st], K=self.K_lin, prec=self.prec, verbose=verbose)
            t_prop = time.time() - t0
            hull = st.hull()
            x_d = self.x_d
            Td = x_d.exp()
            zd = self.ce.eval(x_d)
            ud = plain_from_scaled(x_d, *zd)
            if system == "3D":
                qt = hull[:3]
                q = [qt[0] / Td, qt[1] * Td * Td, qt[2] * Td]
                dqt = hull[3:6] if deriv else None
                dq = [dqt[0] / Td, dqt[1] * Td * Td, dqt[2] * Td] if deriv else None
                Apd, dApd = self._Ap(ud, kappa, q, dq)
                p = [Apd / (Td * Td)] + qt
                dp = [dApd / (Td * Td)] + dqt if deriv else None
            else:
                p, dp = hull[:4], (hull[4:8] if deriv else None)
                Apd = p[0] * Td * Td
            rf = lincentre.RegularFamily(self.S4a if deriv else self.S4, self.ce, kappa, K=self.K_c, prec=self.prec)
            ok, eps, cdet = rf.certify(self.nu_c, prec=self.prec)
            if not ok:
                raise RuntimeError(f"centre certificate failed: {cdet}")
            r1, r2 = rf.eval(x_d, prec=self.prec)
            nrm = (3 * x_d).exp()                                # p~ ~ e^{-3 x}: E normalised to O(1)
            E = det3([r1[:4], r2[:4], p]) * nrm
            dE = (det3([r1[4:8], r2[:4], p]) + det3([r1[:4], r2[4:8], p]) + det3([r1[:4], r2[:4], dp])) * nrm if deriv else None
            return dict(E=E, dE=dE, E_fin=Td * Apd, p=p, kappa=kappa, x_d=x_d, time=time.time() - t0,
                        t_prop=t_prop, width_end=max(float(abs_up(z - z.mid())) for z in hull),
                        centre=cdet, log=log, expansion=ex)

    def krawczyk(self, kappa_c, width, system="3D", m=5):
        """1D Krawczyk test on the box: K = k_c - E(k_c)/dE(k_c) + (1 - dE(box)/dE(k_c)) (box - k_c)
        subset of the interior of the box  =>  exactly one zero of E in the box."""
        with precision(self.prec):
            pt = self.E(kappa_c, 0.0, deriv=True, system=system, m=m)
            bx = self.E(kappa_c, width, deriv=True, system=system, m=m)
            Y = 1 / acb(pt["dE"].real.mid(), pt["dE"].imag.mid())
            Kop = -Y * pt["E"] + (1 - Y * bx["dE"]) * kappa_box(0, width)
            box = kappa_box(0, width)
            ok = bool(box.real.contains_interior(Kop.real)) and bool(box.imag.contains_interior(Kop.imag))
            return ok, Kop, pt, bx
