"""Validated (Lohner-type) integrator for the augmented 7D system z = (u, y),
u = (n, w, v, T), y = d(n, w, v)/dV0  (dT/dV0 = 0 identically and is dropped).

Two solution sets are carried around one reference trajectory m(x) (the
midpoint trajectory, V0 = c):
    point set     z(x; c)    in  m + {A_p r_p}          (tails and rounding only),
    interval set  z(x; V0)   in  m + {A r},  V0 in [c-w, c+w],
each a Lohner set (matrix times a box centred at 0; QR re-orthogonalisation
against wrapping).  A step x -> x - h uses only *point* Taylor data at m, which
keeps every ball tight (interval Taylor coefficients of this system blow up by
~100x per order and are never formed):
 1. Taylor coefficients z_i of the augmented solution through m (order K, tight);
    Phi~(m) = sum z_i (-h)^i + tail, with the x-tail of u from the Banach certificate
    (``tailbound``) and of y from the affine-contraction argument (``lintail``);
    h <= hfrac * nu;
 2. J, an interval enclosure of the Jacobian of the step map over the hull of the
    larger of the two sets: the point fundamental-matrix series Y_K(s) plus a
    Groenwall bound of its defect along the certified tube  sum z_i s^i + tail + rho_R
    that contains every solution started in either hull (``jacobian_step``);
 3. mean-value form: Phi~(m + d) in Phi~(m) + J d for both sets; the tails and the
    rounding radii of the new midpoint go into the sets.
Everything is ball arithmetic: every bound, tail and radius on the certified path
is an exact arb upper bound (floats appear only in control flow, rounded up, and
in logs); the reference trajectory itself is only the midpoint of the point set
(its true value is enclosed by the point set).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from flint import arb, arb_mat

from . import lintail, recursion
from . import shootsys as ss
from .arbseries import Series, abs_upper
from .tailbound import certify_tail, norm_inf

DU, DY = 4, 3
DZ = DU + DY
BLOCKS = [(0, 1), (1, 3), (3, 4), (4, 5), (5, 7)]     # block lower-triangular structure of P~


def _colvec(v):
    return arb_mat([[x] for x in v])


def _rad_ball(x):
    return arb(0, abs_upper(x))


def _float_up(x):
    """A double f >= x (x an arb), i.e. float(x) rounded toward +inf: for control flow only; the
    certificate itself keeps the arb value."""
    f = float(x.upper())
    while not (arb(f) >= x):
        f = math.nextafter(f, math.inf)
    return f


def _amax(vals):
    """Rigorous maximum (as an exact arb upper bound) of nonnegative arb values."""
    out = arb(0)
    for v in vals:
        out = out.max(v)
    return out


class Remainder:
    """Lohner set {A r : |r_j| <= rho_j} in R^m."""

    def __init__(self, m, radii=None):
        self.m = m
        self.A = arb_mat(m, m, [arb(1) if i == j else arb(0) for i in range(m) for j in range(m)])
        self.r = [arb(0, 0 if radii is None else radii[j]) for j in range(m)]

    @classmethod
    def from_segment(cls, tau, w, box):
        """Lohner set for the segment {tau delta : |delta| <= w} + box (``box``: radii), with the
        tangent tau (list of arb) as the first basis direction (QR of [tau | I]); r is computed
        rigorously as A^{-1} (tau [-w, w] + box)."""
        m = len(tau)
        R = cls(m)
        M = np.zeros((m, m))
        M[:, 0] = [float(t) for t in tau]
        M[:, 1:] = np.eye(m)[:, :m - 1] if abs(M[0, 0]) < 0.5 * np.linalg.norm(M[:, 0]) else np.eye(m)[:, 1:]
        Q, _ = np.linalg.qr(M)
        R.A = arb_mat([[arb(float(Q[i, j])) for j in range(m)] for i in range(m)])
        c = sum((tau[i] * R.A[i, 0] for i in range(m)), arb(0))       # tau . q_1 = +-|tau|
        eps = [tau[i] - c * R.A[i, 0] for i in range(m)]               # tau = c q_1 + eps (eps tiny)
        trans = [arb(0, abs_upper(e) * w) + arb(0, b) for e, b in zip(eps, box)]
        v = R.A.inv() * _colvec(trans)
        R.r = [_rad_ball(v[i, 0] + (arb(0, abs_upper(c) * w) if i == 0 else 0)) for i in range(m)]
        return R

    def hull(self):
        v = self.A * _colvec(self.r)
        return [_rad_ball(v[i, 0]) for i in range(self.m)]

    def radii(self):
        """Exact arb radii of the hull (the certified widths; ``widths`` is their float version, logs only)."""
        return [b.rad() for b in self.hull()]

    def widths(self):
        return [float(b.rad()) for b in self.hull()]

    def propagate(self, Jmid, extra, pert):
        """New set containing Jmid (A r) + pert + extra (``pert``, ``extra``: balls centred at 0;
        Jmid a matrix with rounding radii only)."""
        C = Jmid * self.A
        Cm = np.array([[float(C[i, j].mid()) for j in range(self.m)] for i in range(self.m)])
        size = np.array([np.linalg.norm(Cm[:, j]) * float(self.r[j].rad()) for j in range(self.m)])
        perm = np.argsort(-size)
        Q, _ = np.linalg.qr(Cm[:, perm])
        Anew = arb_mat([[arb(float(Q[i, j])) for j in range(self.m)] for i in range(self.m)])
        Ainv = Anew.inv()
        v = (Ainv * C) * _colvec(self.r) + Ainv * _colvec([p + e for p, e in zip(pert, extra)])
        self.A = Anew
        self.r = [_rad_ball(v[i, 0]) for i in range(self.m)]


@dataclass
class State:
    x: arb
    m: list                        # 7 arb (exact midpoints): (n, w, v, T, dn, dw, dv) at V0 = c
    w: float                       # half-width of the V0 interval
    pt: Remainder                  # point set (V0 = c)
    iv: Remainder                  # interval set (all V0 in [c-w, c+w])
    log: list = field(default_factory=list)

    def u_point(self):
        H = self.pt.hull()
        return [self.m[i] + H[i] for i in range(DU)]

    def y_point(self):
        H = self.pt.hull()
        return [self.m[DU + i] + H[DU + i] for i in range(DY)]

    def u_interval(self):
        H = self.iv.hull()
        return [self.m[i] + H[i] for i in range(DU)]

    def y_interval(self):
        """Enclosure of d(n, w, v)/dV0 at x for every V0 in [c-w, c+w]."""
        H = self.iv.hull()
        return [self.m[DU + i] + H[DU + i] for i in range(DY)]


class Integrator:
    def __init__(self, K=28, hmax=0.02, hfrac=0.75, nsub=8, verbose=False):
        from .variational import augment
        self.sys4 = ss.shoot_system()
        self.sys7 = augment(self.sys4, skip=(3,))
        self.eqs4 = ss.regular_level_equations(DU)
        self.eqs7 = ss.regular_level_equations(DZ)
        self.K, self.hmax, self.hfrac, self.nsub, self.verbose = K, hmax, hfrac, nsub, verbose
        self.hess = ss.Hessian(self.sys7)

    def jacobian_step(self, co, tails, h, hull_radii):
        """Interval enclosure J of d Phi~_{-h}/dz over the hull (weighted radius rk) of the sets (7x7).

        Y_K(s) = sum_{i<=K} Y_i s^i is the (tight) fundamental-matrix series through m and
        z_K(s) = sum_{i<=K} z_i s^i the truncated reference solution (x-tail <= ``tails``).  Every
        solution z(s) started within the weighted radius rk of m stays, for |s| <= h, in the tube
            Z(s) = z_K(s) + [tail + rho_R],   rho_R > rk e^{L h}   (Groenwall, L >= sup ||Df~|| on the tube;
        fixed point in rho_R), and its fundamental matrix solves Y' = A(s) Y, A(s) = Df~(z(s)).  With
        A_0(s) = Df~(z_K(s)),
            (Y - Y_K)' = A (Y - Y_K) - E,    E := Y_K' - A Y_K = (Y_K' - A_0 Y_K) + (A_0 - A) Y_K,
        so  ||Y(-h; z) - Y_K(-h)|| <= h e^{L h} sup_{|s|<=h} ||E(s)||.  The two parts of E:
         (a) Y_K' - A_0 Y_K = P~(z_K)^{-1} D^(s) with the exact defect polynomial
                 D^(s) = P~(z_K) Y_K' - [DQ~(z_K) - Psi(z_K, f~(z_K))] Y_K,  Psi(z, v)_{rl} = sum_i v_i dP~_{ri}/dz_l
             (P~ Df~ = DQ~ - Psi(z, f~(z)) by differentiating P~ f~ = Q~).  The polynomial formed below,
             D^_K, uses Psi(z_K, z_K') (the series derivative) instead of Psi(z_K, f~(z_K)).  Psi is
             linear in v and f~(z_K) - z_K' = -P~(z_K)^{-1} R(s), R(s) := P~(z_K) z_K' - Q~(z_K) the ODE
             residual of the truncation (an exact polynomial, O(s^K)); hence
                 D^ = D^_K - Psi(z_K, P~(z_K)^{-1} R) Y_K,   sup||D^|| <= Dsup + |dP~| Pinv Rsup Ysup,
             with |dP~| := max_r sum_{i,l} |dP~_{ri}/dz_l| Sc_r/(Sc_i Sc_l) on the tube (so that
             ||Psi(z, v)|| <= |dP~| ||v||) and Dsup, Rsup, Ysup the coefficient sums of D^_K, R, Y_K on
             |s| <= h (no s-interval dependency blow-up; truncation only).
         (b) ||(A_0 - A) Y_K|| <= ||D^2 f~|| ||z_K(s) - z(s)|| ||Y_K|| <= H2 (tail + rho_R) Ysup
             (mean value on the convex tube box, which contains z_K(s) and z(s)).
        Altogether
            ||Y(-h; z) - Y_K(-h)|| <= h e^{L h} [ Pinv (Dsup + Pinv |dP~| Rsup Ysup) + H2 (tail + rho_R) Ysup ].
        All norms are weighted (y scaled by lam, a power of two <= 1/max|y|, so every weight ratio is
        exact); every bound is an exact arb upper bound, floats only for control flow (rounded up).
        ``hull_radii``: per-component arb radii of every set to be propagated (point and interval set);
        the tube is built for the larger weighted radius, so the bound is valid for both sets.
        Returns (Jmid, bound, Sc) with ``bound`` an arb."""
        K = self.K
        ymax = max(abs(float(c)) for c in co[0][DU:])
        lam = 1.0 if ymax <= 1.0 else 2.0 ** -math.frexp(ymax)[1]      # power of two in (1/(2 ymax), 1/ymax]
        Sc = [1.0] * DU + [lam] * DY
        hh = arb(h)

        def wnorm(M):
            return norm_inf(arb_mat([[M[r, c] * (Sc[r] / Sc[c]) for c in range(DZ)] for r in range(DZ)]))

        def sup_poly(s):
            """sum_k |s_k| h^k  >=  sup_{|s| <= h} |s(s)|  for a polynomial Series."""
            tot, hp = arb(0), arb(1)
            for k in range(len(s)):
                tot += abs_upper(s[k]) * hp
                hp *= hh
            return tot

        Yt = ss.variational_coefficients(self.sys7, co, K, BLOCKS)
        zs = recursion.series_from_coefs(co, DZ, cap=None, extra_zero=False)
        Pser = self.sys7.P_series(zs)
        dQs, Psis = self.sys7.dQ_series(zs), self.sys7.psi_series(zs)
        Rres = self.sys7.residual(zs)                    # R(s) = P~(z_K) z_K' - Q~(z_K), exact polynomial
        YK = [[Series([Yt[i][r, c] for i in range(K + 1)]) for c in range(DZ)] for r in range(DZ)]
        dYK = [[Series([Yt[i][r, c] * i for i in range(1, K + 1)]) for c in range(DZ)] for r in range(DZ)]
        Dsup, Ysup, Rsup = arb(0), arb(0), arb(0)
        for r in range(DZ):
            rowD, rowY = arb(0), arb(0)
            for c in range(DZ):
                acc = Series([arb(0)])
                for l in range(DZ):
                    if len(Pser[r][l]) > 1 or Pser[r][l][0] != 0:
                        acc = acc + Pser[r][l] * dYK[l][c]
                    acc = acc - (dQs[r][l] - Psis[r][l]) * YK[l][c]
                rowD += sup_poly(acc) * (Sc[r] / Sc[c])
                rowY += sup_poly(YK[r][c]) * (Sc[r] / Sc[c])
            Dsup, Ysup = Dsup.max(rowD), Ysup.max(rowY)
            Rsup = Rsup.max(sup_poly(Rres[r]) * Sc[r])
        tail_w = tails[0].max(tails[1] * lam)             # weighted norm of the x-tail (z_K vs reference)
        Id = [[arb(1) if i == j else arb(0) for j in range(DZ)] for i in range(DZ)]
        rk = _amax(Sc[j] * hr[j] for hr in hull_radii for j in range(DZ))   # weighted radius covering every set
        dP = self.sys7.dP()

        def dP_norm(args):
            best = arb(0)
            for r in range(DZ):
                row = arb(0)
                for i in range(DZ):
                    for l in range(DZ):
                        if not dP[r][i][l].is_zero():
                            row += abs_upper(ss.eval_box(dP[r][i][l], args)) * (Sc[r] / (Sc[i] * Sc[l]))
                best = best.max(row)
            return best

        def tube_pass(rho_R):
            Lmax, Pinv, H2, dPn = arb(0), arb(0), arb(0), arb(0)
            vr, dvr = None, None                         # ranges of v and v' over the step tube (A4)
            for j in range(self.nsub):                   # exact balls covering [-h, 0]
                sj = (arb(-h) * (2 * j + 1)) / (2 * self.nsub) + arb(0, hh / (2 * self.nsub))
                Z = ss.horner_vec(co, sj)
                Z = [Z[i] + arb(0, tails[0 if i < DU else 1] + rho_R / Sc[i]) for i in range(DZ)]
                f, P = ss.rhs_enclosure(self.sys7, Z, BLOCKS)
                vr = Z[2] if vr is None else vr.union(Z[2])
                dvr = f[2] if dvr is None else dvr.union(f[2])
                A = ss.jacobian_enclosure(self.sys7, Z, f, P, BLOCKS)
                Lmax = Lmax.max(wnorm(A))
                Pinv = Pinv.max(wnorm(ss.block_solve(P, Id, BLOCKS)))
                H2 = H2.max(self.hess.norm(Z, BLOCKS, Sc))
                dPn = dPn.max(dP_norm([arb(0)] + Z))
            eLh = (Lmax * hh).exp()
            g = _float_up(eLh)
            E = Pinv * (Dsup + Pinv * dPn * Rsup * Ysup) + H2 * (rho_R + tail_w) * Ysup
            bound = abs_upper(eLh * hh * E)
            info = dict(L=float(Lmax), Pinv=float(Pinv), H2=float(H2), dPn=float(dPn),
                        v_range=(float(vr.lower()), float(vr.upper())),
                        dv_range=(float(dvr.lower()), float(dvr.upper())),
                        v_pos=bool(vr > 0), v_neg=bool(vr < 0), dv_sign=bool(dvr > 0) or bool(dvr < 0))
            return g, bound, info

        growth, validated = 1.0, False
        for _ in range(6):                           # fixed point on the tube enlargement rho_R (Groenwall)
            rho_R = abs_upper(rk * growth)
            g, bound, info = tube_pass(rho_R)
            if g > 4.0:
                raise RuntimeError("L h too large for the tube enlargement; reduce the step")
            if g < growth or rk == 0:                # rk e^{L h} <= rk g < rk growth <= rho_R: invariant tube
                validated = True
                break
            growth = 1.5 * g
        if not validated:
            raise RuntimeError("tube enlargement did not converge; reduce the step")
        # refinement (valid because the tube above is): a solution started within rk deviates from the
        # reference by at most (sup|Y_K| + bound) rk for |s| <= h (mean value with the fundamental
        # matrix), so the smaller tube rho_R2 > that also contains every such solution
        rho_R2 = abs_upper(rk * (Ysup + bound) * 1.01)
        if rho_R2 < rho_R:
            g2, bound2, info2 = tube_pass(rho_R2)
            if bound2 < bound:
                bound, info, rho_R = bound2, info2, rho_R2
        Jmid = ss.horner_mat(Yt, arb(-h))
        self.last = dict(Dsup=float(Dsup), Rsup=float(Rsup), bound=float(bound), rho_R=float(rho_R),
                         Ysup=float(Ysup), rk=float(rk), **info)
        return Jmid, bound, Sc

    def step(self, st, x_end):
        K = self.K
        co = ss.taylor_coefficients(self.sys7, st.m, K, blocks=BLOCKS)
        co4 = [c[:DU] for c in co]
        D4, E4 = recursion.structure_matrices(self.sys4, self.eqs4, co4, dm=1)
        cert_u = certify_tail(self.sys4, self.eqs4, co4, D4, E4)
        if not cert_u.ok:
            raise RuntimeError(f"u tail certificate failed at x={st.x}")
        cert_y = lintail.linear_tail_certificate(self.sys7, self.eqs7, co, DU, cert_u)
        nu = float(cert_u.nu)
        h = min(self.hmax, self.hfrac * nu, float(st.x - x_end))
        hull_radii = [S.radii() for S in (st.pt, st.iv)]          # exact arb radii before the step
        for _ in range(12):
            tails = (cert_u.tail_bound(arb(h)), cert_y.tail_bound(arb(h)))      # exact arb upper bounds
            try:
                J = self.jacobian_step(co, tails, h, hull_radii)
                break
            except (RuntimeError, ZeroDivisionError):
                h *= 0.5
        else:
            raise RuntimeError(f"Jacobian enclosure failed at x={st.x}")
        Jmid, bound, Sc = J
        mn = ss.horner_vec(co, arb(-h))
        extra = [arb(0, mn[i].rad() + tails[0 if i < DU else 1]) for i in range(DZ)]
        for S, d in zip((st.pt, st.iv), hull_radii):
            du = _amax(d[:DU])                            # u-rows of the perturbation see (d_u, d_T) only
            dw = _amax(Sc[j] * d[j] for j in range(DZ))   # weighted norm of the whole deviation
            pert = [arb(0, bound * du) for i in range(DY)]
            pert += [arb(0, bound * d[DU - 1])]           # T-row: (e^{-h} - Y_K(-h)_TT) d_T, |.| <= bound d_T
            pert += [arb(0, bound * dw / Sc[DU + i]) for i in range(DY)]
            S.propagate(Jmid, extra, pert)
        st.m = [arb(c.mid()) for c in mn]
        st.x = st.x - arb(h)
        st.log.append(dict(x=float(st.x), h=h, nu=nu, eps_u=float(cert_u.eps), eps_y=float(cert_y.eps),
                           tails=(float(tails[0]), float(tails[1])), J=dict(self.last), pt=st.pt.widths(),
                           iv=st.iv.widths()))
        if self.verbose:
            p, v = st.pt.widths(), st.iv.widths()
            print(f"x={float(st.x):+.5f} h={h:.4f} nu={nu:.3g} L={self.last['L']:.3g} Jbound={self.last['bound']:.1e} "
                  f"pt<={max(p):.1e} iv_u<={max(v[:DU]):.1e} iv_y<={max(v[DU:]):.1e}")

    def integrate(self, st, x_end):
        while float(st.x) > x_end + 1e-15:
            self.step(st, x_end)
        return st
