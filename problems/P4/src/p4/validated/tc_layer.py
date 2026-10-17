"""Theorem C, S4-4: the near-sonic layer engine for the commutator strip
certificate (tau >= TSTAR on the band xi <= 0.036).

The commuted energy condition of notes/s4-commutator.md section 3 is
    S = Re B1 - (tau/c) Im B1 - (f' + rho f)/2 + nu^2/4 < 0,
and on the band it scales as S ell^2/tau^2 = F(ell tau, eps) + o(1).  The
tau-quadratic phase of tc_commutator expands in 1/kappa at fixed ell and loses
the (ell tau)-coupled cancellations.  This module evaluates S with the
tau -> infinity normalization built in: every building block is a ratio of
kappa-polynomials whose coefficients are Laurent Taylor models in ell (tc_tm),
evaluated on boxes in the layer variables

    z = ell kappa,  w = 1/kappa,  zeta = 1/z,

in a "pivot form" that is bounded (and, for denominators, bounded away from
0) uniformly in z in [0, infinity] and ell in [0, ell_max]; the final quantity

    Sigma := S ell^2 / tau^2

is assembled from the normalized blocks with every explicit kappa-power split
exactly into (eps, tau) parts.  Boxes are (ell-strip) x (eps-cell) x
(s = xi tau)-box x (tau)-box, with pseudo-boxes reaching tau = infinity
(w -> 0) and, on the closing piece, ell -> 0.  See notes/s4-4-layer-engine.md.

Point evaluation (``point_*``) is the float-check of the design: it evaluates
S by the direct chain formula (independent of the layer algebra) at a point
(xi, kappa) with arb balls.
"""
from __future__ import annotations

import time

from flint import acb, arb

from .arbseries import precision

TSTAR = 224.0                # the layer engine covers tau >= TSTAR (= 16 b1)
CNUM, CDEN = -1, 32          # c = -1/32 (tc_commutator.CNUM/CDEN)
EPS_MAX = 6.0


def cval():
    return arb(CNUM) / CDEN


# =====================================================================
# point evaluation (design check; not a certificate)
# =====================================================================
def point_pairs(at, ex, xi, prec=320):
    """build_pairs at the single point x = -xi (constant Taylor models), from
    the certified sonic series; returns the pair dict of tc_commutator."""
    from . import tc_commutator as tc
    from .tc_tm import LC, TM, AtomKP, sonic_components
    with precision(prec):
        comps = sonic_components(ex, xi)
        vals = [c.tm.eval(xi, xi) for c in comps]
        cc = tuple(LC(TM.const(v, xi), 0) for v in vals)
        ak = AtomKP(at, cc, sd_shift=0)
        return tc.build_pairs(ak)


def _kpval(kp, slot, k):
    acc = acb(0)
    for x in reversed(kp):
        lc = x[slot]
        v = lc.tm.c[0]
        for cc in lc.tm.c[1:]:
            if not cc.is_zero():
                raise ValueError("point pairs must have constant Taylor models")
        acc = acc * k + v
    return acc


def point_S_direct(pr, eps, tau, prec=320):
    """S at kappa = eps + i tau from the point pairs by the direct chain
    (A0 = nT/dT, B0 = nB/dB, A1 = A0 - (log B0)', B1 = B0 + A0' - A0 (log B0)',
    nu = Im A1 + tau tt, ...), exactly as tc_commutator.direct_S_boxes but at a
    point.  Returns (S, parts dict)."""
    with precision(prec):
        c = cval()
        k = acb(arb(eps), arb(tau))
        tb = arb(tau)
        nT = [_kpval(pr["nT"], s, k) for s in (0, 1)]
        dT = [_kpval(pr["dT"], s, k) for s in (0, 1)]
        nB = [_kpval(pr["nB"], s, k) for s in (0, 1, 2)]
        dB = [_kpval(pr["dB"], s, k) for s in (0, 1, 2)]
        A0 = nT[0] / dT[0]
        A0p = (nT[1] * dT[0] - nT[0] * dT[1]) / (dT[0] * dT[0])
        B0 = nB[0] / dB[0]
        lnB, ldB = nB[1] / nB[0], dB[1] / dB[0]
        lB0 = lnB - ldB
        lB0p = (nB[2] / nB[0] - lnB * lnB) - (dB[2] / dB[0] - ldB * ldB)
        A1 = A0 - lB0
        A1p = A0p - lB0p
        B1 = B0 + A0p - A0 * lB0
        TN = pr["TN"]
        TNt = [TN[-1][s].tm.c[0] for s in (0, 1)]
        dTt = [pr["dT"][-1][s].tm.c[0] for s in (0, 1)]
        tt = TNt[0] / dTt[0]
        ttp = (TNt[1] * dTt[0] - TNt[0] * dTt[1]) / (dTt[0] * dTt[0])
        nu = A1.imag + tb * tt
        nup = A1p.imag + tb * ttp
        rho = A1.real - c * tt
        cf = c + tb * tb / c
        f = cf * tt - (tb / c) * nu
        fp = cf * ttp - (tb / c) * nup
        K = B1.real - (tb / c) * B1.imag - (fp + rho * f) / 2
        S = K + nu * nu / 4
        return S, dict(A1=A1, B1=B1, nu=nu, nup=nup, rho=rho, f=f, fp=fp,
                       K=K, tt=tt, ttp=ttp)


def point_F(at, ex, xi, eps, s, prec=320, pr=None):
    """F := S xi^2 / tau^2 at tau = s / xi (the layer-normalized condition)."""
    if pr is None:
        pr = point_pairs(at, ex, xi, prec)
    tau = s / xi
    S, _ = point_S_direct(pr, eps, tau, prec)
    with precision(prec):
        return S * arb(xi) ** 2 / arb(tau) ** 2


# =====================================================================
# interval helpers
# =====================================================================
INF = float("inf")


def iv(lo, hi):
    """Ball containing the closed real interval [lo, hi] (finite)."""
    return arb(lo).union(arb(hi))


def pw_hull(L, e):
    """Hull of ell^e over the interval L = (lo, hi), lo >= 0."""
    lo, hi = L
    if e == 0:
        return arb(1)
    if e > 0:
        return iv(lo, hi) ** e if lo > 0 else arb(0).union(arb(hi) ** e)
    if not lo > 0:
        raise ArithmeticError("negative ell-power on a strip touching ell = 0")
    return iv(lo, hi) ** e


# =====================================================================
# pivot forms
# =====================================================================
class PForm:
    """Pivot form of a kappa-polynomial P = sum_j a_j kappa^j whose coefficients
    are Laurent TMs a_j = ell^{o_j} TM_j(ell) (tc_tm.LC):

        P = ell^m z^p [ sum_{j >= p} TM_j ell^{e_j} z^{j-p}
                        + sum_{j < p} TM_j ell^{g_j} w^{p-j} ],

    z = ell kappa, w = 1/kappa, m = min_j (o_j - j) (the layer weight),
    p = max{j : o_j = min_k o_k} (the fixed-kappa leading index),
    e_j = o_j - j - m >= 0 and g_j = o_j - m - p.  The bracket Q is a polynomial
    in (z, w) with ell-TM coefficients: bounded on bounded z, and Q -> TM_p(0)
    at (z, w, ell) -> 0 when the pivot lies on the layer line (e_p = 0).  For
    large z the equivalent form Q / z^{d-p} = sum_{j>=p} TM_j ell^{e_j}
    zeta^{d-j} + zeta^{d-p} sum_{j<p} TM_j ell^{g_j} w^{p-j} (zeta = 1/z) is
    used; it tends to TM_d(0) when the top lies on the layer line (e_d = 0).
    Every value is an identity of the true functions: the exponents are
    integers read off the Laurent offsets, nothing is truncated."""

    def __init__(self, kp, name="", pivot=None):
        self.name = name
        kt = [x.trim() for x in kp]
        nz = [j for j, x in enumerate(kt) if not all(c.is_zero() for c in x.tm.c)]
        if not nz:
            raise ArithmeticError("zero polynomial " + name)
        self.d = max(nz)
        self.idx = [j for j in nz if j <= self.d]
        self.tm = {j: kt[j].tm for j in self.idx}
        self.o = {j: kt[j].off for j in self.idx}
        self.m = min(self.o[j] - j for j in self.idx)
        omin = min(self.o.values())
        self.p = max(j for j in self.idx if self.o[j] == omin) if pivot is None else pivot
        self.e = {j: self.o[j] - j - self.m for j in self.idx if j >= self.p}
        self.g = {j: self.o[j] - self.m - self.p for j in self.idx if j < self.p}
        self.pivot_on_layer = self.e.get(self.p, 1) == 0
        self.top_on_layer = self.e.get(self.d, 1) == 0
        self._cache = {}

    def strip(self, L):
        """Coefficient hulls on the ell-strip L = (lo, hi): A[j] (j >= p),
        B[j] (j < p)."""
        key = (L[0], L[1])
        if key not in self._cache:
            A, B = {}, {}
            for j in self.idx:
                h = self.tm[j].eval(L[0], L[1])
                if j >= self.p:
                    A[j] = h * pw_hull(L, self.e[j])
                else:
                    B[j] = h * pw_hull(L, self.g[j])
            self._cache[key] = (A, B)
        return self._cache[key]

    def small(self, L, z, w):
        """Q = P / (ell^m z^p)."""
        A, B = self.strip(L)
        acc = acb(A.get(self.d, arb(0)))
        for j in range(self.d - 1, self.p - 1, -1):
            acc = acc * z + A.get(j, arb(0))
        if B:
            accw = acb(0)
            for j in range(0, self.p):
                accw = accw * w + B.get(j, arb(0))
            acc = acc + accw * w
        return acc

    def large(self, L, zeta, w):
        """Q / z^{d-p} = P / (ell^m z^d)."""
        A, B = self.strip(L)
        acc = acb(A.get(self.p, arb(0)))
        for j in range(self.p + 1, self.d + 1):
            acc = acc * zeta + A.get(j, arb(0))
        if B:
            accw = acb(0)
            for j in range(0, self.p):
                accw = accw * w + B.get(j, arb(0))
            accw = accw * w
            if self.d > self.p:
                accw = accw * zeta ** (self.d - self.p)
            acc = acc + accw
        return acc

    def info(self):
        return dict(name=self.name, d=self.d, p=self.p, m=self.m,
                    offs=[self.o.get(j) for j in range(self.d + 1)],
                    pivot_on_layer=self.pivot_on_layer,
                    top_on_layer=self.top_on_layer)


def ratio_hull(N, D, k, L, lam, z, zeta, w, use_large, laurent):
    """Hull of lam^k * N/D on the box, from the pivot forms N, D.
    Small-z form:  N/D = ell^{mN-mD} z^{pN-pD} QN/QD;
    large-z form:  N/D = ell^{mN-mD} z^{dN-dD} QN^L/QD^L.
    Laurent pieces have lam = ell, so the prefactor is ell^{k+mN-mD}; tube steps
    have no Laurent offsets (mN - mD = dD - dN) and lam = xi (physical)."""
    if use_large:
        q = N.large(L, zeta, w) / D.large(L, zeta, w)
        zp = N.d - D.d
        zb = zeta
        zp = -zp
    else:
        q = N.small(L, z, w) / D.small(L, z, w)
        zp = N.p - D.p
        zb = z
    if zp:
        if zp < 0:
            raise ArithmeticError("negative z-power in ratio %s/%s" % (N.name, D.name))
        q = q * zb ** zp
    dm = N.m - D.m
    if laurent:
        return q * pw_hull(L, k + dm)
    if N.d != D.d or N.p != N.d or D.p != D.d:
        raise ArithmeticError("tube ratio %s/%s: needs balanced pure w-forms"
                              % (N.name, D.name))
    return q * lam ** k if k else q


# =====================================================================
# layer data: the blocks of S as ratios of kappa-polynomials
# =====================================================================
def _lc_zero_like(x):
    from .tc_tm import lc_zero
    return lc_zero(x.tm.R)


class LayerData:
    """The polynomials of the layer form of S on one piece.

    With A0 = nT/dT, B0 = nB/dB (tc_commutator.build_pairs, raised at the
    sonic point by sonic_raise_inputs), tt = TN_top/dT_top and u2 = nB_top/dB_top:

      V  := A0 + kappa tt = nV/dT,  nV = nT + kappa tt dT  (kappa^6 slot = 0),
      V' = nV'/dT - V dT'/dT,      nV' = nT' + kappa (tt dT)'  (kappa^6 slot = 0),
      Y  := (B0 - u2 kappa^2)/kappa = nY/(kappa dB),  nY = nB - u2 kappa^2 dB
                                                      (kappa^9 slot = 0),
      (log B0)' = nB'/nB - dB'/dB,
      (log B0)'' = nB''/nB - (nB'/nB)^2 - dB''/dB + (dB'/dB)^2.

    The dropped top slots vanish identically by the definitions of tt and u2
    (true functions; their TM enclosures are not used there).  At the sonic
    point (laurent=True) the derivative (tt dT_k)' is raised by one Laurent
    order: tt dT_k = TN_top (dT_k/xi)/(dT_top/xi) is analytic at xi = 0 (the
    raised dT_top has a constant term bounded away from 0, checked here), so
    its derivative has no xi^{-1} term (drop_lead checks containment and ve).
    """

    def __init__(self, pr, laurent, xik=None, label=""):
        """laurent=True: sonic pieces, TM variable ell = xi on [0, hi].
        laurent=False: tube steps, TM variable sigma = x_k - x on [0, h] and
        xi = xik + sigma with xik = -x_k (an arb ball)."""
        from .tc_tm import lc_raise
        self.laurent, self.label = laurent, label
        if laurent:
            self.xik = None
            self.xi_of_ell = lambda l: l
            self.ell_of_xi = lambda x: x
        else:
            self.xik = arb(xik)
            xm = float(self.xik.mid())
            if not float(self.xik.rad()) < 1e-14:
                raise ArithmeticError("tube step origin is not a near-point ball")
            self.xi_of_ell = lambda l: xm + l
            self.ell_of_xi = lambda x: x - xm
        s0 = lambda kp, sl: [x[sl].trim() for x in kp]
        nT0, nT1 = s0(pr["nT"], 0), s0(pr["nT"], 1)
        dT0, dT1 = s0(pr["dT"], 0), s0(pr["dT"], 1)
        nB0, nB1, nB2 = (s0(pr["nB"], k) for k in (0, 1, 2))
        dB0, dB1, dB2 = (s0(pr["dB"], k) for k in (0, 1, 2))
        TNt = [pr["TN"][-1][k].trim() for k in (0, 1)]
        dTt = [dT0[-1], dT1[-1]]
        if not (abs(dTt[0].tm.c[0]).lower() > 0):
            raise ArithmeticError("dT top: constant term not bounded away from 0")
        ti = dTt[0].inv()
        tt = (TNt[0] * ti).trim()
        ttp = ((TNt[1] * dTt[0] - TNt[0] * dTt[1]) * (ti * ti)).trim()
        self.tt, self.ttp = tt, ttp
        # V and V'
        nV = []
        nVp = []
        for j in range(len(nT0) - 1):          # kappa^0..kappa^5; kappa^6 == 0
            a = nT0[j]
            ap = nT1[j]
            if j >= 1:
                a = a + tt * dT0[j - 1]
                der = (ttp * dT0[j - 1] + tt * dT1[j - 1]).trim()
                if laurent:
                    val_off = (tt * dT0[j - 1]).trim().off
                    tgt = max(val_off - 1, 0) if val_off >= 0 else val_off - 1
                    if der.off < tgt:
                        der = lc_raise(der, tgt - der.off)
                ap = ap + der
            nV.append(a.trim())
            nVp.append(ap.trim())
        # Y
        if not (abs(dB0[-1].tm.c[0]).lower() > 0):
            raise ArithmeticError("dB top: constant term not bounded away from 0")
        u2 = (nB0[-1] * dB0[-1].inv()).trim()
        self.u2 = u2
        nY = []
        for j in range(len(nB0) - 1):          # kappa^0..kappa^8; kappa^9 == 0
            a = nB0[j]
            if j >= 2:
                a = a - u2 * dB0[j - 2]
            nY.append(a.trim())
        kdB = [_lc_zero_like(dB0[0])] + dB0
        self.forms = dict(
            nV=PForm(nV, "nV"), dT=PForm(dT0, "dT"), nVp=PForm(nVp, "nVp"),
            dT1=PForm(dT1, "dT1"), nY=PForm(nY, "nY"), kdB=PForm(kdB, "kdB"),
            dB=PForm(dB0, "dB"), dB1=PForm(dB1, "dB1"), dB2=PForm(dB2, "dB2"),
            nB=PForm(nB0, "nB"), nB1=PForm(nB1, "nB1"), nB2=PForm(nB2, "nB2"))
        for nm in ("dT", "kdB", "dB", "nB"):
            f = self.forms[nm]
            if not (f.pivot_on_layer and f.top_on_layer):
                raise ArithmeticError("denominator %s: pivot/top off the layer line %s"
                                      % (nm, f.info()))
        self._rcache = {}

    def xi_hull(self, L):
        """Ball containing xi for every TM variable in L = (lo, hi)."""
        if self.laurent:
            return iv(L[0], L[1])
        return self.xik + iv(L[0], L[1])

    def real_hulls(self, L, lam):
        """(t, tp, u2h) = (lam tt, lam^2 tt', lam u2) on the strip."""
        key = (L[0], L[1])
        if key not in self._rcache:
            out = []
            for lc, k in ((self.tt, 1), (self.ttp, 2), (self.u2, 1)):
                h = lc.tm.eval(L[0], L[1])
                if self.laurent:
                    out.append(h * pw_hull(L, k + lc.off))
                else:
                    if lc.off != 0:
                        raise ArithmeticError("tube LC with a Laurent offset")
                    out.append(h * lam ** k)
            self._rcache[key] = tuple(out)
        return self._rcache[key]

    def form_info(self):
        return {k: f.info() for k, f in self.forms.items()}


ZSPLIT = 4.0          # boxes with |z| >= ZSPLIT use the large-z (zeta) forms


def box_geometry(ld, box):
    """Effective hulls on a box (la, lb, e0, e1, s0, s1, t0, t1) of the TM
    variable ell, eps, s = xi tau and tau (s1, t1 may be inf).  The box stands
    for {(ell, eps, tau): ell in [la, lb], eps in [e0, e1], tau in [t0, t1],
    xi(ell) tau in [s0, s1]}; the returned hulls enclose every such point.
    Returns None if the set is empty."""
    la, lb, e0, e1, s0, s1, t0, t1 = box
    up, dn = 1 + 1e-12, 1 - 1e-12
    xa, xb = ld.xi_of_ell(la), ld.xi_of_ell(lb)
    xlo = max(xa, (s0 / t1) * dn if t1 < INF else 0.0)
    xhi = min(xb, (s1 / t0) * up if s1 < INF else INF)
    if xlo > xhi:
        return None
    tlo = max(t0, (s0 / xhi) * dn if xhi > 0 else t0)
    thi = min(t1, (s1 / xlo) * up if (s1 < INF and xlo > 0) else INF)
    if tlo > thi:
        return None
    slo = max(s0, xlo * tlo * dn)
    shi = min(s1, xhi * thi * up) if thi < INF else s1
    if slo > shi:
        return None
    lla = max(la, ld.ell_of_xi(xlo) * dn if ld.laurent else ld.ell_of_xi(xlo) - 1e-15)
    llb = min(lb, ld.ell_of_xi(xhi) * up if ld.laurent else ld.ell_of_xi(xhi) + 1e-15)
    lla = max(lla, 0.0)
    if lla > llb:
        return None
    L = (lla, llb)
    # every hull below is built by ball arithmetic from the float endpoints
    # (which already carry the 1e-12 slack), so it is outward-rounded
    E = iv(e0, e1)
    Tb = arb(tlo) if thi == tlo else iv(tlo, thi) if thi < INF else None
    if thi < INF:
        w = acb(1) / acb(E, Tb)
        itau = arb(1) / Tb
    else:
        rw = arb(0).union(arb(e1) / arb(tlo) ** 2)
        w = acb(rw, (-arb(1) / arb(tlo)).union(arb(0)))
        itau = arb(0).union(arb(1) / arb(tlo))
    lam = ld.xi_hull(L)
    z = zeta = None
    use_large = False
    if ld.laurent:
        rez = iv(L[0], L[1]) * E
        if shi < INF and slo < ZSPLIT:
            z = acb(rez, iv(slo, shi))
        else:
            use_large = True
            if shi < INF:
                zeta = acb(1) / acb(rez, iv(slo, shi))
            else:
                # z = s (r + i), r = ell eps / s in [0, L1 e1 / slo]
                rmax = (arb(L[1]) * arb(e1) / arb(slo)).upper()
                eta = acb(1) / acb(arb(0).union(arb(rmax)), arb(1))
                zeta = acb(arb(0).union(arb(1) / arb(slo))) * eta
    return dict(L=L, E=E, w=w, itau=itau, lam=lam, z=z, zeta=zeta,
                use_large=use_large, tau=(tlo, thi), s=(slo, shi), xi=(xlo, xhi))


def sigma_box(ld, box, geo=None, terms=False):
    """Ball enclosing Sigma = S xi^2/tau^2 over the box (or None if empty).
    Raises ArithmeticError if a denominator hull contains 0.  terms=True
    returns the dict of the assembled terms instead (diagnostics)."""
    g = box_geometry(ld, box) if geo is None else geo
    if g is None:
        return None
    L, E, w, it, lam = g["L"], g["E"], g["w"], g["itau"], g["lam"]
    z, zeta, ul = g["z"], g["zeta"], g["use_large"]
    F = ld.forms
    R = lambda n, d, k: ratio_hull(F[n], F[d], k, L, lam, z, zeta, w, ul, ld.laurent)
    c = cval()
    t, tp, u2h = ld.real_hulls(L, lam)
    V = R("nV", "dT", 1)
    Vp = R("nVp", "dT", 2) - V * R("dT1", "dT", 1)
    Y = R("nY", "kdB", 2)
    dBp, dBpp = R("dB1", "dB", 1), R("dB2", "dB", 2)
    nBp, nBpp = R("nB1", "nB", 1), R("nB2", "nB", 2)
    Lh = nBp - dBp
    Lp = (nBpp - nBp * nBp) - (dBpp - dBp * dBp)
    W, Wp = V - Lh, Vp - Lp
    G = Y + t * Lh - tp
    H = Vp - V * Lh
    nu, nup = W.imag, Wp.imag
    rho = W.real - (E + c) * t
    cf = c * it * it + 1 / c
    T = dict(reB1=u2h * lam * (E * E * it * it - 1) + E * it * it * G.real
             - it * G.imag + it * it * H.real,
             imB1=-(2 * E * u2h * lam + G.real + E * it * G.imag + it * H.imag) / c,
             fp=-cf * tp / 2 + (it / c) * nup / 2,
             rhof=-rho * (cf * t - (it / c) * nu) / 2,
             nu2=nu * nu * it * it / 4)
    if terms:
        T.update(G=G, H=H, W=W, Wp=Wp, rho=rho, t=t, tp=tp, u2h=u2h)
        return T
    return T["reB1"] + T["imB1"] + T["fp"] + T["rhof"] + T["nu2"]


# =====================================================================
# adaptive cover of one piece
# =====================================================================
def eps_cells_default():
    e = [0.0, 0.03, 0.07, 0.12, 0.2, 0.35, 0.6, 1.0, 1.6, 2.4, 3.4, 4.6, 6.0]
    return list(zip(e, e[1:]))


def _normalize(ld, box):
    """Replace the box bounds by the effective ones (same point set)."""
    g = box_geometry(ld, box)
    if g is None:
        return None, None
    la, lb, e0, e1, s0, s1, t0, t1 = box
    L = g["L"]
    nb = (L[0], L[1], e0, e1, max(s0, g["s"][0]), min(s1, g["s"][1]),
          max(t0, g["tau"][0]), min(t1, g["tau"][1]))
    return nb, g


def _children(ld, box, var):
    la, lb, e0, e1, s0, s1, t0, t1 = box
    if var == "s":
        if s1 == INF:
            m = 4.0 * max(s0, 1e-300)
        elif s0 <= 0:
            m = s1 / 4.0
        elif s1 <= s0 * (1 + 1e-9):
            return None
        else:
            m = (s0 * s1) ** 0.5
        return ((la, lb, e0, e1, s0, m, t0, t1), (la, lb, e0, e1, m, s1, t0, t1))
    if var == "t":
        if t1 == INF:
            m = 4.0 * t0
        elif t1 <= t0 * (1 + 1e-9):
            return None
        else:
            m = (t0 * t1) ** 0.5
        return ((la, lb, e0, e1, s0, s1, t0, m), (la, lb, e0, e1, s0, s1, m, t1))
    if var == "e":
        if e1 - e0 <= 1e-7:
            return None
        m = (e0 + e1) / 2
        return ((la, lb, e0, m, s0, s1, t0, t1), (la, lb, m, e1, s0, s1, t0, t1))
    if var == "l":
        if lb - la <= 1e-12 * max(lb, 1e-300):
            return None
        if ld.laurent:
            m = lb / 16.0 if la <= 0 else (la * lb) ** 0.5
        else:
            m = (la + lb) / 2
        return ((la, m, e0, e1, s0, s1, t0, t1), (m, lb, e0, e1, s0, s1, t0, t1))
    raise ValueError(var)


def _upper(ld, box):
    """(normalized box or None if empty, upper bound of Sigma or +inf)."""
    nb, g = _normalize(ld, box)
    if nb is None:
        return None, -INF
    try:
        S = sigma_box(ld, nb, g)
    except (ArithmeticError, ZeroDivisionError, ValueError):
        return nb, INF
    if S is None:
        return None, -INF
    u = S.upper()
    return nb, (float(u) if u.is_finite() else INF)


def certify_piece(ld, la, lb, eps_cells=None, s_lo=None, s_a=1e-3, s_big=None,
                  s_ratio=2.0, target=0.0, max_evals=400000, max_depth=60, log=None,
                  t_pre=2e4, t_ratio=1.5, tau_min=None):
    """Certify Sigma <= target (< 0 suffices; target < 0 certifies a margin) on
    {ell in [la, lb], eps in [0, EPS_MAX], tau >= tau_min} (tau_min = TSTAR by
    default) by an adaptive cover with boxes in (ell, eps, s = xi tau, tau).
    The root boxes cover eps in [0, EPS_MAX] (cells), s in [xi(la) tau_min, inf]
    (s in [0, inf] on a piece touching xi = 0) and tau in [tau_min, inf];
    splits share endpoints.
    Returns a dict with ok, worst (the largest certified upper bound of Sigma,
    i.e. the certified margin is -worst), the number of certified boxes and
    evaluations, and up to 8 failing boxes."""
    t_start = time.time()
    tmin = TSTAR if tau_min is None else float(tau_min)
    if not tmin > 0:
        raise ValueError("tau_min must be positive")
    if not target <= 0:
        raise ValueError("target must be <= 0")
    eps_cells = eps_cells_default() if eps_cells is None else eps_cells
    if not (eps_cells[0][0] == 0.0 and eps_cells[-1][1] == EPS_MAX and all(
            a[1] == b[0] for a, b in zip(eps_cells, eps_cells[1:]))):
        raise ValueError("eps cells must tile [0, EPS_MAX]")
    xa = ld.xi_of_ell(la)
    if s_lo is None:
        s_lo = xa * tmin * (1 - 1e-12)
    if s_big is None:
        # the outer pseudo-box s >= s_big needs |w| = 1/tau small as well as
        # |zeta| = 1/s: start it at tau >= 1e4 on pieces with xi > 0.1
        s_big = max(1e3, 1e4 * ld.xi_of_ell(lb))
    grid = []
    if s_lo <= 0:
        grid.append((0.0, s_a))
        s = s_a
    else:
        s = s_lo
    while s < s_big:
        s2 = min(s * s_ratio, s_big)
        grid.append((s, s2))
        s = s2
    grid.append((s_big, INF))
    # On a piece touching xi = 0 a root box with s0 < t_pre * xi(lb) has an
    # unbounded tau-range starting below t_pre; its 1/tau terms cancel in
    # truth but not in interval arithmetic, so tau is pre-split geometrically
    # there (a tiling of [TSTAR, inf], hence the cover is unchanged).
    edges = [tmin]
    while edges[-1] < t_pre:
        edges.append(min(edges[-1] * t_ratio, t_pre))
    slabs = list(zip(edges, edges[1:])) + [(edges[-1], INF)]
    xb = ld.xi_of_ell(lb)
    stack = []
    for e0, e1 in eps_cells:
        for s0, s1 in grid:
            if xa <= 0 and s0 < t_pre * xb:
                for ta, tb in slabs:
                    stack.append(((la, lb, e0, e1, s0, s1, ta, tb), 0, None))
            else:
                stack.append(((la, lb, e0, e1, s0, s1, tmin, INF), 0, None))
    nev = ncert = nempty = 0
    worst = -INF
    fails = []
    hist = {}
    while stack:
        box, depth, up = stack.pop()
        if up is None:
            nb, up = _upper(ld, box)
            nev += 1
            if nb is None:
                nempty += 1
                continue
            box = nb
        if up < target or (target == 0 and up < 0):
            ncert += 1
            worst = max(worst, up)
            hist[depth] = hist.get(depth, 0) + 1
            continue
        if depth >= max_depth or nev >= max_evals:
            fails.append(dict(box=list(box), upper=up, depth=depth))
            if len(fails) >= 8:
                break
            continue
        best = None
        for var in ("s", "e", "t", "l"):
            ch = _children(ld, box, var)
            if ch is None:
                continue
            res = []
            for cb in ch:
                nb, u = _upper(ld, cb)
                nev += 1
                res.append((nb, u))
            score = max(u for _, u in res)
            if best is None or score < best[0]:
                best = (score, var, res)
        if best is None:
            fails.append(dict(box=list(box), upper=up, depth=depth, why="unsplittable"))
            if len(fails) >= 8:
                break
            continue
        for nb, u in best[2]:
            if nb is None:
                nempty += 1
                continue
            stack.append((nb, depth + 1, u))
    ok = not fails and not stack
    out = dict(ok=ok, worst=(worst if worst > -INF else None), boxes=ncert,
               evals=nev, empty=nempty, fails=fails[:8],
               depth_hist={str(k): v for k, v in sorted(hist.items())},
               time=time.time() - t_start,
               params=dict(s_a=s_a, s_big=s_big, s_ratio=s_ratio, zsplit=ZSPLIT,
                           tstar=TSTAR, tau_min=tmin, eps_cells=eps_cells, la=la, lb=lb,
                           target=target, s_lo=s_lo, t_pre=t_pre, t_ratio=t_ratio))
    if log:
        log("certify_piece %s: ok=%s worst=%s boxes=%d evals=%d t=%.1fs"
            % (ld.label, ok, out["worst"], ncert, nev, out["time"]))
    return out


# =====================================================================
# piece wrappers (called per piece by scripts/s44_layer.py)
# =====================================================================
def sonic_layer_data(at, ex, hi, prec=320, label=None):
    """LayerData of a sonic piece: Laurent TMs on xi in [0, hi], built exactly as
    tc_commutator.run_sonic builds its small-xi pieces (sd_shift = 1 and
    sonic_raise_inputs, justified by E1/E3/E5)."""
    from . import tc_commutator as tc
    from .tc_tm import AtomKP, sonic_components
    with precision(prec):
        comps = sonic_components(ex, hi)
        ak = AtomKP(at, comps, sd_shift=1)
        pr = tc.sonic_raise_inputs(tc.build_pairs(ak))
        return LayerData(pr, True, label=label or "sonic[0,%g]" % hi)


def tube_layer_data(at, sd, prec=320, label=None):
    """LayerData of a tube step (sigma = x_k - x in [0, h]), built exactly as
    tc_commutator.run_tube builds its pairs (no Laurent structure)."""
    from . import tc_commutator as tc
    from .tc_tm import AtomKP, step_components
    with precision(prec):
        comps = step_components(sd)
        ak = AtomKP(at, comps, sd_shift=0)
        pr = tc.build_pairs(ak)
        return LayerData(pr, False, xik=-sd.x,
                         label=label or "step[x=%.6f,h=%g]" % (float(sd.x), sd.h))


MARGIN_LADDER = (-0.75, -0.5, -0.25)


def run_layer_piece(ld, la, lb, prec=320, ladder=MARGIN_LADDER, ladder_evals=20000,
                    **kw):
    """Certify Sigma < 0 for tau >= TSTAR (or tau >= kw["tau_min"]) on the piece.  The margin ladder
    first tries to certify Sigma <= t for t in `ladder` with a small
    evaluation budget each and stops at the first success; otherwise it
    certifies Sigma < 0 with the full budget.  Returns the certify_piece
    record of the successful (or last) attempt plus the attempt list and the
    form structure (pivot data)."""
    attempts = []
    res = None
    with precision(prec):
        for tg in ladder:
            r = certify_piece(ld, la, lb, target=tg, max_evals=ladder_evals, **kw)
            attempts.append(dict(target=tg, ok=r["ok"], evals=r["evals"],
                                 time=r["time"]))
            if r["ok"]:
                res = r
                break
        if res is None:
            res = certify_piece(ld, la, lb, target=0.0, **kw)
            attempts.append(dict(target=0.0, ok=res["ok"], evals=res["evals"],
                                 time=res["time"]))
    res["attempts"] = attempts
    res["certified_bound"] = res["worst"] if res["ok"] else None
    res["forms"] = ld.form_info()
    res["laurent"] = ld.laurent
    return res
