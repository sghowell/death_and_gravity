"""S2 part-(i) probe, route M3: exclusion of QNMs in the barrier-top window by complex rotation of
the tortoise coordinate and a variable-phase coercivity bound certified with a ground state.

Proof: notes/s2-m3.md.  Pre-registration: notes/s2-preregistration.md (window, criteria, fit).

Contour: x(s) = x_P + e^{i theta} s, theta = pi/4, x_P = r*(r_P).  On it, a QNM gives a decaying
solution of  -u'' + Q u = 0,  Q = e^{2 i theta}(V_l(r(s)) - w^2).  For a phase function phi(s),
eps1 = eta |phi'|, c = cos(phi) - eps1 and
    q = Re(e^{-i phi} Q) - |phi'|/(4 eta),
    Re < e^{-i phi}(-u'' + Q u), u >  >=  int c |u'|^2 + q |u|^2,
and for g = e^{-G} > 0 with G' continuous and piecewise linear,
    int c|u'|^2 + q|u|^2  >=  mu ||u||^2   if   E(s) := q + c'G' - c G'^2 + c G'' >= mu  (all s).
Certified here: E >= mu_cert > 0 on s-cells (mean-value form, Arb), for w^2 in the interval hull of
{w^2 : w in the window}; E is affine in w^2, so this is exact in w^2.  Tails |s| > S by the
horizon bounds of m3_contour.

Replay:  python -m p2.m3 verify      Run:  python -m p2.m3 run [l ...]
"""
from __future__ import annotations

import itertools
import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
from flint import acb, arb, arb_series, fmpq

from .arbutil import Timer, git_commit, precision
from .m3_contour import F_of, Ray, rstar
from .params import Point

PREC = 128
CERT_DIR = Path(__file__).resolve().parents[2] / "certificates" / "m3"
H_STEP = 0.0625               # contour step (exact dyadic)
S_LEFT, S_RIGHT = -64.0, 32.0  # validated contour range (exact)

# design parameters per l (float design study; recorded in every certificate)
DESIGN = {
    "s1": 5.0, "p": 4, "phi_inf": -1.4375, "eta": 0.5,
    "lo_gamma0": 8, "hi_gamma0": 6, "slack": 0.05, "node_ds": 0.01,
}


def fa(x):
    """Exact arb from a Fraction / int / exact-binary float."""
    if isinstance(x, Fraction):
        return arb(fmpq(x.numerator, x.denominator))
    return arb(x)


def theta_unit():
    """e^{i pi/4} exactly as a ball."""
    return acb(1, 1) / arb(2).sqrt()


# ---------------------------------------------------------------------------------------------
# geometric bounds for wide boxes (midpoint-radius squaring is pessimistic for wide balls)
# ---------------------------------------------------------------------------------------------
def dist_bounds(z, p):
    """(lower, upper) of |w - p| over the complex box z (acb), p real (arb ball)."""
    a, b = (z.real - p).lower(), (z.real - p).upper()
    c, d = z.imag.lower(), z.imag.upper()
    mx = arb(0) if (a <= 0 <= b) else min(abs(a), abs(b))
    my = arb(0) if (c <= 0 <= d) else min(abs(c), abs(d))
    Mx = max(abs(a), abs(b))
    My = max(abs(c), abs(d))
    lo_ = arb((mx * mx + my * my).sqrt().lower())
    hi_ = arb((Mx * Mx + My * My).sqrt().upper())
    return lo_, hi_


# ---------------------------------------------------------------------------------------------
# potential
# ---------------------------------------------------------------------------------------------
class Pot:
    def __init__(self, pt: Point, l: int):
        self.pt, self.L = pt, l * (l + 1)

    def Delta(self, r):
        prod = acb(-self.pt.Lam / 3)
        for ri in self.pt.roots:
            prod = prod * (r - ri)
        return prod

    def parts(self, r):
        """V_l(r) and dV_l/dr (acb) for an acb r (ball or box)."""
        pt, L = self.pt, self.L
        lam3 = pt.Lam / 3
        D = self.Delta(r)                                  # factored: no cancellation near roots
        Dp = -4 * lam3 * r ** 3 + 2 * r - 2
        Dpp = -12 * lam3 * r * r + 2
        K = L * r * r + r * Dp - 2 * D                    # V = D K / r^6
        Kp = 2 * L * r + Dp + r * Dpp - 2 * Dp
        N = D * K
        Np = Dp * K + D * Kp
        V = N / r ** 6
        Vp = (r * Np - 6 * N) / r ** 7
        return V, Vp


# ---------------------------------------------------------------------------------------------
# phase profile phi(s) = phi_inf s^4/(s1^4 + s^4) and derivatives (exact rational function)
# ---------------------------------------------------------------------------------------------
def phi_derivs(s, design):
    """(phi, phi', phi'', phi''') at s (arb point or interval) via a Taylor series in t."""
    s1 = arb(design["s1"])
    phinf = arb(design["phi_inf"])
    x = arb_series([s, 1], prec=4)
    x4 = x ** 4
    f = phinf * x4 / (s1 ** 4 + x4)
    c = [f.coeffs()[k] if k < len(f.coeffs()) else arb(0) for k in range(4)]
    return c[0], c[1], 2 * c[2], 6 * c[3]


# ---------------------------------------------------------------------------------------------
# window
# ---------------------------------------------------------------------------------------------
def window(pt: Point, l: int, design=DESIGN):
    wl = (l + arb(1) / 2) * pt.OmegaP
    a1 = wl - design["lo_gamma0"] * pt.gamma0
    a2 = wl + design["hi_gamma0"] * pt.gamma0
    b1 = -pt.beta
    return wl, a1, a2, b1


def _iv(lo_, hi_):
    return (lo_ + hi_) / 2 + arb(0, 1) * ((hi_ - lo_) / 2)


def w2_hull(a1, a2, b1, b2=None):
    """Interval hull of {w^2 : Re w in [a1, a2], Im w in [b1, b2]} with 0 < a1, b1 <= b2 <= 0:
    Re w^2 = a^2 - b^2 in [a1^2 - b1^2, a2^2 - b2^2], Im w^2 = 2ab in [2 a2 b1, 2 a1 b2]
    (endpoint arithmetic; midpoint-radius squaring would be pessimistic)."""
    b2 = arb(0) if b2 is None else b2
    A1, A2, B1, B2 = a1.lower(), a2.upper(), b1.lower(), b2.upper()
    re_lo = (A1 * A1 - B1 * B1).lower()
    re_hi = (A2 * A2 - B2 * B2).upper()
    im_lo = (2 * A2 * B1).lower()
    im_hi = (2 * A1 * B2).upper()
    return acb(_iv(re_lo, re_hi), _iv(im_lo, im_hi))


# ---------------------------------------------------------------------------------------------
# ground-state design (float only; any positive g is admissible)
# ---------------------------------------------------------------------------------------------
def float_contour(pt: Point, s_lo, s_hi, h):
    """RK4 float path r(s) on the uniform grid s_lo, s_lo + h, ..., s_hi (design only)."""
    rs = [complex(float(x.mid())) for x in pt.roots]
    lam3 = float((pt.Lam / 3).mid())
    e = np.exp(1j * np.pi / 4)

    def F(r):
        return -lam3 * (r - rs[0]) * (r - rs[1]) * (r - rs[2]) * (r - rs[3]) / r ** 2

    def run(sgn, n):
        r = complex(float(pt.rP.mid()))
        out = [r]
        for _ in range(n):
            k1 = e * F(r); k2 = e * F(r + sgn * h / 2 * k1); k3 = e * F(r + sgn * h / 2 * k2); k4 = e * F(r + sgn * h * k3)
            r = r + sgn * h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            out.append(r)
        return np.array(out)

    nL, nR = round(-s_lo / h), round(s_hi / h)
    left = run(-1, nL)
    right = run(1, nR)
    s = np.concatenate([-h * np.arange(nL, 0, -1), h * np.arange(0, nR + 1)])
    r = np.concatenate([left[::-1][:-1], right])
    return s, r


def float_V(pt: Point, r, l):
    rs = [complex(float(x.mid())) for x in pt.roots]
    lam3 = float((pt.Lam / 3).mid())
    Q2 = float(pt.Q2.mid())
    L = l * (l + 1)
    D = -lam3 * r ** 4 + r ** 2 - 2 * r + Q2
    Dp = -4 * lam3 * r ** 3 + 2 * r - 2
    del rs
    return D * (L * r * r + r * Dp - 2 * D) / r ** 6


def design_ground_state(pt: Point, l: int, w: complex, design=DESIGN, s_lo=-60.0, s_hi=30.0, h=0.001):
    """Float ground state of -(c u')' + q u at the worst corner w; returns (lambda_min, G' nodes)."""
    from scipy.sparse import diags
    from scipy.sparse.linalg import eigsh
    s, r = float_contour(pt, s_lo, s_hi, h)
    Q = 1j * (float_V(pt, r, l) - w * w)          # e^{2 i pi/4} = i
    s1, p, phinf, eta = design["s1"], design["p"], design["phi_inf"], design["eta"]
    phi = phinf * s ** p / (s1 ** p + s ** p)
    dphi = phinf * p * s ** (p - 1) * s1 ** p / (s1 ** p + s ** p) ** 2
    c = np.cos(phi) - eta * np.abs(dphi)
    q = np.real(np.exp(-1j * phi) * Q) - np.abs(dphi) / (4 * eta)
    ds = s[1] - s[0]
    cm = 0.5 * (c[1:] + c[:-1])
    main = (np.concatenate([cm, [cm[-1]]]) + np.concatenate([[cm[0]], cm])) / ds ** 2 + q
    off = -cm / ds ** 2
    Aop = diags([off, main, off], [-1, 0, 1], format="csc")
    val, vec = eigsh(Aop, k=1, sigma=float(q.min()) - 1.0, which="LM")
    g = np.abs(vec[:, 0])
    lam = float(val[0])
    # G' := y solving the Riccati equation c y' = mu_t - q - c'y + c y^2 (so E = mu_t), with
    # mu_t = (1 - slack) lam, from y = 0 at the maximum of g outward until y returns to 0; y = 0
    # beyond (there q > mu_t).  Any positive g is admissible: this is design only.
    mu_t = (1 - design.get("slack", 0.05)) * lam
    cp = np.gradient(c, ds)
    i_star = int(np.argmax(g))

    def rhs(k, y):
        return (mu_t - q[k]) / c[k] - cp[k] / c[k] * y + y * y

    def sweep(direction):
        ys = {i_star: 0.0}
        y = 0.0
        k = i_star
        hstep = ds * direction
        while 0 < k + direction < len(s) - 1:
            # Heun step on the grid
            k1 = rhs(k, y)
            k2 = rhs(k + direction, y + hstep * k1)
            y_new = y + hstep * (k1 + k2) / 2
            k += direction
            if direction * (k - i_star) > 3 and np.sign(y_new) != np.sign(y) and y != 0:
                ys[k] = 0.0
                return ys, k
            y = y_new
            ys[k] = y
        raise ArithmeticError("Riccati sweep did not return to 0")

    yr, _kR = sweep(+1)
    yl, _kL = sweep(-1)
    ys = {**yl, **yr}
    ks = sorted(ys)
    # nodes: subsample the Riccati solution (spacing ~ design['node_ds'])
    stride = max(1, round(design.get("node_ds", 0.01) / ds))
    sel = ks[::stride]
    if sel[-1] != ks[-1]:
        sel.append(ks[-1])
    if sel[0] != ks[0]:
        sel.insert(0, ks[0])
    out = []
    for k in sel:
        out.append((Fraction(round(s[k] * 2 ** 20), 2 ** 20), Fraction(round(ys[k] * 2 ** 30), 2 ** 30)))
    # endpoints exactly 0 (y returned to 0 there)
    out[0] = (out[0][0], Fraction(0))
    out[-1] = (out[-1][0], Fraction(0))
    # make 0 a node (sign change of phi')
    if not any(x == 0 for x, _ in out):
        k = max(i for i, (x, _) in enumerate(out) if x < 0)
        x0, y0 = out[k]
        x1, y1 = out[k + 1]
        y = y0 + (y1 - y0) * (0 - x0) / (x1 - x0)
        out.insert(k + 1, (Fraction(0), Fraction(round(float(y) * 2 ** 30), 2 ** 30)))
    # remove duplicate abscissae
    ded = []
    for x, y in out:
        if not ded or x > ded[-1][0]:
            ded.append((x, y))
    out = ded
    return lam, out


# ---------------------------------------------------------------------------------------------
# the certificate
# ---------------------------------------------------------------------------------------------
class Checker:
    def __init__(self, pt: Point, l: int, rays, design=DESIGN):
        self.pt, self.l, self.design = pt, l, design
        self.pot = Pot(pt, l)
        self.E1 = theta_unit()
        self.E2 = acb(0, 1)                    # e^{2 i pi/4}
        self.E3 = self.E2 * self.E1
        self.rays = rays                       # {-1: left Ray, +1: right Ray}
        self.eta = arb(design["eta"])

    def ray_for(self, s_mid):
        return self.rays[1] if s_mid >= 0 else self.rays[-1]

    def E_parts(self, s, r, W2, Gp, Gpp, sgn):
        """E and E' (as arbs) at s (point or interval) with r enclosing r(s), w^2 in W2,
        G' = Gp, G'' = Gpp (constant on the cell), sgn = sign(phi') on the cell (0 allowed only
        at a point)."""
        phi, d1, d2, d3 = phi_derivs(s, self.design)
        eta = self.eta
        ad1 = d1 * sgn                          # |phi'|
        cphi, sphi = phi.cos(), phi.sin()
        c = cphi - eta * ad1
        cp = -sphi * d1 - eta * sgn * d2
        cpp = -cphi * d1 * d1 - sphi * d2 - eta * sgn * d3
        V, Vp = self.pot.parts(r)
        Q = self.E2 * (V - W2)
        rot = acb(cphi, -sphi)                  # e^{-i phi}
        q = (rot * Q).real - ad1 / (4 * eta)
        Qp = self.E3 * Vp * F_of(self.pt, r)    # dQ/ds = e^{2it} V'(r) e^{it} F(r)
        qp = (rot * (Qp - acb(0, 1) * d1 * Q)).real - sgn * d2 / (4 * eta)
        E = q + cp * Gp - c * Gp * Gp + c * Gpp
        Ep = qp + cpp * Gp + 2 * cp * Gpp - cp * Gp * Gp - 2 * c * Gp * Gpp
        return E, Ep, c

    def check_cell(self, sa, sb, W2, Gpa, Gpb):
        """Lower bound of E on [sa, sb] (exact rationals), G' linear from Gpa to Gpb."""
        sa_, sb_ = fa(sa), fa(sb)
        w = sb_ - sa_
        sm = (sa_ + sb_) / 2
        cell = sm + arb(0, 1) * (w / 2)
        Gpa, Gpb = fa(Gpa), fa(Gpb)
        Gpp = (Gpb - Gpa) / w
        sgn = 1 if float(sm.mid()) > 0 else -1
        sgn = sgn * (1 if self.design["phi_inf"] > 0 else -1)
        ray = self.ray_for(float(sm.mid()))
        r_m = ray.at(sm)
        R = ray.enclose(sa_, sb_)
        Gp_m = (Gpa + Gpb) / 2
        Gp_cell = Gp_m + arb(0, 1) * (abs(Gpb - Gpa) / 2)
        E_m, _, _c_m = self.E_parts(sm, r_m, W2, Gp_m, Gpp, sgn)
        _, Ep_c, c_c = self.E_parts(cell, R, W2, Gp_cell, Gpp, sgn)
        low = E_m.lower() - abs(Ep_c).upper() * (w / 2)
        return arb(low), arb(c_c.lower())

    def cells_for(self, nodes, S_lo, S_hi, h_out):
        """Cell list [(sa, sb, G'a, G'b)] covering [S_lo, S_hi]: spline cells in the core,
        G' = 0 cells of width <= h_out outside; boundaries aligned with contour steps."""
        cells = []
        xs = [x for x, _ in nodes]
        # left outer region
        def fill(a, b):
            n = max(1, int(np.ceil(float(b - a) / h_out)))
            for k in range(n):
                cells.append((a + (b - a) * k / n, a + (b - a) * (k + 1) / n, Fraction(0), Fraction(0)))
        fill(Fraction(S_lo), xs[0])
        for (x0, y0), (x1, y1) in itertools.pairwise(nodes):
            cells.append((x0, x1, y0, y1))
        fill(xs[-1], Fraction(S_hi))
        # split cells at contour step boundaries (multiples of H_STEP) and at 0
        hstep = Fraction(H_STEP)
        out = []
        for (a, b, ya, yb) in cells:
            pts = [a]
            k = (a // hstep) + 1
            while k * hstep < b:
                pts.append(k * hstep)
                k += 1
            if a < 0 < b and 0 not in pts:
                pts.append(Fraction(0))
            pts = sorted(set(pts)) + [b]
            for u, v in itertools.pairwise(pts):
                yu = ya + (yb - ya) * (u - a) / (b - a)
                yv = ya + (yb - ya) * (v - a) / (b - a)
                out.append((u, v, yu, yv))
        return out

    def run(self, nodes, W2, S_lo, S_hi, target, h_out=Fraction(1, 16), max_depth=10):
        """min over cells of the certified lower bound of E (bisecting cells whose bound is below
        `target`, up to max_depth), the min of c over cells, #evaluations, worst cell."""
        cells = self.cells_for(nodes, S_lo, S_hi, h_out)
        mu, cmin = None, None
        n_eval = 0
        worst = None
        for (a, b, ya, yb) in cells:
            stack = [(a, b, ya, yb, 0)]
            while stack:
                a_, b_, ya_, yb_, d = stack.pop()
                low, clow = self.check_cell(a_, b_, W2, ya_, yb_)
                n_eval += 1
                if low > target or d >= max_depth:
                    if mu is None or low < mu:
                        mu = low
                        worst = (float(a_), float(b_))
                    if cmin is None or clow < cmin:
                        cmin = clow
                    continue
                m = (a_ + b_) / 2
                ym = (ya_ + yb_) / 2
                stack.append((m, b_, ym, yb_, d + 1))
                stack.append((a_, m, ya_, ym, d + 1))
        return mu, cmin, n_eval, worst


# ---------------------------------------------------------------------------------------------
# tails |s| > S: |r(s) - r_h| <= d(s) from Re r*(r(s)) = x_P + s cos(theta)
# ---------------------------------------------------------------------------------------------
def tail_bound(pt: Point, side, S, r_S, d0):
    """Certified sup_{|s| >= |S|} |r(s) - r_h| (or None): side -1 -> r_+, +1 -> r_c."""
    xP = rstar(pt, acb(pt.rP))
    ct = arb(2).sqrt() / 2
    rh = pt.rp if side < 0 else pt.rc
    others = [(Ai, ri) for Ai, ri in zip(pt.Acoef, pt.roots) if not ri.overlaps(rh)]
    Ah = next(Ai for Ai, ri in zip(pt.Acoef, pt.roots) if ri.overlaps(rh))
    disc = acb(rh + arb(0, 1) * d0, arb(0, 1) * d0)
    Smin, Smax = arb(0), arb(0)
    for Ai, ri in others:
        lo_, hi_ = dist_bounds(disc, ri)
        if not (lo_ > 0):
            return None
        t1, t2 = Ai * lo_.log(), Ai * hi_.log()          # A_i log|.| is monotone in the distance
        Smin = Smin + arb(t1.lower()).min(arb(t2.lower()))   # exact Arb endpoints (no float round trip)
        Smax = Smax + arb(t1.upper()).max(arb(t2.upper()))
    if not (abs(r_S - rh) < d0):
        return None
    xS = xP + arb(S) * ct
    if side < 0:      # A_+ > 0: |r - r_+| = exp((Re x - Sigma)/A_+) <= exp((x_S - Sigma_min)/A_+)
        d = ((xS - arb(Smin.lower())) / Ah).exp()
    else:             # A_c < 0: |r - r_c| <= exp((x_S - Sigma_max)/A_c)
        d = ((xS - arb(Smax.upper())) / Ah).exp()
    d = arb(d.upper())
    return d if d < d0 else None


def tail_check(pt: Point, l: int, side, S, d, W2, design=DESIGN):
    """Lower bound of E = q on |s| >= |S| (G' = 0 there): q >= -|V| + min Re(-e^{-i phi} i w^2)
    - |phi'(S)|/(4 eta), with phi between phi(S) and phi_inf."""
    pot = Pot(pt, l)
    rh = pt.rp if side < 0 else pt.rc
    disc = acb(rh + arb(0, 1) * d, arb(0, 1) * d)
    V, _ = pot.parts(disc)
    phS, d1S, _, _ = phi_derivs(arb(S), design)
    phinf = arb(design["phi_inf"])
    lo_, hi_ = (phS, phinf) if phS < phinf else (phinf, phS)
    phi_int = (lo_.lower() + hi_.upper()) / 2 + arb(0, 1) * ((hi_.upper() - lo_.lower()) / 2)
    rot = acb(phi_int.cos(), -phi_int.sin())
    term = (rot * acb(0, 1) * (-W2)).real
    q_low = term.lower() - abs(V).upper() - abs(d1S).upper() / (4 * arb(design["eta"]))
    return arb(q_low)


# ---------------------------------------------------------------------------------------------
# geometry certificate (independent of l and w)
# ---------------------------------------------------------------------------------------------
def real_point(pt: Point, X):
    """r in (r_+, r_c) with r*(r) = X (interval Newton on the increasing real map)."""
    import mpmath as mp
    lo_, hi_ = float((pt.rp).mid()) + 1e-9, float(pt.rc.mid()) - 1e-9
    f = lambda r: float(rstar(pt, acb(arb(r))).mid()) - float(arb(X).mid())
    for _ in range(200):
        m = (lo_ + hi_) / 2
        if f(m) > 0:
            hi_ = m
        else:
            lo_ = m
    x0 = arb(mp.mpf((lo_ + hi_) / 2))
    for rad in ("1e-12", "1e-9", "1e-6"):
        Rb = x0 + arb(0, 1) * arb(rad)
        fpr = sum((Ai / (Rb - ri) for Ai, ri in zip(pt.Acoef, pt.roots)), arb(0))   # = r^2/Delta
        N = x0 - (rstar(pt, acb(x0)) - X) / fpr
        if N.lower() > Rb.lower() and N.upper() < Rb.upper():
            return N
    raise ArithmeticError("real point not certified")


class VRay(Ray):
    """Validated integration along x = X + i Y, Y from 0 to Ymax (start on the real axis)."""

    def __init__(self, pt, X, Ymax, h=0.0625, N=20):
        self.pt, self.E, self.N = pt, acb(0, 1), N
        self.steps = []
        r0 = acb(real_point(pt, X))
        from .m3_contour import Step, apriori_box, taylor_coeffs
        s = arb(0)
        S = arb(Ymax)
        hh0 = arb(h)
        while S - s > 0:
            hh = hh0 if S - s > hh0 else S - s          # noqa: FURB136 (arb: certain comparison)
            B = apriori_box(pt, r0, self.E, hh)
            if B is None:
                raise ArithmeticError("vertical ray: no a-priori box")
            m0 = acb(r0.real.mid(), r0.imag.mid())
            D0 = r0 - m0
            co = taylor_coeffs(pt, m0, self.E, N)
            rem = taylor_coeffs(pt, B, self.E, N + 1)[N + 1]
            F0 = F_of(pt, r0)
            J = F_of(pt, B) / F0
            st = Step(s, hh, co, rem, B, N, D0, J)
            Rh = st.eval_mid(hh) + J * D0
            r0 = st.eval_mid(hh) + (F_of(pt, Rh) / F0) * D0
            self.steps.append(st)
            s = s + hh
        self.S = S


def geometry_certificate(pt: Point):
    """A0: contour enclosure + identity, sector analyticity, funnel tail bounds (all l)."""
    with Timer() as tm:
        E = theta_unit()
        left = Ray(pt, E, arb(S_LEFT), arb(H_STEP))
        right = Ray(pt, E, arb(S_RIGHT), arb(H_STEP))
        xP = rstar(pt, acb(pt.rP))
        ct = arb(2).sqrt() / 2
        ident_ok = True
        for ray in (left, right):
            for st in ray.steps[::8]:
                rr = st.eval(arb(0))
                if not (rstar(pt, rr) - (xP + st.s0 * ct)).contains(0):
                    ident_ok = False
        # sector analyticity: singularities of r(x) only where Re x = Re r*(0) or Re x = 0
        x0re = sum((Ai * abs(ri).log() for Ai, ri in zip(pt.Acoef, pt.roots)), arb(0))
        sumA = sum(pt.Acoef, arb(0))
        lower_left_free = bool(x0re > xP)            # lower-left sector has Re x <= x_P
        seg = {}
        for name, X in (("Re_x0", x0re), ("Re_xinf", arb(0))):
            Ymax = X - xP                            # the upper-right sector at Re x = X
            vr = VRay(pt, X, arb(Ymax.upper()))
            rmin = min((abs(st.B).lower() for st in vr.steps), key=lambda z: float(z))
            rmax = max((abs(st.B).upper() for st in vr.steps), key=lambda z: float(z))
            seg[name] = {"X": X.str(12), "Ymax": arb(Ymax.upper()).str(8), "steps": len(vr.steps),
                         "min_abs_r": arb(rmin).str(6), "max_abs_r": arb(rmax).str(6),
                         "ok": bool(rmin > 0 and rmax < 1e6)}
        # tails
        dL = tail_bound(pt, -1, S_LEFT, left.at(arb(S_LEFT)), arb("0.06"))
        dR = tail_bound(pt, +1, S_RIGHT, right.at(arb(S_RIGHT)), arb("0.6"))
    cert = {"theta": "pi/4", "x_P": xP.str(15), "contour_range": [S_LEFT, S_RIGHT], "step": H_STEP,
            "taylor_order": left.N, "identity_Re_rstar_holds": ident_ok,
            "sum_A_contains_0": bool(sumA.contains(0)),
            "Re_rstar_0": x0re.str(12), "lower_left_sector_free": lower_left_free,
            "upper_right_special_segments": seg,
            "tail_left_sup_dist_to_rplus": dL.str(6) if dL is not None else None,
            "tail_right_sup_dist_to_rc": dR.str(6) if dR is not None else None,
            "timing_s": {"wall": round(tm.wall, 2), "cpu": round(tm.cpu, 2)}}
    cert["A0_pass"] = bool(ident_ok and lower_left_free and all(v["ok"] for v in seg.values())
                           and dL is not None and dR is not None and sumA.contains(0))
    return cert, {-1: left, 1: right}, (dL, dR)


def certify_l(pt: Point, l: int, rays, tails, design=DESIGN, nodes=None, lam_float=None, box=None):
    """Certificate for the window of l (or for an explicit box (a1, a2, b1, b2) of w)."""
    if box is None:
        wl, a1, a2, b1 = window(pt, l, design)
        b2 = arb(0)
    else:
        a1, a2, b1, b2 = box
        wl = (a1 + a2) / 2
    W2 = w2_hull(a1, a2, b1, b2)
    if nodes is None:
        wc = complex(float(a1.mid()), float(b1.mid())) if box is None else \
            complex(float(((a1 + a2) / 2).mid()), float(((b1 + b2) / 2).mid()))
        lam_float, nodes = design_ground_state(pt, l, wc, design)
    target = arb(0.9 * (1 - design.get("slack", 0.05)) * lam_float)
    ch = Checker(pt, l, rays, design)
    with Timer() as tm:
        mu, cmin, n_eval, worst = ch.run(nodes, W2, Fraction(S_LEFT), Fraction(S_RIGHT), target)
        dL, dR = tails
        tL = tail_check(pt, l, -1, S_LEFT, dL, W2, design)
        tR = tail_check(pt, l, +1, S_RIGHT, dR, W2, design)
    mu_all = min([mu, tL, tR], key=lambda z: float(z.lower()))
    muHO = 2 * (l + arb(1) / 2) * pt.OmegaP * pt.gamma0
    # Jost decay on the contour: Im(w e^{i pi/4}) = (Re w + Im w)/sqrt(2) > 0 on the window
    decay_ok = bool(a1 + b1 > 0)
    # tails: c >= cos(phi_inf) - eta |phi'(S)| > 0 there (phi between phi(S) and phi_inf)
    c_tail_ok = True
    for S in (S_LEFT, S_RIGHT):
        ph, d1, _, _ = phi_derivs(arb(S), design)
        lo_ = min(float(ph.lower()), design["phi_inf"])
        c_t = arb(lo_).cos() - arb(design["eta"]) * abs(d1)
        c_tail_ok = c_tail_ok and bool(c_t > 0) and abs(arb(lo_)) < arb(1.57)
    # P2-2 overlap: a_l <= window top (checked by the caller when a_l is available)
    return {"l": l, "omega_l_over_kappa_minus": (wl / pt.km).str(10),
            "window_Re_over_kappa_minus": [(a1 / pt.km).str(8), (a2 / pt.km).str(8)],
            "window_Im_over_kappa_minus": [(b1 / pt.km).str(8), "0"],
            "w2_hull": [W2.real.str(8), W2.imag.str(8)],
            "design": {k: (str(v) if not isinstance(v, (int, float)) else v) for k, v in design.items()},
            "ground_state_nodes": [[str(x), str(y)] for x, y in nodes],
            "float_lambda_min": lam_float,
            "cells_evaluated": n_eval, "worst_cell": worst,
            "mu_cells_lower": mu.str(10), "mu_tail_left_lower": tL.str(6), "mu_tail_right_lower": tR.str(6),
            "mu_cert": arb(mu_all.lower()).str(10),
            "mu_HO": muHO.str(10),
            "m_ratio": (arb(mu_all.lower()) / muHO).str(8),
            "jost_decay_condition": decay_ok,
            "kinetic_coefficient_min_lower": cmin.str(8), "c_nonneg": bool(cmin >= 0 and c_tail_ok),
            "A1_pass": bool(mu_all > 0 and decay_ok and cmin >= 0 and c_tail_ok),
            "timing_s": {"wall": round(tm.wall, 1), "cpu": round(tm.cpu, 1)}}


def control_l(pt: Point, l: int, rays, tails, window_nodes, window_lam, design=DESIGN):
    """A2: boxes of half-width gamma0/50 around the l-th PS fundamental mode (float, S1 solver).
    The rigorous checker must NOT certify them, (a) with the window's ground state and (b) with a
    ground state designed at the box (a design failure also counts as 'not certified')."""
    import mpmath as mp

    from .qnm_float import Geom, Solver
    g = Geom.frozen("PB")
    with mp.workdps(30):
        top = (l + mp.mpf(1) / 2) * g.OmegaP
        ps = Solver(g, K=6).find(l, mp.mpc(top, -g.gamma0))
    e = float(pt.gamma0.mid()) / 50
    box = (arb(float(mp.re(ps))) - e, arb(float(mp.re(ps))) + e,
           arb(float(mp.im(ps))) - e, arb(float(mp.im(ps))) + e)
    out = {"PS_mode_over_kappa_minus": [mp.nstr(mp.re(ps) / g.km, 12), mp.nstr(mp.im(ps) / g.km, 12)],
           "box_half_width_over_gamma0": "1/50"}
    ca = certify_l(pt, l, rays, tails, design, nodes=window_nodes, lam_float=window_lam, box=box)
    out["a_window_ground_state"] = {"certified": ca["A1_pass"], "mu_cert": ca["mu_cert"]}
    try:
        cb = certify_l(pt, l, rays, tails, design, box=box)
        out["b_box_design"] = {"certified": cb["A1_pass"], "mu_cert": cb["mu_cert"],
                               "float_lambda_min": cb["float_lambda_min"]}
    except (ArithmeticError, ValueError) as ex:
        # the float margin at the box is <= 0, so the Riccati design cannot start: the checker is
        # NOT run on this branch (S2 review, S4); branch (a) and A2' below are the real controls
        out["b_box_design"] = {"certified": False, "design_failed": str(ex), "checker_run": False}
    out["A2_pass"] = (not out["a_window_ground_state"]["certified"]) and (not out["b_box_design"]["certified"])
    return out


def near_mode_control(pt: Point, l: int, rays, tails, design=DESIGN):
    """A2' (S2 review, section 5): the same checker with a ground state designed close to the l-th
    PS fundamental mode (Re w = Re w_PS, Im w = -0.88 gamma0, where the float margin is still > 0).
    Boxes (real half-width gamma0/50 about Re w_PS):
      IN    : |Im w - Im w_PS| <= gamma0/50          (contains a QNM)   -> must NOT certify
      STRAD : Im w in [-1.05, -0.95] gamma0           (contains the mode) -> must NOT certify
      OUT   : Im w in [-0.90, -0.86] gamma0           (mode-free)         -> reported (may certify)."""
    import mpmath as mp

    from .qnm_float import Geom, Solver
    g = Geom.frozen("PB")
    with mp.workdps(30):
        top = (l + mp.mpf(1) / 2) * g.OmegaP
        ps = Solver(g, K=6).find(l, mp.mpc(top, -g.gamma0))
    g0 = float(pt.gamma0.mid())
    re_ps, im_ps = float(mp.re(ps)), float(mp.im(ps))
    e = g0 / 50
    lam, nodes = design_ground_state(pt, l, complex(re_ps, -0.88 * g0), design)
    muHO = float((2 * (l + arb(1) / 2) * pt.OmegaP * pt.gamma0).mid())
    boxes = {"IN": (re_ps - e, re_ps + e, im_ps - e, im_ps + e),
             "STRAD": (re_ps - e, re_ps + e, -1.05 * g0, -0.95 * g0),
             "OUT": (re_ps - e, re_ps + e, -0.90 * g0, -0.86 * g0)}
    out = {"PS_mode_over_kappa_minus": [mp.nstr(mp.re(ps) / g.km, 12), mp.nstr(mp.im(ps) / g.km, 12)],
           "design_point_Im_over_gamma0": -0.88, "float_lambda_at_design_over_muHO": lam / muHO,
           "ground_state_nodes": [[str(x), str(y)] for x, y in nodes]}
    for name, (a1, a2, b1, b2) in boxes.items():
        c = certify_l(pt, l, rays, tails, design, nodes=nodes, lam_float=lam,
                      box=(arb(a1), arb(a2), arb(b1), arb(b2)))
        mu = float(arb(c["mu_cert"].split()[0].strip("[")).mid())
        out[name] = {"box": [a1, a2, b1, b2], "certified": c["A1_pass"], "mu_cert": c["mu_cert"],
                     "mu_over_muHO": mu / muHO, "cells": c["cells_evaluated"]}
    out["A2prime_pass"] = (not out["IN"]["certified"]) and (not out["STRAD"]["certified"])
    return out


def _overlap(pt: Point, l: int):
    """a_l (P2-2 threshold) <= window top omega_l + 6 gamma0 (certified comparison)."""
    from . import lg_bounds
    f = lg_bounds.CERT_DIR / f"lg_PB_l{l}.json"
    stored = None
    if f.exists():
        stored = json.loads(f.read_text())
    else:
        f = CERT_DIR / f"lg_PB_l{l}.json"
        if f.exists():
            stored = json.loads(f.read_text())
    if stored is None:
        stored = lg_bounds.certify("PB", l)
        stored["note"] = "computed in S2 for the window overlap (same code and criterion as S1)"
        (CERT_DIR / f"lg_PB_l{l}.json").write_text(json.dumps(stored, indent=1) + "\n")
    a_l = fa(Fraction(stored["W_l_exact"]))
    _, _, a2, _ = window(pt, l)
    return {"a_l_exact": stored["W_l_exact"], "a_l_over_kappa_minus": (a_l / pt.km).str(10),
            "window_top_over_kappa_minus": (a2 / pt.km).str(10), "a_l_le_top": bool(a_l < a2),
            "lg_criterion_certified": stored["criterion"]["certified"]}


def fit(ms):
    """Least squares m(l) = m_inf - C/(l + 1/2) and m(l) = m_inf2 - C2/(l + 1/2)^2."""
    ls = sorted(ms)
    out = {}
    for name, pw in (("C_over_l", 1), ("C2_over_l2", 2)):
        x = np.array([1 / (l + 0.5) ** pw for l in ls])
        y = np.array([ms[l] for l in ls])
        Amat = np.vstack([np.ones_like(x), -x]).T
        (m_inf, C), *_ = np.linalg.lstsq(Amat, y, rcond=None)
        resid = float(np.max(np.abs(Amat @ np.array([m_inf, C]) - y)))
        l0 = (C / m_inf) ** (1 / pw) - 0.5 if C > 0 and m_inf > 0 else 0.0
        out[name] = {"m_inf": float(m_inf), "C": float(C), "max_residual": resid,
                     "projected_l0": max(0, int(np.ceil(l0)))}
    return out


def main_run(ls):
    CERT_DIR.mkdir(parents=True, exist_ok=True)
    pt = Point("PB", 256)
    commit = git_commit()
    with precision(PREC):
        geo, rays, tails = geometry_certificate(pt)
        geo["commit"] = commit
        (CERT_DIR / "m3_PB_geometry.json").write_text(json.dumps(geo, indent=1) + "\n")
        print("A0", geo["A0_pass"], geo["timing_s"], flush=True)
        ms, rows = {}, {}
        for l in ls:
            c = certify_l(pt, l, rays, tails)
            nodes = [(Fraction(x), Fraction(y)) for x, y in c["ground_state_nodes"]]
            with Timer() as tc:
                c["A2_control"] = control_l(pt, l, rays, tails, nodes, c["float_lambda_min"])
            c["A2_control"]["timing_s"] = {"wall": round(tc.wall, 1), "cpu": round(tc.cpu, 1)}
            with Timer() as tn:
                c["A2prime_near_mode"] = near_mode_control(pt, l, rays, tails)
            c["A2prime_near_mode"]["timing_s"] = {"wall": round(tn.wall, 1), "cpu": round(tn.cpu, 1)}
            c["A3_mu_cert_le_float"] = bool(float(arb(c["mu_cert"].split()[0].strip("[")).mid())
                                            <= c["float_lambda_min"])
            c["P2_2_overlap"] = _overlap(pt, l)
            c["commit"] = commit
            (CERT_DIR / f"m3_PB_l{l}.json").write_text(json.dumps(c, indent=1) + "\n")
            ms[l] = float(arb(c["m_ratio"].split()[0].strip("[")).mid()) if "+/-" in c["m_ratio"] else float(c["m_ratio"])
            rows[l] = {"A1": c["A1_pass"], "A2": c["A2_control"]["A2_pass"],
                       "A2prime_near_mode": c["A2prime_near_mode"]["A2prime_pass"],
                       "A2prime_OUT_certified": c["A2prime_near_mode"]["OUT"]["certified"],
                       "A3": c["A3_mu_cert_le_float"],
                       "overlap": c["P2_2_overlap"]["a_l_le_top"], "mu_cert": c["mu_cert"], "m": c["m_ratio"],
                       "float_m": c["float_lambda_min"] / float(arb(c["mu_HO"].split()[0].strip("[")).mid())}
            print(l, rows[l], c["timing_s"], flush=True)
        f = fit(ms)
        ok = geo["A0_pass"] and all(r["A1"] and r["A2"] and r["A3"] and r["overlap"] for r in rows.values())
        l0 = f["C_over_l"]["projected_l0"]
        if not rows.get(100, {}).get("A1", False):
            verdict = "NO-GO-for-now (A1 fails at l = 100)"
        elif l0 > 100:
            verdict = "NO-GO-for-now (projected l0 > 100)"
        elif ok:
            verdict = "GO for S3 (barrier-top window certified at l = 20, 50, 100; projected l0 <= 100)"
        else:
            verdict = "INCONCLUSIVE (see failing items)"
        summ = {"point": "PB", "route": "M3", "rows": rows, "fit": f, "verdict": verdict, "commit": commit,
                "preregistration": "notes/s2-preregistration.md",
                "note": "projected_l0 is the formal output of the pre-registered fit. The fit is mis-specified "
                        "(S2 review, section 4: design slack, a ground-state loss ~ 1/l and a spline loss ~ sqrt(l) "
                        "add up to a non-monotone m(l)); it has no predictive content. See notes/s2-report.md section 3."}
        (CERT_DIR / "m3_PB_summary.json").write_text(json.dumps(summ, indent=1) + "\n")
        print(json.dumps(summ, indent=1))


def verify():
    """A4 replay: geometry, every stored window certificate with its stored ground-state nodes,
    and the stored P2-2 overlap certificate for l = 100."""
    from . import lg_bounds
    pt = Point("PB", 256)
    probs = []
    with precision(PREC):
        geo, rays, tails = geometry_certificate(pt)
        if not geo["A0_pass"]:
            probs.append("A0 failed on replay")
        for f in sorted(CERT_DIR.glob("m3_PB_l*.json")):
            st = json.loads(f.read_text())
            nodes = [(Fraction(x), Fraction(y)) for x, y in st["ground_state_nodes"]]
            c = certify_l(pt, st["l"], rays, tails, nodes=nodes, lam_float=st["float_lambda_min"])
            if c["A1_pass"] != st["A1_pass"] or c["mu_cert"] != st["mu_cert"]:
                probs.append(f"{f.name}: stored {st['A1_pass']}/{st['mu_cert']} vs replay {c['A1_pass']}/{c['mu_cert']}")
        f = CERT_DIR / "lg_PB_l100.json"
        if f.exists():
            st = json.loads(f.read_text())
            p1 = Point("PB")
            lg = lg_bounds.LG(p1, 100)
            res = lg.criterion(fmpq(*[int(x) for x in st["W_l_exact"].split("/")]),
                               fmpq(*[int(x) for x in st["matching_point_exact"].split("/")]))
            if res is None or not res["ok"]:
                probs.append("lg_PB_l100: criterion not certified on replay")
    return probs


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "verify":
        pr = verify()
        print("m3 replay:", "OK" if not pr else pr)
        sys.exit(1 if pr else 0)
    ls = [int(x) for x in sys.argv[2:]] if len(sys.argv) > 2 else [20, 50, 100]
    main_run(ls)
