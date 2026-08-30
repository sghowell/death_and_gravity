"""Theorem B, Stage 3 (a): rigorous winding number of the matching function E along the boundary
of the rectangle R = [0, 15] x [-14, 14] (argument principle; analyticity of E inside R is
established in ``analyticity.py``).

E is Stage 2's centre determinant, here evaluated at the tube end x_d (= -3, the end of the A3
tube; ``notes/s2-theorem-b.md`` section 3 records why the zero set and the winding number do not
depend on x_d):  E(kappa) = e^{3 x_d} det [r_1, r_2, p~(x_d)]_{(A^, n~, v~)}  with the 4D solution
p~ of the sonic data (``linsonic4``) propagated in 4D (``lintaylor``) and the regular family r_1,
r_2 at the centre (``lincentre``).  Every ingredient is a *Taylor model in delta = kappa - kappa_c*
on a segment |delta| <= r of the contour:
  sonic data       ``LinSonicExpansion4.delta_model``  (point run + box-run remainder);
  propagation      ``lintaylor.propagate_tm`` (fundamental matrix) composed by ``tm_apply``;
  regular family   ``RegularFamilyTM`` (the lincentre recursion in delta-polynomial arithmetic,
                   remainder from the box-based run and the box tail certificate);
  E                scalar Taylor-model arithmetic (``TM``) of the 3x3 determinant.
On each segment E(delta) in sum_k c_k delta^k + box(rem); the segment is split into sub-segments on
which the enclosure ball must exclude 0, and the argument increment between consecutive endpoint
values is rigorous because the whole sub-segment image lies in an open half-plane.  Summing the
increments around the contour gives 2 pi (winding number) as an arb ball containing one integer.
Segments are processed independently (multiprocessing, spawn); the contour is walked adaptively
(a segment whose model is not tight enough is halved).
"""
from __future__ import annotations

from flint import acb, acb_mat, arb

from . import lincentre, lintaylor
from .arbseries import abs_upper, precision, to_arb
from .linsys import kappa_box, to_acb


# ---------------------------------------------------------------------------------------------
# scalar Taylor models in delta on |delta| <= r
# ---------------------------------------------------------------------------------------------
class TM:
    """sum_{k<=m} c_k delta^k + R,  |R| <= rem  on |delta| <= r  (c_k acb balls, rem arb)."""

    def __init__(self, coefs, rem, r):
        self.c, self.rem, self.r = list(coefs), arb(rem), to_arb(r)
        self.m = len(coefs) - 1

    def _rp(self, n):
        return [self.r**k for k in range(n)]

    def __add__(self, o):
        return TM([a + b for a, b in zip(self.c, o.c)], self.rem + o.rem, self.r)

    def __neg__(self):
        return TM([-a for a in self.c], self.rem, self.r)

    def __sub__(self, o):
        return self + (-o)

    def __mul__(self, o):
        m, rp = self.m, self._rp(2 * self.m + 1)
        c = [sum((self.c[k] * o.c[j - k] for k in range(j + 1)), acb(0)) for j in range(m + 1)]
        trunc = sum((rp[j] * abs_upper(sum((self.c[k] * o.c[j - k] for k in range(j - m, m + 1)), acb(0)))
                     for j in range(m + 1, 2 * m + 1)), arb(0))
        sa = sum((abs_upper(x) * rp[k] for k, x in enumerate(self.c)), arb(0))
        sb = sum((abs_upper(x) * rp[k] for k, x in enumerate(o.c)), arb(0))
        return TM(c, trunc + sa * o.rem + sb * self.rem + self.rem * o.rem, self.r)

    def scale(self, a):
        return TM([x * a for x in self.c], self.rem * abs_upper(a), self.r)

    def eval(self, delta):
        """Ball enclosing the model at delta (acb, possibly a box with |delta| <= r)."""
        acc = acb(0)
        for k in reversed(range(self.m + 1)):
            acc = acc * delta + self.c[k]
        return acc + acb(arb(0, self.rem), arb(0, self.rem))


def det3_tm(cols, rows=(0, 1, 3)):
    """Taylor model of det[cols] restricted to ``rows`` (cols: 3 lists of 4 TMs)."""
    a = [[cols[c][r] for c in range(3)] for r in rows]
    return (a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
            + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]))


# ---------------------------------------------------------------------------------------------
# the regular family at the centre as a Taylor model in delta
# ---------------------------------------------------------------------------------------------
class RegularFamilyTM:
    """r_i(t; kappa_c + delta), i = 1, 2, as delta-polynomials of degree m (point run, exact
    recursion in delta-polynomial arithmetic) with the (m+1)-th coefficient from the run with the
    box kappa_c + [-w, w] + i[-w, w] as base point (encloses r_i^{(m+1)}/(m+1)! on the box) and the
    tail from ``lincentre.RegularFamily(S, ce, box).certify(nu)``."""

    def __init__(self, S, ce, kappa_c, width, m, K=50, nu=0.06, prec=256):
        self.S, self.ce, self.kc, self.w, self.m, self.K = S, ce, to_acb(kappa_c), float(width), m, K
        with precision(prec):
            Pser, Gser = lincentre.coefficient_series(S, ce)
            self.deg = max(len(x) for M in Pser + Gser for row in M for x in row) - 1
            ent = lambda M, k: acb_mat([[acb(M[r][c][k]) for c in range(S.d)] for r in range(S.d)])   # noqa: E731
            self.P = [ent(Pser[0], k) for k in range(self.deg + 1)]
            self.G0 = [ent(Gser[0], k) for k in range(self.deg + 1)]
            self.G1 = [ent(Gser[1], k) for k in range(self.deg + 1)]
            self.rows = [r for r in range(S.d) if r not in lincentre._rows_vanish_at_centre(S)]
            self.point = [self._recur(self.kc, e, m + 2) for e in range(2)]
            self.boxrun = [self._recur(kappa_box(self.kc, self.w), e, m + 2) for e in range(2)]
            self.rf = lincentre.RegularFamily(S, ce, kappa_box(self.kc, self.w), K=K, prec=prec)
            ok, self.eps, self.cdet = self.rf.certify(nu, prec=prec)
            if not ok:
                raise RuntimeError(f"centre certificate failed over the box: {self.cdet}")
            self.nu = arb(nu)

    def _recur(self, base, e, dm):
        d, K = self.S.d, self.K
        Gc = [self.G0[k] + base * self.G1[k] for k in range(self.deg + 1)]
        solved, free = [0, 1], [2, 3]
        M = [acb_mat([[G[r, c] for c in solved] for r in self.rows]) for G in (Gc[0], self.G1[0])]
        rhs = [acb_mat([[-G[r, free[e]]] for r in self.rows]) for G in (Gc[0], self.G1[0])]
        p0 = []
        for j in range(dm):
            b = (rhs[j] if j < 2 else acb_mat(2, 1)) - (M[1] * p0[j - 1] if j else 0)
            p0.append(M[0].solve(b))
        p = [[acb_mat([[p0[j][0, 0]], [p0[j][1, 0]], [acb(int(e == 0 and j == 0))], [acb(int(e == 1 and j == 0))]])
              for j in range(dm)]]
        for n in range(1, K + 1):
            Mn = (self.P[0] * n - Gc[0])
            pn = []
            for j in range(dm):
                b = acb_mat(d, 1)
                for k in range(1, min(n, self.deg) + 1):
                    b += (Gc[k] - self.P[k] * (n - k)) * p[n - k][j]
                    if j:
                        b += self.G1[k] * p[n - k][j - 1]
                if j:
                    b += self.G1[0] * pn[j - 1]
                pn.append(Mn.solve(b))
            p.append(pn)
        return p

    def models(self, x, r):
        """[r_1, r_2] as lists of 4 TMs at t = e^x, valid for |delta| <= r <= width."""
        assert r <= self.w
        t = to_arb(x).exp()
        out = []
        for e in range(2):
            tail = abs_upper(self.eps[e] * (t / self.nu).abs_upper() ** (self.K + 1))
            tms = []
            for i in range(self.S.d):
                coefs = []
                for k in range(self.m + 1):
                    acc = acb(0)
                    for n in reversed(range(self.K + 1)):
                        acc = acc * t + self.point[e][n][k][i, 0]
                    coefs.append(acc)
                top = sum((abs_upper(self.boxrun[e][n][self.m + 1][i, 0]) * t**n for n in range(self.K + 1)), arb(0))
                tms.append(TM(coefs, top * to_arb(r) ** (self.m + 1) + tail, r))
            out.append(tms)
        return out


# ---------------------------------------------------------------------------------------------
# one contour segment: E as a Taylor model, rigorous argument increment
# ---------------------------------------------------------------------------------------------
V0_EC, W_V0 = "0.1124394013880983", 1e-16          # certified V0* enclosure (A3): c* +/- 1e-16
A_STAR = "-0.2123656467659762832750918714540807905889"
MU_STAR = "8.901323275379966931515526907200000000000"
KAPPA1, KGAUGE = "2.8105525488", "0.3556992037"


def certified_centre(K=30):
    """The certified A2 centre expansion at (a*, mu*) (Stage 2 / A3 balls)."""
    from . import centre
    a = arb(A_STAR) + arb("7.031e-12") + arb(0, 8.87e-16)
    mu = arb(MU_STAR) + arb("2.3757e-9") + arb(0, 4.61e-14)
    ce = centre.centre_expansion(mu * (2 * a).exp(), nhat=(-a).exp(), K=K)
    if not ce.certify().ok:
        raise RuntimeError("centre certificate failed")
    return ce


class Context:
    """Everything kappa-independent: tube (A3, to x_d), 4D system, box background, centre series."""

    def __init__(self, tube, x0=-0.05, K_sonic=40, prec=256):
        from . import linmatch, linscaled
        self.tube, self.x0, self.K_sonic, self.prec = tube, x0, K_sonic, prec
        with precision(prec):
            self.S4 = linscaled.full_system()
            self.bg = linmatch.box_background(V0_EC, W_V0, K=K_sonic + 1)
            self.ce = certified_centre()
            last = tube.steps[-1]
            self.x_d = last.x - arb(last.h)
            T0 = to_arb(x0).exp()
            self.scale = [1 / (T0 * T0), T0, 1 / (T0 * T0), 1 / T0]        # p~ = Lambda_4^{-1} p

    def sonic_model(self, kc, w, m, tol=1e-6):
        """Taylor model of p~(x0; kc + delta), |delta| <= w, or None if not tight enough."""
        from . import linsonic4
        ex = linsonic4.linear_sonic_expansion4(self.bg, kc, width=w, m=m, K=self.K_sonic)
        cert = ex.certify()
        if not (cert.ok and cert.nu > abs(self.x0)):
            return None
        co, rem = ex.delta_model(self.x0, w)
        coefs = [[c[i] * self.scale[i] for i in range(4)] for c in co]
        rem0 = [abs_upper(rem[i] * self.scale[i]) for i in range(4)]
        mag = max(float(abs_upper(c)) for c in coefs[0])
        if max(float(x) for x in rem0) > tol * mag:
            return None
        return coefs, rem0, float(cert.nu)

    def E_model(self, kc, w, m, verbose=False):
        """TM of E on |delta| <= w around kc (None if the sonic model is not tight enough)."""
        with precision(self.prec):
            kc = to_acb(kc)
            sm = self.sonic_model(kc, w, m)
            if sm is None:
                return None, dict(reason="sonic")
            coefs, rem0, nu = sm
            Phi, Rcols, log = lintaylor.propagate_tm(self.tube, self.S4, kc, w, m, prec=self.prec, verbose=verbose)
            cf, rm = lintaylor.tm_apply(Phi, Rcols, coefs, rem0, w)
            rtm = RegularFamilyTM(self.S4, self.ce, kc, w, m, prec=self.prec).models(self.x_d, w)
            ptm = [TM([cf[k][i] for k in range(m + 1)], rm[i], w) for i in range(4)]
            E = det3_tm([rtm[0], rtm[1], ptm]).scale((3 * self.x_d).exp())
            info = dict(nu=nu, t_prop=log[-1]["time"], width=log[-2]["width"], rem_prop=max(float(x) for x in rm),
                        mag=max(float(abs_upper(z)) for z in cf[0]))
            return E, info


def segment_increment(E, w, direction, nsub):
    """Rigorous change of arg E along delta = direction * s, s in [-w, w], split into nsub pieces.
    On piece j the enclosure B_j must satisfy Re(B_j e^{-i theta_j}) > 0 (0 excluded, the continuous
    argument on the piece is the principal argument relative to theta_j).  The increments are
    telescoped: sum_j [b_j - a_j] = b_{n-1} - a_0 + sum_j [b_j - a_{j+1}], and b_j - a_{j+1} (two
    arguments of the same value E(delta_{j+1}) relative to theta_j, theta_{j+1}) equals
    arg(rot_j) - arg(rot_{j+1}) + 2 pi l_j exactly, with l_j the unique integer compatible with the ball
    (radius < pi required).  The result's radius is that of two endpoint arguments only.
    Returns (total arb, min|E| float, endpoint values) or None if some piece fails."""
    dr = to_acb(direction)
    minE, ends, a, b, rots = None, [], [], [], []
    for j in range(nsub):
        sa, sb = -w + 2 * w * j / nsub, -w + 2 * w * (j + 1) / nsub
        piece = dr * arb((sa + sb) / 2, (sb - sa) / 2)
        B = E.eval(piece)
        z = complex(float(B.real.mid()), float(B.imag.mid()))
        if z == 0:
            return None
        rot = acb(z.conjugate().real / abs(z), z.conjugate().imag / abs(z))
        Br = (B * rot).real                     # Re(E e^{-i theta}) > 0 on the whole piece: 0 excluded, arg in (-pi/2, pi/2)
        if not (Br > 0):
            return None
        Ea, Eb = E.eval(dr * arb(sa)), E.eval(dr * arb(sb))
        ra, rb = Ea * rot, Eb * rot
        if not (ra.real > 0 and rb.real > 0):
            return None
        a.append(ra.arg())
        b.append(rb.arg())
        rots.append(rot)
        lo = float(Br.mid()) - float(Br.rad())                    # |E| >= Re(E e^{-i theta}) on the piece
        minE = lo if minE is None else min(minE, lo)
        ends.append((sa, complex(float(Ea.real.mid()), float(Ea.imag.mid()))))
    ends.append((w, complex(float(Eb.real.mid()), float(Eb.imag.mid()))))
    twopi = 2 * arb.pi()
    total = b[-1] - a[0]
    for j in range(nsub - 1):
        diff = b[j] - a[j + 1]                                  # ball around arg(rot_j) - arg(rot_{j+1}) + 2 pi l
        if not (diff.rad() < 3):
            return None
        base = rots[j].arg() - rots[j + 1].arg()
        l = int(round(float((diff - base).mid() / twopi.mid())))
        if not (diff - base - twopi * l).contains(arb(0)):
            return None
        total += base + twopi * l
    return total, minE, ends


# ---------------------------------------------------------------------------------------------
# contour walk (worker side) and the parallel driver
# ---------------------------------------------------------------------------------------------
POLE = -0.099                                   # S1: E has a pole at kappa = -0.0990 (A_p(0) = 1 normalisation)
_CTX = None


def enc(a):
    """Exact dyadic encoding of an arb (mid, rad) for JSON."""
    return [str(v) for me in (a.mid().man_exp(), a.rad().man_exp()) for v in me]


def dec(l):
    from flint import fmpz
    return arb(arb(fmpz(l[0])) * arb(2) ** int(l[1]), arb(fmpz(l[2])) * arb(2) ** int(l[3]))


def _init(tube_path, prec):
    global _CTX
    from . import lintube
    _CTX = Context(lintube.Tube.load(tube_path), prec=prec)


def _wcap(kc, wmax):
    """Dyadic half-width <= wmax and <= dist(kc, pole)/10 (the sonic model near the pole)."""
    import math
    w = min(wmax, abs(kc - POLE) / 10)
    return 2.0 ** math.floor(math.log2(w))


def walk(ctx, corner, direction, s_a, s_b, m, wmax=0.25, nsub=8, verbose=False, tol_inc=2e-3):
    """Adaptive tiling of the contour piece corner + direction*s, s in [s_a, s_b]; per tile the
    rigorous arg increment of E (exactly encoded; a tile is halved until its increment ball has radius
    < tol_inc), tightness diagnostics."""
    import time
    out, pos = [], float(s_a)
    while pos < s_b - 1e-12:
        w = min(_wcap(corner + direction * pos, wmax), (s_b - pos) / 2)
        while True:
            t0 = time.time()
            kc = corner + direction * (pos + w)
            E, info = ctx.E_model(acb(kc.real, kc.imag), w, m)
            res = None
            if E is not None:
                for ns in (nsub, 4 * nsub, 16 * nsub):
                    res = segment_increment(E, w, acb(direction.real, direction.imag), ns)
                    if res is not None:
                        break
            if res is not None and not (res[0].rad() < tol_inc):
                res = None                                              # too loose: halve the tile
            if res is not None:
                inc, minE, ends = res
                rec = dict(s0=pos, s1=pos + 2 * w, kc=[kc.real, kc.imag], w=w, m=m, nsub=ns, inc=enc(inc),
                           inc_float=[float(inc.mid()), float(inc.rad())], minE=minE, time=time.time() - t0,
                           ends=[[s, v.real, v.imag] for s, v in ends], **{k: v for k, v in info.items()})
                out.append(rec)
                if verbose:
                    print(f"  seg kc={kc:.4g} w={w} inc={rec['inc_float']} minE={minE:.3g} t={rec['time']:.0f}s", flush=True)
                pos += 2 * w
                break
            w /= 2
            if w < 1e-5:
                raise RuntimeError(f"no tight model at kc={kc}")
    return out


def _task(args):
    corner, direction, s_a, s_b, m, wmax, nsub = args
    return walk(_CTX, complex(*corner), complex(*direction), s_a, s_b, m, wmax, nsub)


def sides_of(rect):
    a, b, c, d = rect
    return [((a, c), (1, 0), b - a), ((b, c), (0, 1), d - c), ((b, d), (-1, 0), b - a), ((a, d), (0, -1), d - c)]


def winding_number(tube_path, rect=(0.0, 15.0, -14.0, 14.0), m=8, wmax=0.25, chunk=0.5, workers=8, prec=256,
                   nsub=8, out=None, verbose=True):
    """Rigorous winding number of E around the rectangle (counterclockwise): dict with the arb ball
    of the total argument change / 2 pi, the certified integer, and all segment records."""
    import json
    import multiprocessing as mp
    import time
    t0 = time.time()
    tasks = []
    for corner, direction, length in sides_of(rect):
        s = 0.0
        while s < length - 1e-12:
            tasks.append((corner, direction, s, min(s + chunk, length), m, wmax, nsub))
            s = min(s + chunk, length)
    with mp.get_context("spawn").Pool(workers, initializer=_init, initargs=(tube_path, prec)) as pool:
        recs = []
        for i, r in enumerate(pool.imap_unordered(_task, tasks)):
            recs.extend(r)
            if verbose:
                print(f"[{i + 1}/{len(tasks)}] {len(r)} tiles, {time.time() - t0:.0f}s", flush=True)
    with precision(prec):
        total = sum((dec(r["inc"]) for r in recs), arb(0))
        wn = total / (2 * arb.pi())
        N = int(round(float(wn.mid())))
        ok = bool((wn - N).abs_upper() < 0.5)
    result = dict(rect=rect, m=m, wmax=wmax, winding=N, certified=ok, winding_ball=[float(wn.mid()), float(wn.rad())],
                  total_arg=enc(total), n_tiles=len(recs), min_abs_E=min(r["minE"] for r in recs),
                  time=time.time() - t0, segments=sorted(recs, key=lambda r: (r["kc"][1] if r["kc"][0] in (rect[0], rect[1]) else r["kc"][0])))
    if out:
        json.dump(result, open(out, "w"))
    return result
