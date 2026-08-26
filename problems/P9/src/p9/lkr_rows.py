"""Single source of truth for the lambda-kappa-rho relaxation: builds every constraint row once,
parameterized by the arithmetic (floats for the solver, Arb balls for the verifier), so that the
row structure and order are identical by construction.

Rows: (coeffs: list[(var, value)], rhs, kind) with kind 'eq' (sum = rhs) or 'le' (sum <= rhs).
The dense SN whitening block and the second-order cone are described separately (see Build).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from .model import Frozen
from .socp2 import MP_BOX

N_TANGENTS = 4


class FloatArith:
    LN10 = math.log(10.0)
    @staticmethod
    def c(x): return float(x)
    @staticmethod
    def exp(x): return math.exp(x)
    @staticmethod
    def log(x): return math.log(x)
    @staticmethod
    def log10(x): return math.log10(x)
    @staticmethod
    def pow10(y): return 10.0 ** y
    @staticmethod
    def sqrt(x): return math.sqrt(x)
    @staticmethod
    def expm1(x): return math.expm1(x)
    @staticmethod
    def absmax(a, b): return max(abs(a), abs(b))
    @staticmethod
    def up(x): return x      # rigorous upper endpoint (identity for floats)
    @staticmethod
    def dn(x): return x
    @staticmethod
    def mid(x): return float(x)


class ArbArith:
    def __init__(self, prec: int = 160):
        from flint import arb, ctx
        ctx.prec = prec
        self.arb = arb
        self.LN10 = arb(10).log()
    def c(self, x): return self.arb(float(x))
    def exp(self, x): return x.exp()
    def log(self, x): return x.log()
    def log10(self, x): return x.log() / self.LN10
    def pow10(self, y): return (y * self.LN10).exp()
    def sqrt(self, x): return x.sqrt()
    def expm1(self, x): return x.exp() - 1
    def absmax(self, a, b): return (abs(a) + abs(b) + abs(abs(a) - abs(b))) / 2
    def up(self, x):
        from .verify import _endpoint
        return self.arb(_endpoint(x, +1))
    def dn(self, x):
        from .verify import _endpoint
        return self.arb(_endpoint(x, -1))
    def mid(self, x): return float(x.mid().str(20, radius=False))


def _I_pm(x, L, sign, ar):
    a = 1 - sign * L
    if abs(ar.mid(a)) < 1e-9:
        return ar.exp(sign * L * x) * x
    return ar.exp(sign * L * x) * (ar.exp(a * x) - 1) / a


def _J_pm(xa, xb, L, sign, ar):
    a = 1 + sign * L
    if abs(ar.mid(a)) < 1e-9:
        return ar.exp(xa) * (xb - xa)
    return ar.exp(xa) * (ar.exp(a * (xb - xa)) - 1) / a


def theta_bounds(x, L, ar):
    em1 = ar.expm1(x)
    return em1 / _I_pm(x, L, +1, ar), em1 / _I_pm(x, L, -1, ar)


def curvature_bound(xa, xb, L, ar):
    """sup |F''| on [xa, xb], F = ln[D_M/(e^x-1)]; see geometry.kappa_second_derivative_bound."""
    if ar.mid(xa) <= 0.0:
        return L / 6 + L * L * (ar.c(0.5) + xb / 6) ** 2
    q = ar.exp(xa) / (ar.exp(xa) - 1)
    th_lo, th_hi = theta_bounds(xb, L, ar)
    dth = ar.absmax(th_lo - 1, th_hi - 1)
    return q * (L * th_hi + dth * (q * (1 + th_hi) - 1))


@dataclass
class Layout:
    N: int; n: int; nbao: int
    lam: np.ndarray; kappa: np.ndarray; Mp: int; ell: np.ndarray
    Ed: np.ndarray; Er: np.ndarray; Es: np.ndarray; yb: np.ndarray; P: np.ndarray; U: np.ndarray
    wb: np.ndarray; ws: np.ndarray; nvar: int
    enodes: list; epos: dict; idx_dm: list; idx_dh: list


@dataclass
class Build:
    layout: Layout
    eq_rows: list = field(default_factory=list)      # (coeffs, rhs)
    le_rows: list = field(default_factory=list)      # (coeffs, rhs)
    # dense SN block: w_s + 5 Wsn ell + (Wsn 1) Mp = Wsn (m - off); handled as matrices in both paths
    sn_rhs: list = None                              # (m - off) per SN (in ar)
    soc_T: object = None                             # sqrt(T) or None
    var_lo: list = None; var_hi: list = None         # bounds on every variable over the feasible set (in ar)
    c_node: list = None


def layout_for(fr: Frozen) -> Layout:
    spc = fr.spec; x = spc.x
    N, n = spc.n_seg, len(fr.sn.m)
    kinds = fr.bao.kind; nbao = len(kinds)
    idx_dm = [r for r, k in enumerate(kinds) if k == "DM_over_rs"]
    idx_dh = [r for r, k in enumerate(kinds) if k == "DH_over_rs"]
    from .geometry import segment_index
    xb = np.log1p(fr.bao.z); kb = segment_index(x, xb)
    enodes = sorted({int(kb[r]) for r in idx_dh} | {int(kb[r]) + 1 for r in idx_dh})
    epos = {node: p for p, node in enumerate(enodes)}
    pos = 0
    def block(m):
        nonlocal pos
        r = np.arange(pos, pos + m); pos += m; return r
    lam = block(N + 1); kappa = block(N + 1); Mp = int(block(1)[0]); ell = block(n)
    Ed = block(N); Er = block(N); Es = block(N); yb = block(len(idx_dm)); P = block(nbao)
    U = block(len(enodes)); wb = block(nbao); ws = block(n)
    return Layout(N, n, nbao, lam, kappa, Mp, ell, Ed, Er, Es, yb, P, U, wb, ws, pos, enodes, epos, idx_dm, idx_dh)


def sandwich_range(lo, hi, ar, K=N_TANGENTS):
    """Range of the sandwiched variable E over y in [lo, hi]: E >= max_k tangent_k(lo) > 0 (tangents are
    increasing in y) and E <= chord(hi) = 10^hi."""
    E_lo = None
    for k in range(1, K + 1):
        yt = lo + (hi - lo) * ar.c(k) / ar.c(K + 1)
        p = ar.pow10(yt)
        v = p * (1 + (lo - yt) * ar.LN10)
        E_lo = v if E_lo is None else (E_lo + v + abs(E_lo - v)) / 2   # max
    return E_lo, ar.pow10(hi)


def sandwich_rows(iE, lin, lo, hi, ar, K=N_TANGENTS):
    """E >= 10^y (tangents), E <= chord(y); y = sum lin coeffs, y in [lo, hi] (bracket rows added by caller)."""
    rows = []
    for k in range(1, K + 1):
        yt = lo + (hi - lo) * ar.c(k) / ar.c(K + 1)
        p = ar.pow10(yt)
        rows.append(([(iE, ar.c(-1))] + [(v, p * ar.LN10 * cf) for v, cf in lin], p * (yt * ar.LN10 - 1)))
    plo, phi = ar.pow10(lo), ar.pow10(hi)
    s = (phi - plo) / (hi - lo)
    rows.append(([(iE, ar.c(1))] + [(v, -s * cf) for v, cf in lin], plo - s * lo))
    return rows


def build(fr: Frozen, br, T, ar) -> Build:
    """br: Brackets3 (floats, certified outward-rounded). T: float or None."""
    lay = layout_for(fr)
    spc = fr.spec
    N, n, nbao = lay.N, lay.n, lay.nbao
    L = ar.c(spc.L)
    x = [ar.c(v) for v in spc.x]
    hs = [x[k + 1] - x[k] for k in range(N)]
    ex = [ar.exp(v) for v in x]
    c_node = [ar.c(0)] + [ar.log10(ex[i] - 1) for i in range(1, N + 1)]
    a = [(ex[k + 1] - (hs[k] + 1) * ex[k]) / hs[k] for k in range(N)]
    b = [((hs[k] - 1) * ex[k + 1] + ex[k]) / hs[k] for k in range(N)]
    lh = [ar.log10(1 + L * hs[k]) for k in range(N)]
    one = ar.c(1)
    B = Build(lay); B.c_node = c_node
    ulo = ar.c(spc.u_box[0]); uhi = ar.c(spc.u_box[1])
    llo_box, lhi_box = ar.dn(ar.log10(ulo)), ar.up(ar.log10(uhi))   # outward: the box must contain the class
    # ---------------- equalities ----------------
    # BAO: w_b + Wb P = Wb d
    Wb = [[ar.c(v) for v in row] for row in fr.Wb]
    d = [ar.c(v) for v in fr.bao.value]
    for r in range(nbao):
        coeffs = [(lay.wb[r], one)] + [(lay.P[s], Wb[r][s]) for s in range(nbao)]
        B.eq_rows.append((coeffs, sum((Wb[r][s] * d[s] for s in range(nbao)), ar.c(0))))
    # SN dense block (matrices): rhs m - off
    off = [ar.c(5) * ar.log10(1 + ar.c(z)) for z in fr.sn.zHEL]
    B.sn_rhs = [ar.c(m) - off[j] for j, m in enumerate(fr.sn.m)]
    # kappa_0 = lambda_0
    B.eq_rows.append(([(lay.kappa[0], one), (lay.lam[0], -one)], ar.c(0)))
    # segment identities
    B.eq_rows.append(([(lay.Ed[0], one), (lay.Es[0], -b[0])], a[0]))
    for i in range(1, N):
        B.eq_rows.append(([(lay.Ed[i], one), (lay.Er[i], -a[i]), (lay.Es[i], -b[i])], one))
    # D_H rows: P_r = (1-t) U_k + t U_{k+1}
    from .geometry import segment_index
    xb_f = np.log1p(fr.bao.z); kb = segment_index(spc.x, xb_f)
    xb = [ar.log(1 + ar.c(z)) for z in fr.bao.z]
    tb = [(xb[r] - x[int(kb[r])]) / hs[int(kb[r])] for r in range(nbao)]
    for r in lay.idx_dh:
        k = int(kb[r]); t = tb[r]
        B.eq_rows.append(([(lay.P[r], one), (lay.U[lay.epos[k]], -(1 - t)), (lay.U[lay.epos[k + 1]], -t)], ar.c(0)))
    # ---------------- inequalities ----------------
    le = B.le_rows
    lam_lo = [ar.absmax(llo_box, llo_box) for _ in range(N + 1)]  # placeholders replaced below
    lam_lo = [llo_box if ar.mid(ar.c(br.lam_lo[i])) < ar.mid(llo_box) else ar.c(br.lam_lo[i]) for i in range(N + 1)]
    lam_hi = [lhi_box if ar.mid(ar.c(br.lam_hi[i])) > ar.mid(lhi_box) else ar.c(br.lam_hi[i]) for i in range(N + 1)]
    for i in range(N + 1):
        le.append(([(lay.lam[i], -one)], -lam_lo[i]))
        le.append(([(lay.lam[i], one)], lam_hi[i]))
    for k in range(N):   # class
        le.append(([(lay.lam[k + 1], one), (lay.lam[k], -one)], lh[k]))
        le.append(([(lay.lam[k + 1], -one), (lay.lam[k], one)], lh[k]))
    # SN interpolation
    xj_f = np.log1p(fr.sn.zHD); kj = segment_index(spc.x, xj_f)
    e_seg = []
    for k in range(N):
        Bk = curvature_bound(x[k], x[k + 1], L, ar)
        e_seg.append(ar.up(hs[k] * hs[k] / 8 * Bk / ar.LN10))
    for j in range(n):
        k = int(kj[j]); xj = ar.log(1 + ar.c(fr.sn.zHD[j])); t = (xj - x[k]) / hs[k]
        cj = ar.log10(ar.c(fr.sn.zHD[j]))          # log10(e^{x_j} - 1) = log10 z_j
        le.append(([(lay.ell[j], one), (lay.kappa[k], -(1 - t)), (lay.kappa[k + 1], -t)], cj + e_seg[k]))
        le.append(([(lay.ell[j], -one), (lay.kappa[k], (1 - t)), (lay.kappa[k + 1], t)], -cj + e_seg[k]))
    # delta / rho / s brackets and sandwiches
    B._ranges = dict(Ed=[None] * N, Er=[None] * N, Es=[None] * N)
    for i in range(N):
        dc = c_node[i + 1] - c_node[i]
        Jm, Jp = _J_pm(x[i], x[i + 1], L, -1, ar), _J_pm(x[i], x[i + 1], L, +1, ar)
        if i == 0:
            r_lo, r_hi = Jm, Jp
        else:
            Im, Ip = _I_pm(x[i], L, -1, ar), _I_pm(x[i], L, +1, ar)
            r_lo, r_hi = 1 + Jm / Ip, 1 + Jp / Im
        lo_d, hi_d = ar.dn(ar.log10(r_lo)), ar.up(ar.log10(r_hi))       # delta_i = kappa_{i+1}-kappa_i+dc in [lo_d, hi_d]
        # bracket rows on kappa_{i+1} - kappa_i
        le.append(([(lay.kappa[i + 1], -one), (lay.kappa[i], one)], -(lo_d - dc)))
        le.append(([(lay.kappa[i + 1], one), (lay.kappa[i], -one)], hi_d - dc))
        # Ed_i = 10^{delta_i} = 10^{dc} * 10^{kappa_{i+1}-kappa_i}; sandwich in y' = kappa_{i+1}-kappa_i on [lo_d-dc, hi_d-dc]
        sc = ar.pow10(dc)
        for coeffs, rhs in sandwich_rows(lay.Ed[i], [(lay.kappa[i + 1], one), (lay.kappa[i], -one)], lo_d - dc, hi_d - dc, ar):
            coeffs = [(v, cf / sc) if v == lay.Ed[i] else (v, cf) for v, cf in coeffs]
            le.append((coeffs, rhs))
        _el, _eh = sandwich_range(lo_d - dc, hi_d - dc, ar); B._ranges["Ed"][i] = (sc * _el, sc * _eh)
        if i >= 1:
            lo_r, hi_r = ar.c(br.rho_lo[i]), ar.c(br.rho_hi[i])
            # rho_i = lambda_i - kappa_i - c_i in [lo_r, hi_r]
            le.append(([(lay.lam[i], -one), (lay.kappa[i], one)], -(lo_r + c_node[i])))
            le.append(([(lay.lam[i], one), (lay.kappa[i], -one)], hi_r + c_node[i]))
            sc = ar.pow10(c_node[i])     # Er = 10^{rho} = 10^{(lam - kappa)} / 10^{c}
            for coeffs, rhs in sandwich_rows(lay.Er[i], [(lay.lam[i], one), (lay.kappa[i], -one)], lo_r + c_node[i], hi_r + c_node[i], ar):
                coeffs = [(v, cf * sc) if v == lay.Er[i] else (v, cf) for v, cf in coeffs]
                le.append((coeffs, rhs))
            _el, _eh = sandwich_range(lo_r + c_node[i], hi_r + c_node[i], ar); B._ranges["Er"][i] = (_el / sc, _eh / sc)
            lin_s = [(lay.lam[i + 1], one), (lay.kappa[i], -one)]
            lo_s, hi_s = lo_r - lh[i], hi_r + lh[i]; shift = c_node[i]
        else:
            lin_s = [(lay.lam[1], one), (lay.lam[0], -one)]
            lo_s, hi_s = -lh[0], lh[0]; shift = ar.c(0)
        le.append(([(v, -cf) for v, cf in lin_s], -(lo_s + shift)))
        le.append(([(v, cf) for v, cf in lin_s], hi_s + shift))
        sc = ar.pow10(shift)
        for coeffs, rhs in sandwich_rows(lay.Es[i], lin_s, lo_s + shift, hi_s + shift, ar):
            coeffs = [(v, cf * sc) if v == lay.Es[i] else (v, cf) for v, cf in coeffs]
            le.append((coeffs, rhs))
        _el, _eh = sandwich_range(lo_s + shift, hi_s + shift, ar); B._ranges["Es"][i] = (_el / sc, _eh / sc)
    # BAO D_M rows
    yb_lo = [ar.c(v) for v in br.yb_lo]; yb_hi = [ar.c(v) for v in br.yb_hi]
    for p, r in enumerate(lay.idx_dm):
        k = int(kb[r]); t = tb[r]; cbr = ar.log10(ar.c(fr.bao.z[r])); e = e_seg[k]
        le.append(([(lay.yb[p], one), (lay.kappa[k], -(1 - t)), (lay.kappa[k + 1], -t)], cbr + e))
        le.append(([(lay.yb[p], -one), (lay.kappa[k], (1 - t)), (lay.kappa[k + 1], t)], -cbr + e))
        le.append(([(lay.yb[p], -one)], -yb_lo[p]))
        le.append(([(lay.yb[p], one)], yb_hi[p]))
        for coeffs, rhs in sandwich_rows(lay.P[r], [(lay.yb[p], one)], yb_lo[p], yb_hi[p], ar):
            le.append((coeffs, rhs))
    # U nodes
    for node in lay.enodes:
        lo_l, hi_l = lam_lo[node], lam_hi[node]
        for coeffs, rhs in sandwich_rows(lay.U[lay.epos[node]], [(lay.lam[node], one)], lo_l, hi_l, ar):
            le.append((coeffs, rhs))
    # Mp box
    le.append(([(lay.Mp, -one)], -ar.c(MP_BOX[0])))
    le.append(([(lay.Mp, one)], ar.c(MP_BOX[1])))
    B.soc_T = ar.sqrt(ar.c(T)) if T is not None else None
    # ---------------- variable bounds (for residual absorption) ----------------
    lo = [None] * lay.nvar; hi = [None] * lay.nvar
    big_neg = -ar.c(1e6); big = ar.c(1e6)
    for i in range(N + 1):
        lo[lay.lam[i]] = lam_lo[i]; hi[lay.lam[i]] = lam_hi[i]
        # kappa_i = lambda_i - rho_i - c_i, rho in bracket (rho_0 = 0)
        rlo = ar.c(br.rho_lo[i]) if i > 0 else ar.c(0); rhi = ar.c(br.rho_hi[i]) if i > 0 else ar.c(0)
        lo[lay.kappa[i]] = lam_lo[i] - rhi - c_node[i]; hi[lay.kappa[i]] = lam_hi[i] - rlo - c_node[i]
    lo[lay.Mp] = ar.c(MP_BOX[0]); hi[lay.Mp] = ar.c(MP_BOX[1])
    for j in range(n):
        k = int(kj[j]); cj = ar.log10(ar.c(fr.sn.zHD[j]))
        klo = ar.absmax(lo[lay.kappa[k]], lo[lay.kappa[k + 1]]); khi = ar.absmax(hi[lay.kappa[k]], hi[lay.kappa[k + 1]])
        m_ = ar.absmax(klo, khi) + abs(cj) + e_seg[k]
        lo[lay.ell[j]] = -m_; hi[lay.ell[j]] = m_
    for i in range(N):
        blo, bhi = B._ranges["Ed"][i]; lo[lay.Ed[i]] = blo; hi[lay.Ed[i]] = bhi
        blo, bhi = B._ranges["Es"][i]; lo[lay.Es[i]] = blo; hi[lay.Es[i]] = bhi
        if i >= 1:
            blo, bhi = B._ranges["Er"][i]; lo[lay.Er[i]] = blo; hi[lay.Er[i]] = bhi
        else:
            lo[lay.Er[0]] = big_neg; hi[lay.Er[0]] = big   # unused variable (no rows): rho stays exactly 0
    for p in range(len(lay.idx_dm)):
        lo[lay.yb[p]] = yb_lo[p]; hi[lay.yb[p]] = yb_hi[p]
    U_rng = {}
    for node in lay.enodes:
        U_rng[node] = sandwich_range(lam_lo[node], lam_hi[node], ar)
        lo[lay.U[lay.epos[node]]], hi[lay.U[lay.epos[node]]] = U_rng[node]
    for p, r in enumerate(lay.idx_dm):
        lo[lay.P[r]], hi[lay.P[r]] = sandwich_range(yb_lo[p], yb_hi[p], ar)
    for r in lay.idx_dh:
        k = int(kb[r])
        lo[lay.P[r]] = (U_rng[k][0] + U_rng[k + 1][0] - abs(U_rng[k][0] - U_rng[k + 1][0])) / 2   # min
        hi[lay.P[r]] = (U_rng[k][1] + U_rng[k + 1][1] + abs(U_rng[k][1] - U_rng[k + 1][1])) / 2   # max
    sT = B.soc_T if B.soc_T is not None else big
    for r in range(nbao):
        lo[lay.wb[r]] = -sT; hi[lay.wb[r]] = sT
    for j in range(n):
        lo[lay.ws[j]] = -sT; hi[lay.ws[j]] = sT
    B.var_lo, B.var_hi = lo, hi
    return B
