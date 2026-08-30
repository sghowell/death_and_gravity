"""Theorem B, Stage 2: the certified background along which the linear systems are propagated.

Two regions, one data structure (``StepData`` per step x_k -> x_k - h_k):
  * tube region [x_c, x0] (default [-3, -0.05]): the A3 integrator (``tmint``) is run for the
    V0 box [c - w, c + w] containing the certified EC value V0*; per step the *point* Taylor
    coefficients z_K of the scaled background through the reference midpoint m_k (exact) and the
    enlargement eps_z = tail_u(h) + rk (sup|Y_K| + bound) such that every solution of the interval
    set (hence the EC solution) satisfies |z(s) - z_K(s)| <= eps_z on [-h, 0]  (rk = weighted
    radius of the sets before the step, bound = tmint's Groenwall bound of the step Jacobian, both
    captured from the integrator as arb; Y_K the 4D fundamental-matrix series: mean value along
    the tube, exactly the refinement argument of tmint.jacobian_step);
  * centre region [x_d, x_c]: the certified A2 centre series gives z(x_k) as balls; the Taylor
    coefficients through the ball are balls containing the true ones, eps_z = the A1/A2 Banach
    tail of that series at h (no enlargement needed).
Derived per step (kappa-independent): sub-boxes Z_j (nsub exact pieces of [-h, 0]) of the tube,
z' = f(Z_j), L_bg = sup|Df|, Pinv_bg, the residual sup of the truncation, eps_dz = L_bg eps_z +
Pinv_bg Rsup, coefficient sums a = sup|z_K|, sup|z_K'|; and, per linear system, the polynomial
coefficient matrices along z_K (truncated at K, tails), box values and increment majorants.
The tube can be saved/loaded exactly (dyadic mantissa/exponent of every ball).
"""
from __future__ import annotations

import json
import time
from dataclasses import dataclass, field

from flint import arb, arb_mat, fmpz

from . import matching, recursion, tmint
from . import shootsys as ss
from .arbseries import Series, abs_upper, precision, to_arb
from .linstep import norm_inf as _lnorm
from .tailbound import certify_tail, norm_inf

BLOCKS4 = [(0, 1), (1, 3), (3, 4)]
_ID4 = [[arb(int(i == j)) for j in range(4)] for i in range(4)]


@dataclass
class StepData:
    x: arb
    h: float
    co: list                      # z_K: list over i <= K of 4 arb (n, w, v, T)
    eps_z: arb
    region: str = "tube"
    derived: dict = None
    cache: dict = field(default_factory=dict)


class _Hooked(tmint.Integrator):
    def jacobian_step(self, co, tails, h, hull_radii):
        out = super().jacobian_step(co, tails, h, hull_radii)
        self._cap = (out[1], hull_radii, h, tails, out[2])
        return out


def _ysup(sys4, co, K, h):
    Y = ss.variational_coefficients(sys4, co, K, BLOCKS4)
    tot, hp = arb_mat(4, 4), arb(1)
    for Yi in Y:
        tot += arb_mat([[abs_upper(Yi[r, c]) for c in range(4)] for r in range(4)]) * hp
        hp *= h
    return norm_inf(tot)


class Tube:
    def __init__(self, steps, K, info=None):
        self.steps, self.K, self.info = steps, K, info or {}
        self.sys4 = ss.shoot_system()
        self.eqs4 = ss.regular_level_equations(4)

    @classmethod
    def build(cls, c, w, x0=-0.05, x_c=-3.0, K=28, hmax=0.02, prec=384, verbose=False):
        t0 = time.time()
        with precision(prec):
            st = matching.sonic_initial_state(c, w, x0=x0)
            it = _Hooked(K=K, hmax=hmax, verbose=verbose)
            steps = []
            while float(st.x) > x_c + 1e-15:
                x, m = st.x, [arb(v) for v in st.m[:4]]
                it.step(st, x_c)
                bound, hull_radii, h, tails, Sc = it._cap
                rk = tmint._amax(arb(Sc[j]) * hr[j] for hr in hull_radii for j in range(7))
                co = ss.taylor_coefficients(it.sys4, m, K, blocks=BLOCKS4)
                rho = abs_upper(rk * (_ysup(it.sys4, co, K, arb(h)) + bound))
                steps.append(StepData(x, h, co, abs_upper(tails[0] + rho)))
            info = dict(V0=str(c), w=w, x0=x0, x_c=x_c, time=time.time() - t0, prec=prec,
                        iv_width=max(st.iv.widths()[:4]), pt_width=max(st.pt.widths()[:4]))
        return cls(steps, K, info)

    def extend_centre(self, ce, x_d, hmax=0.1, prec=384):
        """Append steps from the last x down to x_d along the certified centre expansion ``ce``."""
        with precision(prec):
            x = self.steps[-1].x - arb(self.steps[-1].h) if self.steps else arb(self.info["x_c"])
            while float(x) > x_d + 1e-15:
                z = list(ce.eval(x)) + [x.exp()]
                co = ss.taylor_coefficients(self.sys4, z, self.K, blocks=BLOCKS4)
                D, E = recursion.structure_matrices(self.sys4, self.eqs4, co, dm=1)
                cert = certify_tail(self.sys4, self.eqs4, co, D, E)
                if not cert.ok:
                    raise RuntimeError(f"centre-region tail certificate failed at x={float(x)}")
                h = min(hmax, 0.75 * float(cert.nu), float(x - x_d))
                self.steps.append(StepData(x, h, co, cert.tail_bound(arb(h)), region="centre"))
                x = x - arb(h)
        return self

    # -- kappa-independent derived data ------------------------------------------------
    def derive(self, sd, nsub=8):
        if sd.derived is not None:
            return sd.derived
        h, K = arb(sd.h), len(sd.co) - 1
        boxes, L, Pinv = [], arb(0), arb(0)
        for j in range(nsub):
            sj = (-h * (2 * j + 1)) / (2 * nsub) + arb(0, h / (2 * nsub))
            Z = [zi + arb(0, sd.eps_z) for zi in ss.horner_vec(sd.co, sj)]
            f, P = ss.rhs_enclosure(self.sys4, Z, BLOCKS4)
            L = L.max(norm_inf(ss.jacobian_enclosure(self.sys4, Z, f, P, BLOCKS4)))
            Pinv = Pinv.max(norm_inf(ss.block_solve(P, _ID4, BLOCKS4)))
            boxes.append((Z, f[:3]))
        zs = recursion.series_from_coefs(sd.co, 4, cap=None, extra_zero=False)
        R = self.sys4.residual(zs)
        Rsup = tmint._amax(sum((abs_upper(Rr[k]) * h**k for k in range(len(Rr))), arb(0)) for Rr in R)
        eps_dz = abs_upper(L * sd.eps_z + Pinv * Rsup)
        a = [sum((abs_upper(sd.co[i][k]) * h**i for i in range(K + 1)), arb(0)) for k in range(4)]
        a += [sum((abs_upper(sd.co[i][k]) * i * h ** (i - 1) for i in range(1, K + 1)), arb(0)) for k in range(3)]
        dzs = [s.deriv() for s in zs[:3]]
        amaj = [Series([abs_upper(c) for c in s.coeffs()]) for s in zs + dzs]
        sd.derived = dict(boxes=boxes, L=L, Pinv=Pinv, Rsup=Rsup, eps_dz=eps_dz, a=a, zs=zs, dzs=dzs,
                          amaj=amaj, eps=[sd.eps_z] * 4 + [eps_dz] * 3)
        return sd.derived

    def system_data(self, sd, S, nsub=8):
        """Per linear system S: (Pser, Gser) truncated s-polynomials (lists over j of d x d Series
        coefficient lists), their tails, box values on the sub-boxes, increment majorants."""
        if S.name in sd.cache:
            return sd.cache[S.name]
        dv, K, h = self.derive(sd, nsub), len(sd.co) - 1, arb(sd.h)
        Pser, Gser = S.series_matrices(dv["zs"], dv["dzs"], cap=K + 1)      # exact orders <= K
        Ptail, Gtail = S.majorant_tails(dv["amaj"], dv["a"], K, h)         # orders > K (majorant)
        P = [([[m.coeffs(K + 1) for m in row] for row in M], _lnorm(T)) for M, T in zip(Pser, Ptail)]
        G = [([[m.coeffs(K + 1) for m in row] for row in M], _lnorm(T)) for M, T in zip(Gser, Gtail)]
        bx = [S.box_matrices(Z, dZ) for Z, dZ in dv["boxes"]]
        Pbox = [[b[0][j] for b in bx] for j in range(len(S.P))]            # [j][sub-box]
        Gbox = [[b[1][j] for b in bx] for j in range(len(S.G))]
        incP, incG = S.increments(dv["a"], dv["eps"])
        sd.cache[S.name] = dict(P=P, G=G, Pbox=Pbox, Gbox=Gbox, incP=incP, incG=incG)
        return sd.cache[S.name]

    # -- exact serialisation --------------------------------------------------------------
    def save(self, path):
        enc = lambda a: [str(v) for me in (a.mid().man_exp(), a.rad().man_exp()) for v in me]
        data = dict(K=self.K, info=self.info, steps=[dict(x=enc(s.x), h=s.h, eps=enc(s.eps_z), region=s.region,
                                                          co=[[enc(v) for v in ci] for ci in s.co]) for s in self.steps])
        json.dump(data, open(path, "w"))

    @classmethod
    def load(cls, path, prec=384):
        data = json.load(open(path))
        with precision(prec):
            dec = lambda l: arb(arb(fmpz(l[0])) * arb(2) ** int(l[1]), arb(fmpz(l[2])) * arb(2) ** int(l[3]))
            steps = [StepData(dec(s["x"]), s["h"], [[dec(v) for v in ci] for ci in s["co"]], dec(s["eps"]), s["region"])
                     for s in data["steps"]]
        return cls(steps, data["K"], data["info"])
