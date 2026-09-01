"""Theorem C, S4-2: certified tube sign certificate promoting the P4-3 monotonicity
diagnostics (ledger P4-3, notes/numerics-report.md section D) to CERTIFIED.

Claims certified on the closed sound cone (-inf, 0]:
  (M1) v_rel' > 0            (v_rel = (1+NV)/(N+V), HM01's |V_z|; increases to 1/sqrt3)
  (M2) rho^' < 0             (rho^ = 4 pi rho t^2 = W e^{-2x}/A, decreasing)
  (M3) 0 < v_rel <= 1/sqrt3  (equality only at x = 0), hence
       w = 1/3 - v_rel^2 > 0 on the open cone and w' = -2 v_rel v_rel' < 0
  (M4) the sonic weight of the reduction is positive: -D > 0 on the open cone
       (D = 3N^2V^2 - N^2 + 4NV - V^2 + 3; D(0) = 0 exactly, D'(0) = D_1 > 0), and
       D' > 0 on the closed cone -- so w'/w = D'/D < 0, GHJS's commutator sign,
       holds globally (S4-3 input).
V_R (= KHA's V) is NOT monotone -- one interior zero at x = -0.2528 -- and stays a
float remark (results/ec_monotonicity.json); no claim is made for it.

Three regions, all by interval arithmetic on (u, u') enclosures (no vanishing
denominators anywhere -- the quantities are regular including x = 0):
  * tube steps ([-8, -0.05], the certified A3+A2 tube): scaled (n, w, v, t) sub-boxes
    with z' from ``shootsys.rhs_enclosure``; A via the constraint, A' = A F_A;
  * sonic interval [-0.05, 0]: the certified A1 series (V0 box) + Cauchy-tail
    derivative enclosures (tc_energy.u_prime_balls); includes the endpoint x = 0;
  * centre tail t = e^x in (0, e^{x_tail}]: the certified A2 centre series, with the
    t-factored forms (v_rel = t g, rho^' = (theta w - w F_A)/A, ...) so the signs are
    certified on the whole open tail.

Run:  PYTHONPATH=problems/P4/src uv run python -m p4.validated.tc_monotone <tube.json>
      (writes results/theorem_c/tc_monotone.json)
"""
from __future__ import annotations

import json
import os
import sys
import time

from flint import arb

from .arbseries import precision
from .tc_energy import sonic_expansion_cached, u_prime_balls

BLOCKS4 = [(0, 1), (1, 3), (3, 4)]


def _mins(acc, kv):
    for k, v in kv.items():
        acc[k] = v if k not in acc else acc[k].min(v)


def margins_scaled(Z, dZ):
    """Sign margins from scaled (n, w, v, t) balls and (n', w', v') = dZ balls.
    All margins must be > 0.  Formulas: t' = t; A - 1 = 2 w t^2 (1 + v^2t^2/3
    + 4nv/3)/(1 - v^2t^2); A' = A F_A, F_A = 1 - A + 2wt^2(1 + v^2t^2/3)/(1 - v^2t^2);
    v_rel = a/b, a = t(1+nv), b = n + v t^2; Ds := t^2 D = -n^2 + t^2(3n^2v^2 + 4nv
    + 3) - v^2 t^4."""
    n, w, v, t = Z
    dn, dw, dv = dZ
    S = 1 - (v * t) ** 2
    A = 1 + 2 * w * t * t * (1 + (v * t) ** 2 / 3 + 4 * n * v / 3) / S
    FA = 1 - A + 2 * w * t * t * (1 + (v * t) ** 2 / 3) / S
    dA = A * FA
    a, b = t * (1 + n * v), n + v * t * t
    da = t * (1 + n * v) + t * (dn * v + n * dv)
    db = dn + dv * t * t + 2 * v * t * t
    vrel = a / b
    m = dict(
        vrel_pos=vrel,                                   # (M3) v_rel > 0
        wG_pos=arb(1) / 3 - vrel * vrel,                 # (M3) 1/3 - v_rel^2 > 0
        dvrel=(da * b - a * db) / (b * b),               # (M1) v_rel' > 0
        neg_drho=-(dw * A - w * dA) / (A * A),           # (M2) -rho^' > 0
        neg_Ds=n * n - t * t * (3 * (n * v) ** 2 + 4 * n * v + 3) + (v * t * t) ** 2,
        dD_pos=(6 * n * v * v - 2 * n / (t * t) + 4 * v) * dn
        + (6 * n * n * v + 4 * n - 2 * v * t * t) * dv
        + (2 * n * n / (t * t) - 2 * (v * t) ** 2),      # dD/dx > 0 (t dD/dt-part: t'=t)
    )                                                    # (M4) -t^2 D > 0, D' > 0
    return m


def run_tube(tube, nsub=4, prec=256):
    from . import shootsys as ss
    mins, fails, t0 = {}, [], time.time()
    with precision(prec):
        for sd in tube.steps:
            h = arb(sd.h)
            for j in range(nsub):
                w2 = sd.h / (2 * nsub)
                sm = (-h * (2 * j + 1)) / (2 * nsub) + arb(0, w2)
                Z = [zi + arb(0, sd.eps_z) for zi in ss.horner_vec(sd.co, sm)]
                f, _P = ss.rhs_enclosure(tube.sys4, Z, BLOCKS4)
                m = margins_scaled(Z, f[:3])
                _mins(mins, m)
                if not all(bool(x > 0) for x in m.values()):
                    fails.append((float(sd.x), j))
    return dict(ok=not fails, fails=fails[:8], time=time.time() - t0,
                mins={k: float(v.lower()) for k, v in mins.items()})


def margins_unscaled(u, up, at_zero=False):
    """Sign margins from (A, N, W, V) and (A', N', W', V') balls (sonic region;
    regular at x = 0).  -D > 0 is certified separately via D' > 0 and D(0) = 0."""
    A, N, W, V = u
    dA, dN, dW, dV = up
    vrel = (1 + N * V) / (N + V)
    dvrel = ((dN * V + N * dV) * (N + V) - (1 + N * V) * (dN + dV)) / (N + V) ** 2
    m = dict(vrel_pos=vrel, dvrel=dvrel,
             neg_drho=-(dW / W - 2 - dA / A))
    if not at_zero:
        m["wG_pos"] = arb(1) / 3 - vrel * vrel
    # D' = grad D . u' (D has no A, W): -D > 0 on [-0.05, 0) then follows from the
    # exact D(u(0)) = 0 (closed forms) and D' > 0; likewise v_rel(0) = 1/sqrt3 exactly,
    # so (M1) gives v_rel < 1/sqrt3 and w = 1/3 - v_rel^2 > 0 on the open cone.
    m["dD_pos"] = (6 * N * V * V - 2 * N + 4 * V) * dN \
        + (6 * N * N * V + 4 * N - 2 * V) * dV
    return m


def run_sonic(nx=6, prec=256):
    """[-0.05, 0] in nx sub-intervals, the last one closing at x = 0 (included)."""
    ex, t0 = sonic_expansion_cached(prec), time.time()
    mins, fails = {}, []
    with precision(prec):
        edges = [-0.05 * (1 - i / nx) for i in range(nx + 1)]
        for a, b in zip(edges, edges[1:]):
            xb = arb((a + b) / 2, (b - a) / 2)
            u = ex.eval(xb)
            up = u_prime_balls(ex, xb, prec)
            m = margins_unscaled(u, up, at_zero=(b == 0.0))
            _mins(mins, m)
            if not all(bool(x > 0) for x in m.values()):
                fails.append((a, b))
    return dict(ok=not fails, fails=fails, time=time.time() - t0,
                mins={k: float(v.lower()) for k, v in mins.items()})


def run_centre(x_tail=-8.0, prec=256):
    """The open tail t = e^x in (0, e^{x_tail}]: certified A2 centre series; the
    t-factored margins (v_rel = t g etc.) certify the signs on the whole tail."""
    from . import centre
    t0 = time.time()
    root = os.path.join(os.path.dirname(__file__), "..", "..", "..")
    aj = json.load(open(os.path.join(root, "results", "a3_midpoint.json")))
    with precision(prec):
        ball = lambda s: arb(s[1:-1].split(" +/-")[0]) + arb(0, float(s[1:-1].split(" +/- ")[1]))
        a, mu = ball(aj["a"]), ball(aj["mu"])
        ce = centre.centre_expansion(mu * (2 * a).exp(), nhat=(-a).exp(), K=30)
        cert = ce.certify()
        if not cert.ok:
            raise RuntimeError("centre certificate failed")
        ub = arb(x_tail).exp()
        tb = arb(ub / 2, (ub / 2).upper())          # encloses (0, e^{x_tail}]
        rho = arb(cert.nu) * arb(1) / 3
        tail = cert.tail_bound(rho.upper())
        dtail = (tail / (rho - arb(ub.upper()))).upper()      # Cauchy: |R'| on |t|<=ub
        sers = ce.series()
        nb, wb, vb = [srs(tb) + arb(0, cert.tail_bound(ub.upper()).upper()) for srs in sers]
        # theta Y = t Y': coefficients k Y_k; tails t |R'|
        thn, thw, thv = [sum((k * srs.coeffs()[k] * tb ** k for k in range(1, len(srs.coeffs()))), arb(0))
                         + arb(0, (arb(ub.upper()) * dtail).upper()) for srs in sers]
        # rho^' = t^2 (theta w/t^2 - w F_A/t^2)/A: the k = 1 coefficient of w vanishes
        # EXACTLY (the centre recursion at order 1 is homogeneous), checked here, so
        # theta w / t^2 is regular; its tail via |R_w'(t)|/t <= 2 tail(rho)/(rho-ub)^2.
        w1 = sers[1].coeffs()[1]
        if not (w1.mid() == 0 and float(w1.rad()) == 0.0):
            raise RuntimeError("w_1 is not exactly zero; rho-hat factoring invalid")
        d2tail = (2 * tail / (rho - arb(ub.upper())) ** 2).upper()
        thw_t2 = sum((k * sers[1].coeffs()[k] * tb ** (k - 2)
                      for k in range(2, len(sers[1].coeffs()))), arb(0)) + arb(0, d2tail)
        S = 1 - (vb * tb) ** 2
        Am1_t2 = 2 * wb * (1 + (vb * tb) ** 2 / 3 + 4 * nb * vb / 3) / S       # (A-1)/t^2
        A = 1 + Am1_t2 * tb * tb
        FA_t2 = -Am1_t2 + 2 * wb * (1 + (vb * tb) ** 2 / 3) / S                # F_A / t^2
        g = (1 + nb * vb) / (nb + vb * tb * tb)
        dg_num = (thn * vb + nb * thv) * (nb + vb * tb * tb) \
            - (1 + nb * vb) * (thn + thv * tb * tb + 2 * vb * tb * tb)
        thg = dg_num / (nb + vb * tb * tb) ** 2
        m = dict(
            vrel_pos=g,                                  # v_rel = t g > 0
            dvrel=g + thg,                               # v_rel' = t (g + theta g) > 0
            neg_drho=-(thw_t2 - wb * FA_t2) / A,         # rho^'/t^2 < 0
            wG_pos=arb(1) / 3 - (g * tb) ** 2,           # 1/3 - v_rel^2 > 0
            neg_Ds=nb * nb - tb * tb * (3 * (nb * vb) ** 2 + 4 * nb * vb + 3)
            + (vb * tb * tb) ** 2,                       # -t^2 D > 0
            # D' > 0 (so w'/w = D'/D < 0): t^2 D' = theta(Ds) - 2 Ds
            #   = 2n^2 - 2n theta(n) + t^2 (6nv+4)(v theta(n) + n theta(v))
            #     - 2 v theta(v) t^4 - 2 v^2 t^4   -> 2 n^2 > 0 at t -> 0
            dD_pos=2 * nb * nb - 2 * nb * thn
            + tb * tb * (6 * nb * vb + 4) * (vb * thn + nb * thv)
            - 2 * vb * thv * tb ** 4 - 2 * (vb * tb * tb) ** 2,
        )
        ok = all(bool(x > 0) for x in m.values())
    return dict(ok=ok, x_tail=x_tail, time=time.time() - t0,
                mins={k: float(v.lower()) for k, v in m.items()})


def main(argv):
    from .lintube import Tube
    t0 = time.time()
    res = {}
    res["sonic"] = run_sonic()
    print("sonic:  ok=%s  mins=%s" % (res["sonic"]["ok"], res["sonic"]["mins"]))
    res["centre"] = run_centre()
    print("centre: ok=%s  mins=%s" % (res["centre"]["ok"], res["centre"]["mins"]))
    path = argv[0] if argv else os.environ.get("P4_TUBE_CACHE", "")
    if path and os.path.exists(path):
        tube = Tube.load(path)
        res["tube_info"] = {k: str(v) for k, v in tube.info.items()}
        res["tube"] = run_tube(tube)
        print("tube:   ok=%s  steps=%d  mins=%s" %
              (res["tube"]["ok"], len(tube.steps), res["tube"]["mins"]))
    else:
        print("tube: skipped (no tube JSON; pass a path or set P4_TUBE_CACHE)")
    res["ok"] = all(res[r]["ok"] for r in ("sonic", "centre", "tube") if r in res)
    res["total_time"] = time.time() - t0
    out = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                       "results", "theorem_c", "tc_monotone.json")
    json.dump(res, open(out, "w"), indent=1)
    print("P4-3 tube sign certificate %s; written %s" % ("PASSES" if res["ok"] else "FAILS", out))
    return res


if __name__ == "__main__":
    main(sys.argv[1:])
