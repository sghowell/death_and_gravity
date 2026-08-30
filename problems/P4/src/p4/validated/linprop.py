"""Theorem B, Stage 2 (a): validated propagation of a linear system  P(kappa) y' = G(kappa) y
along the certified background tube (``lintube.Tube``), for complex kappa (point or box).

Per step the kappa-independent data of the step (``Tube.system_data``) is combined with kappa:
P(s) = sum_j kappa^j P_j(s) etc., and ``linstep.step_bound`` gives Y_K(-h) and the Groenwall
bound; every solution set (``linstep.LohnerSet``) is then mapped by the mean-value form.
``propagate`` returns the sets at the end of the range together with a log; ``fundamental_matrix``
propagates the d unit vectors (columns of Phi).
"""
from __future__ import annotations

import time

from flint import acb, acb_mat, acb_poly, arb

from . import linstep
from .arbseries import abs_upper, precision
from .linsys import to_acb


def _kpowers(kap, n):
    out = [acb(1)]
    for _ in range(n - 1):
        out.append(out[-1] * kap)
    return out


def combine(data, kap, d, K):
    """kappa-combined step data: (Pp, Gp acb_poly matrices, Pbox, Gbox acb_mat lists, incP, incG, dPt, dGt)."""
    kp = _kpowers(kap, 3)
    ka = [abs_upper(k) for k in kp]

    def polys(entries):
        out = [[None] * d for _ in range(d)]
        for r in range(d):
            for c in range(d):
                coefs = [acb(0)] * (K + 1)
                for j, (co, _) in enumerate(entries):
                    for i in range(K + 1):
                        coefs[i] += kp[j] * co[r][c][i]
                out[r][c] = acb_poly(coefs)
        return out

    def boxes(lists):
        n = len(lists[0])
        return [acb_mat([[sum((kp[j] * lists[j][b][r][c] for j in range(len(lists))), acb(0)) for c in range(d)]
                         for r in range(d)]) for b in range(n)]

    def inc(mats):
        return [[sum((ka[j] * mats[j][r][c] for j in range(len(mats))), arb(0)) for c in range(d)] for r in range(d)]
    Pp, Gp = polys(data["P"]), polys(data["G"])
    dPt = sum((ka[j] * data["P"][j][1] for j in range(len(data["P"]))), arb(0))
    dGt = sum((ka[j] * data["G"][j][1] for j in range(len(data["G"]))), arb(0))
    return Pp, Gp, boxes(data["Pbox"]), boxes(data["Gbox"]), inc(data["incP"]), inc(data["incG"]), dPt, dGt


def box_extra(data, kap, kc, d):
    """(P(box), G(box), P(box) - P(kc), G(box) - G(kc)) per sub-box for a kappa box ``kap`` around kc."""
    kb, kp = _kpowers(kap, 3), _kpowers(kc, 3)
    dk = [kb[j] - kp[j] for j in range(3)]

    def boxes(lists, w):
        return [acb_mat([[sum((w[j] * lists[j][b][r][c] for j in range(len(lists))), acb(0)) for c in range(d)]
                         for r in range(d)]) for b in range(len(lists[0]))]
    return boxes(data["Pbox"], kb), boxes(data["Gbox"], kb), boxes(data["Pbox"], dk), boxes(data["Gbox"], dk)


def propagate(tube, S, kappa, sets, K=None, x_stop=None, prec=256, nsub=8, verbose=False):
    """Propagate the LohnerSets ``sets`` (in C^d, d = S.d) from the tube's first step down to
    x_stop (default: the tube's end).  A kappa *box* (acb with nonzero radius) is handled by the
    point propagation at its centre plus the rigorous per-step perturbation bound (linstep).
    Returns (sets, log)."""
    kap = to_acb(kappa)
    is_box = float(max(kap.real.rad(), kap.imag.rad())) > 0
    kc = acb(kap.real.mid(), kap.imag.mid()) if is_box else kap
    Sc = list(S.weights)
    log, t0 = [], time.time()
    with precision(prec):
        for sd in tube.steps:
            if x_stop is not None and float(sd.x) <= x_stop + 1e-15:
                break
            data = tube.system_data(sd, S, nsub)
            Kl = len(sd.co) - 1 if K is None else min(K, len(sd.co) - 1)
            Pp, Gp, Pbox, Gbox, incP, incG, dPt, dGt = combine(data, kc, S.d, Kl)
            bd = box_extra(data, kap, kc, S.d) if is_box else None
            Jm, bound, info = linstep.step_bound(sd, kc, Kl, S.d, Pp, Gp, Pbox, Gbox, incP, incG, dPt, dGt, Sc, bd, S.d0)
            for st in sets:
                st.propagate(Jm, bound, Sc=Sc)
            info.update(x=float(sd.x) - sd.h, h=sd.h, width=max(float(abs_upper(z - z.mid())) for st in sets for z in st.hull()))
            log.append(info)
            if verbose:
                print(f"x={info['x']:+.4f} h={sd.h:.4f} L={info['L']:.3g} Pinv={info['Pinv']:.3g} Dsup={info['Dsup']:.1e} "
                      f"bound={info['bound']:.1e} width={info['width']:.1e}", flush=True)
    log.append(dict(time=time.time() - t0))
    return sets, log


def fundamental_matrix(tube, S, kappa, **kw):
    """Enclosure of Phi(x_end, x_start; kappa) (d x d acb_mat, hull of the column sets) and the log."""
    d = S.d
    sets = [linstep.LohnerSet([acb(int(i == j)) for i in range(d)]) for j in range(d)]
    sets, log = propagate(tube, S, kappa, sets, **kw)
    cols = [st.hull() for st in sets]
    return acb_mat([[cols[c][r] for c in range(d)] for r in range(d)]), log
