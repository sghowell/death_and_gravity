"""S4-6 review: an independent pointwise evaluation of the strip condition S.

Not the author's code path: the background is evaluated here from a centre
expansion (series + certified tail), the pairs come from the reviewed
tc_commutator.build_pairs on constant Taylor models, and S is assembled
directly from the definitions of notes/s4-5-assembly.md section 4.1 (S4-S6):

  A0 = nT/dT, B0 = nB/dB, A1 = A0 - (log B0)', B1 = B0 + A0' - A0 (log B0)',
  tt = TN_top/dT_top, nu = Im A1 + tau tt, rho = Re A1 - c tt,
  f = (c + tau^2/c) tt - (tau/c) nu,
  K = Re B1 - (tau/c) Im B1 - (f' + rho f)/2,  S = K + nu^2/4,  c = -1/32.

Sigma_c = S/(1 + (y tau)^2).  Everything is arb/acb ball arithmetic, so the
result is a rigorous enclosure for the background given; with the midpoint
background of the certified (a*, mu*) balls (a point of the certified box)
it is tight and every certified enclosure must contain it.
"""
from flint import acb, arb, ctx

C = arb(-1) / 32


def midpoint_expansion(prec=640, K=30):
    from p4.validated import centre
    from p4.validated.arbseries import precision
    from p4.validated.ec_centre import certified_centre_balls
    with precision(prec):
        a, mu = certified_centre_balls()
        a, mu = arb(a.mid()), arb(mu.mid())
        ce = centre.centre_expansion(mu * (2 * a).exp(), nhat=(-a).exp(), K=K, prec=prec)
        if not ce.certify(prec=prec).ok:
            raise RuntimeError("midpoint centre certificate failed")
    return ce


def background(ce, y):
    """(A, N, W, V) at y = e^x from the centre expansion (balls + tail)."""
    yb = arb(y)
    tb = ce.cert.tail_bound(yb.abs_upper())
    vals = []
    for i in range(3):
        acc = arb(0)
        for k in range(ce.K, -1, -1):
            acc = acc * yb + ce.balls[k][i]
        vals.append(acc + arb(0, tb))
    n, w, v = vals
    S = 1 - yb * yb * v * v
    Ahat = 2 * w * (1 + yb * yb * v * v / 3 + 4 * n * v / 3) / S
    return 1 + yb * yb * Ahat, n / yb, yb * yb * w, yb * v


def pairs_at(at, ce, y):
    from p4.validated import tc_commutator as tc
    from p4.validated.tc_tm import LC, TM, AtomKP
    A, N, W, V = background(ce, y)
    comps = tuple(LC(TM.const(x, y), 0) for x in (A, N, W, V))
    return tc.build_pairs(AtomKP(at, comps, sd_shift=0))


def _val(kp, slot, k):
    acc = acb(0)
    for jet in reversed(kp):
        lc = jet[slot]
        if lc.off != 0 or len(lc.tm.c) != 1:
            # constant TMs at offset 0 throughout (products and sums of them)
            raise ValueError("non-constant point pair")
        acc = acc * k + lc.tm.c[0]
    return acc


def sigma_c(pr, y, eps, tau):
    """(Sigma_c, S) at kappa = eps + i tau."""
    k = acb(arb(eps), arb(tau))
    t = arb(tau)
    nT = [_val(pr["nT"], s, k) for s in (0, 1)]
    dT = [_val(pr["dT"], s, k) for s in (0, 1)]
    nB = [_val(pr["nB"], s, k) for s in (0, 1, 2)]
    dB = [_val(pr["dB"], s, k) for s in (0, 1, 2)]
    A0 = nT[0] / dT[0]
    A0p = nT[1] / dT[0] - nT[0] * dT[1] / (dT[0] * dT[0])
    B0 = nB[0] / dB[0]
    l1n, l1d = nB[1] / nB[0], dB[1] / dB[0]
    L1 = l1n - l1d                                  # (log B0)'
    L2 = (nB[2] / nB[0] - l1n * l1n) - (dB[2] / dB[0] - l1d * l1d)   # (log B0)''
    A1 = A0 - L1
    A1p = A0p - L2
    B1 = B0 + A0p - A0 * L1
    TNt = [pr["TN"][-1][s].tm.c[0] for s in (0, 1)]
    dTt = [pr["dT"][-1][s].tm.c[0] for s in (0, 1)]
    tt = TNt[0] / dTt[0]
    ttp = TNt[1] / dTt[0] - TNt[0] * dTt[1] / (dTt[0] * dTt[0])
    nu = A1.imag + t * tt
    nup = A1p.imag + t * ttp
    rho = A1.real - C * tt
    cf = C + t * t / C
    f = cf * tt - (t / C) * nu
    fp = cf * ttp - (t / C) * nup
    K = B1.real - (t / C) * B1.imag - (fp + rho * f) / 2
    S = K + nu * nu / 4
    s = arb(y) * t
    return S / (1 + s * s), S


def needed_prec(y):
    import math
    return max(384, int(160 + 14 * math.log2(1.0 / y)))
