"""Theorem B, Stage 3 (c), (d): Krawczyk enclosures of the zeros of E near kappa_1 and near the
gauge value kappa-bar, and the identification of the latter with the pure-gauge mode.

Krawczyk (E: C -> C as a map of R^2, Jacobian = the real form of E'):  on the box B = kappa_c + [-w_B, w_B]^2
    K(B) = kappa_c - Y E(kappa_c) + (1 - Y E'(B)) (B - kappa_c),   Y = 1/E'(kappa_c) (midpoint),
K(B) contained in the interior of B  =>  E has exactly one zero in B, and it lies in K(B).
E and E'(B) come from the Taylor model of ``modecount.Context.E_model`` on |delta| <= w = 2 w_B:
E(kappa_c) = c_0 + R(0), E'(delta) = sum_k k c_k delta^{k-1} + R'(delta) with |R'| <= rem/(w - sqrt2 w_B) on B (|delta| <= sqrt2 w_B there)
(Cauchy's estimate for the derivative of the remainder).  The complex box arithmetic of acb
over-approximates the real 2x2 interval products, so the containment test is sound.
Gauge mode: kappa-bar = 2 - A_0 + 2 W_0/3 (certified sonic data); the pure-gauge perturbation
g = (A', N' + kappa-bar N, W', V') / A'(0) of the background (the residual coordinate freedom
x -> x + eps e^{kappa-bar s}, admissible in the sonic-point gauge iff N_p(0) = 0, i.e. iff kappa =
-N'(0)/N(0) = kappa-bar) is checked as a ball identity: it satisfies the 4D linearised system and
the linearised constraint order by order, and coincides with the 4D sonic series at kappa-bar; it
is regular at the centre (background analytic in t = e^x), so E(kappa-bar) = 0.
"""
from __future__ import annotations

from types import SimpleNamespace

from flint import acb, acb_mat, arb

from . import linsonic, linsonic4
from .arbseries import abs_upper, precision
from .linsys import LinSystem, to_acb


def _cbox(rad):
    return acb(arb(0, rad), arb(0, rad))


def krawczyk(ctx, kappa_c, w, m=3, newton=True, verbose=False):
    """Krawczyk test on B = kappa_c' +/- w/2 (kappa_c' = Newton-refined centre if ``newton``)."""
    with precision(ctx.prec):
        kc = complex(kappa_c)
        E, info = ctx.E_model(acb(kc.real, kc.imag), w, m)
        if E is None:
            raise RuntimeError("E model not tight enough")
        if newton:
            step = complex(E.c[0].mid()) / complex(E.c[1].mid())
            kc = kc - step
            E, info = ctx.E_model(acb(kc.real, kc.imag), w, m)
        wB = w / 2
        B = _cbox(wB)
        E0 = E.c[0] + _cbox(E.rem)                                       # E(kappa_c) = c_0 + R(0)
        dE_B = sum((E.c[k] * k * B ** (k - 1) for k in range(1, E.m + 1)), acb(0))
        dE_B += _cbox(abs_upper(E.rem / (arb(w) - arb(2).sqrt() * wB)))    # |delta| <= sqrt2 w_B on the square box
        Y = 1 / acb(E.c[1].real.mid(), E.c[1].imag.mid())
        Kop = -Y * E0 + (1 - Y * dE_B) * B
        ok = bool(B.real.contains_interior(Kop.real)) and bool(B.imag.contains_interior(Kop.imag))
        zero = acb(kc.real, kc.imag) + Kop
        return dict(ok=ok, kappa_c=kc, w_B=wB, K=Kop, zero=zero, E=E0, dE=E.c[1], dE_box=dE_B, rem=float(E.rem),
                    model=E, info=info)


# ---------------------------------------------------------------------------------------------
# the gauge mode
# ---------------------------------------------------------------------------------------------
def gauge_value(bg):
    """kappa-bar = 2 - A_0 + 2 W_0 / 3 as an arb ball (certified A1 sonic data)."""
    return linsonic.gauge_eigenvalue(bg)


def gauge_generator(bg, kbar, K=40):
    """Coefficients g_n = (A', N' + kbar N, W', V')_n / A'_0 (lists of 4 acb) of the pure-gauge mode."""
    u = bg.series()
    du = [s.deriv() for s in u]
    A1 = du[0][0]
    out = []
    for n in range(K + 1):
        out.append([acb(du[0][n] / A1), (acb(du[1][n]) + to_acb(kbar) * u[1][n]) / A1, acb(du[2][n] / A1), acb(du[3][n] / A1)])
    return out


def gauge_checks(bg, K=40, prec=256):
    """Ball identities: (i) the gauge generator satisfies the 4D linearised system at kappa-bar
    (residual orders 0..K-2 contain 0) and the linearised constraint; (ii) it coincides with the 4D
    sonic series over a box containing kappa-bar (all balls overlap, n <= K)."""
    with precision(prec):
        kb = gauge_value(bg)
        kbar = acb(kb)
        g = gauge_generator(bg, kb, K)
        L = LinSystem()
        polys = linsonic4.polys4(L)
        Pser, Gser = linsonic4.series_parts4(polys, bg.series())
        st = linsonic4._Structure4(Pser, Gser, kbar)
        q = [[acb_mat([[x] for x in gn])] for gn in g]
        res = [st.residual(q, j, K, 1)[0] for j in range(K - 1)]
        res_ok = all(res[j][i, 0].contains(acb(0)) for j in range(K - 1) for i in range(4))
        res_max = max(float(abs_upper(res[j][i, 0])) for j in range(K - 1) for i in range(4))
        fake = SimpleNamespace(L=L, K=K, bg=bg, balls=g, kappa=kbar)
        cr = linsonic4.constraint_residual(fake, K - 1)
        con_ok = all(z.contains(acb(0)) for z in cr)
        w = max(2 * float(kb.rad()), 1e-15)
        ex = linsonic4.linear_sonic_expansion4(bg, acb(kb.mid()), width=w, m=2, K=K)
        cert = ex.certify()
        ovl = all(ex.balls[n][i].overlaps(g[n][i]) for n in range(K + 1) for i in range(4))
        dist = max(float(abs_upper(ex.balls[n][i] - g[n][i])) for n in range(11) for i in range(4))    # n <= 10: tight balls
        return dict(kappa_bar=kb, residual_ok=res_ok, residual_max=res_max, constraint_ok=con_ok,
                    constraint_max=max(float(abs_upper(z)) for z in cr), series_overlap=ovl, series_dist=dist,
                    cert_ok=bool(cert.ok), nu=float(cert.nu), sigma=cert.details.get("sigma"))
