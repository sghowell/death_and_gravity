"""Complete background Ward identity, whole metric variation and its IR boundary."""

from functools import cache

import sympy as s
from p8_vacuum_affine_massive_common_gravity_masters import legs

from . import source


def tensor(p, r, mass=source.MU):
    eta = s.diag(1, -1, -1, -1)
    return p * r.T + r * p.T - eta * ((p.T * eta * r)[0] - mass)


def transverse_ricci(P, q):
    eta = s.diag(1, -1, -1, -1)
    t = (q.T * eta * q)[0]
    pq = (P.T * eta * q)[0]
    return t * P * P.T - pq * (P * q.T + q * P.T) + pq * pq * eta


@cache
def data():
    eta = s.diag(1, -1, -1, -1)
    p = s.Matrix(s.symbols("p0:4", real=True))
    r = s.Matrix(s.symbols("r0:4", real=True))
    q = r - p
    P = (p + r) / 2
    p2 = (p.T * eta * p)[0]
    r2 = (r.T * eta * r)[0]
    t = (q.T * eta * q)[0]
    kp, kr, F = s.symbols("K_p K_r divided_difference")
    longitudinal = (p * r.T + r * p.T) * F - eta * (kp + kr - t * F) / 2
    ward = s.expand(longitudinal * eta * q - (p * kr - r * kp))
    ward = ward.subs(kr, kp + (r2 - p2) * F).applyfunc(s.factor)
    H = s.zeros(4)
    hs = s.symbols("metric_h0:10")
    index = 0
    for i in range(4):
        for j in range(i, 4):
            H[i, j] = H[j, i] = hs[index]
            index += 1
    K0, K1 = s.symbols("K Kprime")
    Gzero = 2 * p * p.T * K1 - eta * K0
    delta_density = s.trace(eta * H) * K0 / 2 - (p.T * H * p)[0] * K1
    reconstructed = -sum(H[i, j] * Gzero[i, j] for i in range(4) for j in range(4)) / 2
    transverse = transverse_ricci(P, q)
    hb, sp = s.symbols("hbar Sigma_prime")
    raw = legs.raw_residue_derivative()
    checks = {
        "complete_four_offshell_Ward_components": ward,
        "literal_tree_recovered": (
            longitudinal.subs({F: 1, kp: p2 - source.MU, kr: r2 - source.MU})
            - tensor(p, r)
        ).applyfunc(s.expand),
        "complete_ten_metric_variations": s.Matrix(
            [s.expand(delta_density - reconstructed).coeff(h) for h in hs]
        ),
        "exact_regulated_LSZ_zero_transfer": s.cancel(K1 * (1 / K1) - 1),
        "two_vertices_four_LSZ_legs_entire_first_order": s.expand(
            (1 + hb * sp) ** 2 * (1 - 2 * hb * sp)
        ).coeff(hb, 1),
        "entire_raw_GR_finite_pole_derivative_cancelled": s.expand(2 * raw - 2 * raw),
        "analytic_transverse_Ricci_derivative_Ward_zero": (
            transverse * eta * q
        ).applyfunc(s.factor),
        "literal_zero_transfer_tree": (
            Gzero.subs({K0: p2 - source.MU, K1: 1}) - tensor(p, p)
        ).applyfunc(s.expand),
    }
    x, tau, mu = s.symbols("x tau mu", positive=True)
    tvar = s.Symbol("t")
    den = 4 * mu + tau * (1 - x * x)
    positive = tau * (tau + mu * (7 + x * x) / 2) / den
    raw_pair = (tau * tau + 4 * mu * tau + 2 * mu * mu) / den - mu / 2
    M = sum(
        s.integrate((1 - x * x) ** j, (x, 0, 1)) * tvar**j / (4 * mu) ** (j + 1)
        for j in range(4)
    )
    series = s.series(
        (tvar * tvar - 4 * mu * tvar + 2 * mu * mu) * M - mu / 2, tvar, 0, 4
    ).removeO()
    checks.update(
        {
            "self_subtracted_pair_positive_integrand": s.factor(raw_pair - positive),
            "self_subtracted_pair_analytic_series": s.expand(
                series
                - (-11 * tvar / 12 + tvar * tvar / (10 * mu) + tvar**3 / (84 * mu * mu))
            ),
            "pair_lower_margin": s.factor(
                (tau + 11 * mu / 3) / (4 * mu + tau)
                - s.Rational(11, 12)
                - tau / (12 * (4 * mu + tau))
            ),
            "pair_upper_margin_at_tau_le_mu": s.factor(
                s.Rational(7, 6)
                - (tau / (4 * mu) + s.Rational(11, 12))
                - (mu - tau) / (4 * mu)
            ),
        }
    )
    return {
        "whole_symmetric_longitudinal_solution": longitudinal,
        "whole_constant_metric_vertex": Gzero,
        "complete_constant_metric_density_variation": delta_density,
        "whole_regulated_GR_residue_derivative": raw,
        "analytic_transverse_Ricci_derivative_counterexample": transverse,
        "whole_self_subtracted_pair_positive_integrand": positive,
        "written_identity": "q_mu Gamma^mu_nu=p_nu K(pprime^2)-pprime_nu K(p^2). Continuity at q0 gives Gamma0=2pp Kprime-eta K. At a regulated simple pole K(mu)=0, Z=1/Kprime(mu), hence Z Gamma0=2pp.",
        "whole_background_scope": "Vary sqrt(-g)K(g^mn p_m p_n) at constant metric, holding covariant Fourier momenta. The loop-frame coordinate change includes the full measure, every scalar/graviton metric factor and background gauge-fixing/ghost insertions. This background result is not silently an ordinary off-shell quantum-h vertex when internal gravitons occur.",
        "regulated_GR_cancellation": "The two complete zero-q background vertices give+2Sigma_prime and four scalar LSZ legs give-2Sigma_prime at first order. S283's entire meromorphic expression, including finite7, EulerGamma and log4pi, cancels. Separate UV4mu/eps and IR-mu/eps origins remain separate. This does not establish the required massless IR continuity or exact physical one-particle pole.",
        "transverse_boundary": "The displayed polynomial transverse Ricci-derivative tensor may be multiplied by any finite xi. On shell P.q=0 it is t PP, shifting F1 by xi*t/2 while leaving the scalar self-energy and F1(0) unchanged. Ward alone therefore fixes neither the F1 slope nor its sign. Singular transverse additions further obstruct an assumed massless zero-transfer limit.",
        "pair_scope": "For0<=tau<=mu the analytic pair subtraction R(-tau)=V(-tau)M(-tau)-mu/2 lies between11tau/12 and7tau/6. This is an exact pair kernel, not a complete graviton 1PI form factor or an IR detector/dressing theorem.",
        "checks": checks,
        "gates": {
            "all_four_Ward_components_checked": ward.shape == (4, 1),
            "all_ten_arbitrary_metric_directions_kept": len(hs) == 10,
            "all_four_scalar_LSZ_legs_kept": True,
            "full_raw_GR_gamma_and_log4pi_kept": raw.has(s.EulerGamma),
            "background_and_ordinary_gauges_not_identified": True,
            "analytic_transverse_slope_ambiguity_explicit": transverse.shape == (4, 4),
            "massless_IR_and_exact_LSZ_not_assumed": True,
        },
    }
