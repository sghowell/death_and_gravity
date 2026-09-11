"""Literal constant-mass operator, readout and homogeneous-source bridge."""

from functools import cache

import sympy as s
from p8_proca_nonlocal_response import source as old
from p8_vacuum_affine_proca_gaussian.bridge import KAPPA, MASS, ZETA

__all__ = ["KAPPA", "MASS", "ZETA", "data", "homogeneous"]


@cache
def data():
    N, a, k, m = s.symbols("N a k m", positive=True)
    n, v = old.n[0], old.v[0]
    q = k * k / a**2
    w2 = N * N * (q + m * m)
    rows, checks = {}, {}

    def change(expr):
        return n * s.diff(expr, N) + v * a * s.diff(expr, a)

    def on(expr):
        expr = s.cancel(expr.subs(N, 1))
        return s.factor(expr.subs(k * k, a * a * m * m * old.z / (1 - old.z)))

    for sector in ("T", "L"):
        # The S6.85 longitudinal coordinate is the scalar potential:
        # A_i = i k_i potential. Its g^2 differs by the FIXED k^2 only.
        g2 = a / N if sector == "T" else a * m * m / (N * (q + m * m))
        old_g2 = a / N if sector == "T" else a**3 * m * m * q / (N * (q + m * m))
        factor = s.S.One if sector == "T" else k * k
        r = on(change(w2) / (2 * w2))
        rg = on(change(g2) / (2 * g2))
        checks[sector + "_whole_metric_family_canonical_factor"] = s.cancel(
            old_g2 - factor * g2
        )
        checks[sector + "_actual_frequency_source"] = old.clean(
            r - old.data(sector)["r"]
        )
        checks[sector + "_actual_normalization_source"] = old.clean(
            rg - old.data(sector)["rg"]
        )
        rows[sector] = {
            "g_squared": g2,
            "frequency_squared": w2,
            "delta_log_frequency": r,
            "delta_log_g": rg,
        }
        for name, sign, derivative in (
            ("energy", -s.S.One, lambda x: s.diff(x, N)),
            ("pressure", s.Rational(1, 3), lambda x: a * s.diff(x, a)),
        ):
            raw = (derivative(g2) / g2, -derivative(g2 * w2) / (g2 * w2))
            normal = 3 * v if name == "energy" else n + 3 * v
            target = old.readout(sector, name)
            actual = {}
            for key, weight in zip(("A", "B"), raw):
                value = sign * on(weight)
                contact = sign * on(change(weight)) - normal * value
                checks[sector + "_" + name + "_weight_" + key] = old.clean(
                    value - target[key]
                )
                checks[sector + "_" + name + "_contact_" + key] = old.clean(
                    contact - target["delta_" + key]
                )
                actual[key] = value
                actual["delta_" + key] = contact
            rows[sector][name] = actual
    gL = rows["L"]["g_squared"]
    checks["zero_momentum_cartesian_longitudinal_normalization"] = s.factor(
        gL.subs(k, 0) - a / N
    )
    checks["canonical_mass_same"] = MASS**2 * ZETA - 1
    return {
        "actual_lapse_scale_operator_and_readout_vertices": rows,
        "state": "The S6.176 source-aware conditional state uses the same S6.55 all-order Cauchy covariance at t0=-1/2, canonical mass1000 and physical metric. Sources vanish on an initial neighborhood; there is no independent state variation.",
        "measure": "Only the three physical Proca canonical modes are quantized. The metric-independent constant longitudinal k map is used only at k>0; the regular Cartesian k=0 limit keeps the third polarization.",
        "prescription": "The S6.176 ordinary covariant dimensional prescription at mu=m is exactly the S6.84/S6.85 physical vector prescription. All mass-source jets vanish as functions; the scalar profile in the OLD total clock response is not included.",
        "checks": checks,
    }


@cache
def homogeneous():
    t = s.symbols("t", real=True)
    mean = s.Function("S0")(t)
    W0, W1, metric, measure = s.symbols(
        "W0 W1 inverse_temporal_metric volume", real=True
    )
    dW0, dW1, dS, dm, dvol = s.symbols("dW0 dW1 dS dm dvol", real=True)
    F, dF = s.symbols("field_strength field_strength_variation", real=True)
    density = measure * (
        -ZETA * F * F / 4 + metric * (W0 - mean) ** 2 / 2 - W1 * W1 / 2
    )
    variation = sum(
        s.diff(density, x) * dx
        for x, dx in ((W0, dW0), (W1, dW1), (F, dF), (metric, dm), (measure, dvol))
    )
    variation += s.diff(density, mean) * dS
    point = {W0: mean, W1: 0, F: 0}
    # The temporal one-form may be nonzero; its exterior derivative vanishes.
    u = s.Function("u")(t)
    N = s.Function("N")(t)
    a = s.Function("a")(t)
    X = s.diff(u, t) ** 2 / N**2
    Box = (
        s.diff(u, t, 2) / N**2
        + 3 * s.diff(a, t) * s.diff(u, t) / (a * N**2)
        - s.diff(N, t) * s.diff(u, t) / N**3
    )
    Z = s.diff(u, t) ** 2 * (s.diff(u, t, 2) - s.diff(N, t) * s.diff(u, t) / N) / N**4
    return {
        "arbitrary_homogeneous_invariants": {"X": X, "Box_u": Box, "u_Hess_u_u": Z},
        "source": "For every sufficiently small smooth homogeneous physical metric/clock history inside the same analytic domain, the complete source is S=S0(t)dt. Thus dS=0 and the compatible unique classical mean is Wbar=S, even when S0 is nonzero.",
        "exact_sector_equality": "On these histories F(Wbar)=0 and Wbar-S=0. Every first metric and light variation of the complete mean action vanishes before restricting its test variation. Therefore the full conditional mean stress equals the ordinary connected Proca stress on the SAME varied metric, with no source/contact contribution.",
        "preparation": "All histories agree with the CD clock and zero source near the original Cauchy slice. A homogeneous smooth coherent shift defines the local sourced algebra state; no finite spatial energy or global Fock implementability is inferred.",
        "checks": {
            "closed_source_spatial_derivatives": s.diff(mean, s.Symbol("x", real=True)),
            "complete_mean_density_zero": density.subs(point),
            "every_first_mean_variation_zero": s.expand(variation.subs(point)),
            "homogeneous_X_on_clock": X.subs({u: t, N: 1}).doit() - 1,
        },
    }
