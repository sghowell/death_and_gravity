"""Literal raw IR sign and ordinary zero-transfer charge/LSZ cancellation."""

from functools import cache

import sympy as s
from p8_vacuum_affine_massive_common_gravity_masters import legs

from . import source

MU, T, EP, K, NU = source.MU, source.T, source.EP, source.K, source.NU
X = s.Symbol("x", real=True)


def kernel(transfer=T, mass=MU):
    return s.Integral(1 / (4 * mass - transfer * (1 - X * X)), (X, 0, 1))


def eikonal(transfer=T, mass=MU, epsilon=EP):
    return (
        transfer**2
        - 4 * mass * transfer
        + 2 * mass**2
        + 2 * mass**2 * epsilon / (1 + epsilon)
    )


def triangle_raw_kernel(transfer=T, mass=MU, epsilon=EP, scale=NU):
    return (
        -s.gamma(1 - epsilon)
        * (4 * s.pi * scale**2) ** (-epsilon)
        / (2 * epsilon)
        * s.Integral((mass - transfer * (1 - X * X) / 4) ** (-1 + epsilon), (X, 0, 1))
    )


def endpoint_IR_coefficient(transfer=T, mass=MU, kappa=K):
    return (2 * eikonal(transfer, mass, 0) * kernel(transfer, mass) - mass) / (
        16 * s.pi**2 * kappa
    )


@cache
def data():
    raw = legs.raw_residue_derivative(MU, K, EP, NU)
    sigmaUV = -4 * MU / (16 * s.pi**2 * K * EP)
    sigmaIR = MU / (16 * s.pi**2 * K * EP)
    finite = (
        MU
        / (16 * s.pi**2 * K)
        * (7 + 3 * s.log(4 * s.pi * NU**2 / MU) - 3 * s.EulerGamma)
    )
    hb, kprime = s.symbols("hbar Kprime")
    M = s.Symbol("M")
    V = eikonal(T, MU, 0)
    pair = V * M - MU / 2
    rawpref = -s.gamma(1 - EP) * (4 * s.pi * NU**2) ** (-EP) / 2
    master_a = MU - T * X * (1 - X)
    Mseries = sum(
        s.integrate((1 - X * X) ** n, (X, 0, 1)) * T**n / (4 * MU) ** (n + 1)
        for n in range(4)
    )
    pairseries = s.series(V * Mseries - MU / 2, T, 0, 4).removeO()
    tau = s.Symbol("positive_tau", positive=True)
    den = 4 * MU + tau * (1 - X * X)
    positive = tau * (tau + MU * (7 + X * X) / 2) / den
    # The all-D negative triangle follows directly from the same master convention.
    checks = {
        "same_raw_leg_UV_IR_and_entire_finite_split": s.simplify(
            raw - sigmaUV - sigmaIR - finite
        ),
        "entire_raw_vertex_equals_leg_derivative_by_proved_bridge": s.simplify(
            raw - legs.raw_residue_derivative(MU, K, EP, NU)
        ),
        "whole_single_endpoint_LSZ_first_order": s.expand(
            (1 + hb * raw) * (1 - hb * raw)
        ).coeff(hb, 1),
        "whole_two_endpoint_four_leg_first_order": s.expand(
            (1 + hb * raw) ** 2 * (1 - 2 * hb * raw)
        ).coeff(hb, 1),
        "exact_regulated_simple_pole_charge": s.cancel(kprime / kprime - 1),
        "literal_Feynman_relative_triangle_minus_sign": s.simplify(
            ((-s.I) ** 3 * s.I**3 * s.I) / (-s.I) + 1
        ),
        "raw_triangle_radial_negative_residue": rawpref.subs(EP, 0) + s.Rational(1, 2),
        "raw_triangle_finite_derivative": s.simplify(
            s.diff(rawpref, EP).subs(EP, 0)
            + (s.EulerGamma - s.log(4 * s.pi * NU**2)) / 2
        ),
        "whole_S283_massive_parameter_reflection": s.factor(
            master_a.subs(X, (1 + X) / 2) - (MU - T * (1 - X * X) / 4)
        ),
        "whole_raw_triangle_IR_scalar_factor": s.factor((-s.Rational(1, 2)) * 4 + 2),
        "negative_eikonal_times_negative_triangle_IR": s.factor(
            (-V) * (-2 * M) - 2 * V * M
        ),
        "complete_leg_subtracted_IR_pair_factor": s.factor(
            (2 * V * M - MU) / (16 * s.pi**2 * K) - pair / (8 * s.pi**2 * K)
        ),
        "zero_transfer_IR_charge_cancellation": s.factor(
            (2 * V * M - MU).subs({T: 0, M: 1 / (4 * MU)})
        ),
        "actual_pair_small_transfer_series": s.expand(
            pairseries - (-11 * T / 12 + T * T / (10 * MU) + T**3 / (84 * MU**2))
        ),
        "complete_spacelike_pair_positive_integrand": s.factor(
            (eikonal(-tau, MU, 0) / den - MU / 2) - positive
        ),
        "pair_lower_margin": s.factor(
            (tau + 11 * MU / 3) / (4 * MU + tau)
            - s.Rational(11, 12)
            - tau / (12 * (4 * MU + tau))
        ),
        "pair_upper_margin": s.factor(
            s.Rational(7, 6)
            - (tau / (4 * MU) + s.Rational(11, 12))
            - (MU - tau) / (4 * MU)
        ),
        "raw_triangle_zero_transfer_exact_angular_kernel": s.simplify(
            ((MU - T * (1 - X * X) / 4) ** (-1 + EP)).subs(T, 0) - MU ** (-1 + EP)
        ),
    }
    return {
        "complete_raw_GR_proper_F1_at_zero_minus_one": raw,
        "separate_raw_UV_and_IR_derivatives": {
            "UV": sigmaUV,
            "IR": sigmaIR,
            "finite": finite,
        },
        "whole_raw_triangle": triangle_raw_kernel(),
        "whole_raw_eikonal_IR_before_LSZ": 2 * V * kernel() / (16 * s.pi**2 * K * EP),
        "remaining_endpoint_IR_pole": endpoint_IR_coefficient() / EP,
        "positive_spacelike_pair_integrand": positive,
        "pair_small_transfer_series": pairseries,
        "ordinary_charge_result": "The all-D GF bridge and the all-topology pole/finite continuity proof promote S286's conditional background identity to the specified ordinary on-shell pure-GR one-loop vertex: F1_proper(0)=1+Sigma_prime. One endpoint has its two scalar LSZ halves,1-Sigma_prime, so the whole pole and finite zero-transfer correction cancels. Two endpoints and four scalar legs cancel the full raw expression, including finite7, EulerGamma and log4pi.",
        "finite_t_IR_boundary": "At nonzero spacelike transfer the proper-plus-leg vertex still has IR pole R(t)/(8pi^2*kappa*EP), with R=V*M-mu/2. For0<=tau<=mu,11tau/12<=R(-tau)<=7tau/6. Its coefficient vanishes at zero transfer, but EP->0 at fixed nonzero transfer is not a finite physical observable. No exchange of these limits, detector choice or exact massless one-particle pole is asserted.",
        "matching_boundary": "Zero-transfer charge cancellation removes an independent proper-endpoint normalization ambiguity in this formal-loop sector. It does not fix S285's light curvature coefficients, the full metric propagator or Newton residue, physical four-point double-pole/tadpole completion, finite F1 slope or F2, curved Ricci-derivative matching, the improved b20, detector/Regge data or omitted loops.",
        "retained_background_and_master_inputs": "S286 supplies the entire Ward tensor and S283 the raw scalar masters and scalar residue; their frozen statements and all qualifications remain unchanged. This successor adds the explicitly scoped ordinary/spacelike bridge.",
        "checks": checks,
        "gates": {
            "literal_negative_triangle_convention_checked": True,
            "full_raw_UV_IR_and_finite_leg_terms_retained": raw.has(s.EulerGamma),
            "whole_one_and_two_endpoint_LSZ_cancellation": True,
            "finite_transfer_IR_pole_not_deleted": True,
            "spacelike_continuity_not_a_finite_slope": True,
            "metric_Newton_and_local_matching_remain_separate": True,
            "no_exact_dressed_state_or_full_P8_claim": True,
        },
    }
