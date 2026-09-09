"""Independent closed center invariant calculation, without labelled geometry."""

from functools import cache

import sympy as sp


@cache
def data():
    k = sp.Symbol("positive_wave_scale", positive=True)
    v, p, chi, pc, T, PT, S, PS, Wx, Wy, Wz, Px, Py, Pz = sp.symbols(
        "v p chi pc T PT S PS Wx Wy Wz Px Py Pz", real=True
    )
    fields = (v, p, chi, pc, T, PT, S, PS, Wx, Wy, Wz, Px, Py, Pz)
    ell = sp.Rational(1, 10)
    Je = sp.Rational(243, 160)
    # The linear York correction is tracefree and fixes pi_xx=ell*chi/2.
    pi = sp.Matrix(
        [
            [ell * chi / 2, 0, 0],
            [0, (p - ell * chi) / 4 + PT / 2, PS / 2],
            [0, PS / 2, (p - ell * chi) / 4 - PT / 2],
        ]
    )
    dp1 = p / 3
    shear = (pi * pi).trace() - pi.trace() ** 2 / 3
    dc1 = pc - 3 * ell * v
    dc2 = sp.Rational(15, 2) * ell * v * v + ell * (T * T + S * S) / 2 - 3 * v * pc
    R1 = 4 * k * k * v
    # These are the ZERO-OUTPUT quadratic curvature coefficients, after IBP.
    R2 = -10 * k * k * v * v - k * k * (T * T + S * S) / 2
    gradient = k * k * chi * chi
    gauss = k * k * Px * Px
    electric = 10**6 * (Px * Px + Py * Py + Pz * Pz)
    magnetic = sp.Rational(2, 10**6) * k * k * (Wy * Wy + Wz * Wz)
    vector = Wx * Wx + Wy * Wy + Wz * Wz
    F1 = -dc1 / 20 + R1 / 4
    n1 = F1 / (2 * Je)
    F2 = (
        -dc2 / 20
        + R2 / 4
        - sp.Rational(9, 8) * dp1 * dp1
        - dc1 * dc1 / 4
        - gauss / 4
        + 3 * shear
        + electric / 4
        + magnetic / 8
        + 3 * vector / 4
        + 3 * gradient / 4
    )
    L2 = 3 * dc1 / 40 - 3 * R1 / 8
    n2 = (F2 + sp.Rational(13821, 800) * n1 * n1 + L2 * n1) / (2 * Je)
    common = (
        ell * dc2
        + dc1 * dc1 / 2
        - 3 * ell * dc1 * n1
        - sp.Rational(3, 2) * ell * ell * n2
        + 3 * ell * ell * n1 * n1
    )
    rho = sp.factor(common + gradient / 2)
    pressure = sp.factor(common - gradient / 6)
    missing = -sp.Rational(3, 2) * ell * ell * n2
    matrices = {
        key: sp.hessian(expr, fields).applyfunc(sp.factor)
        for key, expr in (
            ("rho", rho),
            ("pressure", pressure),
            ("missing_if_lapse_truncated", missing),
        )
    }
    # Direct center map to Y=(k*b,k*chi,p_tilde_b,p_chi), including BOTH
    # canonical matter boundary shifts; the gamma moving shift is zero here.
    C = sp.Matrix(
        [
            [0, 0, 1 / (2 * k * k), 0],
            [-2 * k, 3 * ell / k, 0, 0],
            [0, 1 / k, 0, 0],
            [0, 0, 3 * ell / (2 * k * k), 1],
        ]
    )
    packet = (C.T * matrices["rho"][:4, :4] * C).applyfunc(sp.factor)
    leading = packet.applyfunc(lambda val: sp.limit(val, k, sp.oo))
    negative = sp.Matrix([1340, 1, 0, 0])
    return {
        "fields": fields,
        "rho": rho,
        "pressure": pressure,
        "n1": sp.factor(n1),
        "n2": sp.factor(n2),
        "matrices": matrices,
        "center_packet_map": C,
        "actual_packet_density_Hessian": packet,
        "principal_packet_density_Hessian": leading,
        "negative_principal_configuration_test": sp.factor(
            (negative.T * leading * negative)[0]
        ),
        "checks": {
            "independent_linear_York_trace": sp.factor(pi.trace() - p / 2),
            "independent_linear_York_longitudinal_momentum": sp.factor(
                pi[0, 0] - ell * chi / 2
            ),
            "independent_center_retuned_lapse_Hessian": -2 * Je + sp.Rational(243, 80),
            "independent_tensor_momentum_density_Hessian": matrices["rho"][5, 5]
            + sp.Rational(2, 135),
            "independent_transverse_Proca_momentum_density_Hessian": matrices["rho"][
                12, 12
            ]
            + sp.Rational(200000, 81),
            "physical_scalar_principal_density_negative_control": sp.factor(
                (negative.T * leading * negative)[0] + sp.Rational(134, 135)
            ),
        },
    }
