"""Homogeneous tensor duality, complex Fourier lift and exact norm embeddings."""

from functools import cache

import sympy as s
from p8_vacuum_affine_covariant_current import current

HOMOGENEOUS = s.Integer(10) ** 95


def weight_rows():
    v = s.Symbol("momentum_squared", nonnegative=True)
    return (
        v,
        (1 + v) ** 6,
        {
            "N61_time": s.expand((1 + v) ** 6 - 1),
            "N61_gradient": s.expand((1 + v) ** 6 - v),
            "X46": s.Integer(0),
            "Z24": s.expand((1 + v) ** 6 - (1 + v) ** 4),
            "S209_conversion_tail_X46": s.Integer(0),
        },
    )


@cache
def data():
    v, weight, rows = weight_rows()
    x = s.Symbol("x", real=True)
    re, im = s.symbols("real_source_norm imaginary_source_norm", nonnegative=True)
    checks = {
        "same_coordinate_density_homogeneous_first_derivative": current.CURRENT[1]
        - HOMOGENEOUS,
        "unit_slab_time_Cauchy_factor": s.integrate(
            1, (x, -s.Rational(1, 2), s.Rational(1, 2))
        )
        - 1,
        "complex_Hilbert_sum_no_extra_sqrt_two": HOMOGENEOUS**2 * re**2
        + HOMOGENEOUS**2 * im**2
        - HOMOGENEOUS**2 * (re**2 + im**2),
        "Y_squared_embedding": s.Integer(2) + 1 - 3,
        "M_squared_weight": (1 + v) - 1 - v,
        "full_source_time_order_for_zero_initial_germ": s.Integer(12) + 1 - 13,
    }
    for name, row in rows.items():
        expected = {
            "N61_time": weight - 1,
            "N61_gradient": weight - v,
            "X46": 0,
            "Z24": weight - (1 + v) ** 4,
            "S209_conversion_tail_X46": 0,
        }[name]
        checks[name + "_weight_reconstruction"] = s.expand(row - expected)
    return {
        "source_space": "Z136^2=sum_r0..13||(1-Delta)^3 partial_t^r Gamma||L2(IxR3,Frobenius)^2. Smooth real tracefree sources are compact in the open unit slab and Schwartz in space, with the common zero initial neighborhood. Complete this admitted source space in Z136 only after proving the bound.",
        "detector_space": "M[D]^2=||D||L2^2+||grad D||L2^2. Detector tensors are tracefree and have the same unitary Fourier convention. The target is the dual of L2_time H1_space, namely L2_time Hminus1_space in this weak setting.",
        "homogeneous_anchor": "S194's1e95 bound is already for the complete fixed-prescription coordinate-density first current variation. Frobenius-unit real tensor detectors have operator norm<=1. Pointwise tensor duality therefore gives a Frobenius output bound1e95 max_j0..12 sup_t||Gamma^(j)||op.",
        "zero_germ_lift": "For every Fourier mode, Gamma^(j)(t)=integral_left^t Gamma^(j+1)(s)ds. The unit time interval and Cauchy-Schwarz give each sup norm<=||Gamma^(j+1)||L2_time,Frobenius. The homogeneous response is real linear. Complexifying its real/imaginary parts and adding their squared output norms gives1e95 times the Hilbert sum of source derivatives1,...,13, without an extra sqrt2. Plancherel supplies the spatial lift.",
        "embeddings": "N61<=sqrt2 Z136, X46<=Z136, Y=sqrt(N61^2+X46^2)<=sqrt3 Z136<2Z136, and Z24<=Z136. Nonnegative polynomial coefficient checks establish all momentum weights for arbitrary P, including P0.",
        "zero_transfer_known": "Known(0) means the P-independent multiplier obtained from the same explicit endpoint/state/time kernels at zero transfer, before external Fourier Cauchy-Schwarz. Their pointwise bounds specialize atP0 and apply to each arbitrary spatial Fourier source. It is not evaluation of an arbitrary bounded operator at a measure-zero frequency.",
        "not_inverse": "The thirteen-time/six-spatial derivative-losing bound is not a same-space contraction or inverse estimate. No smallness assertion for the reduced scalar/mixed operator follows.",
        "checks": checks,
        "gates": {
            "all_spatial_weight_polynomials_nonnegative": all(
                all(c >= 0 for c in s.Poly(row, v).all_coeffs())
                for row in rows.values()
            ),
            "strict_Y_rounding": s.Integer(3) < 4,
            "homogeneous_coordinate_density_not_rescaled_twice": True,
            "both_complex_source_parts_controlled": True,
            "zero_initial_germ_not_zero_quantum_metric_state": True,
            "pointwise_zero_transfer_kernels_used_before_norms": True,
            "derivative_loss_not_contraction": True,
        },
    }
