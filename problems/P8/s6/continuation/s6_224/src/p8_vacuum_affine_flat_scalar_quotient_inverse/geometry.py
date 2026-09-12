"""Full flat scalar projectors, ordered gauge quotient and fixed Hessian."""

from functools import cache

import sympy as s
from p8_proca_rank_one_inverse import spectral as trace
from p8_vacuum_affine_flat_spatial_conversion import tensor as spatial
from p8_vacuum_affine_isolated_shear_resolvent import normalization as shear
from p8_vacuum_affine_local_tensor_response import prescription

lam, q = s.symbols("Laplace_frequency spatial_momentum_squared", positive=True)
p = lam * lam + q
B = s.Matrix([[lam * lam + 2 * q / 3, -lam * lam / 3], [-q, -lam * lam]])
M = s.Matrix([[q / 3, lam * lam + 2 * q / 3, -lam / 3], [q, -q, -lam]])
T = s.Matrix([[lam, 0, 0], [0, 1, 0], [q, 0, lam]])
f0, f2 = s.symbols("trace_factor shear_factor", nonzero=True)


@cache
def data():
    basis = [2 * s.eye(3), -2 * s.diag(1, 0, 0)]
    columns = s.Matrix.hstack(*(s.Matrix(list(Q)) for Q in basis))
    mapping = {spatial.P[0]: s.sqrt(q), spatial.P[1]: 0, spatial.P[2]: 0}
    projected = {
        spin: (
            columns.T * spatial.numerator(spin, -p).subs(mapping) * columns
        ).applyfunc(s.expand)
        for spin in (0, 2)
    }
    scalar = B[0:1, :]
    weyl = B[1:2, :]
    weight = s.diag(f0, s.Rational(8, 3) * f2)
    full3 = M.subs(lam, -lam).T * weight * M
    full2 = B.T * weight * B
    gauge = s.Matrix([lam, 0, q])
    quotient = s.zeros(2, 3)
    quotient[:, 1:3] = B
    charted = s.zeros(3)
    charted[1:3, 1:3] = full2
    fixed = prescription.data()["literal_fixed_finite_density_before_64_pi_squared"]
    names = {str(v): v for v in fixed.free_symbols}
    C = s.expand(fixed).coeff(names["Weyl_squared"])
    R2 = s.expand(fixed).coeff(names["R_old"], 2)
    metric_gram = columns.T * columns
    checks = {
        "full_S200_spin0_projector_not_COM_replacement": projected[0]
        - 12 * scalar.T * scalar,
        "full_S200_spin2_projector_not_COM_replacement": projected[2]
        - s.Rational(8, 3) * weyl.T * weyl,
        "fixed_Weyl_Hessian_matches_S223": C
        - shear.data()["actual_fixed_fourth_coefficient"],
        "fixed_scalar_curvature_Hessian_matches_S86": 72 * R2
        + trace.data()["H_at_zero"],
        "ordered_source_gauge_kernel": M * gauge,
        "ordered_detector_gauge_kernel": M.subs(lam, -lam) * gauge.subs(lam, -lam),
        "complete_source_gauge_quotient": M * T - quotient,
        "complete_detector_gauge_quotient": M.subs(lam, -lam) * T.subs(lam, -lam)
        - quotient,
        "full_three_source_ordered_quotient_block": (
            T.subs(lam, -lam).T * full3 * T - charted
        ).applyfunc(s.expand),
        "quotient_determinant": s.factor(B.det() + lam * lam * p),
        "complete_gauge_chart_determinant": T.det() - lam * lam,
        "physical_spatial_metric_Gram": metric_gram - s.Matrix([[12, -4], [-4, 4]]),
        "physical_metric_lower_margin": (metric_gram - 2 * s.eye(2)).det() - 4,
        "physical_metric_upper_margin": (14 * s.eye(2) - metric_gram).det() - 4,
        "constant_coefficient_quotient_formal_transpose": full2 - full2.T,
        "full_ordered_three_source_adjoint_parity": full3 - full3.subs(lam, -lam).T,
    }
    return {
        "curvature_channel_map_B": B,
        "full_three_source_channel_map": M,
        "ordered_gauge_chart": T,
        "full_three_source_reference": full3.applyfunc(s.expand),
        "two_channel_quotient_reference": full2.applyfunc(s.expand),
        "fixed_finite_channel_coefficients": {"trace": 72 * R2, "shear": C},
        "metric_Gram": metric_gram,
        "normalization": "For Q=2zeta I-2c Pi and p=lambda^2+q, p^2 Pi0=12 S_D S_G and p^2 Pi2=(8/3)W_D W_G. The full reference is B^T diag(Ftrace(p),(8/3)F2(p)) B, with both original finite coefficients.",
        "ordered_gauge": "In(n,zeta,b), S=q n/3+(D^2+2q/3)zeta-D b/3 and W=q(n-zeta)-D b. Source n=D eta,b=D c+q eta and the opposite detector derivatives give the same B quotient and a zero gauge block.",
        "scope": "The full three-source block is not invertible. The flat two-channel spatial quotient is the target. Its q0 continuous extension is not a solution of the literal homogeneous lapse/shift constraint. Curved one-current contacts and matter/tree terms remain outside this reference.",
        "checks": checks,
        "gates": {
            "both_full_spatial_spin_sectors_retained": projected[0] != s.zeros(2)
            and projected[2] != s.zeros(2),
            "full_three_source_has_a_gauge_kernel": full3.det().expand() == 0,
            "quotient_generically_invertible": s.factor(full2.det()) != 0,
            "gauge_chart_has_no_inverse_transfer": all(
                not entry.has(1 / q) for entry in T
            ),
            "full_three_source_not_symmetric_at_fixed_Laplace_frequency": s.factor(
                full3[0, 2] - full3[2, 0]
            )
            != 0,
            "physical_metric_norm_lower_diagonal_positive": (
                metric_gram - 2 * s.eye(2)
            )[0, 0]
            > 0,
            "physical_metric_norm_upper_diagonal_positive": (
                14 * s.eye(2) - metric_gram
            )[0, 0]
            > 0,
        },
    }
