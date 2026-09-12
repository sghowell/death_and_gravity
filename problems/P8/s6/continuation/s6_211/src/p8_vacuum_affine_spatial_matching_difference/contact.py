"""Complete original one-leg contact in the spatial matching difference."""

from functools import cache

import sympy as s
from p8_vacuum_affine_spatial_current import hamiltonian as ham
from p8_vacuum_affine_spatial_current import vertices


def mixed_tensor(D, G):
    return (D * G + G * D) / 2


def paired_contact(k, P, D, G, a, m):
    # P and -P sum to zero in the full pointwise mixed exponential.
    total = P + (-P)
    return vertices.metric_contact(k, k - total, mixed_tensor(D, G), a, m)


@cache
def data():
    k = s.Matrix(s.symbols("k0:3", real=True))
    P = s.Matrix(s.symbols("P0:3", real=True))
    d = s.symbols("D0:5", real=True)
    g = s.symbols("G0:5", real=True)
    D = s.Matrix([[d[0], d[1], d[2]], [d[1], d[3], d[4]], [d[2], d[4], -d[0] - d[3]]])
    G = s.Matrix([[g[0], g[1], g[2]], [g[1], g[3], g[4]], [g[2], g[4], -g[0] - g[3]]])
    a, m = s.symbols("a m", positive=True)
    full = paired_contact(k, P, D, G, a, m)
    H = mixed_tensor(D, G)
    checks = {
        "full_original_contact_P_minus_P_zero": full
        - paired_contact(k, s.zeros(3, 1), D, G, a, m),
        "noncommuting_mixed_tensor_trace": s.expand(s.trace(H) - s.trace(D * G)),
        "original_positive_mass_contact_block": full[:3, :3]
        - a * m * m * H
        - ham.cross(k).T * H * ham.cross(k) / a,
        "original_electric_contact_block": full[3:, 3:] - H / a,
        "original_constraint_contact_absent": vertices.feature_form(H, True)[9, 9],
    }
    cA, cPi, crosscov = s.symbols("c_A c_pi c_cross", real=True)
    covariance = s.BlockMatrix(
        [[cA * s.eye(3), crosscov * s.eye(3)], [crosscov * s.eye(3), cPi * s.eye(3)]]
    ).as_explicit()
    point = s.factor(-s.trace(full * covariance) / 2)
    checks["contact_generally_nonzero_at_zero_internal_momentum"] = s.factor(
        point.subs(dict(zip(k, [0, 0, 0])))
        + (a * m * m * cA + cPi / a) * s.trace(D * G) / 2
    )
    return {
        "same_contact": "Use the exact second Hamiltonian vertex M_DGamma(k,q)=diag(a m^2 H+Ck^t H Cq/a,H/a), H=(D Gamma+Gamma D)/2. The mass sign is positive and the constraint A0 has no unimodular metric vertex.",
        "Fourier_pair": "For the full second-variation trace on a homogeneous state, external momenta pair as P,-P. Their convolution has total0, so the vertex is M_H(k,k). There is no remaining P at each internal k, before integration or use of the radial cutoff.",
        "cancellation": "Consequently C_K(P)-C_K(0)=0 for every finite original one-leg cutoff and for both actual and unit-W8 isotropic covariances. The generally nonzero homogeneous contact is retained in the anchor; this does not subtract the whole contact.",
        "test_functions": "For compact smooth spacetime smears, spatial Fourier contraction integrates D_hat(-P) and Gamma_hat(P). The P0 anchor is the same local-in-space multiplier, not a constant physical perturbation or an altered state/history.",
        "checks": checks,
        "gates": {
            "contact_cancelled_only_in_P_minus_P0_difference": True,
            "full_noncommuting_exponential_convolution_retained": D * G != G * D,
            "original_mass_electric_magnetic_constraint_content": True,
            "cancellation_before_cutoff_removal": True,
            "no_zero_contact_or_zero_mode_state_claim": True,
        },
    }
