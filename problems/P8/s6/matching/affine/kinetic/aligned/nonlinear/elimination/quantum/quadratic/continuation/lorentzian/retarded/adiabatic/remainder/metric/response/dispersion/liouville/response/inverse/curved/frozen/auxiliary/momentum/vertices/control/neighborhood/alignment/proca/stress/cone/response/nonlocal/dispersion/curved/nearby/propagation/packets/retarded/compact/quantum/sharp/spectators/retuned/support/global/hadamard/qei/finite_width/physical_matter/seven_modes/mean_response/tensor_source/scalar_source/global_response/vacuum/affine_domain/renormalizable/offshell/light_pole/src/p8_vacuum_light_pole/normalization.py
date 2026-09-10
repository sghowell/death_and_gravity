"""Full finite-regulator trace normalization before continuum continuation."""

from functools import cache

import sympy as sp
from p8_polynomial_vacuum import gaussian


@cache
def data():
    d = gaussian.data()
    phi = sp.Matrix(sp.symbols("phi1 phi2", real=True))
    G, L = sp.symbols("cubic_G quartic_L", real=True)
    a, b, c = sp.symbols("positive_K11 K12 positive_K22", real=True)
    kd, ke, kf = sp.symbols("light_K11 light_K12 light_K22", real=True)
    r, s, t = sp.symbols(
        "free_light_covariance11 free_light_covariance12 free_light_covariance22",
        real=True,
    )
    K = sp.Matrix([[a, b], [b, c]])
    Km = sp.Matrix([[kd, ke], [ke, kf]])
    covariance = sp.Matrix([[r, s], [s, t]])
    J = sp.Matrix([v * v for v in phi])
    actual_W = sp.hessian(d["effective_quartic_action"], tuple(phi)) - Km
    local = sp.diag(
        *[L * phi[i] ** 2 / 2 - G**2 * (K.inv() * J)[i] / 2 for i in range(2)]
    )
    mixed = -(G**2) * sp.diag(*phi) * K.inv() * sp.diag(*phi)
    halftrace = sp.trace(covariance * actual_W) / 2
    pure_tadpole = L * sum(covariance[i, i] * phi[i] ** 2 for i in range(2)) / 4
    heavy_stationary_tadpole = (
        -(G**2) * sum(covariance[i, i] * (K.inv() * J)[i] for i in range(2)) / 4
    )
    mixed_trace = (
        -(G**2) * sp.trace(covariance * sp.diag(*phi) * K.inv() * sp.diag(*phi)) / 2
    )
    return {
        "actual_light_quadratic_Hessian_insertion": actual_W,
        "local_insertion": local,
        "mixed_insertion": mixed,
        "half_trace": halftrace,
        "mixed_trace": mixed_trace,
        "continuum_boundary": "For translation-invariant kinetic kernels, the first two trace terms are momentum-independent local light mass terms. Only the mixed trace supplies the displayed nonlocal two-point bubble. All local terms cancel in the same on-shell affine subtraction.",
        "checks": {
            "actual_full_reduced_Hessian_insertion": (
                actual_W - local - mixed
            ).applyfunc(sp.factor),
            "one_loop_half_trace_includes_every_local_and_mixed_term": sp.factor(
                halftrace - pure_tadpole - heavy_stationary_tadpole - mixed_trace
            ),
            "mixed_term_not_replaced_by_field_independent_heavy_logdet": sp.factor(
                sp.diff(mixed_trace, G, 2) / 2
                + sp.trace(covariance * sp.diag(*phi) * K.inv() * sp.diag(*phi)) / 2
            ),
        },
    }
