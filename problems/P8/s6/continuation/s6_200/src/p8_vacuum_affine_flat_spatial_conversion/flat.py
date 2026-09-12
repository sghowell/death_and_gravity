"""Actual flat sine-kernel extraction and nonlocal conversion diagnostic."""

from functools import cache

import sympy as s
from p8_vacuum_affine_flat_tensor_cut import dispersion

from . import tensor as ten


@cache
def data():
    e, w, q, sigma = s.symbols("E w q sigma", positive=True)
    t = s.Symbol("t", real=True)
    f = s.Function("f")(t)
    partial = f / e**2 - s.diff(f, t, 2) / e**4 + s.diff(f, t, 4) / e**6
    checks = {
        "complete_sine_kernel_equation": s.expand(
            s.diff(partial, t, 2) + e**2 * partial - f - s.diff(f, t, 6) / e**6
        ),
        "six_time_derivative_remainder_transform": s.factor(
            1 / (e**2 - w)
            - 1 / e**2
            - w / e**4
            - w * w / e**6
            - w**3 / (e**6 * (e**2 - w))
        ),
    }
    D = s.diag(1, -1, 0)
    vec = s.Matrix(list(D))
    sub = {ten.P[0]: 0, ten.P[1]: 0, ten.P[2]: s.sqrt(q)}
    checks["complete_transverse_tensor_channel"] = s.expand(
        (vec.T * ten.numerator(2) * vec)[0].subs(sub) - 2 * ten.Z**2
    )
    checks["complete_transverse_scalar_channel"] = s.expand(
        (vec.T * ten.numerator(0) * vec)[0].subs(sub)
    )
    J0 = 1 / (sigma**3 * (sigma + q))
    J1 = 1 / (sigma**3 * (sigma + q) ** 2)
    J2 = 1 / (sigma**3 * (sigma + q) ** 3)
    checks["TT_first_time_coefficient"] = s.factor(
        3 * q * q * J0
        - q**3 * J1
        - q * q * (3 * sigma + 2 * q) / (sigma**3 * (sigma + q) ** 2)
    )
    checks["TT_second_time_coefficient"] = s.factor(
        -3 * q * J0
        + 3 * q * q * J1
        - q**3 * J2
        + q * (3 * sigma**2 + 3 * sigma * q + q * q) / (sigma**3 * (sigma + q) ** 3)
    )
    return {
        "physical_flat_kernel": "The complete spatial current spectral sine kernel is theta(t) integral sum_i rho_i(s)Pi_i(s,p) sin(sqrt(s+q)t)/sqrt(s+q) ds. Pi_i=Q_i/s^2; the S199 normalization already includes the full physical modes, Wick exchange and i/4.",
        "exact_retarded_extraction": "For y=G_E,ret Gamma with E^2=s+q, y''+E^2 y=Gamma and common zero initial jets give y=Gamma/E^2-Gamma''/E^4+Gamma''''/E^6-G_E,ret Gamma^(6)/E^6. The final transform is exactly w^3/[E^6(E^2-w)].",
        "no_curved_odd_term_discard": "The flat real tensor spectral sine kernel has only these even endpoint derivatives. This follows after the full flat tensor cut is reconstructed, not by dropping odd terms in the curved complex-amplitude S198 identity.",
        "nonlocality": "For a unit-norm spatial TT direction transverse to p, A0_TT=-q^3 integral rho2(s)/[s^3(s+q)]ds. Around q=0 its q^(3+n) coefficient is(-1)^(n+1) moment(2,3+n), nonzero for every n>=0. Therefore the finite equal-time conversion is not a finite spatial differential counterterm.",
        "positive_TT_moments": tuple(dispersion.moment(2, 3 + n) for n in range(4)),
        "regulator_boundary": "The spectral integrals here and their tails define continuum flat quantities. They do not identify the original contact one-mode band with the memory two-mode band, or assert that a finite hard-band subtraction has the same physical local polynomial.",
        "scope": "A full spatial flat nonlocal-representative conversion, not the curved CD endpoint/contact matching, fixed finite local action, full inverse, interacting background, cutoff or original V/G/B.",
        "checks": checks,
        "gates": {
            "complete_even_time_order_four": s.diff(partial, t, 4).has(s.diff(f, t, 8)),
            "all_nonzero_displayed_TT_moments": all(
                dispersion.moment(2, 3 + n) > 0 for n in range(4)
            ),
            "same_full_tensor_nonzero_transfer": ten.Q != 0,
        },
    }
