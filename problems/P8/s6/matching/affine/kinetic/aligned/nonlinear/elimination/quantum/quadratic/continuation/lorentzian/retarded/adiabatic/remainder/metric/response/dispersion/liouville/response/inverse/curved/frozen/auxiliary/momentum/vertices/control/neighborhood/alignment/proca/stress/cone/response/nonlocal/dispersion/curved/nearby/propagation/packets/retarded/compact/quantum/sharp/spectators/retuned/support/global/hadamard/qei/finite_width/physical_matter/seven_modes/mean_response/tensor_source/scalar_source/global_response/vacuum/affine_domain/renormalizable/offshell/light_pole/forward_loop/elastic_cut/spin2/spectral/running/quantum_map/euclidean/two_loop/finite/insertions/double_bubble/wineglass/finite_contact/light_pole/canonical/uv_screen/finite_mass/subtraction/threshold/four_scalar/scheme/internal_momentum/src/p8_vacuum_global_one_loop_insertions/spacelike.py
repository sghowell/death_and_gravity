"""Global spacelike logarithmic envelope and exact concavity identities."""

from functools import cache

import sympy as sp

from . import halfplane


def envelope(mass, yukawa_upper, Q_lower, t):
    base = halfplane.enclosure(mass, yukawa_upper, Q_lower)
    t = halfplane.exact(t)
    if t < 0:
        raise ValueError("Require Euclidean momentum squared t>=0")
    m, Y, Q = (base[k] for k in ("fermion_mass", "Y_upper", "Q_lower"))
    d, K = t + 1, 4 * m * m - 1
    C = 2 * halfplane.N * Y / Q
    return {
        "t": t,
        "d": d,
        "K": K,
        "minus_fermion_remainder_upper": C * d * sp.log(1 + d / K),
        "insertion_decay_envelope": C * sp.log(1 + d / K) / d,
        "coarse_uniform_insertion_upper": C / K,
        "scope": "All t>=0, one-loop only. Logarithmic decay of the inserted propagator is not a complete outer-integral convergence proof.",
    }


@cache
def data():
    m, A, d = sp.symbols("m A d", positive=True)
    s = sp.symbols("s")
    K = 4 * m * m - 1
    b = A / (m * m - A)
    h = d * b
    R = -(K + d) * sp.log(1 + h) + K * h
    original = -(4 * m * m - s) * sp.log(1 - A * s / (m * m))
    second_original = sp.diff(original, s, 2).subs(s, 1 - d)
    z = sp.symbols("nonnegative_z", positive=True)
    H = z - sp.log(1 + z)
    checks = {
        "spacelike_remainder_second_derivative": sp.factor(
            sp.diff(R, d, 2) - second_original
        ),
        "spacelike_remainder_zero_value": R.subs(d, 0),
        "spacelike_remainder_zero_slope": sp.diff(R, d).subs(d, 0),
        "reference_gap_factorization": sp.factor(
            m * m - A * (1 - d) - (m * m - A) * (1 + d * b)
        ),
        "maximum_parameter_ratio_margin": sp.factor(
            1 / K - b - m * m * (1 - 4 * A) / (K * (m * m - A))
        ),
        "remainder_log_envelope_difference": sp.expand(
            -R - d * sp.log(1 + h) + K * (h - sp.log(1 + h))
        ),
        "log_below_argument_derivative": sp.factor(sp.diff(H, z) - z / (1 + z)),
        "log_below_argument_anchor": H.subs(z, 0),
        "log_positive_derivative": sp.diff(sp.log(1 + z), z) - 1 / (1 + z),
        "large_t_envelope_prefactor": sp.factor(
            (sp.log(1 + d / K) / d) * d - sp.log(1 + d / K)
        ),
    }
    return {
        "symbols": {"m": m, "A": A, "d": d},
        "spacelike_subtracted_parameter_kernel": R,
        "positive_gap": m * m - A,
        "b": b,
        "K": K,
        "sign": "For t>=0, d=t+1: f_R(-t)<0 when Y>0. Concavity follows from the half-plane resolvent decomposition.",
        "global_envelope": "0<=-f_R(-t)<=(2NY/Q)(t+1) log(1+(t+1)/(4m^2-1)).",
        "inserted_propagator_envelope": "|f_R(-t)|/(t+1)^2 <= (2NY/Q) log(1+(t+1)/K)/(t+1) <= (2NY/Q)/K.",
        "scope": "All spacelike momenta in the fixed one-loop functional, not an all-orders positive spectrum or a globally positive inverse at arbitrary energy.",
        "checks": checks,
    }
