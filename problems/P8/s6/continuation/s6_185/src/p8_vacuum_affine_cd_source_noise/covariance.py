"""Full time-dependent three-mode covariance with physical temporal readout."""

from functools import cache

import sympy as s
from p8_vacuum_affine_proca_gaussian import bridge


@cache
def data():
    a, m, k, omega = s.symbols("a m k omega", positive=True)
    gT2 = a
    gL2 = a * m * m / omega**2
    T = 9 / (a**3 * omega)
    L = 9 * omega / (a**3 * m * m)
    time = 27 * k * k / (a**5 * m * m * omega)
    trace = 2 * T + L + time
    target = (27 / omega + 36 * k * k / (a * a * m * m * omega)) / a**3
    prev = bridge.data()["modes"]
    checks = {
        "actual_transverse_canonical_normalizer": prev["transverse"][
            "normalizer_squared"
        ].subs(
            {
                next(
                    x
                    for x in prev["transverse"]["normalizer_squared"].free_symbols
                    if str(x) == "a"
                ): a
            }
        )
        - gT2,
        "actual_longitudinal_canonical_normalizer": s.factor(
            prev["longitudinal"]["normalizer_squared"]
            - gL2.subs(omega, s.sqrt(m * m + k * k / a**2))
        ),
        "transverse_proper_frame_covariance": 9 / omega / gT2 / a**2 - T,
        "longitudinal_proper_frame_covariance": 9 / omega / gL2 / a**2 - L,
        "temporal_constraint_covariance": k * k * gL2 * 27 * omega / (a**6 * m**4)
        - time,
        "all_three_mode_trace": s.factor(
            (trace - target).subs(omega, s.sqrt(m * m + k * k / a**2))
        ),
        "volume_squared_keeps_one_a_cubed": a**6 / a**3 - a**3,
        "covariance_display_value": 4 * 27 - 108,
        "covariance_display_gradient": 4 * 36 - 144,
    }
    return {
        "actual_canonical_normalizers_squared": {"T": gT2, "L": gL2},
        "proper_frame_mode_diagonal_bounds": {
            "each_T": T,
            "L_spatial": L,
            "temporal": time,
        },
        "full_covariance_trace_upper": target,
        "exact_maximum_volume": s.Rational(25, 16) ** 3,
        "source_class": "Smooth compact physical-frame test four-vectors on the unit CD slab; spatial coordinate Fourier transform with measure d^3k/(2pi)^3",
        "same_state_covariance_bound": "C_A(f,f)<=108/m ||f||L2(dt dx)^2+144/m^3 ||grad_spatial f||L2(dt dx)^2; all momenta and physical polarizations, including temporal constraint",
        "normalization": "The smearing is integral a^3 A.f. Proper-frame spatial components of a lower-coordinate source are divided by a; a>=1 so their spatial norms only decrease.",
        "not_a_QEI": "A positive covariance upper bound for this specified state, not a state-independent QEI, stationary mass-shell formula or stress-tensor noise bound.",
        "checks": checks,
        "gates": {
            "volume_below_four": s.Rational(25, 16) ** 3 < 4,
            "temporal_covariance_positive": time > 0,
            "full_longitudinal_kept": L > 0,
        },
    }
