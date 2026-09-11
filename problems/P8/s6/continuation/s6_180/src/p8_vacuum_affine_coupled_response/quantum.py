"""Uncancelled ordinary-Proca current, exact time covariance and local Ward channel."""

from functools import cache

import sympy as s
from p8_proca_nonlocal_response import source
from p8_proca_prepared_inverse import coordinates as old_coordinates

u = s.symbols("time", real=True)
H = s.Function("H")(u)
rho, P = (s.Function(name)(u) for name in ("scaled_rho", "scaled_pressure"))
eta, w = (s.Function(name)(u) for name in ("prepared_time_shift", "invariant_scale"))
R, Q = (
    s.Function(name)(u)
    for name in ("pure_scale_energy_response", "pure_scale_pressure_response")
)
S = rho + P


def D(value):
    return s.diff(value, u) + 3 * H * value


def ward(value):
    expr = s.expand(value.doit())
    expr = expr.subs(s.diff(R, u), -3 * H * (R + Q) - 3 * S * s.diff(w, u))
    expr = expr.subs(
        s.diff(rho, u, 2),
        -3 * s.diff(H, u) * S - 3 * H * (s.diff(rho, u) + s.diff(P, u)),
    )
    return s.factor(expr.subs(s.diff(rho, u), -3 * H * S))


@cache
def data():
    n = s.diff(eta, u)
    v = w + H * eta
    drho = R + s.diff(rho, u) * eta
    dP = Q + s.diff(P, u) * eta
    physical = s.Matrix([-drho - 3 * rho * v, 3 * dP + 3 * P * (n + 3 * v)])
    transformed = s.Matrix([-D(physical[0]) + H * physical[1], physical[1]])
    local_eta = -3 * P * s.diff(w, u) - 3 * s.diff(H, u) * P * eta
    local_w = 3 * P * s.diff(eta, u) + 3 * (s.diff(P, u) + 3 * H * P) * eta + 9 * P * w
    full = s.Matrix([local_eta, 3 * Q + local_w])
    density = (
        -s.Rational(3, 2) * s.diff(H, u) * P * eta**2
        - 3 * P * eta * s.diff(w, u)
        + s.Rational(9, 2) * P * w * w
    )
    euler = lambda f: s.diff(density, f) - D(s.diff(density, s.diff(f, u)))
    checks = {
        "actual_uncancelled_time_channel_from_full_varied_Ward": ward(
            transformed[0] - local_eta
        ),
        "actual_uncancelled_scale_channel": s.expand(transformed[1] - 3 * Q - local_w),
        "local_cross_density_eta": s.factor(euler(eta) - local_eta),
        "local_cross_density_scale": s.factor(euler(w) - local_w),
        "time_channel_has_no_eta_derivative": s.diff(local_eta, s.diff(eta, u)),
        "time_channel_no_fourth_derivative": s.diff(local_eta, s.diff(eta, u, 4)),
        "time_scale_channel_at_most_first": s.diff(local_eta, s.diff(w, u, 2)),
        "only_scale_nonlocal_channel": s.diff(full[0], R) + s.diff(full[0], Q),
        "new_pure_scale_minus_old_physical_pressure_is_local": s.factor(
            full[1].subs({eta: 0, s.diff(eta, u): 0}) - 3 * Q - 9 * P * w
        ),
        "prepared_source_map_roundtrip": s.expand(v - H * eta - w),
    }
    # The old total contains a DIFFERENT already-fixed scalar profile;
    # its extra local terms are not present in the current parent.
    old_eta = (
        S * s.diff(eta, u, 2)
        + (s.diff(S, u) + 3 * H * S) * n
        - 3 * s.diff(H, u) * S * eta
        - 3 * S * s.diff(w, u)
    )
    old_scale = 3 * Q + 3 * s.diff(P, u) * eta + 3 * S * n
    return {
        "normalization": "rho,P are 64*pi^2 times the actual unscaled physical conditional vector stress. Q_vector=64*pi^2 delta_e and gamma=1/(64*pi^2*kappa); the classical normalized tree has kappa divided out.",
        "exact_time_change": "For a prepared diffeomorphism psi=t+epsilon eta, N=psi',a=a0(psi),f_N=f0(psi)/sqrt(psi'),p_N=sqrt(psi')p0(psi). The same Cauchy data and covariant finite prescription give delta rho=rho' eta and delta P=P' eta.",
        "physical_sources": s.ImmutableMatrix([n, v]),
        "uncancelled_physical_current_variation": s.ImmutableMatrix(physical),
        "local_time_channel": local_eta,
        "local_scale_channel_addition": local_w,
        "complete_adapted_vector_response": s.ImmutableMatrix(full),
        "local_cross_density": density,
        "difference_from_OLD_profile_cancelled_time_channel": s.expand(
            local_eta - old_eta
        ),
        "difference_from_OLD_profile_cancelled_scale_channel": s.expand(
            full[1] - old_scale
        ),
        "checks": checks,
    }


@cache
def mode_covariance():
    # Re-evaluate the physically matched mode identities, not the old
    # profile-total current formula. This routine contains no profile.
    checks = old_coordinates.mode_checks()
    result = {
        "literal_current_time_covariance_" + key: value for key, value in checks.items()
    }
    e = s.Function("time_diffeomorphism")(source.u)
    mapping = {field: s.diff(e, source.u, j + 1) for j, field in enumerate(source.n)}
    mapping.update(
        {field: s.diff(source.H * e, source.u, j) for j, field in enumerate(source.v)}
    )
    for sector in ("T", "L"):
        rate = source.data(sector)["rg"].subs(mapping, simultaneous=True)
        expected = source.data(sector)["rate"] * e - s.diff(e, source.u) / 2
        result[sector + "_full_moving_coordinate_normalization"] = s.factor(
            rate - expected
        )
    return {
        "source_policy": "Only S6.87's exact ordinary mode covariance/readout identities are recomputed. Its scalar-profile current formulas and full-parent inverse are not copied.",
        "checks": result,
    }
