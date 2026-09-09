"""Blockwise leading seven-mode response and an exact-frame geometric corollary."""

from functools import cache

import sympy as sp
from p8_proca_mean_response import band, mean
from p8_proca_mean_response import tails as proca_bounds
from p8_tensor_mean_source import bounds as tensor_bounds

from . import bounds, calibration

eta_s, eta_t, eta_p = sp.symbols("eta_scalar eta_tensor eta_proca", nonnegative=True)


@cache
def data():
    s = bounds.data()
    t = tensor_bounds.data()
    p = proca_bounds.data()
    coeffs = {
        "lapse": (
            sp.Integer(10**64),
            t["global_lapse_upper_per_eta"],
            p["global_lapse_upper_per_eta"],
        ),
        "actual_frame_log_scale": (
            sp.Integer(2 * 10**64),
            t["global_actual_frame_log_scale_upper_per_eta"],
            p["global_actual_frame_log_scale_upper_per_eta"],
        ),
        "matter_field": (
            sp.Integer(10**64),
            t["global_matter_field_upper_per_eta"],
            p["global_matter_field_upper_per_eta"],
        ),
    }
    totals = {
        key: sum(a * b for a, b in zip(row, (eta_s, eta_t, eta_p)))
        for key, row in coeffs.items()
    }
    upper_sub = {
        eta_s: calibration.DEFAULT_AMPLITUDE,
        eta_t: band.DEFAULT_AMPLITUDE,
        eta_p: band.DEFAULT_AMPLITUDE,
    }
    lapse_max = totals["lapse"].subs(upper_sub)
    scale_max = totals["actual_frame_log_scale"].subs(upper_sub)
    # Positive Proca and tensor scale jets, and positive Proca lapse jet,
    # are inherited. N>=1/2 bounds the negative tensor/scalar lapse jets.
    bounce_lower = (
        4 - 3 * lapse_max - sp.Rational(3, 2) * lapse_max**2 - 240 * eta_s - 30 * eta_t
    )
    u = mean.u
    tau = 1 + u * u
    density = (
        10**64 * eta_s / tau**4
        + (
            t["global_scalar_density_upper_per_eta"] * eta_t
            + p["global_scalar_density_upper_per_eta"] * eta_p
        )
        / tau**12
    )
    return {
        "independent_amplitudes": (eta_s, eta_t, eta_p),
        "exact_amplitude_box": upper_sub,
        "rounded_global_response_coefficients": coeffs,
        "global_response_upper": totals,
        "complete_physical_matter_density_and_pressure_absolute_upper": density,
        "joint_maximum_lapse_bound": lapse_max,
        "joint_maximum_exact_frame_log_scale_bound": scale_max,
        "center_exact_log_scale_second_derivative_lower": bounce_lower,
        "bounds": {
            "scalar_phase_bound_is_below_rounded_joint_value": s[
                "weighted_mean_phase_upper_per_eta"
            ]
            < 10**64,
            "scalar_lapse_bound_is_below_rounded_joint_value": s[
                "global_lapse_upper_per_eta"
            ]
            < 10**64,
            "scalar_actual_frame_bound_is_below_rounded_joint_value": s[
                "global_actual_frame_scale_upper_per_eta"
            ]
            < 2 * 10**64,
            "scalar_field_bound_is_below_rounded_joint_value": s[
                "global_matter_field_upper_per_eta"
            ]
            < 10**64,
            "joint_lapse_is_below_one_e_minus_ten": lapse_max < sp.Rational(1, 10**10),
            "joint_exact_frame_scale_is_below_one_e_minus_ten": scale_max
            < sp.Rational(1, 10**10),
            "joint_lapse_is_globally_above_one_half": 1 - lapse_max > sp.Rational(1, 2),
            "joint_center_exact_scale_acceleration_exceeds_three": bounce_lower.subs(
                upper_sub
            )
            > 3,
        },
    }


def calibrated(
    scalar=calibration.DEFAULT_AMPLITUDE,
    tensor=band.DEFAULT_AMPLITUDE,
    proca=band.DEFAULT_AMPLITUDE,
):
    values = {
        eta_s: calibration.exact_amplitude(scalar),
        eta_t: band.exact_amplitude(tensor),
        eta_p: band.exact_amplitude(proca),
    }
    return {
        key: sp.factor(value.subs(values))
        for key, value in data()["global_response_upper"].items()
    }


@cache
def identities():
    bg = mean.data()
    xi = sp.symbols("xi_s xi_t xi_p", real=True)
    dp = sp.symbols("dp_s dp_t dp_p", real=True)
    Fn = sp.symbols("Fn_s Fn_t Fn_p", real=True)
    Fx = sp.symbols("Fx_s Fx_t Fx_p", real=True)
    Fp = sp.symbols("Fp_s Fp_t Fp_p", real=True)
    Fc = sp.symbols("Fc_s Fc_t Fc_p", real=True)

    def flow(x, d, nforce, xforce, pforce, cforce):
        lapse = (bg["alpha"] * d - 3 * bg["ell"] * bg["beta"] * x + nforce) / (
            2 * bg["J"]
        )
        return sp.Matrix(
            [
                lapse,
                -d / 2 + bg["alpha"] * lapse / 3 + xforce,
                -3 * bg["H"] * d
                + bg["ell"] * bg["beta"] * lapse
                - 3 * bg["ell"] ** 2 * x
                + pforce,
                bg["beta"] * lapse - 3 * bg["ell"] * x + cforce,
                -3 * bg["ell"] ** 2 * (x + lapse / (2 * bg["h"])),
            ]
        )

    separate = sum((flow(*row) for row in zip(xi, dp, Fn, Fx, Fp, Fc)), sp.zeros(5, 1))
    combined = flow(*(sum(row) for row in (xi, dp, Fn, Fx, Fp, Fc)))
    # Scalar-tensor quadratic cross contractions vanish in every direction.
    kx, ky, kz = sp.symbols("kx ky kz", real=True)
    kv = sp.Matrix([kx, ky, kz])
    q = (kv.T * kv)[0]
    Pi = sp.eye(3) - kv * kv.T / q
    e = sp.symbols("e0:6", real=True)
    E = sp.Matrix([[e[0], e[1], e[2]], [e[1], e[3], e[4]], [e[2], e[4], e[5]]])
    TT = Pi * E * Pi - sp.trace(Pi * E) * Pi / 2
    return {
        "physical_mean_and_induced_matter_observable_superpose_once": (
            combined - separate
        ).applyfunc(sp.factor),
        "TT_has_no_isotropic_scalar_contraction": sp.factor(sp.trace(TT)),
        "TT_has_no_scalar_momentum_contraction": (TT * kv).applyfunc(sp.factor),
        "TT_has_no_longitudinal_York_shear_contraction": sp.factor(
            sp.trace(TT * (kv * kv.T / q - sp.eye(3) / 3))
        ),
    }
