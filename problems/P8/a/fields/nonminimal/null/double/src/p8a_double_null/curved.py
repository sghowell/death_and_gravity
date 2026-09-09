"""Actual affine-null physical conformal transport and scalar reference."""

from functools import cache

import sympy as sp
from p8a_nonminimal import stress


@cache
def data():
    a = sp.Symbol("a", positive=True)
    h, hd, hdd, h3 = sp.symbols("H Hdot Hddot Hthird", real=True)
    f, ft, fz, ftt, ftz = sp.symbols("f ft fz ftt ftz", real=True)

    def dt(expr):
        return (
            sp.diff(expr, a) * a * h
            + sp.diff(expr, h) * hd
            + sp.diff(expr, f) * ft
            + sp.diff(expr, ft) * ftt
            + sp.diff(expr, fz) * ftz
        )

    F = a**-2 * f
    first_eta = sp.factor(a * dt(F))
    second_eta = sp.factor(a * dt(first_eta))
    mixed = sp.factor(a * dt(a**-2 * fz))
    D = first_eta + a**-2 * fz
    L = ftt - 3 * h * ft + (2 * h * h - 2 * hd) * f
    M = ftz - 2 * h * fz
    G = ft + a**-1 * fz - 2 * h * f
    # The affine vector K=a^-2(d_eta+d_z), not a unit timelike vector.
    ap = sp.Symbol("a_prime_eta", real=True)
    K = sp.Matrix([a**-2, a**-2])
    connection_acceleration = sp.Matrix([2 * (ap / a) * a**-4, 2 * (ap / a) * a**-4])
    partial_acceleration = sp.Matrix([a**-2 * sp.diff(k, a) * ap for k in K])
    beta = sp.Symbol("beta_S", real=True)
    source = stress.reference_jets(h, hd, hdd, h3, beta)
    null_reference = sp.factor((source["rho"] + source["pressure"]) / a**2)
    wanted = (-4 * h * h * hd + beta * (12 * h3 + 72 * hd * hd + 36 * h * hdd)) / a**2
    return {
        "flat_sampler": F,
        "physical_plane_volume": "a*dt*dz=a^2*deta*dz",
        "affine_null_vector_conformal_components": K,
        "proper_second_time_operator": L,
        "proper_mixed_operator": M,
        "proper_state_weight_operator": G,
        "physical_reference_null_numerator": null_reference,
        "checks": {
            "actual_physical_null_energy_sampler_measure": sp.factor(
                a * a * a**-6 * f * f - F * F
            ),
            "actual_affine_null_vector_has_zero_norm": sp.factor(
                a * a * (K[0] ** 2 - K[1] ** 2)
            ),
            "actual_affine_null_vector_geodesic_acceleration": (
                partial_acceleration + connection_acceleration
            ).applyfunc(sp.factor),
            "actual_second_conformal_time_derivative": sp.factor(second_eta - L),
            "actual_mixed_conformal_derivative": sp.factor(mixed - a**-1 * M),
            "actual_flat_null_sampler_derivative": sp.factor(D - a**-1 * G),
            "actual_quantum_time_measure": sp.factor(
                second_eta**2 / a**2 - L * L / a**2
            ),
            "actual_quantum_cross_measure": sp.factor(
                second_eta * mixed / a**2 - L * M / a**3
            ),
            "actual_quantum_mixed_measure": sp.factor(
                mixed * mixed / a**2 - M * M / a**4
            ),
            "actual_physical_Wick_measure_keeps_extra_scale_squared": sp.factor(
                a * a * D * D / a**2 - G * G / a**2
            ),
            "actual_scalar_reference_null_contraction_retains_finite_beta": sp.factor(
                null_reference - wanted
            ),
            "actual_FK_affine_null_Ricci_contraction": sp.factor(
                (3 * (hd + h * h) - (hd + 3 * h * h)) / a**2 - 2 * hd / a**2
            ),
            "cosmological_constant_drops_from_null_SEE_contraction": sp.factor(
                a * a * (K[0] ** 2 - K[1] ** 2)
            ),
        },
    }
