"""Physical Proca stress, spatial-only lapse chain and exact conservation."""

from functools import cache

import sympy as sp
from p8_proca_physical_matter import center
from p8_proca_seven_modes import model

a, m, k, q = model.a, model.m, model.kvec, model.q
H = model.H
h = sp.Symbol("positive_clock_h", positive=True)
N = sp.Symbol("positive_unchanged_physical_lapse", positive=True)


@cache
def data():
    d = model.data()
    G = sp.diag(d["Proca_coordinate_Hessian"], d["Proca_momentum_Hessian"])
    R = G / a**3
    P = -sp.diff(G, a) / (3 * a * a)
    M = d["Proca_cartesian_generator"]
    Z = sp.zeros(3)
    pieces = {
        "electric": sp.diag(Z, m * m * sp.eye(3) / a**4),
        "longitudinal_temporal": sp.diag(Z, k * k.T / a**6),
        "spatial_mass": sp.diag(sp.eye(3) / a**2, Z),
        "magnetic": sp.diag((q * sp.eye(3) - k * k.T) / (m * m * a**4), Z),
    }
    # Only the spatial metric is rescaled. N itself is unchanged.
    e = ((h - 1 + N**-2) / h) ** (-sp.Rational(1, 4))
    physical = N * G.subs(a, e * a)
    force = sp.diff(physical, N).subs(N, 1) / a**3
    false_force = sp.diff(e * physical, N).subs(N, 1) / a**3
    expected = R - 3 * P / (2 * h)
    aligned = {a: 1, m: 1000, k[1]: 0, k[2]: 0, k[0]: sp.Symbol("k", real=True), h: 1}
    anchor_force = model.clean(force.subs(aligned, simultaneous=True))
    density = center.data()["matrices"]["rho"][8:14, 8:14]
    assert len(density.free_symbols) == 1
    center_k = next(iter(density.free_symbols))
    anchor_force = anchor_force.subs(aligned[k[0]], center_k)
    return {
        "physical_density_Hessian": R,
        "physical_isotropic_pressure_Hessian": P,
        "four_positive_density_pieces": pieces,
        "actual_lapse_force_Hessian": force,
        "actual_center_force_Hessian": anchor_force,
        "N_is_unchanged_by_the_spatial_hat_chart": True,
        "false_four_dimensional_conformal_extra_source": R / (2 * h),
        "checks": {
            "false_four_dimensional_rescaling_has_nonzero_extra_density_source": model.clean(
                false_force - force - R / (2 * h)
            ),
            "physical_density_is_sum_of_four_positive_pieces": model.clean(
                R - sum(pieces.values(), sp.zeros(6))
            ),
            "physical_pressure_has_actual_four_piece_weights": model.clean(
                P
                - pieces["electric"] / 3
                - pieces["longitudinal_temporal"]
                + pieces["spatial_mass"] / 3
                - pieces["magnetic"] / 3
            ),
            "physical_Proca_density_is_Hamiltonian_per_physical_volume": model.clean(
                a**3 * R - G
            ),
            "physical_pressure_is_full_spatial_scale_variation": model.clean(
                3 * a * a * P + sp.diff(G, a)
            ),
            "complete_Proca_stress_conservation_on_actual_time_dependent_generator": model.clean(
                a * H * sp.diff(R, a) + M.T * R + R * M + 3 * H * (R + P)
            ),
            "actual_spatial_only_lapse_force_keeps_pressure": model.clean(
                force - expected
            ),
            "all_six_center_induced_matter_channels_match_actual_stress_reconstruction": model.clean(
                density + sp.Rational(2, 405) * anchor_force
            ),
            "physical_density_and_pressure_have_symmetric_real_Hessians": model.clean(
                R - R.T + P - P.T
            ),
        },
    }
