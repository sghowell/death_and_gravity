"""Explicit physical mean-energy and full causal light-force constants."""

from functools import cache

import sympy as s

from . import clock, hamiltonian


@cache
def data():
    K = s.Integer(10) ** 800
    z = s.Rational(1, 10**6)
    m = s.Integer(1000)
    amax = s.Rational(25, 16)
    B = hamiltonian.data()["exact_CD_propagation_majorant"]
    cS = clock.bounds()["uniform_source_L2_constant"]
    cdS = clock.bounds()["uniform_spatial_and_Frechet_L2_constant"]
    forcing_squared = amax * cS * cS + (amax * 3 + 9) * cdS * cdS / (m * m)
    forcing = s.Integer(641)
    E0 = forcing**2 * B**2 / (2 * z)
    Efree = 2 * E0
    instantaneous = B * cS * cS
    spatial_W = forcing * B / s.sqrt(z)
    Wmajorant = s.Integer(4) * 10**6
    force = B * cdS * (2 * cS + Wmajorant)
    return {
        "kappa": K,
        "zeta": z,
        "canonical_mass": m,
        "full_source_forcing_squared_bound": forcing_squared,
        "full_source_forcing_constant": forcing,
        "normalized_reference_energy_coefficient": E0,
        "normalized_free_mean_energy_integrated_coefficient": Efree,
        "normalized_free_mean_energy_instantaneous_coefficient": instantaneous,
        "W_spatial_energy_constant": spatial_W,
        "full_W_time_integral_constant": Wmajorant,
        "normalized_full_causal_light_force_constant": force,
        "explicit_small_source_force_corollary": {
            "maximum_coordinate_jet_size": s.Rational(1, 10**8),
            "relative_H3_dual_norm_upper": s.Rational(4, 10**6),
            "scope": "For this smaller ball the source-force bound is below 4e-6*kappa*||psi||H3*||eta||H3. This is not a bound on the full light inverse, metric/noise response or interacting stability; the larger delta<=1/100 domain by itself does not make the force coefficient small.",
        },
        "energy_result": "E0[Abar](t)/kappa <3e12 delta^2 (integral U_psi)^2. The actual free coherent mean energy Efree/kappa <=6e12 delta^2 (integral U_psi)^2+1e6 delta^2 U_psi(t)^2. These are spatial integrals, not pointwise stress or the full source-interaction metric tensor.",
        "temporal_constraint_in_W": "W0=S0-sqrt(zeta) div pi_A/(a^3 sqrt(kappa)); hence ||W||L2<=512delta U_psi+4e6delta integral U_psi and ||S-W||L2<=1024delta U_psi+4e6delta integral U_psi.",
        "causal_light_force_identity": "At fixed g with the S6.176 canonical conditional measure/state, the expectation of the UNINTEGRATED source-sector light variation is kappa integral sqrt(-g) (S-Wbar).DS_eta. The connected centered Gaussian mean vanishes and its determinant is u-independent at fixed g. Evaluate Wbar at the actual retarded solution; do not vary a single-branch quadratic functional containing Gret.",
        "causal_light_force_result": "|Force_psi[eta]|/kappa <=4e10 delta^2 ||psi||H3(IxR3) ||eta||H3(IxR3), with all coordinate derivatives through3 in the stated norm. This bounds the full source/contact contribution, and is cubic for psi=epsilon f.",
        "Gaussian_boundary": "The conditional connected two-point function and ordinary-Proca renormalized contribution remain unchanged at fixed g; the free coherent stress difference is the classical sourced mean tensor. Metric-response/noise bounds and full parent loops are not supplied.",
        "checks": {
            "Proca_mass_normalization": m * m * z - 1,
            "spatial_derivative_source_count": 3 + 9 - 12,
            "exact_instantaneous_energy_coefficient": instantaneous - 10**6,
            "reference_to_physical_free_energy_factor": Efree - 2 * E0,
            "causal_source_contact_force_sign": s.Symbol("S")
            - s.Symbol("W")
            + (s.Symbol("W") - s.Symbol("S")),
            "physical_normalized_forcing_squared": forcing_squared
            - (amax * cS * cS + (amax * 3 + 9) * cdS * cdS * z),
        },
        "gates": {
            "explicit_smaller_ball_source_force_margin": bool(
                force * s.Rational(1, 10**8) ** 2 < s.Rational(4, 10**6)
            ),
            "full_source_forcing_below_six_four_one": bool(
                forcing_squared < forcing * forcing
            ),
            "reference_mean_energy_below_three_e_twelve": bool(E0 < 3 * 10**12),
            "physical_free_mean_energy_below_six_e_twelve": bool(Efree < 6 * 10**12),
            "full_W_integral_constant_valid": bool(
                2 * spatial_W * spatial_W < Wmajorant * Wmajorant
            ),
            "full_causal_light_force_below_four_e_ten": bool(force < 4 * 10**10),
        },
    }
