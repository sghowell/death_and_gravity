"""Complete one-loop local vacuum, tadpole and scalar two-point ledger."""

from functools import cache

import sympy as sp
from p8_vacuum_light_pole import kernel as scalar_parent


@cache
def data():
    L, G, M, g, m, Y, ell = sp.symbols(
        "MS_quartic MS_cubic MS_heavy_mass_squared MS_cubic_squared reference_fermion_mass MS_Yukawa_squared scale_logarithm",
        positive=True,
    )
    x, s, phi = sp.symbols("unit_parameter continued_s light_field", real=True)
    Ibar, T = sp.symbols("MS_dimensional_reference regulated_light_tadpole", real=True)
    Q, N = 16 * sp.pi**2, sp.Integer(6)
    A = x * (1 - x)
    Delta = x * M + 1 - x - A * s
    Delta1 = Delta.subs(s, 1)
    bubble = sp.Integral(ell - sp.log(Delta), (x, 0, 1)) / Q
    tadpole = -(ell + 1) / Q
    fsc = L * tadpole / 2 - g * bubble
    r_integrand = g * A / (Q * Delta1)
    r = sp.Integral(r_integrand, (x, 0, 1))
    fscUV = -(L / 2 + g) * Ibar / Q
    fFUV = 2 * N * Y * (6 * m * m - s) * Ibar / Q
    JH = -G * T / 2
    stationary_correction = -G * JH * phi**2 / (2 * M)
    reduced_mass = (L / 2 - g / (2 * M)) * T
    sc_vac = (
        -ell - sp.Rational(3, 2) + M * M * (sp.log(M) - ell - sp.Rational(3, 2))
    ) / (4 * Q)
    ferm_vac = 63 * m**4 / Q
    fsc1 = fsc.subs(s, 1)
    fF1, fp = sp.symbols(
        "finite_fermion_mass_anchor finite_fermion_slope_anchor", real=True
    )
    F = sp.Function("fermion_MS_kernel")
    P = sp.Function("scalar_MS_Pi")
    kappa = 1 + r - fp
    refmass = 1 + P(1) - fF1
    fR = F(s) - fF1 - (s - 1) * fp
    PiR = P(s) - P(1) - (s - 1) * r
    Gamma = (refmass - s - P(s) + F(s)) / kappa
    old = scalar_parent.data()
    z = sp.symbols("positive_mass_squared", positive=True)
    checks = {
        "on_shell_parameter_denominator": sp.factor(Delta1 - x * M - (1 - x) ** 2),
        "parameter_denominator_minimum_factorization": sp.factor(
            Delta1 - 1 - x * (M - 2 + x)
        ),
        "parameter_denominator_maximum_factorization": sp.factor(
            M - Delta1 - (1 - x) * (M - 1 + x)
        ),
        "same_frozen_scalar_parameter_denominator": sp.factor(
            Delta1.subs({M: old["heavy_mass_squared"], x: old["Feynman_parameter"]})
            - old["on_shell_Delta"]
        ),
        "same_frozen_scalar_momentum_kernel": sp.factor(
            (-sp.log(Delta)).subs(
                {
                    M: old["heavy_mass_squared"],
                    x: old["Feynman_parameter"],
                    s: old["Minkowski_invariant"],
                }
            )
            - old["self_energy_momentum_kernel"]
        ),
        "scalar_Euclidean_derivative_sign": sp.factor(
            sp.diff(sp.log(Delta), s) + A / Delta
        ),
        "finite_tadpole_H_one_point_cancellation": G * T / 2 + JH,
        "stationary_H_tadpole_and_counterterm_cancel": sp.factor(
            (reduced_mass + sp.diff(stationary_correction, phi, 2)).subs(G**2, g)
            - L * T / 2
        ),
        "whole_scalar_two_point_UV_reference": L * (-Ibar / Q) / 2
        - g * Ibar / Q
        - fscUV,
        "combined_UV_second_derivative_zero": sp.diff(fscUV + fFUV, s, 2),
        "whole_UV_slope_is_fermion_wavefunction_pole": sp.diff(fscUV + fFUV, s)
        + 2 * N * Y * Ibar / Q,
        "two_scalar_vacuum_normalization": sp.diff(
            z * z * (sp.log(z) - ell - sp.Rational(3, 2)) / (4 * Q), z
        )
        - z * (sp.log(z) - ell - 1) / (2 * Q),
        "massless_gauge_vacuum_limit": sp.limit(
            z * z * (sp.log(z) - sp.Rational(5, 6)), z, 0, dir="+"
        ),
        "vacuum_energy_fixed_all_species": -(sc_vac + ferm_vac) + sc_vac + ferm_vac,
        "full_canonical_Phi_kernel_reference_identity": sp.factor(
            Gamma - (1 - s) - (fR - PiR) / kappa
        ),
        "mass_anchor_identity": sp.factor(Gamma.subs(s, 1).subs(F(1), fF1)),
        "residue_anchor_identity": sp.factor(
            sp.diff(Gamma, s)
            .subs(s, 1)
            .subs({sp.diff(P(s), s).subs(s, 1): r, sp.diff(F(s), s).subs(s, 1): fp})
            + 1
        ),
    }
    return {
        "parameters": {"L": L, "G": G, "g": g, "M": M, "m": m, "Y": Y, "ell": ell},
        "scalar_parameter_denominator": Delta,
        "on_shell_scalar_parameter_denominator": Delta1,
        "finite_MS_light_tadpole": tadpole,
        "finite_MS_mixed_bubble": bubble,
        "complete_scalar_MS_Euclidean_increment": fsc,
        "positive_scalar_Pi_slope": r,
        "complete_scalar_MS_mass_anchor": fsc1,
        "entire_scalar_two_point_UV_reference": fscUV,
        "entire_fermion_two_point_UV_reference": fFUV,
        "heavy_one_point_counterterm": JH,
        "two_scalar_MS_vacuum_constant": sc_vac,
        "all_fourteen_fermion_MS_vacuum_constant": ferm_vac,
        "total_vacuum_energy_reference": -sc_vac - ferm_vac,
        "finite_scalar_reference_mass_parameter": refmass,
        "finite_canonical_kinetic_normalization": kappa,
        "complete_canonical_Phi_kernel": Gamma,
        "scope": "One-loop scalar and fermion local terms, with H one-point function and vacuum energy fixed to zero and the physical scalar pole fixed to one. The reference propagator inside perturbative scalar graphs has mass one; the large local mass counterterm is paired with its subgraph, not used as a negative free mass. Gauge/ghost one-loop vacuum terms are scaleless in this regulator, not a claim about nonperturbative gauge vacuum energy. No global quantum potential is asserted.",
        "checks": checks,
    }
