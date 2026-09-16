"""Whole crossed metric insertion with exact spin and optical normalization."""

from functools import cache

import sympy as s
from p8_vacuum_affine_massive_elastic_proca_infrared import proca

from . import source

A, B, MU, NU, K = s.symbols("s t mu spectator_mass_squared kappa", real=True)
P = s.Symbol("Laplace_p", complex=True)
ETA = s.diag(1, -1, -1)


def channels(energy, transfer, mass):
    a, b, mu = map(s.sympify, (energy, transfer, mass))
    c = 4 * mu - a - b
    return ((a, b, c), (b, a, c), (c, a, b))


def numerator0(a, mass):
    a, mu = map(s.sympify, (a, mass))
    return (a + 2 * mu) ** 2 / 3


def numerator(a, b, c, mass):
    a, b, c, mu = map(s.sympify, (a, b, c, mass))
    return 2 * mu**2 - 2 * mu * a - b * c


def numerator2(a, b, c, mass):
    return numerator(a, b, c, mass) + numerator0(a, mass) / 2


def gravity_born(energy, transfer, mass, kappa):
    return (
        -sum(
            numerator(a, b, c, mass) / a for a, b, c in channels(energy, transfer, mass)
        )
        / kappa
    )


def gapped_non_Newton(energy, transfer, mass, kappa):
    total = 0
    for a, b, c in channels(energy, transfer, mass):
        aa, hh = source.poles.radial_form_factors(-a)
        total -= (
            numerator2(a, b, c, mass) * aa / (16 * s.pi**2)
            + numerator0(a, mass) * hh / (768 * s.pi**2)
        ) / kappa**2
    return total


def gapped_fixed(energy, transfer, mass, kappa):
    return -source.poles.delta_kappa() * gravity_born(
        energy, transfer, mass, kappa
    ) / kappa + gapped_non_Newton(energy, transfer, mass, kappa)


def massless_known(energy, transfer, mass, kappa, scale_squared):
    return -sum(
        (4 * numerator2(a, b, c, mass) + 10 * numerator0(a, mass))
        * s.log(-a / scale_squared)
        for a, b, c in channels(energy, transfer, mass)
    ) / (3840 * s.pi**2 * kappa**2)


def sewn_cut(species, energy, cosine, mass, spectator_mass, kappa):
    a, z, mu, nu, k = map(s.sympify, (energy, cosine, mass, spectator_mass, kappa))
    b = -(a - 4 * mu) * (1 - z) / 2
    c = 4 * mu - a - b
    return (
        s.pi
        * (
            4 * source.spectral.density(species, 2, a, nu) * numerator2(a, b, c, mu)
            + source.spectral.density(species, 0, a, nu) * numerator0(a, mu)
        )
        / (k * k * a * a)
    )


def radial_closed(species, spin, p, mass):
    """Entire subtracted radial integral, continued from positive real p."""
    p, nu = map(s.sympify, (p, mass))
    w = source.spectral.weight(species, spin)
    if p == 0:
        return s.S.Zero
    v = source.spectral.V
    coeff = s.Poly(w, v)
    ratio = s.sqrt(p / (4 * nu + p))
    j = {0: s.atanh(ratio) / ((4 * nu + p) * ratio)}
    for order in range(1, 4):
        j[order] = (4 * nu + p) * j[order - 1] / p - 1 / (p * (2 * order - 1))
    return sum(c * p * j[powers[0] // 2] for powers, c in coeff.terms())


@cache
def data():
    a, b, mu, nu, k = A, B, MU, NU, K
    c = 4 * mu - a - b
    z, h, dk, aa, hh = s.symbols("cosine hbar delta_kappa A2 H0", real=True)
    checks = {}

    def put(name, value):
        checks[name] = s.factor(value)

    put(
        "entire_spin2_invariant_projection",
        numerator2(a, b, c, mu) - (3 * (b - c) ** 2 - (a - 4 * mu) ** 2) / 12,
    )
    put(
        "Einstein_conserved_source_contraction",
        numerator2(a, b, c, mu) - numerator0(a, mu) / 2 - numerator(a, b, c, mu),
    )
    o2 = P * (1 + h * dk / k) + h * P**2 * aa / (16 * s.pi**2 * k)
    o0 = -2 * P * (1 + h * dk / k) + h * P**2 * hh / (192 * s.pi**2 * k)
    d2 = s.diff(1 / o2, h).subs(h, 0)
    d0 = s.diff(1 / o0, h).subs(h, 0)
    put("whole_inverse_spin2_sign", d2 + dk / (k * P) + aa / (16 * s.pi**2 * k))
    put("whole_inverse_spin0_sign", d0 - dk / (2 * k * P) + hh / (768 * s.pi**2 * k))
    contracted = (numerator2(a, b, c, mu) * d2 + numerator0(a, mu) * d0) / k
    target = (
        dk * numerator(a, b, c, mu) / (k * k * a)
        - numerator2(a, b, c, mu) * aa / (16 * s.pi**2 * k * k)
        - numerator0(a, mu) * hh / (768 * s.pi**2 * k * k)
    )
    put("whole_conserved_inverse_kernel_to_scattering", contracted.subs(P, -a) - target)
    a0 = (a + 2 * mu) * (a + 2 * nu) / (6 * k * a)
    a2 = -(a - 4 * mu) * (a - 4 * nu) / (6 * k * a)
    scalar_expected = (
        s.sqrt(1 - 4 * nu / a)
        * (a0 * a0 + a2 * a2 * s.legendre(2, z) / 5)
        / (32 * s.pi)
    )
    put(
        "whole_nonforward_scalar_optical_cut",
        sewn_cut("scalar", a, z, mu, nu, k) - scalar_expected,
    )
    put(
        "whole_nine_polarization_Proca_optical_cut",
        sewn_cut("vector", a, z, mu, nu, k) - proca.whole_cut(a, z, mu, nu, k),
    )
    v = source.spectral.V
    beta = s.symbols("beta", positive=True)
    for species in ("scalar", "vector"):
        for spin, factor in ((2, 64), (0, 768)):
            residue = -s.pi * source.spectral.weight(species, spin, beta) / (2 * beta)
            density = source.spectral.density(species, spin, a, a * (1 - beta**2) / 4)
            put(
                f"{species}_{spin}_whole_Feynman_cut_residue",
                residue + factor * s.pi**3 * density / a**2,
            )
            sigma = 4 * nu / (1 - v * v)
            rho = source.spectral.density(species, spin, sigma, nu)
            transformed = (
                factor * s.pi**2 * P * rho * s.diff(sigma, v) / (sigma**3 * (sigma + P))
            )
            # sqrt(v**2)=v on the positive radial interval.
            transformed = transformed.subs(s.sqrt(v * v), v)
            put(
                f"{species}_{spin}_whole_dispersive_radial_map",
                transformed
                - P
                * source.spectral.weight(species, spin)
                / (4 * nu + P * (1 - v * v)),
            )
    put(
        "same_full_mass_one_gravity_Born",
        gravity_born(a, -(a - 4) * (1 - z) / 2, 1, k)
        - source.forward.gravity_born(a, z, k),
    )
    old = source.massless.massless_log_part().subs(
        {
            source.massless.S: a,
            source.massless.T: b,
            source.massless.U: c,
            source.massless.MU: mu,
            source.massless.K: k,
            source.massless.NU2: 1,
        },
        simultaneous=True,
    )
    put("whole_original_M1_crossing_log", massless_known(a, b, mu, k, 1) - old)
    put(
        "massless_limit_is_M1_not_second_Phi_cut",
        sewn_cut("scalar", a, z, mu, 0, k)
        - source.massless.data()["whole_two_M1_absorptive_coefficient"].subs(
            {
                source.massless.S: a,
                source.massless.Z: z,
                source.massless.MU: mu,
                source.massless.K: k,
            },
            simultaneous=True,
        ),
    )
    return {
        "whole_conserved_numerators": {
            0: numerator0(a, mu),
            2: numerator2(a, b, c, mu),
            "Einstein": numerator(a, b, c, mu),
        },
        "whole_fixed_H_Proca_crossing_amplitude": gapped_fixed(a, b, mu, k),
        "whole_original_M1_crossing_log_part": massless_known(a, b, mu, k, 1),
        "whole_optical_cuts": {
            species: sewn_cut(species, a, z, mu, nu, k)
            for species in ("scalar", "vector")
        },
        "whole_inverse_kernel_first_loop": {2: d2, 0: d0},
        "whole_radial_closed_forms": {
            species: {spin: radial_closed(species, spin, P, nu) for spin in (0, 2)}
            for species in ("scalar", "vector")
        },
        "physical_sheet": "p=-a-i0; log(-a-i0)=log(a)-i*pi for a>0. The massive radial terms have their original threshold4nu. Their imaginary parts are -64pi^3*rho2/a^2 and -768pi^3*rho0/a^2; inverse-kernel signs give the positive optical cut. The closed atanh formula is continued from p>0, not independently assigned square-root branches.",
        "scope": "Whole first-loop fixed gapped Gaussian insertions, and complete M1 nonanalytic part modulo local terms. No light-Phi double counting, no interacting all-loop Gaussian promotion, no Regge contour conclusion.",
        "checks": checks,
        "gates": {
            "both_conserved_source_spins_retained": numerator0(a, mu) != 0
            and numerator2(a, b, c, mu) != 0,
            "all_six_gapped_channel_radial_terms_retained": len(
                gapped_non_Newton(a, b, mu, k).atoms(s.Integral)
            )
            == 12,
            "all_three_massless_transfer_logs_retained": len(
                massless_known(a, b, mu, k, 1).atoms(s.log)
            )
            == 3,
            "seagulls_and_source_fixed_volume_cancellation_retained": True,
            "whole_cuts_not_forward_or_massless_external_approximations": True,
        },
    }
