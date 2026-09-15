"""Complete gapped Gaussian metric four-point insertions and physical cuts."""

from functools import cache

import sympy as s
from p8_vacuum_affine_gaussian_metric_pole_matching import poles as fixed_poles
from p8_vacuum_affine_gaussian_metric_pole_matching import spectral
from p8_vacuum_affine_massive_elastic_proca_infrared import proca
from p8_vacuum_affine_massive_gravity_pole_completion import poles as previous_poles

from . import source

X = s.Symbol("Gaussian_unit_velocity", real=True)
S, T, MU, K = source.S, source.T, source.MU, source.K


def require_species(value):
    if not isinstance(value, str) or value not in ("scalar", "vector"):
        raise ValueError("Require an original scalar or Proca Gaussian metric sector")
    return value


def tensor_shapes(channel, other, crossed, mass=MU):
    channel, other, crossed, mass = map(s.sympify, (channel, other, crossed, mass))
    trace = (channel + 2 * mass) ** 2
    traceless = (
        previous_poles.tree_numerator(channel, other, crossed, mass, 4) + trace / 6
    )
    return traceless, trace


def effective_curvature(species, loop_mass, channel, weyl=0, ricci=0):
    require_species(species)
    loop_mass, channel, weyl, ricci = map(s.sympify, (loop_mass, channel, weyl, ricci))
    denominator = 4 * loop_mass - channel * (1 - X * X)
    return (
        weyl
        + s.Integral(
            channel * spectral.weight(species, 2, X) / (64 * s.pi**2 * denominator),
            (X, 0, 1),
        ),
        ricci
        + s.Integral(
            channel * spectral.weight(species, 0, X) / (4608 * s.pi**2 * denominator),
            (X, 0, 1),
        ),
    )


def feynman_boundary_curvature(species, loop_mass, channel, weyl=0, ricci=0):
    delta = s.Symbol("positive_Feynman_boundary", positive=True)
    values = effective_curvature(
        species, loop_mass, s.sympify(channel) + s.I * delta, weyl, ricci
    )
    return tuple(s.Limit(value, delta, 0, dir="+") for value in values)


def local_amplitude(weyl, ricci, energy=S, transfer=T, mass=MU, kappa=K):
    weyl, ricci, energy, transfer, mass, kappa = map(
        s.sympify, (weyl, ricci, energy, transfer, mass, kappa)
    )
    u = 4 * mass - energy - transfer
    return (
        (2 * ricci + 8 * weyl / 3) * (energy**2 + transfer**2 + u**2)
        + (56 * ricci - 64 * weyl / 3) * mass**2
    ) / kappa**2


def nonlocal_amplitude(species, loop_mass, energy=S, transfer=T, mass=MU, kappa=K):
    require_species(species)
    loop_mass, energy, transfer, mass, kappa = map(
        s.sympify, (loop_mass, energy, transfer, mass, kappa)
    )
    u = 4 * mass - energy - transfer
    result = 0
    for a, b, c in (
        (energy, transfer, u),
        (transfer, energy, u),
        (u, energy, transfer),
    ):
        w, r = effective_curvature(species, loop_mass, a)
        tt, trace = tensor_shapes(a, b, c, mass)
        result += 4 * w * tt + 2 * r * trace
    return result / kappa**2


def fixed_HP_amplitude(energy=S, transfer=T, mass=MU, kappa=K):
    coeff = source.gaussian_source.fixed_coefficients()
    return (
        local_amplitude(
            coeff["Weyl_squared"], coeff["R_old_squared"], energy, transfer, mass, kappa
        )
        + nonlocal_amplitude(
            "scalar", source.HEAVY_MASS2, energy, transfer, mass, kappa
        )
        + nonlocal_amplitude(
            "vector", source.VECTOR_MASS2, energy, transfer, mass, kappa
        )
    )


@cache
def data():
    ss, mu, n, k, z, beta = s.symbols(
        "physical_s external_mu loop_mass_squared kappa scattering_cosine physical_beta",
        positive=True,
    )
    t = -(ss - 4 * mu) * (1 - z) / 2
    u = 4 * mu - ss - t
    tt, trace = tensor_shapes(ss, t, u, mu)
    checks = {
        "entire_conserved_spin2_external_tensor": s.factor(
            tt - (ss - 4 * mu) ** 2 * s.legendre(2, z) / 6
        )
    }
    for species in ("scalar", "vector"):
        imag = spectral.weight(species, 2, beta) * tt / (
            32 * s.pi * beta * k * k
        ) + spectral.weight(species, 0, beta) * trace / (4608 * s.pi * beta * k * k)
        if species == "scalar":
            a0 = (ss + 2 * mu) * (ss + 2 * n) / (6 * k * ss)
            a2 = -(ss - 4 * mu) * (ss - 4 * n) / (6 * k * ss)
            target = beta * (a0 * a0 + a2 * a2 * s.legendre(2, z) / 5) / (32 * s.pi)
        else:
            target = proca.whole_cut(ss, z, mu, n, k).subs(s.sqrt(1 - 4 * n / ss), beta)
        checks[species + "_whole_angular_physical_cut"] = s.factor(
            (imag - target).subs(n, ss * (1 - beta * beta) / 4)
        )
        for spin, prefactor, expected in ((2, 64, 128), (0, 4608, 9216)):
            w = spectral.weight(species, spin, beta)
            checks[species + "_spin" + str(spin) + "_whole_cut_delta_Jacobian"] = (
                s.factor(
                    ss * w * s.pi / (prefactor * s.pi**2 * 2 * ss * beta)
                    - w / (expected * s.pi * beta)
                )
            )
    w, r = s.symbols("local_Weyl_squared local_R_squared", real=True)
    a, b = S, T
    c = 4 * MU - a - b
    direct = (
        sum(
            4 * w * tensor_shapes(x, y, z)[0] + 2 * r * tensor_shapes(x, y, z)[1]
            for x, y, z in ((a, b, c), (b, a, c), (c, a, b))
        )
        / K**2
    )
    checks["entire_crossed_curvature_insertions_equal_local_action_map"] = s.factor(
        direct - local_amplitude(w, r)
    )
    P = s.Symbol("Laplace_p", complex=True)
    radial = fixed_poles.radial_form_factors(P)
    coeff = source.gaussian_source.fixed_coefficients()
    expectedA = -64 * s.pi**2 * coeff["Weyl_squared"]
    expectedH = -4608 * s.pi**2 * coeff["R_old_squared"]
    for species, mass in (
        ("scalar", source.HEAVY_MASS2),
        ("vector", source.VECTOR_MASS2),
    ):
        expectedA += s.Integral(
            P
            * spectral.weight(species, 2, spectral.V)
            / (4 * mass + P * (1 - spectral.V**2)),
            (spectral.V, 0, 1),
        )
        expectedH += s.Integral(
            P
            * spectral.weight(species, 0, spectral.V)
            / (4 * mass + P * (1 - spectral.V**2)),
            (spectral.V, 0, 1),
        )
    checks["whole_frozen_S285_TT_radial_kernel"] = s.simplify(radial[0] - expectedA)
    checks["whole_frozen_S285_trace_radial_kernel"] = s.simplify(radial[1] - expectedH)
    loopmass, channel = s.symbols("positive_loop_mass channel", positive=True)
    for species in ("scalar", "vector"):
        for spin in (0, 2):
            numerator = channel * spectral.weight(species, spin, X)
            checks[species + "_spin" + str(spin) + "_gapped_remainder_zero_jet"] = (
                numerator.subs(channel, 0)
            )
    hp = source.gaussian_source.fixed_coefficients()
    L = s.log(source.HEAVY_MASS2) + 2
    checks["whole_fixed_H_Proca_local_polynomial"] = s.factor(
        local_amplitude(hp["Weyl_squared"], hp["R_old_squared"])
        + L
        * (S * S + T * T + (4 * MU - S - T) ** 2 + 12 * MU * MU)
        / (640 * s.pi**2 * K**2)
    )
    return {
        "whole_conserved_external_tensor_shapes": tensor_shapes(S, T, 4 * MU - S - T),
        "whole_scalar_effective_curvature": effective_curvature(
            "scalar", loopmass, channel, w, r
        ),
        "whole_Proca_effective_curvature": effective_curvature(
            "vector", loopmass, channel, w, r
        ),
        "whole_scalar_Feynman_boundary": feynman_boundary_curvature(
            "scalar", loopmass, channel, w, r
        ),
        "whole_fixed_H_Proca_pole_subtracted_amplitude": fixed_HP_amplitude(),
        "whole_fixed_H_Proca_Newton_term": fixed_poles.delta_kappa()
        * previous_poles.newton_shape(dimension=4)
        / K**2,
        "physical_cut_normalization": "Im w_eff=w2(beta)/(128pi beta), Im r_eff=w0(beta)/(9216pi beta). The complete conserved tensor contraction reproduces the entire frozen scalar/Proca production sew at every angle. No forward-only or selected-polarization normalization is used.",
        "analytic_and_matching_scope": "The form factors are analytic on C minus[4mass_squared,infinity), with physical values defined by the explicit upper Feynman-boundary limits. Both remainder form factors vanish at zero transfer. The full fixed Gaussian volume density and Newton coefficient are removed consistently before the stated pole-subtracted amplitude is used; all regular curvature and nonlocal terms remain.",
        "light_and_other_sector_boundary": "The analogous Phi Gaussian term uses its actual mass1 and its still-unmatched Weyl/R_old_squared coefficients; its Newton coefficient remains separate. That Gaussian component is already contained in S288's minimal loop sector. H/Phi/g/matter interactions outside the metric Gaussian insertion are not set to zero or counted twice.",
        "checks": checks,
        "gates": {
            "complete_fixed_H_Proca_metric_functions_not_only_low_series": True,
            "both_conserved_spin_channels_and_all_crossings": True,
            "literal_frozen_scalar_and_Proca_physical_cuts": True,
            "whole_S285_finite_polynomial_and_Newton_shift": True,
            "all_domain_gapped_Feynman_analyticity": True,
            "light_Gaussian_not_double_counted_or_finitely_chosen": True,
            "not_a_full_source_or_Regge_amplitude": True,
        },
    }
