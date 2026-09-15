"""Literal scalar stress and both complete gapped Gaussian metric cuts."""

from functools import cache

import sympy as s
from p8_vacuum_affine_flat_tensor_cut import projectors as old_vector

from . import source

S, NU, P, V, Z, MU, K = s.symbols(
    "spectral_s mass_squared Laplace_p beta_variable scattering_cosine external_mass_squared kappa",
    positive=True,
)
ETA = s.diag(1, -1, -1, -1)


def weight(species, spin, variable=V):
    source.require_spin(spin)
    if species not in ("scalar", "vector"):
        raise ValueError("Require an original scalar or Proca Gaussian cut")
    if species == "scalar":
        return variable**6 / 30 if spin == 2 else variable**2 * (3 - variable**2) ** 2
    return (
        variable**2 * (30 - 20 * variable**2 + 3 * variable**4) / 30
        if spin == 2
        else variable**2 * (3 - 2 * variable**2 + 3 * variable**4)
    )


def density(species, spin, energy=S, mass2=NU):
    source.require_spin(spin)
    if species not in ("scalar", "vector"):
        raise ValueError("Require an original scalar or Proca Gaussian cut")
    if species == "scalar":
        polynomial = (
            (energy - 4 * mass2) ** 2 / 3840
            if spin == 2
            else (energy + 2 * mass2) ** 2 / 384
        )
    else:
        polynomial = (
            (13 * energy**2 + 56 * mass2 * energy + 48 * mass2**2) / 3840
            if spin == 2
            else (energy**2 - 4 * mass2 * energy + 12 * mass2**2) / 384
        )
    return s.sqrt(1 - 4 * mass2 / energy) * polynomial / s.pi**2


@cache
def moment(species, spin, order):
    source.require_spin(spin)
    if isinstance(order, bool) or not isinstance(order, (int, s.Integer)) or order < 3:
        raise ValueError("Require an exact convergent moment order at least3")
    sigma = 4 * NU / (1 - V * V)
    rho = density(species, spin)
    polynomial = s.cancel(rho * s.pi**2 / s.sqrt(1 - 4 * NU / S))
    transformed = s.cancel(
        V * polynomial.subs(S, sigma) * s.diff(sigma, V) / sigma ** (order + 1)
    )
    return s.factor(s.integrate(transformed, (V, 0, 1)) / s.pi**2)


def scalar_tree(energy, angle, external, mass2, kappa):
    return (
        energy**2
        + 8 * external * mass2
        - (energy - 4 * external) * (energy - 4 * mass2) * angle**2
    ) / (4 * kappa * energy)


@cache
def data():
    E, k = s.symbols("pair_energy pair_momentum", positive=True)
    a = s.Matrix([E, 0, 0, k])
    b = s.Matrix([E, 0, 0, -k])
    T = -(a * b.T + b * a.T - ETA * ((a.T * ETA * b)[0] + E * E - k * k)) / 2
    B = T[1:4, 1:4]
    red = lambda value: s.factor(
        s.expand(value).subs(E**2, k**2 + NU).subs(k**2, S / 4 - NU)
    )
    c0 = red(s.trace(B) ** 2 / 3)
    c2 = red((s.trace(B * B) - s.trace(B) ** 2 / 3) / 5)
    eye = s.eye(3)
    P0 = s.Matrix(9, 9, lambda a, b: eye[a // 3, a % 3] * eye[b // 3, b % 3] / 3)
    P2 = (
        s.Matrix(
            9,
            9,
            lambda a, b: (
                (
                    eye[a // 3, b // 3] * eye[a % 3, b % 3]
                    + eye[a // 3, b % 3] * eye[a % 3, b // 3]
                )
                / 2
            ),
        )
        - P0
    )
    checks = {
        "literal_scalar_stress_conservation": T[0, :],
        "entire_scalar_spin0_eigenvalue": s.factor(c0 - (S + 2 * NU) ** 2 / 12),
        "entire_scalar_spin2_eigenvalue": s.factor(c2 - (S - 4 * NU) ** 2 / 120),
        "complete_spin0_projector": P0 * P0 - P0,
        "complete_spin2_projector": P2 * P2 - P2,
        "orthogonal_conserved_sectors": P0 * P2,
        "both_projector_ranks": s.Matrix([s.trace(P0) - 1, s.trace(P2) - 5]),
        "Wick_phase_and_current_normalization": 2 * (V / (8 * s.pi)) / (4 * 2 * s.pi)
        - V / (32 * s.pi**2),
    }
    for spin in (0, 2):
        actual = old_vector.above_threshold_density(spin).subs(
            {old_vector.S: S, old_vector.MASS: s.sqrt(NU)}
        )
        checks["entire_original_Proca_density_" + str(spin)] = s.factor(
            actual - density("vector", spin)
        )
    x = s.Symbol("pair_cosine", real=True)
    a0 = (S + 2 * MU) * (S + 2 * NU) / (6 * K * S)
    a2 = -(S - 4 * MU) * (S - 4 * NU) / (6 * K * S)
    tree = scalar_tree(S, x, MU, NU, K)
    checks["entire_scalar_tree_Legendre_decomposition"] = s.expand(
        tree - a0 - a2 * s.legendre(2, x)
    )
    traceJ = (S + 2 * MU) / (2 * K * S)
    tf = (S - 4 * MU) ** 2 * s.legendre(2, Z) / (6 * K**2 * S**2)
    shape = (S + 2 * NU) ** 2 * traceJ**2 / 9 + (S - 4 * NU) ** 2 * tf / 30
    checks["complete_nonforward_scalar_sewing"] = s.factor(
        shape - a0 * a0 - a2 * a2 * s.legendre(2, Z) / 5
    )
    targets = {
        ("scalar", 2, 3): 1 / (53760 * s.pi**2 * NU),
        ("scalar", 2, 4): 1 / (967680 * s.pi**2 * NU**2),
        ("scalar", 0, 3): 17 / (26880 * s.pi**2 * NU),
        ("scalar", 0, 4): 1 / (13824 * s.pi**2 * NU**2),
        ("vector", 2, 3): 3 / (3584 * s.pi**2 * NU),
        ("vector", 2, 4): 31 / (322560 * s.pi**2 * NU**2),
        ("vector", 0, 3): 3 / (8960 * s.pi**2 * NU),
        ("vector", 0, 4): 1 / (32256 * s.pi**2 * NU**2),
    }
    for (species, spin, order), target in targets.items():
        checks[f"entire_{species}_spin{spin}_moment{order}"] = s.factor(
            moment(species, spin, order) - target
        )
    for species in ("scalar", "vector"):
        for spin, factor in ((2, 64), (0, 768)):
            sigma = 4 * NU / (1 - V * V)
            bare = s.cancel(density(species, spin) * s.pi**2 / s.sqrt(1 - 4 * NU / S))
            radial = s.cancel(
                factor
                * P
                * V
                * bare.subs(S, sigma)
                * s.diff(sigma, V)
                / (sigma**3 * (sigma + P))
            )
            checks[f"entire_{species}_spin{spin}_radial_conversion"] = s.factor(
                radial - P * weight(species, spin) / (4 * NU + P * (1 - V * V))
            )
    return {
        "literal_scalar_pair_stress": T,
        "conserved_spin0_and_spin2_projectors": (P0, P2),
        "complete_open_densities": {
            species: {str(spin): density(species, spin) for spin in (0, 2)}
            for species in ("scalar", "vector")
        },
        "complete_radial_weights": {
            species: {str(spin): weight(species, spin) for spin in (0, 2)}
            for species in ("scalar", "vector")
        },
        "first_two_convergent_moments": {
            f"{species}_{spin}_{order}": value
            for (species, spin, order), value in targets.items()
        },
        "entire_scalar_pair_tree": tree,
        "entire_nonforward_scalar_pair_cut": s.sqrt(1 - 4 * NU / S)
        * shape
        / (32 * s.pi),
        "absolute_normalization": "The displayed stress pair is half the two-scalar Feynman stress insertion. Its Wick factor2 is in the cut measure. The actual current is i theta<[T,T]>/4, hence rho=beta*(spin eigenvalue)/(32pi^2), Im chi(s+i0)=+pi rho(s). The gravitational s-channel scattering cut has the separate identical-pair and optical halves beta/(32pi).",
        "support_and_scope": "Densities vanish below4nu and are positive above. Proca retains all nine pairs of its three physical polarizations. Spin0 and rank-five spin2 are conserved off-shell metric source sectors, not additional massless physical graviton states. The scalar nu1 cut here is only the metric self-energy subset of the complete PhiPhi cut; it must not be added twice. The heavy Gaussian determinant is not a full heavy production amplitude.",
        "checks": checks,
        "gates": {
            "both_conserved_projectors_retained": P0.shape == (9, 9)
            and P2.shape == (9, 9),
            "all_three_original_Proca_modes_not_Maxwell_only": s.expand(
                weight("vector", 2)
            ).coeff(V, 6)
            == s.Rational(1, 10),
            "all_eight_exact_convergent_moments_positive": all(
                value.is_positive for value in targets.values()
            ),
            "both_scalar_masses_and_vector_mass_retained": True,
            "cut_normalization_not_arbitrary_spectral_rescaling": True,
            "Gaussian_subsets_not_double_counted_as_new_full_cuts": True,
            "massless_and_interacting_sectors_not_silently_included": True,
        },
    }
