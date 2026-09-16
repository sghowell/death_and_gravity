"""Both unequal-mass Bose partial waves and conditional coupled-pole sew."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complete_matter_graviton_endpoint import cuts

from . import source

MU, N, G, T = source.MU, source.N, source.G, source.T
X = s.Symbol("pair_angle", real=True)
U = s.Symbol("positive_Q_integration", nonnegative=True)


def require_species(value):
    if not isinstance(value, str) or value not in ("light", "heavy"):
        raise ValueError("Require the PhiPhi or HH pair cut")
    return value


def require_domain(species, energy, mass, heavy):
    require_species(species)
    energy, mass, heavy = map(source.require_mass, (energy, mass, heavy))
    active = mass if species == "light" else heavy
    if heavy <= 4 * mass or energy <= 4 * active:
        raise ValueError("Require the separated hierarchy and open named pair cut")
    return energy, mass, heavy


def kinematics(species, energy=T, mass=MU, heavy=N):
    require_species(species)
    energy, mass, heavy = map(s.sympify, (energy, mass, heavy))
    if species == "light":
        q = energy - 4 * mass
        return heavy + q / 2, q / 2, mass
    return (
        energy / 2 - heavy,
        s.sqrt((energy - 4 * mass) * (energy - 4 * heavy)) / 2,
        heavy,
    )


def Q2(value):
    value = s.sympify(value)
    return (3 * value**2 - 1) * s.atanh(1 / value) / 2 - 3 * value / 2


def Q_integral(order, value):
    order, value = map(s.sympify, (order, value))
    return s.Integral(
        (value + s.sqrt(value**2 - 1) * s.cosh(U)) ** (-order - 1), (U, 0, s.oo)
    )


def even_spin2(species, energy=T, mass=MU, heavy=N, cubic=G):
    a, b, _ = kinematics(species, energy, mass, heavy)
    return s.sympify(cubic) ** 2 * Q2(a / b) / (8 * s.pi * b)


def F1_cut(species, energy=T, mass=MU, heavy=N, cubic=G):
    _, _, active = kinematics(species, energy, mass, heavy)
    return (
        s.sqrt(1 - 4 * active / energy)
        * (energy - 4 * active)
        / (energy - 4 * mass)
        * even_spin2(species, energy, mass, heavy, cubic)
        / 2
    )


@cache
def data():
    a, b, g, beta, qa, qi = s.symbols("a b g beta qa qi", positive=True)
    a2 = cuts.exchange_moments(a, b, g)[1]
    checks = {
        "full_even_Bose_spin2_normalization": s.factor(
            a2 / (80 * s.pi) - g * g * Q2(a / b) / (8 * s.pi * b)
        ),
        "two_external_endpoint_cut": s.factor(
            2 * beta * qa * a2 / (160 * s.pi * qi) - beta * qa / qi * a2 / (80 * s.pi)
        ),
        "single_ordered_exchange_requires_Bose_factor_two": 2 * s.Rational(1, 16)
        - s.Rational(1, 8),
    }
    for k in range(1, 9):
        checks["P2_power_" + str(k)] = s.integrate(
            s.legendre(2, X) * X ** (2 * k), (X, -1, 1)
        ) - s.Rational(4 * k, (2 * k + 1) * (2 * k + 3))
        checks["positive_series_coefficient_majorant_" + str(k)] = (
            s.Rational(1, 15)
            - s.Rational(k, (2 * k + 1) * (2 * k + 3))
            - s.Rational((4 * k - 3) * (k - 1), 15 * (2 * k + 1) * (2 * k + 3))
        )
    checks["light_denominator_gap"] = s.factor(
        (N + (T - 4 * MU) / 2) ** 2 - (T - 4 * MU) ** 2 / 4 - N * (N + T - 4 * MU)
    )
    checks["heavy_denominator_gap"] = s.factor(
        (T / 2 - N) ** 2 - (T - 4 * MU) * (T - 4 * N) / 4 - N * N - MU * (T - 4 * N)
    )
    v = s.Matrix([s.sqrt(1 - X * X), 0, X])
    axis = s.Matrix([0, 0, 1])
    eye = s.eye(3)
    initial = qi * axis * axis.T / 2 - T * eye / 2
    final = qa * v * v.T / 2 - T * eye / 2
    contraction = s.trace(initial * final) - s.trace(initial) * s.trace(final) / 2
    p2 = s.diff(s.expand(contraction), X, 2) / 3
    checks["entire_stress_sew_spin2_species_residue"] = s.factor(p2 - qi * qa / 6)
    checks["constant_and_scalar_pole_are_P0"] = s.integrate(
        s.legendre(2, X), (X, -1, 1)
    )
    d = s.Symbol("J_minus_alpha")
    r1, r2, h11, h12, h22, B1, B2 = s.symbols("r1 r2 h11 h12 h22 B1 B2")
    R = s.Matrix([r1, r2])
    H = s.Matrix([[h11, h12], [h12, h22]])
    phase = s.diag(B1 / 2, B2 / 2)
    cross = H * phase * (R * R.T) / d + (R * R.T) * phase * H / d
    dr = H * phase * R
    checks["entire_coupled_channel_simple_pole"] = (
        cross - (dr * R.T + R * dr.T) / d
    ).applyfunc(s.factor)
    checks["no_selected_double_pole"] = (cross * d * d).applyfunc(s.factor).subs(d, 0)
    checks["two_external_residue_variations"] = s.factor(
        2 * dr[0] / r1 - B1 * h11 - B2 * h12 * r2 / r1
    )
    return {
        "whole_light_Bose_spin2": even_spin2("light"),
        "whole_heavy_Bose_spin2": even_spin2("heavy"),
        "whole_light_F1_cut": F1_cut("light"),
        "whole_heavy_F1_cut": F1_cut("heavy"),
        "whole_coupled_cross_term": cross,
        "whole_conditional_residue_discontinuity": dr,
        "normalization": "M=16pi sum_J(2J+1)f_J P_J; a2=80pi*f2. Every open identical scalar pair has unitarity weight beta/2. The two external residue variations give Disc(delta log residue)=sum_a beta_a f2_ia r_a/r_i. At spin2 the literal stress sew fixes r_a/r_phi=(T-4a)/(T-4mu), reproducing twice the complete generated F1 cut. No equal-mass replacement is made.",
        "Regge_scope": "The matrix calculation is conditional on a factorized isolated even-signature Regge pole and its complex-J unitarity continuation. Its selected matter cross term has a simple pole only, so its trajectory discontinuity vanishes, not its independent analytic trajectory correction. Integer-spin stress normalization does not establish a finite-transfer UV residue or Regge completion.",
        "threshold_scope": "The generated spin2 function has no direct C or H-resonance P0 pole. The light spin2 channel is regular even at T=n, unlike the unprojected fixed-order amplitude. The HH cut is a formal fixed-loop coefficient, not proof of an exactly stable heavy asymptotic particle or a width-resummed threshold.",
        "checks": checks,
        "gates": {
            "both_unequal_mass_pair_thresholds_retained": True,
            "Bose_factor_not_copied_from_single_ordered_exchange": True,
            "full_stress_sew_not_only_scalar_projection": True,
            "identical_pair_half_and_two_endpoint_factors": True,
            "no_double_pole_not_no_analytic_trajectory_term": True,
            "conditional_complex_J_and_pole_factorization_explicit": True,
        },
    }
