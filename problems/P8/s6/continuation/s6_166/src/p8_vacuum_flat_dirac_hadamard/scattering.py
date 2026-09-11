"""Uniform one-particle scattering errors and the source convention dictionary."""

from functools import cache

import sympy as s
from p8_vacuum_flat_dirac_production import calibration as parent

from . import symbols


def norm_caps(amplitude, timescale, time_magnitude):
    d, tau, T = map(symbols.rational, (amplitude, timescale, time_magnitude))
    if d < 0 or tau <= 0 or T < tau:
        raise ValueError("Require Delta>=0, tau>0 and |T|>=tau")
    wave = d * tau**8 / (56 * T**7)
    return {
        "one_particle_Moller_operator_norm_error": wave,
        "transported_asymptotic_projector_norm_error": 2 * wave,
    }


def asymptotic_masses(mean_mass, amplitude):
    m, d = map(symbols.rational, (mean_mass, amplitude))
    if m <= 0 or d < 0 or d >= m:
        raise ValueError("Require positive mean mass and 0<=Delta<mean mass")
    return m - d, m + d


def matrix_norm_squared(matrix):
    return s.simplify(sum(s.conjugate(v) * v for v in matrix))


@cache
def data():
    zero = s.zeros(2)
    pauli = (s.Matrix([[0, 1], [1, 0]]), s.Matrix([[0, -s.I], [s.I, 0]]), s.diag(1, -1))
    B = s.diag(1, 1, -1, -1)
    physical = [B] + [
        zero.row_join(sigma).col_join((-sigma).row_join(zero)) for sigma in pauli
    ]
    source = [-s.I * physical[0]] + [s.I * g for g in physical[1:]]
    metric = s.diag(-1, 1, 1, 1)
    checks = {}
    for a in range(4):
        for b in range(a, 4):
            checks[f"source_Clifford_relation_{a}_{b}"] = matrix_norm_squared(
                source[a] * source[b]
                + source[b] * source[a]
                - 2 * metric[a, b] * s.eye(4)
            )
        checks[f"source_spinor_hermitian_form_{a}"] = matrix_norm_squared(
            source[a].H * B + B * source[a]
        )
    checks["positive_time_spinor_inner_product"] = matrix_norm_squared(
        s.I * B * source[0] - s.eye(4)
    )
    checks["source_time_Dirac_coefficient"] = matrix_norm_squared(
        -source[0] - s.I * physical[0]
    )
    for j in range(1, 4):
        checks[f"source_spatial_Dirac_coefficient_{j}"] = matrix_norm_squared(
            source[j] - s.I * physical[j]
        )
    momenta = s.symbols("p1:4", real=True)
    mass = s.Symbol("M", real=True)
    H = sum((B * physical[j + 1] * momenta[j] for j in range(3)), s.zeros(4)) + B * mass
    checks["physical_Hamiltonian_squared"] = matrix_norm_squared(
        H * H - (sum(p * p for p in momenta) + mass * mass) * s.eye(4)
    )
    # The paper writes partial_t-i H_source=0; its H_source is -H_physical.
    Hsource = (
        sum((-source[0] * source[j + 1] * momenta[j] for j in range(3)), s.zeros(4))
        - s.I * source[0] * mass
    )
    checks["opposite_source_Hamiltonian_label"] = matrix_norm_squared(Hsource + H)
    p = parent.data()
    d, tau, m = (
        p["profile_amplitude_rational_upper"],
        p["same_profile_time_scale"],
        p["quadratic_subsystem_mean_mass"],
    )
    caps = norm_caps(d, tau, 1)
    checks["all_momentum_mass_tail_integral_coefficient"] = s.Rational(
        1, 8 * 7
    ) - s.Rational(1, 56)
    checks["projector_two_factor_bound"] = 2 * s.Rational(1, 56) - s.Rational(1, 28)
    checks["same_fixed_mass_and_switching_scale"] = d * tau**8 - s.Rational(3, 10**603)
    return {
        "physical_gamma": physical,
        "source_frame_gamma": source,
        "source_spinor_form": B,
        "physical_Hamiltonian": H,
        "source_Hamiltonian": Hsource,
        "source_mass_dictionary": "Set the paper's real mass to -M(t). Its (-,+,+,+) frame gamma matrices are Gamma_0=-i gamma_+^0 and Gamma_j=i gamma_+^j. Its slash plus mass equals the physical i gamma_+^mu partial_mu-M. H_source=-H_physical, so the spectral +/- labels exchange. This is a convention dictionary, not a negative physical mass or a different state.",
        "operator_norm_proof": "W(T)=U(0,T)exp(-i H_asym T). The derivative is i U(0,T)(H(T)-H_asym)exp(-i H_asym T), of norm <=|M(T)-M_asym|. Its integrable tail makes W and W* norm Cauchy; the limit is unitary on the one-particle L2 space. Transported spectral projectors differ by at most twice the wave-operator error.",
        "same_positive_physical_asymptotic_masses": asymptotic_masses(m, d),
        "actual_unit_time_norm_enclosures": caps,
        "norm_scope": "Uniform one-particle L2 operator norm, not a second-quantized global Fock implementer, a local Sobolev/Hadamard remainder seminorm or a numerical stress bound.",
        "checks": checks,
        "gates": {
            "same_asymptotic_physical_gap_positive": bool(m - d > 0),
            "actual_unit_time_in_tail_domain": bool(1 >= tau),
            "actual_one_particle_wave_error_below_one_e_minus_604": bool(
                caps["one_particle_Moller_operator_norm_error"] < s.Rational(1, 10**604)
            ),
            "actual_projector_error_below_two_e_minus_604": bool(
                caps["transported_asymptotic_projector_norm_error"]
                < s.Rational(2, 10**604)
            ),
        },
    }
