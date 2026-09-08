"""Explicit prepared retarded estimate, not a full nonlinear EFT theorem."""
from functools import cache

import sympy as sp
from p8_affine_retuned.bounds import exact

ZETA_MAX = sp.Rational(1, 20000)
D_MIN = 1/ZETA_MAX-15


@cache
def identities():
    W = sp.Symbol("W", positive=True)
    Wp, Wpp = sp.symbols("W_prime W_second", real=True)
    J, Jp, Jpp, v, vp, F = sp.symbols("J J_prime J_second v v_prime force", real=True)
    Q = J/W

    def derivative(value):
        return sp.diff(value, W)*Wp+sp.diff(value, Wp)*Wpp+sp.diff(value, J)*Jp+sp.diff(value, Jp)*Jpp

    expected_second = Jpp/W-2*Jp*Wp/W**2-J*Wpp/W**2+2*J*Wp**2/W**3
    energy_prime = 2*vp*(F-W*v)+Wp*v**2+2*W*v*vp
    return {"quasistatic_second_derivative": sp.factor(derivative(derivative(Q))-expected_second),
            "exact_forced_energy_identity": sp.expand(energy_prime-2*vp*F-Wp*v**2),
            "zero_prepared_value": Q.subs(J, 0),
            "zero_prepared_derivative": derivative(Q).subs({J: 0, Jp: 0}),
            "second_derivative": expected_second}


@cache
def proof_checks():
    minimum = D_MIN
    amplification = minimum/(minimum-39)
    error_first = sp.Rational(21, 20)/141*(1+sp.Rational(868, 1)/minimum+sp.Rational(12168, 1)/minimum**2)
    error_total = 1/(100*ZETA_MAX*minimum)+16/minimum
    force_constant = (1+4*3+5*15+2*119)*sp.Rational(5, 4)
    spatial_constant = sp.Rational(101, 100)*force_constant/80
    temporal_constant = 1+sp.Rational(21, 20)*force_constant*sp.Rational(101, 100)*(1+sp.Rational(3, 141))
    return {"minimum_positive_squared_frequency": minimum == 19985,
            "square_root_floor_above_141": bool(minimum > 141**2),
            "energy_exponent_below_one": bool(0 < 39/minimum < 1),
            "energy_amplification_below_21_over_20": bool(amplification < sp.Rational(21, 20)),
            "quasistatic_residual_W_inverse_square": 2*78+712 == 868,
            "quasistatic_residual_W_inverse_cube": 2*78**2 == 12168,
            "first_error_below_one_over_100D": bool(error_first < sp.Rational(1, 100)),
            "total_error_below_zeta_A_over_80": bool(error_total < sp.Rational(1, 80)),
            "compact_scale_factor_upper": (1+sp.Rational(1, 2)**2)**2 == sp.Rational(25, 16),
            "square_root_one_plus_zeta_below_101_over_100": bool(1+ZETA_MAX < sp.Rational(101, 100)**2),
            "zeroth_force_jet_below_second_envelope": bool(1+2*3 < 326),
            "first_force_jet_below_second_envelope": bool(1+3*3+2*15 < 326),
            "second_force_jet_envelope": 1+4*3+5*15+2*119 == 326,
            "canonical_source_envelope": force_constant == sp.Rational(815, 2),
            "physical_spatial_error_below_six_zeta_E": bool(spatial_constant < 6),
            "physical_temporal_error_below_500_zeta_E": bool(temporal_constant < 500)}


def bound_values(coupling, left, right, comoving_squared, source_jet):
    zeta, start, end, k2, envelope = [exact(value, name) for value, name in
                                    ((coupling, "zeta"), (left, "left"), (right, "right"),
                                     (comoving_squared, "k_com_squared"), (source_jet, "E"))]
    if not bool(0 < zeta <= ZETA_MAX):
        raise ValueError("Require 0<zeta<=1/20000")
    if not bool(-sp.Rational(1, 2) <= start < end <= sp.Rational(1, 2)):
        raise ValueError("Require a nonempty interval inside [-1/2,1/2]")
    if not bool(0 <= k2 <= 1) or envelope.is_nonnegative is not True:
        raise ValueError("Require 0<=k_com_squared<=1 and nonnegative source jet envelope")
    homogeneous = k2.is_zero is True
    A = sp.Rational(815, 2)*sp.sqrt(zeta*k2)*envelope
    return {"zeta": zeta, "left": start, "right": end, "k_com_squared": k2, "source_jet_E": envelope,
            "branch": "homogeneous_algebraic_temporal_source" if homogeneous else "prepared_leading_longitudinal_response",
            "squared_frequency_floor": 1/zeta-15,
            "canonical_force_jet_A": A,
            "canonical_error_against_zeta_J": zeta*A/80,
            "spatial_readout_error_against_zeta_sqrt_q_times_Sprime_plus_2rhoS": 0 if homogeneous else 6*zeta*envelope,
            "temporal_readout_error_against_S": 0 if homogeneous else 500*zeta*envelope,
            "requires_zero_heavy_data_and_three_prepared_source_jets": True,
            "requires_three_time_derivative_source_norm_not_just_spatial_bound": True,
            "full_nonlinear_or_quantum_remainder_claim": False}


@cache
def checks():
    return {name: value for name, value in identities().items() if name != "second_derivative"}
