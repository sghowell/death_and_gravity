"""Physical-polarization cut in the leading matched local operator theory."""

from fractions import Fraction
from functools import cache

import sympy as sp

DEFAULT_INVARIANT = sp.Rational(3, 2)


def rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction, sp.Rational)):
        raise TypeError("Require a finite exact rational channel invariant")
    return sp.Rational(value)


def point(invariant=DEFAULT_INVARIANT):
    s = rational(invariant)
    if not 0 < s < 4:
        raise ValueError("Require the overlapping-cut interval 0<s<4")
    d = data()
    return {
        "s": s,
        "u": 4 - s,
        "s_channel_phase_space_nonzero": True,
        "normalized_forward_boundary_jump": sp.expand(
            d["normalized_forward_boundary_jump"].subs(d["s"], s)
        ),
        "scope": "Leading local-operator gauge bubble only. Initial massive scalar momenta are analytically continued for 0<s<4; this is not a physical scalar scattering cross section there.",
    }


@cache
def data():
    E = sp.symbols("positive_gauge_energy", positive=True)
    C = sp.symbols("real_matched_operator_coefficient", real=True)
    s = sp.symbols("real_channel_invariant", real=True)
    nu = sp.symbols("positive_subtraction_scale", positive=True)
    eta = sp.diag(1, -1, -1, -1)
    k = sp.Matrix([E, 0, 0, E])
    l = sp.Matrix([E, 0, 0, -E])
    dot = lambda v, w: (v.T * eta * w)[0]
    pol = [sp.Matrix([0, 1, 0, 0]), sp.Matrix([0, 0, 1, 0])]
    T = dot(k, l) * eta - l * k.T
    amplitudes = [
        sp.expand(8 * C * (dot(k, l) * dot(v, w) - dot(k, w) * dot(l, v)))
        for v in pol
        for w in pol
    ]
    physical_sum = sp.expand(sum(v * v for v in amplitudes))
    lorentz_sum = sp.trace(T * eta * T.T * eta)
    colors = 8
    full_sum = 32 * colors * C * C * s * s
    phase_space = (
        sp.integrate(sp.S.One, (sp.Symbol("cosine"), -1, 1))
        * 2
        * sp.pi
        / (32 * sp.pi**2)
    )
    absorptive = sp.factor(full_sum * phase_space / 4)
    coefficient = colors * C * C / sp.pi**2
    upper = -coefficient * (
        s * s * (sp.log(s / nu**2) - sp.I * sp.pi)
        + (4 - s) ** 2 * (sp.log((4 - s) / nu**2) + sp.I * sp.pi)
    )
    lower = -coefficient * (
        s * s * (sp.log(s / nu**2) + sp.I * sp.pi)
        + (4 - s) ** 2 * (sp.log((4 - s) / nu**2) - sp.I * sp.pi)
    )
    normalized_jump = sp.factor((upper - lower) / coefficient)
    checks = {
        "first_external_gauge_null": dot(k, k),
        "second_external_gauge_null": dot(l, l),
        "first_Ward_identity": k.T * eta * T,
        "second_Ward_identity": T * eta * l,
        "independent_Lorentz_polarization_sum": sp.expand(
            64 * C * C * lorentz_sum - physical_sum
        ),
        "explicit_transverse_polarization_sum": sp.expand(
            physical_sum - 32 * C * C * (4 * E * E) ** 2
        ),
        "massless_phase_space": phase_space - 1 / (8 * sp.pi),
        "optical_and_identical_particle_factors": sp.factor(
            absorptive - colors * C * C * s * s / sp.pi
        ),
        "forward_two_cut_jump": sp.expand(
            normalized_jump - 16 * sp.I * sp.pi * (s - 2)
        ),
        "crossing_reverses_jump": sp.expand(
            normalized_jump.subs(s, 4 - s) + normalized_jump
        ),
        "accidental_center_jump_zero": normalized_jump.subs(s, 2),
        "center_derivative_jump_nonzero_value": sp.diff(normalized_jump, s).subs(s, 2)
        - 16 * sp.I * sp.pi,
        "second_derivative_jump_zero": sp.diff(normalized_jump, s, 2),
        "small_transfer_t_squared_log_limit": sp.limit(
            s * s * sp.log(s), s, 0, dir="+"
        ),
    }
    return {
        "s": s,
        "coefficient": C,
        "gauge_group_dimension": colors,
        "quartic_vertex_normalization": 8 * C,
        "transverse_polarization_amplitudes": amplitudes,
        "summed_color_and_polarization_amplitude_squared": full_sum,
        "integrated_massless_two_body_phase_space": phase_space,
        "s_channel_Im_amplitude": absorptive,
        "forward_log_coefficient": coefficient,
        "upper_physical_boundary": upper,
        "lower_physical_boundary": lower,
        "normalized_forward_boundary_jump": normalized_jump,
        "crossed_cut_intervals": {
            "s_channel": "[0,+infinity)",
            "u_channel": "(-infinity,4] at t=0",
        },
        "original_scalar_first_sheet_disc": "|s-2|<=1",
        "first_full_parent_loop_order_for_this_cut": 3,
        "scope": "Exact leading local Phi^2 F^2 operator calculation: two one-loop matched vertices sewn through one massless gauge loop. Opposite i0 signs are retained. The second derivative of the jump vanishes at this operator order; this does NOT restore the original gapped dispersion hypotheses. No impossibility of choosing a local logarithm branch, finite-mass momentum remainder, confinement spectrum or original-model three-loop result is claimed.",
        "checks": checks,
    }


def bad_cases():
    bad = (
        True,
        False,
        1.0,
        sp.Float(1),
        "1",
        None,
        sp.oo,
        -sp.oo,
        sp.zoo,
        sp.I,
        sp.nan,
        sp.Symbol("s"),
    )
    return [(f"inexact_invariant_{i}", point, (v,)) for i, v in enumerate(bad)] + [
        (f"outside_overlap_{i}", point, (v,)) for i, v in enumerate((-1, 0, 4, 5))
    ]
