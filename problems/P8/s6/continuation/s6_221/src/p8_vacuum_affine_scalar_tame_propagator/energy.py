"""Exact two-chart energy identity and polynomial-momentum phase propagator."""

from functools import cache

import sympy as s

from . import charts as c
from .majorants import ENTRY, LOWER, UPPER

CONVERSION = s.Integer(10) ** 12
ENERGY_RATE = s.Integer(10) ** 28
LOG_PROPAGATOR = s.Integer(10) ** 29
SPATIAL_LOSS = 12


def symmetric(name):
    a, b, d = s.symbols(name + "0:3", real=True)
    return s.Matrix([[a, b], [b, d]])


@cache
def data():
    y = s.Matrix(s.symbols("y0:2", real=True))
    dy = s.Matrix(s.symbols("dy0:2", real=True))
    K, K1, G, G1 = [symmetric(name) for name in ("K", "Kd", "G", "Gd")]
    R = s.Matrix(2, 2, s.symbols("R0:4", real=True))
    omega, H, q = s.symbols("omega H q", real=True)
    skew = s.Matrix([[0, omega], [-omega, 0]])
    Kacc = -(K1 + 3 * H * K + skew) * dy - (q * G + R) * y
    derivative = (dy.T * Kacc)[0] + (dy.T * K1 * dy)[0] / 2 + q * (y.T * G * dy)[0]
    derivative += q * (y.T * (G1 - 2 * H * G) * y)[0] / 2
    target = -(dy.T * K1 * dy)[0] / 2 - 3 * H * (dy.T * K * dy)[0] - (dy.T * R * y)[0]
    target += q * (y.T * (G1 - 2 * H * G) * y)[0] / 2
    pv = s.Symbol("pv", real=True)
    lc = (
        c.th * pv
        - c.w * c.ps
        + 3 * c.th * c.l * c.sigma
        - (2 * c.E * c.q + 3 * c.T) * c.v
    )
    n = lc / (2 * c.J)
    Hphase = c.ham().subs(c.b, -pv / (2 * c.q))
    pv_dot = -3 * c.H * pv - s.diff(Hphase, c.v)
    b_dot = -pv_dot / (2 * c.q) - c.H * pv / c.q
    wanted_b = (
        -(1 + 9 * c.A / (2 * c.q)) * c.v
        - (c.E + 3 * c.T / (2 * c.q)) * n
        + c.H * pv / (2 * c.q)
    )
    central = c.central()
    yy = s.Matrix([c.b, c.sigma])
    ydot = s.Matrix([c.bd, c.sd])
    direct = central["inverse"] * (central["M"] * ydot + central["source"] * yy)
    viaK = (
        central["M"].T.inv()
        * central["K"]
        * (ydot + central["M"].inv() * central["source"] * yy)
    )
    norm = 2 * ENTRY
    growth = (norm + 2 * norm + norm + 4 * norm) / LOWER + 12
    small = s.Rational(1, 10**6)
    conversion = {
        "forward": 100 * (500 + 2),
        "outer_inverse": 5000 * 30000,
        "central_inverse": 2 * s.Integer(10) ** 7 * 30000,
    }
    source_bound = (2 + 300 * q) + 45 + 20 * q + 3
    qstar = s.Rational(152, 25)
    bounce = {
        c.E: s.Rational(-1, 2),
        c.l: s.Rational(1, 10),
        c.J: s.Rational(243, 160),
        c.A: 0,
        c.T: 0,
    }
    checks = {
        "general_complete_energy_identity_without_symmetric_R": s.expand(
            derivative - target
        ),
        "skew_velocity_work_zero": s.expand((dy.T * skew * dy)[0]),
        "central_b_velocity_from_full_weighted_canonical_equation": s.factor(
            b_dot - wanted_b
        ),
        "central_inverse_reconstruction_from_K": (direct - viaK).applyfunc(s.factor),
        "bare_finite_q_auxiliary_degeneracy_control": s.factor(
            central["XX"].det().subs(bounce).subs(c.q, qstar)
        ),
        "two_switch_total_lambda_exponent": 2 + 2 * 2 - 6,
        "phase_Sobolev_derivative_loss": 2 * 6 - SPATIAL_LOSS,
        "low_band_lambda": 1 + 100**2 - 10001,
    }
    return {
        "general_energy_derivative": target,
        "energy_growth_majorant": growth,
        "safe_root_energy_rate": ENERGY_RATE,
        "phase_energy_conversion_bounds": conversion,
        "common_conversion_constant": CONVERSION,
        "chart_transport_entry_sum_majorant": source_bound,
        "high_transfer_propagator": "1e72 exp(1e28 |t-s|)(1+|P|^2)^6 for |P|>=100; at most two chart switches",
        "all_transfer_propagator_log_constant": LOG_PROPAGATOR,
        "spatial_derivative_loss": SPATIAL_LOSS,
        "finite_q_degenerate_auxiliary_control": qstar,
        "energy_argument": "For both real and imaginary Fourier parts use E=(ydot^T K ydot+q y^T G y)/2. The exact skew term does no work. Complete coefficient bounds and q>=1 give |E'|<1e28 E and the safely enlarged root-energy Gronwall rate1e28.",
        "conversion_argument": "For Z=(v,sigma,pv,ps), lambda=1+|P|^2, both charts obey sqrt(E)<=1e12 lambda |Z| and |Z|<=1e12 lambda sqrt(E). Initial/final conversions and at most two switches give1e72 lambda^6. Low |P|<=100 uses the regular first-order bound exp[1000(10001)^2]. With log10<3, the common constant is exp(1e29).",
        "Duhamel_Sobolev_statement": "For any real r, the classical coefficient-sector phase system on the unit slab satisfies ||Z||Linf H^r <= exp(1e29)(||Z(initial)||H^(r+12)+||f||L1 H^(r+12)). The propagator is continuous on this scale, by pointwise ODE continuity and dominated convergence. The additive f is mathematical phase forcing, not an unproved map from physical quantum sources.",
        "strict_scope": "A polynomial-momentum classical comparison with12 lost spatial derivatives is established. The constant exp(1e29) is not small relative to kappa, the estimate is not same-space or a quantum contraction, and neither finite nonlinear stability nor the full quantum constraint inverse follows. A measure-zero P=0 assignment is not the literal homogeneous constraint.",
        "checks": checks,
        "gates": {
            "full_energy_growth_bound": growth < ENERGY_RATE,
            "all_phase_conversions_below_common": all(
                v < CONVERSION for v in conversion.values()
            ),
            "outer_forward_velocity_margin": 2 * 200 + s.Rational(1, 20) + 1 + 20 < 500,
            "central_forward_velocity_margin": 1
            + 9 * small / 2
            + (1 + 3 * small / 2) * 200
            + 2
            + 1
            + 20
            < 500,
            "full_central_transport_bound": all(
                v >= 0 for v in s.Poly(1000 * (1 + q) - source_bound, q).all_coeffs()
            ),
            "central_inverse_arithmetic": 2 * s.Integer(10) ** 7 > UPPER + 10**7 + 3,
            "energy_lower_to_vector_conversion": 2 * 4 / LOWER < 30000**2,
            "low_band_exponent": s.Integer(10) ** 3 * 10001**2 < s.Integer(10) ** 12,
            "log10_below_three_from_exp_series": 1
            + 3
            + s.Rational(9, 2)
            + s.Rational(27, 6)
            > 10,
            "full_log_constant_absorption": 72 * 3 + ENERGY_RATE < LOG_PROPAGATOR,
            "auxiliary_singularity_is_below_chart_threshold": 0 < qstar < 4096,
            "polynomial_not_uniform_same_space_inverse": True,
            "quantum_constraint_and_derivative_loss_still_open": True,
        },
    }
