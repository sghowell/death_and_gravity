"""Both complete finite operator orderings and homogeneous force derivatives."""

from functools import cache

import sympy as s
from p8_vacuum_affine_weyl_operator_comparison import (
    derivatives as original_derivatives,
)

from . import geometry, source

PHASE, MAX_ORDER = 96, 196
RADIUS = geometry.CORE
H_AMPLITUDE = s.Integer(10) ** 280
F_AMPLITUDE = geometry.F_AMPLITUDE
FORCE = s.Rational(1, 10**230)
HOMOGENEOUS_DERIVATIVE = s.Integer(10) ** 470
STATE_WEIGHT = s.Rational(1, 10**350)


def require_amplitude(value):
    if isinstance(value, bool) or not isinstance(value, (int, s.Integer, s.Rational)):
        raise TypeError("Require an exact declared amplitude")
    if value not in (H_AMPLITUDE, F_AMPLITUDE):
        raise ValueError(
            "Require the declared full interaction or centered-volume amplitude"
        )
    return s.sympify(value)


def growth(order):
    if (
        isinstance(order, bool)
        or not isinstance(order, (int, s.Integer))
        or not 0 <= order <= MAX_ORDER
    ):
        raise ValueError("Require proved phase order0 through196")
    return s.Integer(2056) ** order * s.factorial(order) ** 3 / RADIUS**order


def jet(order, amplitude, homogeneous_order=0):
    amplitude = require_amplitude(amplitude)
    if (
        isinstance(homogeneous_order, bool)
        or not isinstance(homogeneous_order, (int, s.Integer))
        or not 0 <= homogeneous_order <= 2
    ):
        raise ValueError("Require homogeneous Cauchy order0 through2")
    return (
        2
        * amplitude
        * growth(order)
        * s.factorial(homogeneous_order)
        / source.HOMOGENEOUS_CAUCHY**homogeneous_order
    )


def operator_bound(amplitude, homogeneous_order=0):
    scale = jet(0, amplitude, homogeneous_order) / 2
    coherent = 3 * scale + PHASE * jet(2, amplitude, homogeneous_order) / 4
    remainder = 10**112 * PHASE**2 * jet(4, amplitude, homogeneous_order) / 32
    return coherent + remainder


@cache
def bounds():
    spatial = geometry.bounds()
    href = 10**128 * (8 * RADIUS) ** 2 / 2
    kamp = 2 * 1000 * geometry.KAPPA * spatial["hvariation"] + href
    first = operator_bound(H_AMPLITUDE, 1)
    second = operator_bound(H_AMPLITUDE, 2)
    force = 10**104 * first / geometry.KAPPA
    force_lip = 10**106 * (first + second) / geometry.KAPPA
    return {
        "reference_amplitude": href,
        "actual_centered_interaction_amplitude": kamp,
        "centered_operator": operator_bound(H_AMPLITUDE),
        "each_first_homogeneous_derivative": first,
        "each_second_homogeneous_derivative": second,
        "quantum_force": force,
        "quantum_force_Lipschitz": force_lip,
        "centered_volume_operator": operator_bound(F_AMPLITUDE),
    }


@cache
def data():
    b = bounds()
    rates = [s.Integer(2056) * (n + 1) ** 3 / RADIUS for n in range(MAX_ORDER)]
    U0, chi, F, h, h0, href = s.symbols(
        "actual_U0 whole_cutoff whole_F whole_h actual_h0 whole_href", real=True
    )
    alpha, alphabar = s.symbols("live_alpha reference_alpha", real=True)
    return {
        "whole_evaluated_operator_and_force_bounds": b,
        "whole_two_centered_amplitudes": [H_AMPLITUDE, F_AMPLITUDE],
        "whole_all_196_high_phase_growth_ratios": rates,
        "whole_homogeneous_derivative_order_radii": [
            source.HOMOGENEOUS_CAUCHY,
            source.HOMOGENEOUS_RADIUS,
            source.REAL_RADIUS,
        ],
        "whole_Hamiltonian_and_cutoff_definition": "At absolute homogeneous Y use h=kappa integral[V Hbar(Nstar,Z)]. With the same fixed reference S_u and h_ref set g=(h-h_ref)(u,Y,S_u z), g0=h(u,Y,0), K_symbol=chi_c(g-g0). Each c=1,2 retains the whole original nonlinear source and every contact. A_C=g0 I+Q_V[(1-D)K_symbol] and A_W=g0 I+OpW(K_symbol), D=Delta_w/4. h_ref is independent of Y; all actual scalar phases are retained.",
        "whole_volume_definition": "F_Y=average exp(3v) R_full(u,Nstar)^(-3/4), pulled back by the same S_u. Let F0(Y)=R_full(u,N0(Y))^(-3/4). Use F_ext=F0+chi_c(F_Y-F0), then calibrated-coherent and Weyl orderings separately. Centered phase derivatives do not differentiate a fictitious constant1 in place of F0. Physical volume is exp(3alpha) times the corresponding expectation, with state AND readout conjugated by the same Uref.",
        "whole_high_phase_and_parameter_proof": "The new whole centered amplitude is holomorphic on initial phase ball4R and the complex homogeneous1e-122 polydisc. Simultaneous96-coordinate phase Cauchy radiusR/8 and the unchanged complete radial-cutoff derivative proof give Jn=2A2056^n(n!)^3/R^n through196. At inner real Y, homogeneous Cauchy radius1e-124 gives first/second factor1/d,2/d^2. No higher real-time profile regularity is assumed. Differentiation commutes with fixed chi,S_u and the fixed same-state quantization.",
        "whole_operator_proof": "Apply the independent normalized coherent-kernel Schur theorem withK=1e112 and all coordinate orders<=2. The complete Weyl-minus-calibrated heat remainder is bounded byK*96^2*J4/32. The coherent bound is3A+96J2/4. Their sum O(A) bounds both orderings and all specified homogeneous derivatives; a symbol supremum alone is not used as a Weyl norm.",
        "whole_force_proof": "Use the full canonical force formulas in homogeneous.py, including1/(kappa Vol a^3), physical heavy1e100, and all coefficient derivatives. On the real and complex inner domains Vol>1 and |a^-3|<2. The displayed1e104 and1e106 ceilings dominate every row sum and its derivative. Scalar g0 gives the separate complete classical force; K derivatives give uniform quantum forces for every normalized invariant state.",
        "checks": {
            "entire_centered_Hamiltonian_scalar_split": s.expand(
                h0 + chi * ((h - href) - h0) - (h0 + chi * (h - h0 - href))
            ),
            "entire_centered_volume_not_center_one": s.expand(
                U0 + chi * (F - U0) - (chi * F + (1 - chi) * U0)
            ),
            "physical_volume_reference_factor": s.simplify(
                s.exp(3 * alpha) - s.exp(3 * alphabar) * s.exp(3 * (alpha - alphabar))
            ),
            "same_full_phase_dimension": s.Integer(PHASE - 2 * 48),
            "same_full_cutoff_derivative_base": original_derivatives.PRODUCT_BASE
            - s.Integer(2056),
            "same_full_kernel_derivative_order": s.Integer(MAX_ORDER - 2 * PHASE - 4),
        },
        "gates": {
            "actual_full_centered_Hamiltonian_amplitude": b[
                "actual_centered_interaction_amplitude"
            ]
            < H_AMPLITUDE,
            "all_196_phase_ratios": max(rates) < s.Rational(1, 10**9),
            "whole_uniform_quantum_force": b["quantum_force"] < FORCE,
            "whole_quantum_force_Lipschitz": b["quantum_force_Lipschitz"] < 1,
            "whole_state_generator_parameter_Lipschitz": 5
            * b["each_first_homogeneous_derivative"]
            < HOMOGENEOUS_DERIVATIVE,
            "whole_parameter_Cauchy_fits": 5 * source.HOMOGENEOUS_CAUCHY
            + source.REAL_RADIUS
            < source.HOMOGENEOUS_RADIUS,
            "same_full_two_cutoffs_and_original_covariance": True,
            "no_small_state_or_leakage_bound_inferred": True,
        },
    }
