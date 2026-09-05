"""Exact tensor operator reduction in the unchanged physical metric variable.

beta=cC/M² is constant.  No scalar constraint reduction or full-loop
effective action is claimed.  The fourth-order operator is used only for
its exact algebraic residual, not as a well-posed low-energy prescription.
"""

from functools import cache
from itertools import pairwise

import sympy as sp

H, H1, H2, H3, H4 = sp.symbols("H Hdot Hddot H3 H4", real=True)
Q = sp.Symbol("q", positive=True)
BETA = sp.Symbol("beta", real=True)
A_SCALE = sp.Symbol("a", positive=True)
GAMMA, VELOCITY = sp.symbols("gamma velocity")
X = sp.Symbol("x", real=True)
H_JETS = (H, H1, H2, H3, H4)
B = H2-2*H*H1
F = 8*B
C = -16*H1*Q
J = 2*(H**2-H1)
W = Q-3*H1/2-9*H**2/4
G = C-3*H*F/2


def coefficient_derivative(expression):
    return sp.expand(sum(sp.diff(expression, jet)*following for jet, following in pairwise(H_JETS))
                     - 2*H*Q*sp.diff(expression, Q)+A_SCALE*H*sp.diff(expression, A_SCALE))


def phase_derivative(expression, friction=None, frequency=None):
    """Cosmic derivative with the selected second-order phase flow."""
    friction = 3*H+BETA*F if friction is None else friction
    frequency = Q+BETA*C if frequency is None else frequency
    return sp.expand(coefficient_derivative(expression)+VELOCITY*sp.diff(expression, GAMMA)
                     - (friction*VELOCITY+frequency*GAMMA)*sp.diff(expression, VELOCITY))


@cache
def operator_checks():
    """Independent fourth-jet identity, with no equation of motion imposed."""
    jets = sp.symbols("g0:5")

    def dt(expression):
        return sp.expand(coefficient_derivative(expression)
                         + sum(sp.diff(expression, jets[n])*jets[n+1] for n in range(4)))

    def conformal_flat_wave(expression):
        return A_SCALE**2*(dt(dt(expression))+H*dt(expression)+Q*expression)

    def l0(expression):
        return dt(dt(expression))+3*H*dt(expression)+Q*expression

    fourth = sp.cancel(conformal_flat_wave(conformal_flat_wave(jets[0]))/A_SCALE**4)
    factorized = l0(l0(jets[0]))+J*l0(jets[0])+4*H1*Q*jets[0]-2*B*jets[1]
    on_shell = {jets[2]: -3*H*jets[1]-Q*jets[0],
                jets[3]: (9*H**2-3*H1-Q)*jets[1]+5*H*Q*jets[0]}
    on_shell[jets[4]] = sp.expand(dt(on_shell[jets[3]]).subs(on_shell, simultaneous=True))
    return {"full_off_shell_operator_factorization": sp.expand(fourth-factorized),
            "fourth_operator_on_baseline_equation": sp.expand(fourth.subs(on_shell, simultaneous=True)
                                                               - 4*H1*Q*jets[0]+2*B*jets[1])}


@cache
def field_map_checks():
    """Off-shell action cancellation and branch-restricted physical map.

    gamma=y+2 beta E0[y]-8 beta H ydot is the full first-order off-shell
    redefinition.  gamma=y-8 beta H ydot is valid only modulo O(beta²)
    on the selected perturbative branch.  Neither renames the matter metric.
    """
    e0, velocity = sp.symbols("E0 v")
    action_variation_per_c_a3 = -e0*(2*e0-8*H*velocity)/4
    weyl_action_per_c_a3 = (e0-2*H*velocity)**2/2
    action_check = sp.expand(action_variation_per_c_a3+weyl_action_per_c_a3-2*H**2*velocity**2)
    y_friction = 3*H+32*BETA*H*H1
    y_frequency = Q-16*BETA*H**2*Q

    def dt(expression):
        return phase_derivative(expression, y_friction, y_frequency)

    physical = GAMMA-8*BETA*H*VELOCITY
    transformed = sp.expand(dt(dt(physical))+(3*H+BETA*F)*dt(physical)+(Q+BETA*C)*physical)
    return {"off_shell_action_redefinition": action_check,
            "branch_map_order_zero": transformed.coeff(BETA, 0),
            "branch_map_order_one": transformed.coeff(BETA, 1)}


@cache
def canonical_checks():
    y, p, acceleration = sp.symbols("Y P Yddot")
    physical = A_SCALE**(-sp.Rational(3, 2))*y

    def dt(expression):
        return coefficient_derivative(expression)+p*sp.diff(expression, y)+acceleration*sp.diff(expression, p)

    transformed = sp.expand(A_SCALE**sp.Rational(3, 2)
                            *(dt(dt(physical))+(3*H+BETA*F)*dt(physical)+(Q+BETA*C)*physical))
    return {"old_canonical_tensor_equation": sp.expand(transformed-acceleration-BETA*F*p-(W+BETA*G)*y)}


@cache
def residual_coefficients():
    """Coefficients of Bop L1 gamma on the exact chosen reduced flow.

    Efull=L0-4 beta a^-4 D², Lred=L0+beta L1,
    Efull[gamma_red]=4 beta² (L0+J)L1 gamma_red exactly.
    Values are grouped by extra beta power, then gamma/velocity.
    """
    l1 = F*VELOCITY+C*GAMMA
    result = sp.expand(phase_derivative(phase_derivative(l1))+3*H*phase_derivative(l1)+(Q+J)*l1)
    return {power: {"gamma": sp.factor(result.coeff(BETA, power).coeff(GAMMA)),
                    "velocity": sp.factor(result.coeff(BETA, power).coeff(VELOCITY))}
            for power in range(3)}


@cache
def residual_checks():
    dt = phase_derivative
    # a^-4 D² = (dt+3H)(dt+2H)(dt+H)dt plus the q terms;
    # compute it from two conformal-wave applications instead.
    def wave(expression):
        return A_SCALE**2*(dt(dt(expression))+H*dt(expression)+Q*expression)

    fourth = sp.cancel(wave(wave(GAMMA))/A_SCALE**4)
    e0 = dt(dt(GAMMA))+3*H*dt(GAMMA)+Q*GAMMA
    actual = sp.expand(e0-4*BETA*fourth)
    expected = 4*BETA**2*sum(BETA**n*(record["gamma"]*GAMMA+record["velocity"]*VELOCITY)
                              for n, record in residual_coefficients().items())
    return {"exact_fourth_order_residual_on_reduced_solution": sp.expand(actual-expected),
            "residual_order_zero": actual.coeff(BETA, 0),
            "residual_order_one": actual.coeff(BETA, 1)}


def compact_derivative(expression, dimension):
    return sp.expand((1-X**2)*sp.diff(expression, X)-dimension*X*expression)


def compact_hubble_jets(order=4):
    return [4*(-1)**n*sp.factorial(n)*sp.chebyshevt(n+1, X) for n in range(order+1)]


def compact_coefficients():
    return {"Hdot": 4-8*X**2, "B": -56*X+96*X**3,
            "friction_F": -448*X+768*X**3,
            "speed_shift_over_beta": -64+128*X**2,
            "canonical_mass_G": (-64+128*X**2)*Q+2688*X**2-4608*X**4,
            "canonical_G_non_q_part": 2688*X**2-4608*X**4}


def compact_checks():
    jets = compact_hubble_jets()
    compact = compact_coefficients()
    replace = dict(zip(H_JETS, jets))
    values = {"hubble_Chebyshev_derivative_"+str(n): compact_derivative(jets[n], n+1)-jets[n+1]
              for n in range(len(jets)-1)}
    values.update({"compact_Hdot": H1.subs(replace)-compact["Hdot"],
                   "compact_B": B.subs(replace)-compact["B"],
                   "compact_F": F.subs(replace)-compact["friction_F"],
                   "compact_speed_shift": (-16*H1).subs(replace)-compact["speed_shift_over_beta"],
                   "compact_G": G.subs(replace)-compact["canonical_mass_G"]})
    return {key: sp.expand(value) for key, value in values.items()}


def controls():
    compact = compact_coefficients()
    kappa, beta = sp.symbols("kappa beta")
    return {"flat_background_reduced_correction": (F+C).subs({H: 0, H1: 0, H2: 0}),
            "physical_speed_shift_at_bounce_per_beta_over_tau2": compact["speed_shift_over_beta"].subs(X, 0),
            "physical_speed_shift_at_u_sqrt3_per_beta_over_ell2": compact["speed_shift_over_beta"].subs(X, sp.sqrt(3)/2),
            "speed_shift_zero_at_u_one": compact["speed_shift_over_beta"].subs(X, 1/sp.sqrt(2)),
            "naive_y_speed_used_as_physical_at_bounce_error": sp.Integer(64),
            "omit_canonical_volume_mass_error_at_x_half": (3*H*F/2).subs(dict(zip(H_JETS, compact_hubble_jets()))).subs(X, sp.Rational(1, 2)),
            "wrong_real_TT_normalization_speed_error_at_bounce": sp.Integer(-64),
            "drop_off_shell_field_map_E0_term_error": sp.Symbol("E0")**2/2,
            "constant_frequency_speed_error_accumulates": sp.sin(kappa*(1+beta))-sp.sin(kappa)}
