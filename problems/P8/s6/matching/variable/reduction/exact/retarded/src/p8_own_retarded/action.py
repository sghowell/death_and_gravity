"""Independent physical-time normalization and original-g Euler readout.

The retarded own-f map is defined by an equation and initial conditions.
This module does not substitute a causal inverse into a single-copy action
and silently assume the resulting kernel is variational or symmetric.
"""
from fractions import Fraction
from functools import cache

import sympy as sp
from p8_exact_stationary import stationary as parent

X = sp.Symbol("inner_x", real=True)
M2, TAU, DELTA = sp.symbols("M_squared tau delta", positive=True)
C = sp.Symbol("c", positive=True)


@cache
def normalization():
    q, hidden, k, source = (sp.Function(name)(X) for name in ("physical_q", "hidden_Q", "kinetic_k", "interaction_s"))
    scale = TAU*sp.sqrt(DELTA)
    qx, hx = sp.diff(q, X), sp.diff(hidden, X)
    b, lapse, beta1 = sp.symbols("b N beta1", positive=True)
    Kf = M2*b**3/lapse
    physical_density = M2*(qx/scale)**2/4+Kf*(hx/scale)**2/4-beta1*b*(hidden-q)**2/2
    inner_density = sp.expand(scale*physical_density)
    normalized = qx**2/4+k*hx**2/4-source*(hidden-q)**2/4
    dictionary = {k: b**3/lapse, source: 2*TAU**2*DELTA*beta1*b/M2}
    own_euler = -2*(sp.diff(normalized, hidden)-sp.diff(sp.diff(normalized, hx), X))
    physical_euler = -2*(sp.diff(normalized, q)-sp.diff(sp.diff(normalized, qx), X))
    momentum = sp.Function("momentum_p")(X)
    matrix = sp.Matrix([[0, 1/k], [-source, 0]])
    forcing = sp.Matrix([0, source*q])
    return {"q": q, "hidden": hidden, "k": k, "s": source, "scale": scale,
            "b": b, "N": lapse, "beta1": beta1, "K_f": Kf,
            "physical_density": physical_density, "inner_density": inner_density,
            "normalized_density": normalized, "dictionary": dictionary,
            "own_Euler": sp.expand(own_euler), "physical_Euler": sp.expand(physical_euler),
            "first_order_matrix": matrix, "first_order_forcing": forcing,
            "momentum": momentum, "state": sp.Matrix([hidden, momentum]),
            "physical_Euler_conversion": M2/(TAU**2*DELTA)}


@cache
def readout():
    eta = sp.Symbol("pulse_amplitude", positive=True)
    hidden0 = sp.Symbol("hidden_at_observation", positive=True)
    cmax = sp.Rational(1251, 625)
    s0 = 128/C
    response_lower = sp.Rational(63, 1024)*eta
    normalized_lower = sp.factor((s0.subs(C, cmax))*response_lower)
    return {"eta": eta, "hidden0": hidden0, "c_max": cmax,
            "s_center": s0, "normalized_original_g_Euler_at_zero_germ": -s0*hidden0,
            "response_lower": response_lower, "normalized_Euler_magnitude_lower": normalized_lower,
            "margin_over_three_amplitudes": sp.factor(normalized_lower-3*eta),
            "physical_Euler_magnitude_lower": M2*normalized_lower/(TAU**2*DELTA),
            "original_g_Euler_convention": "E_g=-2 delta S_T/delta q; literal action derivative has half this magnitude",
            "is_a_conserved_matter_source_or_on_shell_physical_g_claim": False}


@cache
def identities():
    d, r = normalization(), readout()
    q, h, k, s = (d[name] for name in ("q", "hidden", "k", "s"))
    frozen = parent.tensors()
    old = parent.own_equations()
    # The primary action is already a literal four-coordinate audited action.
    # Independently replace its physical derivatives by inner derivatives.
    replacement = {parent.M**2: M2, old["b"]: d["b"], old["N"]: d["N"],
                   old["beta1"]: d["beta1"], frozen["q"]: q, frozen["Q_tensor"]: h,
                   sp.diff(frozen["q"], parent.TIME): sp.diff(q, X)/d["scale"],
                   sp.diff(frozen["Q_tensor"], parent.TIME): sp.diff(h, X)/d["scale"]}
    p = parent.profile()
    s_center_from_profile = 2*DELTA*p["B1"].subs({parent.V: 0, parent.C: C})*2
    values = {"literal_physical_action_bridge": frozen["density"].xreplace(replacement)-d["physical_density"],
              "physical_to_inner_action_measure": d["inner_density"]-M2*d["normalized_density"].subs(d["dictionary"])/d["scale"],
              "full_variable_coefficient_own_Euler": d["own_Euler"]-sp.diff(k*sp.diff(h, X), X)-s*(h-q),
              "original_physical_g_Euler": d["physical_Euler"]-sp.diff(q, X, 2)-s*(q-h),
              "momentum_coordinate_first_equation": k*(d["first_order_matrix"]*d["state"]+d["first_order_forcing"])[0]-d["momentum"],
              "momentum_coordinate_second_equation": (d["first_order_matrix"]*d["state"]+d["first_order_forcing"])[1]+s*(h-q),
              "first_order_trace_zero": sp.trace(d["first_order_matrix"]),
              "exact_center_s": s_center_from_profile.subs(DELTA, C-2)-r["s_center"],
              "literal_delta_max": r["c_max"]-2-sp.Rational(1, 625),
              "zero_germ_readout": d["physical_Euler"].subs({q: 0}).doit().subs({s: r["s_center"], h: r["hidden0"]})
                                      -r["normalized_original_g_Euler_at_zero_germ"]}
    return {name: sp.factor(sp.expand(sp.together(value).as_numer_denom()[0])) for name, value in values.items()}


def independent_volterra_constants():
    """Separately written rational interval inequalities, no Green engine import."""
    lo, hi, strength, length = Fraction(19, 5), Fraction(21, 5), Fraction(65), Fraction(1, 4)
    flux_lower = 1-strength*length**2/(2*lo)
    derivative_lower = flux_lower/hi
    ratio_lower = (1-strength*length**2/(6*lo))/hi
    area, separation = Fraction(1, 16), Fraction(1, 8)
    return {"positive_flux": flux_lower, "R_x": derivative_lower, "R_over_duration": ratio_lower,
            "R_x_margin": derivative_lower-Fraction(1, 10),
            "R_margin": ratio_lower-Fraction(3, 16),
            "Q0_lower_per_amplitude": 42*Fraction(3, 16)*separation*area,
            "Qx0_lower_per_amplitude": 42*Fraction(1, 10)*area,
            "proof_boundary": "constant inequalities in a variable-coefficient first-zero/Volterra argument, not a constant-frequency approximation"}
