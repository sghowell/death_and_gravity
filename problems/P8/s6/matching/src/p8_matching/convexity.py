"""Stationary elimination, pinned second-jet sign and a mixed-gap control."""

from functools import cache

import sympy as sp
from p8.signs import even_rational_positive
from p8_s5 import nonlinear_d as model
from p8_uv import vacuum


@cache
def stationary_checks():
    X, z0, V0 = sp.symbols("X z0 V0", real=True)
    h1, h2, j1, j2 = sp.symbols("h1 h2 j1 j2", real=True)
    k11, k12, k22 = sp.symbols("k11 k12 k22", real=True)
    K = sp.Matrix([[k11, k12], [k12, k22]])
    h, j = sp.Matrix([h1, h2]), sp.Matrix([j1, j2])
    lagrangian = z0*X/2-V0+X*(j.T*h)[0]-(h.T*K*h)[0]/2
    solution = sp.solve([sp.diff(lagrangian, hi) for hi in h], list(h))
    effective = sp.factor(lagrangian.subs(solution))
    curvature = sp.factor(sp.diff(effective, X, 2))
    positive_form = j1**2/k11+(k11*j2-k12*j1)**2/(k11*K.det())
    residuals = {"two_heavy_envelope": sp.cancel(curvature-(j.T*K.inv()*j)[0]),
                 "two_heavy_positive_form": sp.cancel(curvature-positive_form)}
    # A non-Gaussian coupling to X exercises the X-dependent stationary Hessian.
    mu, b, source, h0 = sp.symbols("mu b j h", real=True)
    Z = z0+2*source*h0+b*h0**2
    nonlinear = Z*X/2-V0-mu**2*h0**2/2
    root = sp.solve(sp.diff(nonlinear, h0), h0)[0]
    Hessian = -sp.diff(nonlinear, h0, 2)
    nonlin_P = sp.factor(nonlinear.subs(h0, root))
    nonlin_Pxx = sp.factor(sp.diff(nonlin_P, X, 2))
    current = sp.diff(Z, h0).subs(h0, root)/2
    residuals.update({"nonlinear_stationarity": sp.cancel(sp.diff(nonlinear, h0).subs(h0, root)),
                      "nonlinear_envelope": sp.cancel(nonlin_Pxx-current**2/Hessian),
                      "nonlinear_expected": sp.cancel(nonlin_Pxx-source**2*mu**4/(mu**2-b*X)**3)})
    # The S6.1 positive vacuum mediator is a member of this affine-X class.
    Z0, g = sp.symbols("Z0 g", real=True)
    simple = Z0*(1+g*h0)*X/2-mu**2*h0**2/2
    simple_root = sp.solve(sp.diff(simple, h0), h0)[0]
    simple_P = sp.factor(simple.subs(h0, simple_root))
    lam = sp.diff(simple_P, X, 2)/(2*Z0**2)
    residuals["vacuum_mediator_lambda"] = sp.cancel(lam-g**2/(8*mu**2))
    old_lam = vacuum.mediator**2/(8*vacuum.heavy_mass**2)
    residuals["prior_vacuum_bridge"] = sp.cancel(lam-old_lam.subs({vacuum.mediator: g, vacuum.heavy_mass: mu}))
    return {"residuals": residuals, "two_heavy_P": effective, "two_heavy_Pxx": curvature,
            "two_heavy_positive_form": positive_form, "nonlinear_root": root,
            "nonlinear_Hessian": Hessian, "nonlinear_Pxx": nonlin_Pxx,
            "vacuum_mediator_lambda": sp.factor(lam)}


@cache
def target_checks():
    u, X = model.u, model.X
    F = model.functions()["F"]
    Fxx = sp.factor(sp.diff(F, X, 2))
    numerator = 45*u**8+122*u**6+96*u**4-12*u**2+5
    square = 45*u**8+122*u**6+96*(u**2-sp.Rational(1, 16))**2+sp.Rational(37, 8)
    M, tau = sp.symbols("M tau", positive=True)
    physical_Fxx = sp.diff(M**2*F/tau**2, X, 2)
    return {"Fxx": Fxx, "positive_numerator_decomposition": square,
            "minus_Fxx_positive": even_rational_positive(-Fxx, u),
            "bounce_jets": [sp.diff(F, X, n).subs({u: 0, X: 1}) for n in range(3)],
            "physical_remainder_second_jet_magnitude_lower_bound": 10*M**2/tau**2,
            "residuals": {"Fxx_formula": sp.cancel(Fxx+2*numerator/(1+u**2)**8),
                          "positive_square": sp.expand(numerator-square),
                          "physical_bounce_jet": sp.cancel(physical_Fxx.subs(u, 0)+10*M**2/tau**2)}}


@cache
def mixing_control():
    omega, k2, kappa, c = sp.symbols("omega k2 kappa c", real=True)
    pi0, h0 = sp.symbols("pi h", real=True)
    # Direct E-L Fourier equations for L=1/2(dot pi^2-k2*pi^2+dot h^2
    # -(k2+kappa)*h^2)+c*h*dot pi, with exp(-i*omega*t).
    equations = [(k2-omega**2)*pi0-sp.I*c*omega*h0,
                 (k2+kappa-omega**2)*h0+sp.I*c*omega*pi0]
    matrix = sp.Matrix([[sp.diff(eq, field) for field in (pi0, h0)] for eq in equations])
    determinant = sp.expand(matrix.det())
    expected = (omega**2-k2)*(omega**2-k2-kappa)-c**2*omega**2
    velocity_pi, velocity_h = sp.symbols("vpi vh", real=True)
    kinetic = (velocity_pi**2+velocity_h**2)/2+c*h0*velocity_pi
    gap = kappa+c**2
    speed = kappa/gap
    y = sp.Symbol("y", real=True)
    polynomial = determinant.subs(omega**2, y)
    slope_residual = sp.expand(polynomial.subs(y, speed*k2)).coeff(k2, 1)
    roots = sp.solve(polynomial.subs({kappa: -1, c: sp.sqrt(2), k2: sp.Rational(1, 10)}), y)
    return {"determinant": determinant, "stationary_Hessian": kappa,
            "physical_gap_squared": gap, "low_k_speed_squared": speed,
            "negative_control": {"kappa": -1, "mixing_squared": 2, "gap_squared": 1,
                                 "k_squared": "1/10", "roots_omega_squared": list(map(str, roots)),
                                 "lower_root_negative": bool(roots[0] < 0),
                                 "warning": "IR-unstable control, not a healthy bounce or complete UV model"},
            "residuals": {"direct_equation_determinant": sp.expand(determinant-expected),
                          "zero_k_factor": sp.expand(polynomial.subs(k2, 0)-y*(y-gap)),
                          "low_k_slope": sp.cancel(slope_residual),
                          "bare_velocity_Hessian": sp.hessian(kinetic, (velocity_pi, velocity_h))-sp.eye(2)}}


def derive():
    stationary, target, mixing = stationary_checks(), target_checks(), mixing_control()
    exact = {**stationary["residuals"], **target["residuals"],
             **{key: value for key, value in mixing["residuals"].items() if key != "bare_velocity_Hessian"}}
    exact.update({f"bare_velocity_{i}{j}": mixing["residuals"]["bare_velocity_Hessian"][i, j]
                  for i in range(2) for j in range(2)})
    if any(sp.cancel(value) != 0 for value in exact.values()):
        raise ValueError("Stationary matching or mixing control failed")
    return {"stationary_examples": {key: str(value) for key, value in stationary.items() if key != "residuals"},
            "stationary_positive_domain": "K positive definite; in 2D k11>0 and det K>0; nonlinear mu^2-bX>0",
            "target_Fxx": str(target["Fxx"]),
            "target_positive_numerator_decomposition": str(target["positive_numerator_decomposition"]),
            "target_minus_Fxx_sign": target["minus_Fxx_positive"],
            "target_bounce_F_Fx_Fxx": list(map(str, target["bounce_jets"])),
            "physical_second_jet_remainder_minimum": str(target["physical_remainder_second_jet_magnitude_lower_bound"]),
            "relative_second_jet_remainder_minimum": "1",
            "mixed_gap_control": {key: value if isinstance(value, dict) else str(value)
                                  for key, value in mixing.items() if key != "residuals"},
            "exact_residuals": {key: "0" for key in exact},
            "scope": "Stable affine-X stationary matching in a fixed coefficient convention only; not DHOST positivity"}
