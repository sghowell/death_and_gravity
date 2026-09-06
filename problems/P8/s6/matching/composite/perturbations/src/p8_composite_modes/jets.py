"""Analytic physical-time bounce jets; the parent scale and CD duration differ."""

from functools import cache

import sympy as sp
from p8_composite import bounce, reconstruction

from . import model as m
from . import shift


@cache
def derive():
    old = reconstruction.derive()
    y, rho, h = (old[key] for key in ("y", "rho", "h"))
    accel = sp.Symbol("A", real=True)
    initial = {y: 1, rho: sp.Rational(1, 2), h: 0}
    velocity = {y: old["yprime"].subs(initial), rho: old["rhoprime"].subs(initial), h: accel}

    def derivative_at_zero(value):
        return sp.simplify(sum(sp.diff(value, variable).subs(initial)*speed
                               for variable, speed in velocity.items()))

    ratio = old["Nf"]/old["Ng"]
    yp = velocity[y]
    cp = derivative_at_zero(ratio)
    u, y2, c2 = sp.symbols("u y_second c_second", real=True)
    yu, cu = 1+yp*u+y2*u**2/2, 1+cp*u+c2*u**2/2
    mu = yu*(yu-1)*(yu-cu)/(1+yu)
    mu2 = sp.simplify(sp.diff(mu, u, 2).subs(u, 0)/2)
    assignment = {m.alpha: 1, m.beta: 1, m.y: 1, m.c: 1,
                  m.rho: m.M4/2, m.p: m.M4/2,
                  **dict(zip(m.betas, (0, 0, 1, 0, 0), strict=True))}
    xi0 = sp.simplify(shift.derive()["Xi"].subs(assignment, simultaneous=True)/m.M4)
    mtau = sp.Symbol("m_tau", positive=True)
    return {"A": accel, "u": u, "m_tau": mtau,
            "X_prime": derivative_at_zero(old["X"]), "Y_prime": derivative_at_zero(old["Y"]),
            "y_prime": yp, "rho_prime": velocity[rho], "c_prime": cp,
            "log_c_over_y_prime": sp.simplify(cp-yp),
            "mu_quadratic": mu2, "Xi_initial": xi0,
            "physical_speed_quadratic": sp.simplify(mu2/xi0),
            "free_A": bounce.time_jets()["H_prime"],
            "CD_A": sp.diff(4*u/(mtau**2+u**2), u).subs(u, 0)}


def checks():
    d = derive()
    root = sp.sqrt(sp.Rational(7, 3))
    aa, ratio = d["A"], d["m_tau"]
    free = {aa: sp.Rational(7, 36)}
    frozen = {aa: 4}
    return {"g_null_equation_jet": d["X_prime"]+2,
            "f_null_equation_jet": d["Y_prime"]+2,
            "ratio_actual_clock": d["y_prime"]+root,
            "density_at_bounce": d["rho_prime"],
            "lapse_ratio_actual_clock": sp.simplify(d["c_prime"]+(5+12*aa)/(3*root)),
            "relative_lapse_ratio_actual_clock": sp.simplify(d["log_c_over_y_prime"]-2*(1-6*aa)/(3*root)),
            "generic_acceleration_mu": sp.simplify(d["mu_quadratic"]-(sp.Rational(1, 3)-2*aa)),
            "positive_finite_initial_Xi": d["Xi_initial"]-6,
            "generic_physical_speed": sp.simplify(d["physical_speed_quadratic"]-(1-6*aa)/18),
            "free_actual_acceleration": d["free_A"]-sp.Rational(7, 36),
            "free_negative_mu": d["mu_quadratic"].subs(free)+sp.Rational(1, 18),
            "free_negative_speed": d["physical_speed_quadratic"].subs(free)+sp.Rational(1, 108),
            "frozen_CD_negative_mu": d["mu_quadratic"].subs(frozen)+sp.Rational(23, 3),
            "frozen_CD_negative_speed": d["physical_speed_quadratic"].subs(frozen)+sp.Rational(23, 18),
            "independent_interaction_to_duration_scale": d["CD_A"]-4/ratio**2,
            "independent_scale_sign_threshold": sp.factor(d["mu_quadratic"].subs(aa, d["CD_A"])
                                                            -(ratio**2-24)/(3*ratio**2))}


def scale_checks():
    mass, interaction, duration = sp.symbols("M m tau", positive=True)
    hbar, nbar, mubar, xibar = sp.symbols("hbar nbar mubar Xibar", real=True)
    return {"Einstein_interaction_energy_units": sp.cancel(mass**2*(interaction*hbar)**2/(mass**2*interaction**2)-hbar**2),
            "canonical_scalar_energy_units": sp.cancel((mass*interaction*sp.sqrt(nbar))**2/(mass**2*interaction**2)-nbar),
            "TT_and_shift_stiffness_units": sp.cancel((mass**2*interaction**2*mubar)/(mass**2*interaction**2*xibar)-mubar/xibar),
            "CD_acceleration_in_interaction_clock": sp.cancel((4/duration**2)/interaction**2-4/(interaction*duration)**2)}


def negative_controls():
    d = derive()
    return {"free_frozen_negative_speed_coefficient": d["physical_speed_quadratic"].subs(d["A"], d["free_A"]),
            "CD_frozen_negative_speed_coefficient": d["physical_speed_quadratic"].subs(d["A"], 4),
            "independent_scale_reverses_leading_sign": d["physical_speed_quadratic"].subs(d["A"], sp.Rational(4, 25)),
            "threshold_does_not_determine_higher_jets": d["mu_quadratic"].subs(d["A"], sp.Rational(1, 6))}
